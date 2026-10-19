"""validate2.py -- independent end-to-end checks of the coset projection.

(1) Bosonisation at the eigenvalue level: for every neutral fermion word (psi^+- letters) of doubled weight <= 8 the NS
    cylinder vacuum eigenvalue computed by the registered fermionic evaluator (word rewritten in psi_0, psi_1 letters,
    metric (s, -s)) must equal the eigenvalue of its bosonic image (X_3 letters, metric 1) at X_3-momentum 0 (the NS
    vacuum of the complex fermion is the charge-0 Fock vacuum).  Also split-independence of the recursion.
(2) End to end: for each tower charge I at the fibre, the eigenvalue on |P_0, P_1> x (NS vacuum), computed by gm.vev_ns on
    the ORIGINAL two-superfield density, must equal sum_k q^k vev_uv(D_k)(P_u, P_v) of the projected density with
    P_u = P_0, P_v = -P_1, q = -(2/s) P_1 (the H~ momentum of that state; X_3 momentum 0).
Usage: validate2.py <t>
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
import gm
from super_engine import gen_basis_super, canonicalize
import n2frame as N2
import coset
import glob, json
from towers import cur_from_json


def load_charges(tw, sv):
    stag = str(sv).replace('/', '-')
    ch = {}
    for fn in glob.glob(f'towers_{tw}_s{stag}_w*.json'):
        if 'perturb' in fn:
            continue
        d = json.load(open(fn))
        for w2, r in d['results'].items():
            ch[int(w2)] = [cur_from_json(k) for k in r['kernel']]
    return ch
from coset import PP, PM

t = F(sys.argv[1])
s = 2 * (t - 1) / (t + 1)
Bz = coset.Bosoniser(s)
eta = N2.eta_of(s)


def pm_to_01(cur):
    """psi^+ = psi_0 + psi_1, psi^- = psi_0 - psi_1 (letters (1,a,2),(1,a,3) -> (1,a,0),(1,a,1))."""
    out = {}
    for mono, c in cur.items():
        words = [((), F(1))]
        for (p, m, j) in mono:
            if p == 0:
                words = [(w + ((p, m, j),), cc) for w, cc in words]
            else:
                new = []
                for w, cc in words:
                    new.append((w + ((1, m, 0),), cc))
                    new.append((w + ((1, m, 1),), cc * (1 if j == PP else -1)))
                words = new
        for w, cc in words:
            cm, sg = canonicalize(w)
            if cm is None:
                continue
            out[cm] = out.get(cm, 0) + sg * cc * c
    return gm.clean(out)


# (1)
words = []
for w2 in range(2, 9):
    for mono in gen_basis_super(w2, 4):
        if all(l[0] == 1 and l[2] in (PP, PM) for l in mono) and coset.charge(mono) == 0:
            words.append(mono)
ok1 = True
nq = 0
for w in words:
    ferm = pm_to_01({w: F(1)})
    ef = gm.vev_ns(ferm, eta, [0, 0])
    img = Bz.image_word(w)
    eb = gm.vev_ns(img, (1, 1, 1), [0, 0, 0])
    if sp.simplify(ef - eb) != 0:
        ok1 = False
        print('  vev MISMATCH', w, ef, eb)
    if len(w) >= 4:
        nq += 1
        # split independence: image via the last psi^+ / last psi^- pair instead of the first
        ip = max(k for k, l in enumerate(w) if l[2] == PP)
        jm = max(k for k, l in enumerate(w) if l[2] == PM)
        B = (w[ip], w[jm]) if ip < jm else (w[jm], w[ip])
        Mp = tuple(l for k, l in enumerate(w) if k not in (ip, jm))
        cm, sg = canonicalize(B + Mp)
        NP = coset.ope_pm_mono(B, Mp, 0, s)
        corr = dict(NP)
        corr[cm] = corr.get(cm, 0) - sg
        corr = gm.clean(corr)
        res = gm.normal_prod(Bz.image_word(B), Bz.image_word(Mp), (1, 1, 1))
        for m, c in corr.items():
            for mm, cc in Bz.image_word(m).items():
                res[mm] = res.get(mm, 0) - c * cc
        res = gm.scale(gm.clean(res), sg)
        if gm.add((1, res), (-1, img)):
            ok1 = False
            print('  split-dependence', w)
print(f"(1) bosonisation vev check on {len(words)} neutral words (w2 <= 8; {nq} of length >= 4, split-independence too): {'ALL OK' if ok1 else 'FAILURES'}")

# (2)
P0, P1 = sp.symbols('P0 P1')
Pu, Pv = sp.symbols('P_u P_v')
q = -2 * P1 / sp.Rational(s)
eta_uv = (F(s), 4 - F(s))
for tower in ('Z', 'K'):
    ch = load_charges(tower, s)
    for w2 in sorted(ch):
        if w2 < 4:
            continue
        for i, I in enumerate(ch[w2]):
            Ef = gm.vev_ns(I, eta, [P0, P1])
            Dq, info = coset.coset_project(I, s, Bz)
            Eb = 0
            for k, cur in Dq.items():
                Eb += q ** k * gm.vev_ns(cur, eta_uv, [Pu, Pv]).subs({Pu: P0, Pv: -P1})
            diff = sp.simplify(sp.expand(Ef - Eb))
            print(f"(2) tower {tower} w2 = {w2} charge {i}: full-frame eigenvalue == projected eigenvalue at q = -2P_1/s: {diff == 0}"
                  + ('' if diff == 0 else f'   diff = {diff}'))
