"""b2b.py -- B2b: do Solution 1's charges (lifted to the three-boson algebra) commute with the N = 2 KdV towers?

Usage: b2b.py <tower Z|K> <t>
For s, r in {2, 3, 5} (I_2 of Z only): Res(L_s(z) I_r(w)) mod d in the bosonised (u, v, H~) algebra, metric (s, 4-s, (s-4)/s).
Also: [L_s, L_r] (Solution 1's charges among themselves), [I_r, I_r'] in the fermionic frame and after bosonisation
(consistency of the bosonisation on commutators), and the negative control (lift with kappa -> kappa + 1/7, Q(kappa) arithmetic).
"""
import sys
sys.dont_write_bytecode = True
import os as _os_kit; sys.path.insert(0, _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', 'b2'))  # KIT PATCH: was an absolute path to a private copy of b2/ (all *.py byte-identical to the archived b2/)
import json
import time
from fractions import Fraction as F

import sympy as sp

import gm
from gm import X, add, scale, clean
from super2_screening import boson_vertex, basis
from sparse_linalg import SparsePivots, _to_int_vec
from super_engine import gen_basis_super, d_mono, mono_weight2
import n2frame as N2
import coset
from towers import cur_from_json

tower, t = sys.argv[1], F(sys.argv[2])
s = 2 * (t - 1) / (t + 1)
stag = str(s).replace('/', '-')
ttag = str(t).replace('/', '-')
eta3 = (F(s), 4 - F(s), (F(s) - 4) / F(s))
eta_uv = (F(s), 4 - F(s))
t0 = time.time()


def P(msg):
    print(msg, flush=True)


# ------------------------------------------------------------------ bosonic total-derivative test, nf fields
_PIV = {}


def bos_basis(w2, nf):
    if w2 < 0:
        return []
    if w2 == 0:
        return [()]
    return [m for m in gen_basis_super(w2, nf) if all(l[0] == 0 for l in m)]


def is_d_exact(cur, nf):
    if not cur:
        return True
    w2 = mono_weight2(next(iter(cur)))
    key = (w2, nf)
    if key not in _PIV:
        index = {m: i for i, m in enumerate(bos_basis(w2, nf))}
        piv = SparsePivots()
        for u in bos_basis(w2 - 2, nf):
            iv, _ = _to_int_vec(d_mono(u), index)
            if iv:
                piv.insert(iv)
        _PIV[key] = (index, piv)
    index, piv = _PIV[key]
    iv, _ = _to_int_vec(cur, index)
    return not piv.reduce(iv)


def commute(A, B, eta):
    """[oint A, oint B] = 0 mod d ?  (bosonic currents, nf = len(eta)); returns (bool, #monomials of the residue)."""
    R = gm.residue(A, B, eta)
    return (is_d_exact(R, len(eta)) if R else True), len(R)


# ------------------------------------------------------------------ Q(kappa) arithmetic for the control: numbers (a, b) = a + b kappa
k2 = 1 - F(s) / 4            # kappa^2 (rational)


def qmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d * k2, a * d + b * c)


def qpow(x, n):
    r = (F(1), F(0))
    for _ in range(n):
        r = qmul(r, x)
    return r


# ------------------------------------------------------------------ load
def load_charges(tw):
    ch = {}
    for fn in (f'inputs/towers_{tw}_s{stag}_w2_3_4_5_6_7_8.json', f'inputs/towers_{tw}_s{stag}_w10_12.json'):
        d = json.load(open(fn))
        for w2, r in d['results'].items():
            ch[int(w2)] = [cur_from_json(k) for k in r['kernel']]
    return ch


ch = load_charges(tower)
sol1 = json.load(open(f'inputs/sol1_t{ttag}.json'))
I3uv = cur_from_json(sol1['I3_uv'])
Vuv = [boson_vertex((F(1, 2), F(-1, 2))), boson_vertex((F(1, 2), F(1, 2))), boson_vertex((F(-1), F(0)))]
gens6 = gm.joint_kernel([{m: F(1)} for m in basis(6, 2) if all(l[0] == 0 for l in m)], Vuv, eta_uv)
gens12 = gm.joint_kernel([{m: F(1)} for m in basis(12, 2) if all(l[0] == 0 for l in m)], Vuv, eta_uv)
assert len(gens6) == 1
L2 = scale(gens6[0], F(1) / gens6[0][(X(1, 0), X(1, 0), X(1, 1))])
L3 = I3uv
# L5: the generator proportional to the certified I_5 eigenvalue (Q^6-normalised)
PX, Ppi = sp.symbols('P_X pi'); Pu, Pv = sp.symbols('P_u P_v')
p_ = sp.sqrt(sp.Rational(s)); rho = p_ / 2 - 1 / p_; kappa = sp.sqrt(1 - sp.Rational(s) / 4)
P2, Q2 = sp.symbols('P2 Q2')
E5 = sp.expand(sp.sympify(sol1['E5_alt']).subs({P2: PX ** 2, Q2: Ppi ** 2 - rho ** 2}))
mom = {Pu: -p_ * Ppi, Pv: -2 * kappa * PX}
L5 = None
for g in gens12:
    Eg = sp.expand(gm.vev_ns(g, eta_uv, [Pu, Pv]).subs(mom))
    c = sp.Symbol('c')
    sol = sp.solve(sp.Poly(sp.expand(Eg - c * E5), PX, Ppi).coeffs(), c, dict=True)
    if sol and sol[0][c] != 0:
        cc = sp.nsimplify(1 / sol[0][c])
        L5 = scale(g, F(int(cc.p), int(cc.q)))
assert L5 is not None
P(f"tower {tower}, t = {t}, s = {s}: L_2 ({len(L2)} monomials), L_3 ({len(L3)}), L_5 ({len(L5)}; (u,v) commutant dim {len(gens12)} at w2 = 12)")
Ls = {2: L2, 3: L3, 5: L5}

# tower charges: fermionic frame and bosonised
Bz = coset.Bosoniser(s)
spins = {2: 6, 3: 8, 5: 12}
Iferm, Ibos = {}, {}
for r, w2 in spins.items():
    if not ch.get(w2):
        continue
    cur = ch[w2][0]
    key = (X(1, 0),) * (w2 // 2)
    cur = scale(cur, F(1) / cur[key]) if key in cur else cur
    Iferm[r] = cur
    pm = coset.to_pm(cur)
    parts = coset.split_charge(pm)
    for q_, part in parts.items():
        if q_ != 0:
            assert gm.is_total_derivative(part, 4), (r, q_)
    bos = Bz.image(parts.get(0, {}))
    Ibos[r] = coset.rotate(bos, s)
    P(f"  I_{r} of {tower}: {len(cur)} monomials; charge sectors {sorted(parts)} (non-neutral parts total derivatives: asserted); "
      f"bosonised+rotated {len(Ibos[r])} monomials")

out = dict(tower=tower, t=str(t), s=str(s), pairs={}, controls={}, consistency={})
P("\n== main test: [oint L_s, oint I_r] = 0 mod d in the three-boson algebra")
for sL, L in Ls.items():
    for r, I in Ibos.items():
        t1 = time.time()
        ok, n = commute(L, I, eta3)
        out['pairs'][f'({sL},{r})'] = dict(commute=ok, residue_monomials=n, secs=round(time.time() - t1, 1))
        P(f"   [L_{sL}, I_{r}] = 0 mod d: {ok}   (residue {n} monomials, {time.time() - t1:.1f}s)")

P("\n== consistency: Solution 1's charges among themselves (u,v frame) and the tower's charges among themselves")
for (a, b) in ((2, 3), (2, 5), (3, 5)):
    ok, n = commute(Ls[a], Ls[b], eta_uv)
    out['consistency'][f'[L_{a},L_{b}]'] = ok
    P(f"   [L_{a}, L_{b}] = 0 mod d: {ok}")
rs = sorted(Ibos)
for i, a in enumerate(rs):
    for b in rs[i + 1:]:
        Rf = gm.residue(Iferm[a], Iferm[b], N2.eta_of(s))
        okf = gm.is_total_derivative(Rf, 2) if Rf else True
        okb, n = commute(Ibos[a], Ibos[b], eta3)
        out['consistency'][f'[I_{a},I_{b}] fermionic'] = okf
        out['consistency'][f'[I_{a},I_{b}] bosonised'] = okb
        P(f"   [I_{a}, I_{b}] = 0 mod d: fermionic frame {okf}; bosonised frame {okb}")

P("\n== negative control: lift with kappa -> kappa + 1/7  (X -> -v/(2 kappa') instead of -v/(2 kappa): factor (kappa/kappa')^(n_v))")
# kappa/kappa' = kappa (kappa - 1/7)/(kappa^2 - 1/49) = a + b kappa
den = k2 - F(1, 49)
ratio = (k2 / den, F(-1, 7) / den)
for sL, L in Ls.items():
    LA, LB = {}, {}
    for m, c in L.items():
        nv = sum(1 for l in m if l[2] == 1)
        a, b = qpow(ratio, nv)
        if a * c:
            LA[m] = LA.get(m, 0) + a * c
        if b * c:
            LB[m] = LB.get(m, 0) + b * c
    LA, LB = clean(LA), clean(LB)
    for r, I in Ibos.items():
        okA, nA = commute(LA, I, eta3) if LA else (True, 0)
        okB, nB = commute(LB, I, eta3) if LB else (True, 0)
        out['controls'][f'({sL},{r})'] = dict(commute=(okA and okB), partA=okA, partB=okB)
        P(f"   control [L_{sL}', I_{r}] = 0 mod d: {okA and okB}   (rational part {okA}, kappa part {okB})")
# the control lift must also break the self-commutation of Solution 1's charges? (not required; reported)
json.dump(out, open(f'b2b_{tower}_t{ttag}.json', 'w'), indent=1)
P(f"\ndone in {time.time() - t0:.1f}s")
