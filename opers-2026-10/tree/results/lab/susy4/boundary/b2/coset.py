"""coset.py -- B2 Level II machinery: bosonise the fermions, rotate to (u, v, H~), project onto the H~-oscillator vacuum.

Frame (n2frame): fields (X_0, psi_0), (X_1, psi_1), metric (s, -s).  Null fermions psi^+- = psi_0 +- psi_1,
<psi^+ psi^-> = 2s/(z-w), <psi^+- psi^+-> = 0.  Bosonisation psi^+- = sqrt(2s) e^{+-X_3}, <X_3 X_3> = +log (metric 1):
    psi^+(z) psi^-(w) = 2s (z-w)^-1 :exp(X_3(z) - X_3(w)): = 2s sum_n (z-w)^(n-1) S_n(w),  S_n = Schur polynomials in d^k X_3/k!.
Point-split bilinears (the (z-w)^0 term of d_z^a psi^+(z) d_w^b psi^-(w)):
    :d^a psi^+ d^b psi^-: = 2s sum_{j=0}^b C(b,j) (-1)^(b-j) (a+b-j)! d^j S_{a+b-j+1}.
Longer neutral Wick words are reduced recursively through the vertex-algebra identity
    :B M':_Wick = NP(B, M') - (contraction corrections),  image(NP(B, M')) = NP_X3(image B, image M'),
with NP the point-split (-1)-product (bilinear B = first psi^+ letter with first psi^- letter).  The map is a vertex-algebra
homomorphism from the charge-0 fermion algebra to the Heisenberg algebra of X_3; validate() checks it on OPE residues,
normal products and derivatives of random neutral words.

J = (psi_0 psi_1 - 2 dX_1)/s -> a.dY with a = (0, -2/s, -1) in Y = (X_0, X_1, X_3) (level1.py); H~ = a.Y, v = -X_1 - 2X_3,
u = X_0; X_1 = s(v - 2H~)/(4-s), X_3 = (sH~ - 2v)/(4-s); metric (u, v, H~) = (s, 4-s, (s-4)/s).
Projection: P_H (oint sum A_i B_i) P_H = oint sum A_i <q|B_i|q>, the H~-part evaluated by the cylinder Wick recipe
(super_cyl_vev kernel 'B', metric g_H per contraction, unpaired dH~ -> q, higher derivatives -> 0): a (u, v) density with
coefficients polynomial in q.
"""
import sys
sys.dont_write_bytecode = True
import random
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial

import sympy as sp

import gm
from gm import X, PSI, add, scale, clean
from super_engine import _koszul, canonicalize, d_mono, is_odd
from ope_engine import _compositions
import super_cyl_vev as SCV

PP, PM, X3 = 2, 3, 2   # fermion field indices psi^+ (2), psi^- (3); the boson X_3 carries index 2 (no clash: parity 0)


# ------------------------------------------------------------------ A. (psi_0, psi_1) -> (psi^+, psi^-)

def to_pm(cur):
    """psi_0 = (psi^+ + psi^-)/2, psi_1 = (psi^+ - psi^-)/2, letter by letter in word order, then canonicalised."""
    out = {}
    for mono, c in cur.items():
        words = [((), 1)]
        for (p, m, j) in mono:
            if p == 0:
                words = [(w + ((p, m, j),), sg) for w, sg in words]
            else:
                new = []
                for w, sg in words:
                    new.append((w + ((1, m, PP),), sg * F(1, 2)))
                    new.append((w + ((1, m, PM),), sg * (F(1, 2) if j == 0 else F(-1, 2))))
                words = new
        for w, sg in words:
            cm, s_ = canonicalize(w)
            if cm is None:
                continue
            out[cm] = out.get(cm, 0) + s_ * sg * c
    return clean(out)


def charge(mono):
    return sum(1 for l in mono if l[0] == 1 and l[2] == PP) - sum(1 for l in mono if l[0] == 1 and l[2] == PM)


def split_charge(cur):
    out = {}
    for m, c in cur.items():
        out.setdefault(charge(m), {})[m] = c
    return out


# ------------------------------------------------------------------ B. OPE products in the psi^+- basis

def _pair_pm(lp, lq, s):
    """pairing of two fermion letters in the null basis: psi^+ with psi^- only, value 2s (-1)^a (a+b)!/(z-w)^(a+b+1)."""
    p1, a, i = lp
    p2, b, j = lq
    if p1 != 1 or p2 != 1 or {i, j} != {PP, PM}:
        return None
    return (2 * s * (-1) ** a * factorial(a + b), a + b + 1)


def ope_pm_mono(P, Q, n, s):
    """coefficient of (z-w)^(-n) in P(z) Q(w) for fermion words in the psi^+- basis (n <= 0 includes the
    contraction-free Taylor term)."""
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(0, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, D = 1, 0
                for a_, b_ in zip(psel, qsel):
                    pr = _pair_pm(P[a_], Q[b_], s)
                    if pr is None:
                        num = 0
                        break
                    num *= pr[0]
                    D += pr[1]
                if not num or D < n:
                    continue
                if k:
                    num *= _koszul(P, Q, psel, qsel)
                restP = [P[a_] for a_ in range(np_) if a_ not in psel]
                restQ = [Q[b_] for b_ in range(nq) if b_ not in set(qsel)]
                r = len(restP)
                if r == 0:
                    if D == n:
                        cm, sg = canonicalize(tuple(restQ))
                        if cm is not None:
                            out[cm] = out.get(cm, 0) + sg * num
                    continue
                for ts in _compositions(D - n, r):
                    den = 1
                    sh = []
                    for (p, m, j), t in zip(restP, ts):
                        den *= factorial(t)
                        sh.append((p, m + t, j))
                    cm, sg = canonicalize(tuple(sh) + tuple(restQ))
                    if cm is None:
                        continue
                    out[cm] = out.get(cm, 0) + F(sg * num, den)
    return {m: c for m, c in out.items() if c}


def ope_pm(A, B, n, s):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            for mo, c in ope_pm_mono(ma, mb, n, s).items():
                out[mo] = out.get(mo, 0) + c * ca * cb
    return clean(out)


# ------------------------------------------------------------------ C. bosonisation

@lru_cache(maxsize=None)
def schur(n):
    """S_n = coefficient of t^n in exp(sum_k t^k y_k), y_k = d^k X_3 / k!, as a current in X_3 letters."""
    if n == 0:
        return {(): F(1)}
    t = sp.Symbol('t')
    ys = sp.symbols(f'y1:{n + 1}')
    e = sp.series(sp.exp(sum(t ** k * ys[k - 1] for k in range(1, n + 1))), t, 0, n + 1).removeO()
    poly = sp.Poly(sp.expand(e.coeff(t, n)), *ys)
    out = {}
    for mon, c in poly.terms():
        letters = []
        coef = F(int(c.p), int(c.q))
        for k, e_ in enumerate(mon, start=1):
            for _ in range(e_):
                letters.append((0, k, X3))
                coef /= factorial(k)
        out[tuple(sorted(letters))] = out.get(tuple(sorted(letters)), 0) + coef
    return clean(out)


def dcur(cur):
    out = {}
    for m, c in cur.items():
        for mm, cc in d_mono(m).items():
            out[mm] = out.get(mm, 0) + c * cc
    return clean(out)


def bilinear_image(a, b, s):
    """:d^a psi^+ d^b psi^-:  ->  2s sum_{j<=b} C(b,j)(-1)^(b-j)(a+b-j)! d^j S_{a+b-j+1}."""
    out = {}
    for j in range(b + 1):
        term = schur(a + b - j + 1)
        for _ in range(j):
            term = dcur(term)
        coef = 2 * s * comb(b, j) * (-1) ** (b - j) * factorial(a + b - j)
        for m, c in term.items():
            out[m] = out.get(m, 0) + coef * c
    return clean(out)


class Bosoniser:
    def __init__(self, s):
        self.s = F(s)
        self.memo = {}

    def image_word(self, word):
        """image of the Wick monomial `word` (canonical fermion word in psi^+- letters, neutral) in X_3 letters."""
        if word in self.memo:
            return self.memo[word]
        s = self.s
        if not word:
            res = {(): F(1)}
        elif len(word) == 2:
            (p1, a, i), (p2, b, j) = word
            if (i, j) == (PP, PM):
                res = bilinear_image(a, b, s)
            else:
                assert (i, j) == (PM, PP)
                res = scale(bilinear_image(b, a, s), -1)
        else:
            assert charge(word) == 0
            ip = next(k for k, l in enumerate(word) if l[2] == PP)
            jm = next(k for k, l in enumerate(word) if l[2] == PM)
            B = (word[ip], word[jm]) if ip < jm else (word[jm], word[ip])
            Mp = tuple(l for k, l in enumerate(word) if k not in (ip, jm))
            cm, sg = canonicalize(B + Mp)
            assert cm == word, (cm, word)
            NP = ope_pm_mono(B, Mp, 0, s)              # point-split product (contraction-free term + corrections)
            corr = dict(NP)
            corr[cm] = corr.get(cm, 0) - sg               # subtract the concatenation term
            corr = clean(corr)
            imB, imM = self.image_word(B), self.image_word(Mp)
            res = gm.normal_prod(imB, imM, (1, 1, 1))
            for m, c in corr.items():
                for mm, cc in self.image_word(m).items():
                    res[mm] = res.get(mm, 0) - c * cc
            res = scale(clean(res), sg)
        self.memo[word] = res
        return res

    def image(self, cur):
        """linear extension to currents in (X_0, X_1, psi^+-) letters; bosons pass through (concatenation)."""
        out = {}
        for mono, c in cur.items():
            bos = tuple(l for l in mono if l[0] == 0)
            fer = tuple(l for l in mono if l[0] == 1)
            assert charge(fer) == 0, mono
            cm0, sg0 = canonicalize(bos + fer)
            assert cm0 == mono and sg0 == 1
            for mm, cc in self.image_word(fer).items():
                key = tuple(sorted(bos + mm))
                out[key] = out.get(key, 0) + c * cc
        return clean(out)


def validate(s, trials=60, seed=1, maxw2=8):
    """vertex-algebra homomorphism checks on random neutral words."""
    rng = random.Random(seed)
    Bz = Bosoniser(s)
    from super_engine import gen_basis_super
    words = []
    for w2 in range(2, maxw2 + 1):
        for mono in gen_basis_super(w2, 4):
            if all(l[0] == 1 and l[2] in (PP, PM) for l in mono) and charge(mono) == 0:
                words.append(mono)
    eta3 = (1, 1, 1)
    ok = True
    n = 0
    for _ in range(trials):
        A, B = rng.choice(words), rng.choice(words)
        if gm.mono_weight2(A) + gm.mono_weight2(B) > 12:
            continue
        n += 1
        imA, imB = Bz.image_word(A), Bz.image_word(B)
        for nn in (1, 0, 2):
            lhs = Bz.image(ope_pm({A: F(1)}, {B: F(1)}, nn, F(s)))
            rhs = gm.ope_coeff(imA, imB, nn, eta3)
            if add((1, lhs), (-1, rhs)):
                ok = False
                print('   MISMATCH ope', nn, A, B)
        if add((1, Bz.image(dcur({A: F(1)}))), (-1, dcur(imA))):
            ok = False
            print('   MISMATCH d', A)
    print(f"bosonisation validation at s = {s}: {n} random word pairs, residue/normal-product/double-pole/derivative "
          f"homomorphism: {'ALL OK' if ok else 'FAILURES'}")
    return ok


# ------------------------------------------------------------------ D. rotation to (u, v, H~)

def rotate(cur, s):
    """(X_0, X_1, X_3) -> (u, v, H~): X_1 = s(v - 2H~)/(4-s), X_3 = (sH~ - 2v)/(4-s); output indices u 0, v 1, H~ 2."""
    s = F(s)
    out = {}
    for mono, c in cur.items():
        words = [((), c)]
        for (p, m, j) in mono:
            assert p == 0
            if j == 0:
                words = [(w + ((0, m, 0),), cc) for w, cc in words]
            elif j == 1:
                new = []
                for w, cc in words:
                    new.append((w + ((0, m, 1),), cc * s / (4 - s)))
                    new.append((w + ((0, m, 2),), cc * (-2 * s) / (4 - s)))
                words = new
            else:
                new = []
                for w, cc in words:
                    new.append((w + ((0, m, 2),), cc * s / (4 - s)))
                    new.append((w + ((0, m, 1),), cc * F(-2) / (4 - s)))
                words = new
        for w, cc in words:
            key = tuple(sorted(w))
            out[key] = out.get(key, 0) + cc
    return clean(out)


def metric_uvH(s):
    s = F(s)
    return (s, 4 - s, (s - 4) / s)


# ------------------------------------------------------------------ E. projection onto the H~-oscillator vacuum

def project(cur, s):
    """returns {q-power: (u,v)-current (Fractions)}; the H~ letters of each monomial replaced by their cylinder vev at
    H~-momentum q (metric g_H)."""
    s = F(s)
    gH = (s - 4) / s
    out = {}
    for mono, c in cur.items():
        uv = tuple(l for l in mono if l[2] != 2)
        hh = tuple((m, 0) for (p, m, j) in mono if j == 2)
        if not hh:
            out.setdefault(0, {})[uv] = out.setdefault(0, {}).get(uv, 0) + c
            continue
        part = SCV._boson_part(hh, 1, 'B')          # {(e,): value} at unit metric
        n = len(hh)
        for (e,), val in part.items():
            npair = (n - e) // 2
            coef = c * val * gH ** npair
            d = out.setdefault(e, {})
            d[uv] = d.get(uv, 0) + coef
    return {k: clean(v) for k, v in out.items() if clean(v)}


def coset_project(cur, s, Bz=None):
    """full pipeline for a density in (X_0, X_1, psi_0, psi_1) letters: to_pm -> neutral part -> bosonise -> rotate ->
    project.  Returns (q-expansion dict, info)."""
    Bz = Bz or Bosoniser(s)
    pm = to_pm(cur)
    parts = split_charge(pm)
    info = {ch: (len(v), gm.is_total_derivative(v, 4)) for ch, v in parts.items()}
    neutral = parts.get(0, {})
    bos = Bz.image(neutral)
    rot = rotate(bos, s)
    return project(rot, s), info


def project_vertex_prefactor(W, s, Bz=None):
    """the prefactor current of a J-neutral vertex (bosonised, rotated, projected): {q-power: (u,v)-current}."""
    return coset_project(W, s, Bz)[0]


if __name__ == '__main__':
    sv = F(sys.argv[1]) if len(sys.argv) > 1 else F(2, 3)
    validate(sv)
    # J-direction check: J -> a.dY, a = (0, -2/s, -1)
    import n2frame as N2
    C = N2.currents(sv)
    Bz = Bosoniser(sv)
    Jpm = to_pm(C['J'])
    Jb = Bz.image(Jpm)
    print('J bosonised:', Jb, ' expected {dX_1: -2/s, dX_3: -1}')
    Jr = rotate(Jb, sv)
    print('J rotated:', Jr, ' (should be dH~ only)')
    S = N2.kl_vertices(sv)
    bottom, top = N2.w0_vertex(sv, 'G')
    print('W0 prefactor projected:', project_vertex_prefactor(top.W, sv, Bz))
