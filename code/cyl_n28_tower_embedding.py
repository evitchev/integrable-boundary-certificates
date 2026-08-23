"""Do the N = 28 self-dual embeddings of the cylindrical families persist
up the tower?  (Follow-up to cyl_hidden_kdv_scan.py, Part C.)

At (N, s) = (28, -9) [Solution 1] and (28, 9) [Solution 2] the scan
found I_3 and the unique I_5 conserved to first order under
e^{+-sqrt2 X} -- the self-dual pair whose joint kernel is L(1,0)_X, so
"X enters only through its stress tensor".  This certificate tests
the promoted Q(t)-symbolic charges I_7, I_9, I_11 (results/
qt_promotion.json, coordinates A(N) + s B(N) on the even mixed basis)
specialized to those points, together with I_5 recomputed there:

 1. Reconstruction is pole-free at N = 28, each specialized density is
    a NONZERO CHARGE CLASS (not merely a nonzero density -- a total
    derivative would pass every later test trivially), and each
    commutes with I_3 in-engine (transcription/specialization guard).
 2. e^{sqrt2 X} conservation, EXACT over Q(sqrt2): residue in the
    mixed representation, twisted derivative d + sqrt2 dX, the
    Q(sqrt2) question doubled to a Q question (v = v_r + sqrt2 v_q;
    sqrt2 acts as (v_r, v_q) -> (2 v_q, v_r)), membership by
    fraction-free sparse elimination (SparsePivots).
 3. Controls: I_3 at the same points is screened (agrees with the
    scan); the same I_7 at a generic curve point is NOT; the N = 28
    I_7 is NOT screened at the non-root a = sqrt3.
 4. MERGER: Solution 1 at (28,-9) and Solution 2 at (28,9) have the
    same I_3 (exactly) and the same normalized I_5..I_11 classes --
    the two sheets' families coincide through the known range: one
    tower, a finite-range merger/embedding through I_11 (no all-spin
    claim).  The common I_3 is exactly 4[7 :T_X^2: + 6 T_X T_Y + :T_Y^2:].

Fail-closed: exits nonzero unless every check passes.
"""
import json
import sys
from fractions import Fraction as F
from itertools import combinations
from math import factorial

import sympy as sp

from invariant_engine import curve_point, gen_inv_basis
from mixed_engine import (MixedSector, canon_mixed, cyl_P4, d_mixed_mono,
                          gen_mixed_basis, residue_mixed_cur)
from sparse_linalg import SparsePivots, lcm

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


# ---------------------------------------------------------------- basis
def gen_mixed_full(s):
    def xparts(total):
        res = [[]] if total == 0 else []

        def rec(mx, remaining, acc):
            if remaining == 0:
                res.append(list(acc))
                return
            for v in range(min(mx, remaining), 0, -1):
                acc.append(v)
                rec(v, remaining - v, acc)
                acc.pop()
        rec(total, total, [])
        return res
    if s == 0:
        return [canon_mixed((), ())]
    out = []
    for xspin in range(0, s + 1):
        for xs in xparts(xspin):
            rem = s - xspin
            if rem == 0:
                if xs:
                    out.append(canon_mixed(tuple(xs), ()))
                continue
            for yp in gen_inv_basis(rem):
                out.append(canon_mixed(tuple(xs), yp))
    return sorted(set(out))


def apply_d(cur):
    out = {}
    for m, c in cur.items():
        for mm, k in d_mixed_mono(m).items():
            out[mm] = out.get(mm, 0) + c * k
    return {m: c for m, c in out.items() if c}


# -------------------------------------------------- residue over Q(sqrt(r))
def residue_x_quadratic(P, r):
    """Res P(z) e^{aX}(w) with a = sqrt(r): returns {mono: (u, v)} meaning
    u + v*sqrt(r), Fractions."""
    out = {}
    for (xs, yp), cP in P.items():
        n = len(xs)
        for k in range(1, n + 1):
            for sel in combinations(range(n), k):
                num = F(1)
                D = 0
                for i in sel:
                    m = xs[i]
                    num *= (-1) ** (m - 1) * factorial(m - 1)
                    D += m
                rest_xs = tuple(xs[i] for i in range(n) if i not in sel)
                rest = canon_mixed(rest_xs, yp)
                T = D - 1
                cur = {rest: F(1)}
                for _ in range(T):
                    cur = apply_d(cur)
                # a^k = r^{k//2} (k even) or r^{(k-1)//2} sqrt(r) (k odd)
                pw = F(r) ** (k // 2)
                for m, c in cur.items():
                    val = cP * num * c / factorial(T) * pw
                    u, v = out.get(m, (F(0), F(0)))
                    if k % 2 == 0:
                        out[m] = (u + val, v)
                    else:
                        out[m] = (u, v + val)
    return {m: uv for m, uv in out.items() if uv[0] or uv[1]}


def is_screened(P, w, r):
    """Exact Q(sqrt r) test: is Res(P e^{sqrt(r) X}) a twisted total
    derivative d J + sqrt(r) dX J?"""
    s = w - 1
    tgt = []
    for t in range(0, s + 1):
        tgt += gen_mixed_full(t)
    lower = gen_mixed_full(s - 1)
    index = {m: i for i, m in enumerate(tgt)}
    n = len(tgt)
    piv = SparsePivots()
    for u in lower:
        du = d_mixed_mono(u)                      # rational part of D(u)
        xu = canon_mixed(u[0] + (1,), u[1])        # sqrt(r) part: dX * u
        col = {}
        for m, c in du.items():
            col[index[m]] = col.get(index[m], 0) + int(c)
        col[n + index[xu]] = col.get(n + index[xu], 0) + 1
        piv.insert({k: c for k, c in col.items() if c})
        # sqrt(r) * D(u): (v_r, v_q) -> (r v_q, v_r)
        col2 = {}
        col2[index[xu]] = col2.get(index[xu], 0) + r
        for m, c in du.items():
            col2[n + index[m]] = col2.get(n + index[m], 0) + int(c)
        piv.insert({k: c for k, c in col2.items() if c})
    res = residue_x_quadratic(P, r)
    den = 1
    for u, v in res.values():
        den = lcm(den, u.denominator)
        den = lcm(den, v.denominator)
    vec = {}
    for m, (u, v) in res.items():
        if u:
            vec[index[m]] = int(u * den)
        if v:
            vec[n + index[m]] = int(v * den)
    return not piv.reduce(vec)


# ------------------------------------------- promoted charges at a point
from pathlib import Path

QT_PATH = Path(__file__).resolve().parent.parent / "results" / "qt_promotion.json"


def promoted_density(sol, sigma, N, s):
    data = json.load(open(QT_PATH))[f"cyl{sol}_sigma{sigma}"]
    if data.get("symbolic_verified") is not True:
        raise ValueError(f"cyl{sol}_sigma{sigma}: promotion not symbolically "
                         "verified (fail-closed)")
    basis = gen_mixed_basis(sigma + 1)
    Nsym = sp.Symbol('N')
    Nq = sp.Rational(N.numerator, N.denominator)
    cur = {}
    for i, (A, B) in data["coeffs_A_B_of_N"].items():
        # bind the bare 'N' in the stored strings to a Symbol (bare
        # sympify would read it as sympy's N() function)
        Av = sp.sympify(A, locals={'N': Nsym}).subs(Nsym, Nq)
        Bv = sp.sympify(B, locals={'N': Nsym}).subs(Nsym, Nq)
        if Av in (sp.zoo, sp.nan) or Bv in (sp.zoo, sp.nan):
            raise ZeroDivisionError(f"pole at N={N} in coordinate {i}")
        val = sp.cancel(Av + sp.Rational(s.numerator, s.denominator) * Bv)
        if not val.is_Rational:
            raise ValueError(f"coordinate {i} not rational at the point: {val}")
        c = F(int(val.p), int(val.q))
        if c:
            cur[basis[int(i)]] = c
    return cur


POINTS = [(1, F(28), F(-9)), (2, F(28), F(9))]
from mixed_engine import genuine_kernel_mixed

classes = {}      # (sol, sigma) -> reduced class vector, for the merger test
for sol, N, s in POINTS:
    P4 = cyl_P4(sol, N, s)
    check(f"C0: Solution {sol} at ({N},{s}): I_3 screened at a = sqrt2 "
          f"(agrees with cyl_hidden_kdv_scan)", is_screened(P4, 4, 2))
    sec4 = MixedSector(4)
    classes[(sol, 3)] = sec4.reduce_mod_d(sec4.vec(P4))
    # I_5 from the ad-kernel at the point (unique), tested here too so the
    # marker's I_3..I_11 claim is self-contained
    ker5 = genuine_kernel_mixed(P4, gen_mixed_basis(6), N)
    check(f"C1: Solution {sol} at ({N},{s}): I_5 unique", len(ker5) == 1)
    if len(ker5) == 1:
        sec6 = MixedSector(6)
        v5 = sec6.reduce_mod_d(sec6.vec(ker5[0]))
        check(f"C1: Solution {sol} sigma=5: nonzero class", any(v5))
        check(f"C3: Solution {sol} sigma=5 at ({N},{s}): I_5 screened by "
              f"e^{{sqrt2 X}}", is_screened(ker5[0], 6, 2))
        j5 = next(i for i, x in enumerate(v5) if x)
        classes[(sol, 5)] = [x / v5[j5] for x in v5]
    for sigma in (7, 9, 11):
        w = sigma + 1
        try:
            P = promoted_density(sol, sigma, N, s)
        except (ZeroDivisionError, ValueError) as exc:
            check(f"C1: Solution {sol} sigma={sigma}: reconstruction at N=28 "
                  f"({exc})", False)
            continue
        sec = MixedSector(w)
        v = sec.reduce_mod_d(sec.vec(P))
        check(f"C1: Solution {sol} sigma={sigma} at ({N},{s}): reconstructed "
              f"density is a NONZERO CHARGE CLASS ({len(P)} monomials)", any(v))
        if not any(v):
            continue
        j = next(i for i, x in enumerate(v) if x)
        classes[(sol, sigma)] = [x / v[j] for x in v]
        comm = residue_mixed_cur(P4, P, N)
        tgt = MixedSector(w + 3)
        check(f"C2: Solution {sol} sigma={sigma}: [I_3, I_{sigma}] = 0 mod d "
              f"in-engine at the point", not any(tgt.reduce_mod_d(tgt.vec(comm))))
        check(f"C3: Solution {sol} sigma={sigma} at ({N},{s}): I_{sigma} "
              f"screened by e^{{sqrt2 X}} -- embedding persists",
              is_screened(P, w, 2))

# ---- the two N = 28 families coincide (deck-conjugate points, one tower)
for sigma in (3, 5, 7, 9, 11):
    a1, a2 = classes.get((1, sigma)), classes.get((2, sigma))
    ok = a1 is not None and a2 is not None
    if ok and sigma == 3:
        ok = a1 == a2                       # I_3 coincide exactly (same normalization)
    elif ok:
        ok = a1 == a2                       # normalized class vectors coincide
    check(f"C6: MERGER -- Solution 1 at (28,-9) and Solution 2 at (28,9) have "
          f"the SAME I_{sigma} class", ok)

# exact T-composite form of the common I_3 at N = 28
P4c = cyl_P4(2, F(28), F(9))
from mixed_engine import _M
check("C7: the common N=28 I_3 density is exactly "
      "7(dX)^4 + 6(dX)^2(dY.dY) + (dY.dY)^2 - 14(d2X)^2 - 2(d2Y.d2Y) "
      "= 4[7 :T_X^2: + 6 T_X T_Y + :T_Y^2:] with no further terms",
      P4c == {_M((1, 1, 1, 1)): F(7), _M((1, 1), [(1, 1)]): F(6),
              _M((), [(1, 1), (1, 1)]): F(1), _M((2, 2)): F(-14),
              _M((), [(2, 2)]): F(-2)})

# ---- controls on the N = 28 I_7 of Solution 2
P7 = promoted_density(2, 7, F(28), F(9))
check("C4 control: the N=28 I_7 is NOT screened at the non-root a = sqrt3",
      not is_screened(P7, 8, 3))
Ng, sg = curve_point(F(2))
P7g = promoted_density(2, 7, Ng, sg)
check(f"C5 control: Solution 2 I_7 at a generic point (N={Ng}, s={sg}) is "
      f"NOT screened at a = sqrt2", not is_screened(P7g, 8, 2))
check(f"C5 control: ... and its I_3 there is not screened either",
      not is_screened(cyl_P4(2, Ng, sg), 4, 2))

print()
if failures:
    print(f"N=28 TOWER-EMBEDDING CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("N=28 TOWER EMBEDDING CERTIFIED (I_3, I_5, I_7, I_9, I_11 of Solutions 1 "
      "and 2 -- which coincide at (28,-9)/(28,9) -- are nonzero classes conserved "
      "under e^{+-sqrt2 X}; finite-range merger/embedding through I_11; "
      "controls negative)")
