"""Fermionic screening analysis for the super-KdV tower (Phase 2).

V_alpha = :psi e^{alpha X}: (one N=1 superfield).  oint V_alpha
commutes with oint P iff Res_{z->w} P(z) V_alpha(w) is a *twisted*
total derivative: Res = D_alpha(J) V with D_alpha J = dJ + alpha(dX)J.

Contractions of P(z) against V(w): every subset of P's boson letters
hits the exponential (d^m X -> alpha (-1)^(m-1) (m-1)! (z-w)^-m), and
at most one fermion letter hits the psi prefactor
(d^a psi -> (-1)^a a! (z-w)^-(a+1)), with Koszul sign
(-1)^{# odd letters after it} since psi sits at the end of the word.
The surviving psi (if unconsumed) is a w-letter: it never Taylor-soaks.
Every emitted word has doubled weight w2(P) - 1.

AUDIT-CORRECTED (2026-08-04): the composite pencil uses the
POINT-SPLIT product (super_engine.normal_prod).  Findings at c = 3/2:

  - alpha = +-1 annihilates EVERY member of the pencil
    :TT: + kappa :G dG: — the dual screening pair psi e^{+-X} of
    Chistyakova-Litvinov-Orlov (2110.05870) at the free point: in
    their conventions c-hat = 1 + 2Q^2 with Q = b + 1/b, so c-hat = 1
    forces Q = 0, b = -1/b = +-i, and with X = i Phi the code
    exponents are alpha = -ib = +-1.  Screening therefore does NOT
    select the super-KdV ray here; G-invariance does (kappa = -1/4,
    cf. super_kdv.select_I3).
  - the higher tower is annihilated on the same locus: every
    screening condition of I5 (and of I7 under --deep) is divisible
    by (alpha^2 - 1), with alpha = 0 excluded.
  - the bare-fermion branch {alpha = 0, kappa = 1/2} exists: oint psi
    annihilates the ray :TT: + 1/2 :G dG: only.

CAVEAT (audit finding 4): conditions_sympy eliminates over QQ(alpha)
with pivots that may depend on alpha — a generic localization, so
even substituting claimed points into its OUTPUT is not branch-safe.
Locus statements are therefore "generic stratum" statements, and the
claimed points are certified separately by exact rational-point
membership (twisted_exact_fraction: integer elimination at fixed
alpha) at alpha = +-1 and at the bare alpha = 0 branch;
exceptional-branch completeness would need Groebner/Fitting-ideal
elimination.

--naive reproduces the legacy concatenation-pencil fit (screening
locus 2 alpha^4 - 9 alpha^2 + 2 with kappa = -1/4 forced INSIDE that
slice — an artifact of the naive pencil, retained as data for the
background-charge hypothesis; see results/super_kdv_calibration.md).
"""
import json
import sys
from fractions import Fraction as F
from itertools import combinations
from math import factorial

import sympy as sp

from ope_engine import _compositions
from super_engine import (canonicalize, cur_add, d_cur, d_mono,
                          gen_basis_super, genuine_kernel_super, is_odd,
                          residue_cur)
from sparse_linalg import SparsePivots, _to_int_vec

NF = 1
DX, PSI = (0, 1, 0), (1, 0, 0)
T = {(DX, DX): F(1, 2), (PSI, (1, 1, 0)): F(-1, 2)}
G = {(DX, PSI): F(1)}


def residue_screen_super(P, alpha):
    """Res of P(z) V_alpha(w), the exponential implicit in every word.
    Scalar arithmetic is generic: alpha may be a sympy Symbol or a
    Fraction."""
    out = {}
    for mono, cP in P.items():
        c0 = (sp.Rational(cP.numerator, cP.denominator)
              if isinstance(cP, F) and isinstance(alpha, sp.Basic)
              else cP)
        n = len(mono)
        bos = [i for i in range(n) if not is_odd(mono[i])]
        fer = [i for i in range(n) if is_odd(mono[i])]
        for kb in range(len(bos) + 1):
            for bsel in combinations(bos, kb):
                for fi in (None,) + tuple(fer):
                    num, D = c0, 0
                    for i in bsel:
                        m = mono[i][1]
                        num = num * alpha * ((-1) ** (m - 1)
                                             * factorial(m - 1))
                        D += m
                    keep_psi = True
                    if fi is not None:
                        a = mono[fi][1]
                        sign = ((-1) ** sum(1 for k in fer if k > fi))
                        num = num * sign * ((-1) ** a * factorial(a))
                        D += a + 1
                        keep_psi = False
                    if D == 0:
                        continue
                    sel = set(bsel) | ({fi} if fi is not None else set())
                    rest = [mono[i] for i in range(n) if i not in sel]
                    r = len(rest)
                    if r == 0:
                        if D == 1:
                            word = (PSI,) if keep_psi else ()
                            cm, sg = canonicalize(word)
                            out[cm] = out.get(cm, 0) + sg * num
                        continue
                    for ts in _compositions(D - 1, r):
                        den = 1
                        shifted = []
                        for (p, m, j), t in zip(rest, ts):
                            den *= factorial(t)
                            shifted.append((p, m + t, j))
                        word = tuple(shifted) + ((PSI,) if keep_psi
                                                 else ())
                        cm, sg = canonicalize(word)
                        if cm is None:
                            continue
                        out[cm] = out.get(cm, 0) + sg * num / den
    if isinstance(alpha, sp.Basic):
        return {m: sp.expand(c) for m, c in out.items()
                if sp.expand(c) != 0}
    return {m: c for m, c in out.items() if c}


def twisted_d_super(mono, alpha):
    """D_alpha J = dJ + alpha (dX) J on a super monomial."""
    out = {m: (sp.Integer(c) if isinstance(alpha, sp.Basic) else F(c))
           for m, c in d_mono(mono).items()}
    cm, sg = canonicalize(mono + (DX,))
    if cm is not None:
        out[cm] = out.get(cm, 0) + sg * alpha
    return {m: c for m, c in out.items() if c}


def conditions_sympy(P, alpha):
    """Polynomial alpha-conditions for Res(P V) in Im D_alpha, via the
    coefficient-column-pivot elimination of screening_fit.py."""
    res = residue_screen_super(P, alpha)
    w2r = max(sum(2 * m + (1 if p else 0) for p, m, _ in mono)
              for mono in res if mono) if res else 0
    tgt = gen_basis_super(w2r, NF) + ([()] if w2r == 0 else [])
    lower = gen_basis_super(w2r - 2, NF) + ([()] if w2r == 2 else [])
    index = {m: i for i, m in enumerate(tgt)}
    rows = [[sp.Integer(0)] * (len(lower) + 1) for _ in tgt]
    for cidx, u in enumerate(lower):
        for m, c in twisted_d_super(u, alpha).items():
            rows[index[m]][cidx] = c
    for m, c in res.items():
        rows[index[m]][len(lower)] = c
    ncols = len(lower)
    used = set()
    for c in range(ncols):
        piv, best = None, None
        for r in range(len(rows)):
            if r in used:
                continue
            e = sp.cancel(rows[r][c])
            rows[r][c] = e
            if e != 0:
                if not e.free_symbols:
                    piv = r
                    break
                if best is None:
                    best = r
        if piv is None:
            piv = best
        if piv is None:
            continue
        used.add(piv)
        pe = rows[piv][c]
        for r in range(len(rows)):
            if r == piv:
                continue
            f = sp.cancel(rows[r][c])
            if f != 0:
                rows[r] = [sp.cancel(x - f * y / pe)
                           for x, y in zip(rows[r], rows[piv])]
    conds = []
    for r in range(len(rows)):
        if r in used:
            continue
        if all(sp.cancel(rows[r][c]) == 0 for c in range(ncols)):
            e = sp.cancel(rows[r][ncols])
            if e != 0:
                conds.append(sp.factor(sp.numer(sp.together(e))))
    return conds


def twisted_exact_fraction(P, aval):
    """Exact membership Res(P V_a) in Im D_a at a rational point."""
    res = residue_screen_super(P, aval)
    if not res:
        return True
    w2r = max(sum(2 * m + (1 if p else 0) for p, m, _ in mono)
              for mono in res if mono)
    tgt = gen_basis_super(w2r, NF)
    index = {m: i for i, m in enumerate(tgt)}
    piv = SparsePivots()
    for u in gen_basis_super(w2r - 2, NF):
        iv, _ = _to_int_vec(twisted_d_super(u, aval), index)
        if iv:
            piv.insert(iv)
    iv, _ = _to_int_vec(res, index)
    return not piv.reduce(iv)



def pencil_sympy(kappa, prod):
    """:TT: + kappa :G dG: (with `prod` as the composite product),
    sympy coefficients."""
    P = {}
    for m, c in prod(T, T).items():
        P[m] = sp.Rational(c.numerator, c.denominator)
    for m, c in prod(G, d_cur(G)).items():
        P[m] = sp.expand(P.get(m, 0)
                         + kappa * sp.Rational(c.numerator,
                                               c.denominator))
    return {m: c for m, c in P.items() if c != 0}


if __name__ == "__main__":
    from super_engine import normal_prod
    from super_kdv import derive_tower, naive_mult, select_I3

    deep = "--deep" in sys.argv
    naive = "--naive" in sys.argv
    prod = naive_mult if naive else normal_prod
    tag = "_naive" if naive else ""
    out = {"product": "naive-concatenation" if naive else "point-split",
           "stage0": {}, "pencil_conditions": [], "branches": {},
           "selection": None, "locus": None, "verified": {}}
    a, k = sp.symbols('alpha kappa')
    failed = False

    # stage 0: G and T must impose no alpha-conditions.
    c_G = conditions_sympy(G, a)
    c_T = conditions_sympy(T, a)
    out["stage0"] = {"G": [str(c) for c in c_G],
                     "T": [str(c) for c in c_T]}
    print(f"stage 0: alpha-conditions from G: {c_G or 'none'} "
          f"(expect none: top components are G-invariant for all "
          f"alpha)")
    print(f"         alpha-conditions from T: {c_T or 'none'}")
    failed = bool(c_G or c_T)

    locus = None
    if not failed:
        conds = conditions_sympy(pencil_sympy(k, prod), a)
        out["pencil_conditions"] = [str(sp.factor(c)) for c in conds]
        print("stage 1: pencil conditions on (alpha, kappa):")
        for c in conds:
            print(f"           {sp.factor(c)}")
        sols = sp.solve(conds, [a, k], dict=True)
        bare = sorted({str(s[k]) for s in sols
                       if s.get(a) == 0 and k in s})
        out["branches"]["bare_psi_kappa"] = bare
        if naive:
            locus = sp.Poly(2 * a**4 - 9 * a**2 + 2, a)
            kap = {sp.simplify(s[k]) for s in sols
                   if s.get(a, 1) != 0 and k in s}
            ok_k = kap == {sp.Rational(-1, 4)}
            ok_div = all(sp.div(sp.Poly(
                sp.expand(c.subs(k, sp.Rational(-1, 4))), a),
                locus)[1].is_zero for c in conds)
            out["branches"]["nonzero_alpha_kappa"] = [str(x)
                                                     for x in kap]
            print(f"         nonzero-exponential screenings force "
                  f"kappa = -1/4 within the naive slice "
                  f"(ok={ok_k and ok_div}); locus "
                  f"{locus.as_expr()} = 0")
            failed = failed or not (ok_k and ok_div)
        else:
            locus = sp.Poly(a**2 - 1, a)
            ok_pm1 = all(sp.expand(c.subs(a, v)) == 0
                         for c in conds for v in (1, -1))
            # branch-safe certificate (audit round 2): the generic
            # elimination may divide by alpha-dependent pivots, so
            # substitution into its output is not a direct check at a
            # point.  Certify by EXACT rational-point membership:
            # both pencil generators screen at alpha = +-1, hence
            # every kappa by linearity.
            ok_exact = all(twisted_exact_fraction(gen, F(v))
                           for gen in (prod(T, T), prod(G, d_cur(G)))
                           for v in (1, -1))
            out["branches"]["alpha_pm1_screens_every_kappa"] = ok_pm1
            out["branches"]["generators_exact_at_pm1"] = ok_exact
            print(f"         alpha = +-1 satisfies every condition "
                  f"for ALL kappa: {ok_pm1}; exact rational-point "
                  f"certificates on both generators: {ok_exact}")
            failed = failed or not (ok_pm1 and ok_exact)
        # exact certificate for the bare-fermion branch: oint psi
        # annihilates :TT: + 1/2 :G dG: (plain d-exactness at alpha=0)
        bare_ray = dict(prod(T, T))
        cur_add(bare_ray, prod(G, d_cur(G)), F(1, 2))
        ok_bare = twisted_exact_fraction(
            {m: c for m, c in bare_ray.items() if c}, F(0))
        out["branches"]["bare_psi_exact_at_kappa_half"] = ok_bare
        out["locus"] = str(locus.as_expr())
        print(f"         bare-fermion branch (alpha = 0) at kappa = "
              f"{bare}; exact certificate at kappa = 1/2: {ok_bare}")
        failed = failed or not ok_bare

    if not failed:
        kappa_sel, _ = select_I3(prod)
        out["selection"] = (
            f"G-invariance selects kappa = {kappa_sel}"
            + ("" if naive else " (screening does not select in the "
                                "point-split pencil)"))
        print(f"         selection: {out['selection']}")
        failed = failed or kappa_sel != F(-1, 4)

    # stage 2: the locus annihilates the higher tower (generic-stratum
    # statement; see module docstring caveat).
    if not failed:
        tower = derive_tower(prod, weights=(12, 16) if deep else (12,))
        for w2 in (12, 16) if deep else (12,):
            cs = conditions_sympy(tower[w2], a)
            div = all(sp.div(sp.Poly(sp.expand(c), a),
                             locus)[1].is_zero for c in cs)
            a0 = any(sp.expand(c).subs(a, 0) != 0 for c in cs)
            # branch-safe exact certificates at the rational locus
            # points (corrected tower only: the naive locus is
            # irrational, its divisibility statement stays
            # generic-stratum)
            exact = (None if naive else
                     all(twisted_exact_fraction(tower[w2], F(v))
                         for v in (1, -1)))
            name = f"I{w2 // 2 - 1}"
            out["verified"][name] = {"n_conditions": len(cs),
                                     "all_divisible_by_locus": div,
                                     "alpha0_excluded": a0,
                                     "exact_at_pm1": exact}
            print(f"stage 2: {name}: {len(cs)} conditions, all "
                  f"divisible by locus: {div}, alpha=0 excluded: {a0}"
                  + ("" if exact is None else
                     f", exact at alpha=+-1: {exact}"))
            failed = failed or not (div and a0 and exact is not False)

    path = (f"../results/super_screening_fit{tag}.json" if not failed
            else f"../results/super_screening_fit{tag}.partial.json")
    json.dump(out, open(path, "w"), indent=1)
    print("SUPER-SCREENING FIT COMPLETE" if not failed
          else "SUPER-SCREENING FIT FAILED")
    sys.exit(1 if failed else 0)
