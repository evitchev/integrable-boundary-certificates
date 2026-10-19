"""Phase 3: super-KdV with background charge (one N=1 superfield).

Deformation (coefficient 2 forced by the algebra, verified in check 0):

    T_lam = 1/2 (dX)^2 + lam d2X - 1/2 psi dpsi
    G_lam = psi dX + 2 lam dpsi,      Res(G_lam, G_lam) = 2 T_lam.

Note oint G_lam = oint G and oint T_lam = oint T as charges (the
deformations are total derivatives), so the free-field CHARGE algebra
is lambda-independent; what moves with lambda is the composite pencil

    P(kappa, lam) = :T_lam T_lam: + kappa :G_lam dG_lam:

(point-split products) and hence the selected I3(lam) ray and its
screening locus.  Everything here is exact sympy in (alpha, kappa,
lam) on doubled-weight <= 9 sectors, with every headline claim
ASSERTED (fail-closed), not just printed (audit round 3, finding 3).

Computations:
 0. algebra check: Res(G_lam, G_lam) = 2 T_lam exactly, symbolic lam;
 1. selection: Res(G, P(kappa, lam)) in Im d => kappa(lam), asserted
    to equal -1/4 identically;
 2. HYPOTHESIS ADJUDICATION (results/super_kdv_calibration.md sec. 4):
    the naive-pencil ray (dX)^4 - 3 (dX)^2 psi dpsi lies on the
    lambda-family; the solution set is asserted to be EXACTLY
    {kappa = -1/4, s = 1/4, lam = +-sqrt(10)/4}, and both solutions
    are back-substituted into the original membership equations;
 3. screening locus of I3(lam): asserted to be exactly
    F- * F+ = (alpha^2 - 2 lam alpha - 1)(alpha^2 + 2 lam alpha - 1),
    and to reproduce the naive quartic 2 alpha^4 - 9 alpha^2 + 2 at
    lam = +-sqrt(10)/4;
 4. screening-selection structure (audit round 3, finding 1): the
    single pair F- DIVIDES EVERY free-kappa pencil condition, so one
    screening family is kappa-blind at every lam; kappa = -1/4
    adjoins the sign-reflected pair F+.  It is the DUAL S+- pair
    (CLO's psi e^{+-b Phi}, arXiv:2110.05870) that selects the
    super-KdV ray for lam != 0; at lam = 0 the pairs coincide and
    screening cannot select at all.  Certified branch-safely at
    lam = 3/4 (rational roots: F- = {2, -1/2}, F+ = {1/2, -2}) with
    positive AND negative controls via exact integer elimination.

Outputs results/super_lambda.json; fail-closed with .partial.
"""
import json
import sys
from fractions import Fraction as F

import sympy as sp

from super_engine import (cur_add, d_cur, d_mono, gen_basis_super,
                          mono_weight2, normal_prod, residue_cur)
from super_screening_fit import conditions_sympy, twisted_exact_fraction

NF = 1
DX, D2X, PSI, DPSI = (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 1, 0)
lam, a, kap, s_ = sp.symbols('lam alpha kappa s')
G = {(DX, PSI): sp.Integer(1)}
NAIVE_I3 = {(DX, DX, DX, DX): sp.Integer(1),
            (DX, DX, PSI, DPSI): sp.Integer(-3)}
QUARTIC = 2 * a**4 - 9 * a**2 + 2
F_MINUS = a**2 - 2 * a * lam - 1
F_PLUS = a**2 + 2 * a * lam - 1


def T_l(l):
    return {(DX, DX): sp.Rational(1, 2), (D2X,): l,
            (PSI, DPSI): sp.Rational(-1, 2)}


def G_l(l):
    return {(DX, PSI): sp.Integer(1), (DPSI,): 2 * l}


def clean(cur):
    return {m: sp.expand(c) for m, c in cur.items()
            if sp.expand(c) != 0}


def membership_eqs(cur, unknowns_prefix):
    """Equations for cur in Im d: cur = sum u_i d(b_i) on its sector.
    Returns (equations, derivative-coefficient unknowns)."""
    w2 = mono_weight2(next(iter(cur)))
    tgt = gen_basis_super(w2, NF)
    lower = gen_basis_super(w2 - 2, NF)
    us = sp.symbols(f'{unknowns_prefix}0:{len(lower)}')
    acc = {m: -sp.expand(c) for m, c in cur.items()}
    for u, b in zip(us, lower):
        for m, c in d_mono(b).items():
            acc[m] = acc.get(m, 0) + u * c
    eqs = [sp.expand(acc.get(m, 0)) for m in tgt]
    return eqs, list(us)


def to_fraction_cur(cur):
    return {m: F(int(sp.numer(c)), int(sp.denom(c)))
            for m, c in cur.items() if c != 0}


def main():
    out = {"checks": {}, "kappa_of_lambda": None, "hypothesis": {},
           "locus": {}, "selection": {}, "anomalies": []}

    # 0. algebra: Res(G_lam, G_lam) = 2 T_lam exactly, symbolic lam
    S = clean(residue_cur(G_l(lam), G_l(lam)))
    ok0 = S == clean({m: 2 * c for m, c in T_l(lam).items()})
    out["checks"]["GG_eq_2T_symbolic"] = ok0
    print(f"check 0: Res(G_lam, G_lam) = 2 T_lam exactly "
          f"(symbolic lam): {ok0}", flush=True)
    if not ok0:
        out["anomalies"].append("deformed algebra check failed")
        finish(out, failed=True)

    TTl = clean(normal_prod(T_l(lam), T_l(lam)))
    GdGl = clean(normal_prod(G_l(lam), d_cur(G_l(lam))))

    # 1. kappa(lam) from G-invariance — asserted equal to -1/4
    cur = dict(clean(residue_cur(G, TTl)))
    for m, c in clean(residue_cur(G, GdGl)).items():
        cur[m] = sp.expand(cur.get(m, 0) + kap * c)
    eqs1, us1 = membership_eqs(clean(cur), 'u')
    sols = sp.solve(eqs1, [kap] + us1, dict=True)
    kl = ({sp.cancel(s[kap]) for s in sols if kap in s}
          if sols else set())
    if kl != {sp.Rational(-1, 4)}:
        out["anomalies"].append(f"kappa(lam) != -1/4: {kl}")
        finish(out, failed=True)
    out["kappa_of_lambda"] = "-1/4"
    print("check 1: G-invariance selects kappa(lam) = -1/4 "
          "(asserted, identically in lam)", flush=True)

    # 2. hypothesis: P(kappa, lam) = s * I3_naive mod Im d —
    #    solution set asserted EXACTLY, then back-substituted
    cur = dict(TTl)
    for m, c in GdGl.items():
        cur[m] = sp.expand(cur.get(m, 0) + kap * c)
    for m, c in NAIVE_I3.items():
        cur[m] = sp.expand(cur.get(m, 0) - s_ * c)
    eqs2, us2 = membership_eqs(clean(cur), 'v')
    sols = [s for s in sp.solve(eqs2, [kap, s_, lam] + us2, dict=True)
            if s.get(s_) != 0]
    got = {(sp.simplify(s[kap]), sp.simplify(s[s_]),
            sp.simplify(s[lam])) for s in sols}
    expect = {(sp.Rational(-1, 4), sp.Rational(1, 4),
               sp.sqrt(10) / 4),
              (sp.Rational(-1, 4), sp.Rational(1, 4),
               -sp.sqrt(10) / 4)}
    ok_set = got == expect
    ok_back = all(sp.simplify(e.subs(s)) == 0
                  for s in sols for e in eqs2)
    out["hypothesis"] = {
        "solutions": sorted(str(g) for g in got),
        "solution_set_exact": ok_set,
        "back_substitution_ok": ok_back}
    print(f"check 2: naive ray on the lambda-family: solutions "
          f"(kappa, s, lam) = {sorted(str(g) for g in got)}",
          flush=True)
    print(f"         exact expected set: {ok_set}; back-substitution "
          f"into original equations: {ok_back}", flush=True)
    if not (ok_set and ok_back):
        out["anomalies"].append("hypothesis solution set/back-subst")
        finish(out, failed=True)

    # 3. screening locus of I3(lam) — asserted F- * F+
    P = dict(TTl)
    for m, c in GdGl.items():
        P[m] = sp.expand(P.get(m, 0) - sp.Rational(1, 4) * c)
    conds = conditions_sympy(clean(P), a)
    out["locus"]["conditions"] = [str(sp.factor(c)) for c in conds]
    gexpr = sp.Integer(0)
    for c in conds:
        gexpr = sp.gcd(gexpr, sp.expand(c))
    q3, r3 = sp.div(sp.expand(gexpr), sp.expand(F_MINUS * F_PLUS), a)
    ok_locus = (r3 == 0 and sp.degree(gexpr, a) == 4)
    out["locus"]["gcd"] = str(sp.factor(gexpr))
    out["locus"]["is_Fminus_times_Fplus"] = ok_locus
    print(f"check 3: locus of I3(lam) = {sp.factor(gexpr)}; equals "
          f"F-*F+: {ok_locus}", flush=True)
    ok_q = True
    for sgn in (1, -1):
        gl = sp.expand(gexpr.subs(lam, sgn * sp.sqrt(10) / 4))
        lc = sp.Poly(gl, a).LC()
        okq = sp.expand(2 * gl / lc - QUARTIC) == 0
        ok_q = ok_q and okq
        print(f"         at lam = {sgn}*sqrt(10)/4: matches naive "
              f"quartic: {okq}", flush=True)
    out["locus"]["matches_naive_quartic_at_lambda_star"] = ok_q
    if not (ok_locus and ok_q):
        out["anomalies"].append("locus assertion failed")

    # 4. screening-selection structure (audit round 3, finding 1)
    Pk = dict(TTl)
    for m, c in GdGl.items():
        Pk[m] = sp.expand(Pk.get(m, 0) + kap * c)
    conds_k = conditions_sympy(clean(Pk), a)
    out["selection"]["conditions"] = [str(sp.factor(c))
                                      for c in conds_k]
    print("check 4: pencil conditions with free (kappa, lam):",
          flush=True)
    for c in conds_k:
        print(f"         {sp.factor(c)}", flush=True)
    # (a) F- divides EVERY condition: the single pair is kappa-blind
    ok_fm = all(sp.div(sp.expand(c), F_MINUS, a)[1] == 0
                for c in conds_k)
    # (b) at kappa = -1/4 the locus is exactly the dual pair F- * F+
    g14 = sp.Integer(0)
    for c in conds_k:
        g14 = sp.gcd(g14, sp.expand(c.subs(kap, sp.Rational(-1, 4))))
    ok_pair = (sp.div(sp.expand(g14),
                      sp.expand(F_MINUS * F_PLUS), a)[1] == 0
               and sp.degree(g14, a) == 4)
    # (c) branch-safe certificate at lam = 3/4 (rational roots),
    #     positive AND negative controls, exact integer elimination
    lamv = sp.Rational(3, 4)
    TT34 = clean(normal_prod(T_l(lamv), T_l(lamv)))
    GdG34 = clean(normal_prod(G_l(lamv), d_cur(G_l(lamv))))

    def pencil34(kv):
        P34 = dict(to_fraction_cur(TT34))
        for m, c in to_fraction_cur(GdG34).items():
            P34[m] = P34.get(m, F(0)) + kv * c
        return {m: c for m, c in P34.items() if c}

    fm_roots = (F(2), F(-1, 2))      # F- roots at lam = 3/4
    fp_roots = (F(1, 2), F(-2))      # F+ roots at lam = 3/4
    table = {}
    ok_table = True
    for kv in (F(0), F(1), F(-1, 4)):
        row = {}
        for al in fm_roots + fp_roots:
            screens = twisted_exact_fraction(pencil34(kv), al)
            expected = (al in fm_roots) or kv == F(-1, 4)
            row[str(al)] = screens
            ok_table = ok_table and screens == expected
        table[str(kv)] = row
    out["selection"]["F_minus_divides_every_condition"] = ok_fm
    out["selection"]["kappa_m14_locus_is_dual_pair"] = ok_pair
    out["selection"]["branch_safe_table_lam_3_4"] = table
    out["selection"]["branch_safe_table_ok"] = ok_table
    out["selection"]["statement"] = (
        "the single pair F- is kappa-blind at every lam; kappa = -1/4 "
        "adjoins the sign-reflected pair F+; the DUAL S+- pair "
        "selects the super-KdV ray for lam != 0; at lam = 0 the "
        "pairs coincide and screening cannot select")
    print(f"         F- divides every condition (single pair "
          f"kappa-blind at all lam): {ok_fm}", flush=True)
    print(f"         kappa=-1/4 locus is the dual pair F-*F+: "
          f"{ok_pair}", flush=True)
    print(f"         branch-safe lam=3/4 table (pos+neg controls): "
          f"{ok_table}  {table}", flush=True)
    print(f"         => {out['selection']['statement']}", flush=True)
    if not (ok_fm and ok_pair and ok_table):
        out["anomalies"].append("selection-structure assertion failed")

    verdict = ("CONFIRMED: naive I3 ray = I3(lambda*) member "
               "(ray-level; higher stored charges pinned by "
               "centralizer equality + multiplicity one in range)"
               if not out["anomalies"] else "FAILED")
    out["hypothesis"]["verdict"] = verdict
    print(f"\nhypothesis verdict: {verdict}", flush=True)
    finish(out, failed=bool(out["anomalies"]))


def finish(out, failed):
    path = ("../results/super_lambda.json" if not failed
            else "../results/super_lambda.partial.json")
    json.dump(out, open(path, "w"), indent=1)
    print("SUPER-LAMBDA PHASE-3 COMPLETE" if not failed
          else "SUPER-LAMBDA PHASE-3 FAILED", flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
