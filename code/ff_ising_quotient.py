"""FF quotient experiment I: the A2(2) spin-5 charge dies in the
Ising quotient -- on the physical t = +7 sheet only (the certified
one-sheet resolution below; the t = -7 charge survives as a
nonzero class).

The certified dictionary (a2_infinitude_certificate.py) gives, in
the convention <dX dX> = (z-w)^-2 with T = (dX)^2/2 + Q d^2X,
Q = a/2 - 1/a, c = 1 - 12Q^2 = 13 - 3a^2 - 12/a^2:

  - Vitchev's parameters: c = N, s^2 = (N-25)(N-1), with the curve
    N = (t^2-25)/(t^2-1), s = -24t/(t^2-1).  At t = +-7:
    N = 24/48 = 1/2 -- the ISING central charge -- and s = -+7/2.
  - The two dual Coulomb-gas points with c = 1/2 are a^2 = 3/2
    (sheet coordinate 3a^2 - 12/a^2 = -7/2, i.e. t = +7) and
    a^2 = 8/3 (+7/2, i.e. t = -7).
  - Vitchev's I5 density is the T-composite pencil
    P6 = 8 :T^3: + K :(dT)^2:, with the certified sheet
    coefficients K/8 = (12a^4-5a^2+3)/(8a^2) (sheet 12) and
    (3a^4-20a^2+192)/(32a^2) (sheet 21).

The Feigin--Frenkel picture predicts: the free-field (Fock) charge
exists at every t, but at a minimal-model point its image in the
IRREDUCIBLE quotient can die.  The Ising vacuum module has its
first nontrivial singular vector at level (3-1)(4-1) = 6 --
exactly the weight of the I5 density.  This certificate tests, in
exact arithmetic:

 1. dictionary coherence: c(a) = 13 - 3a^2 - 12/a^2 re-derived
    from ope_pole(T,T,4) = c/2; curve_point(+-7) = (1/2, -+7/2);
    the two Ising a-points hit c = 1/2 with matching sheet signs;
 2. the level-6 vacuum NULL FIELD derived in-engine: the primary
    weight-6 T-composite (all T(z)N(w) poles 3..8 vanish) on the
    4-dim composite basis {:T^3:, :(dT)^2:, :(d^2T)T:, d^4T} --
    dimension 1 at both Ising points, dimension 0 at the control
    a = 3 (c = -46/3, not a level-6-null central charge);
 3. THE VERDICT (computed, not pre-asserted): for each Ising point
    and each sheet, whether the P6 charge class equals the null
    charge class modulo total derivatives -- equality means the
    spin-5 charge is a null charge, i.e. VANISHES in the Ising
    quotient; a discriminating control (:T^3: alone must NOT be
    proportional to the null class) guards against a degenerate
    test.

Scope: t = +-7 (Ising) only.  The tricritical experiment (t = +-9,
c = 7/10) needs the level-12 vacuum null and is filed separately.
"""
import sys

import sympy as sp

from invariant_engine import curve_point
from ope_engine import d_mono, gen_basis
from ope_poles import ope_pole
from super_engine import d_cur as _sd_cur
from super_engine import normal_prod

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


a = sp.Symbol('a', nonzero=True)
Q = a / 2 - 1 / a
X1 = ((0, 1, 0),)
X2 = ((0, 2, 0),)
T_SYM = {tuple(sorted(X1 + X1)): sp.Rational(1, 2), X2: Q}


def dsc(cur):
    return {m: sp.expand(c) for m, c in _sd_cur(cur).items()
            if sp.expand(c) != 0}


def at_a(cur, aval):
    out = {}
    for m, c in cur.items():
        v = sp.simplify(c.subs(a, aval))
        if v != 0:
            out[m] = v
    return out


def bos(cur):
    return {tuple((m, j) for _, m, j in mono): c
            for mono, c in cur.items()}


# ---- 1. dictionary coherence ---------------------------------------
c_sym = sp.factor(2 * ope_pole(T_SYM, T_SYM, 4).get((), 0))
check("c(a) = 13 - 3a^2 - 12/a^2 from ope_pole(T,T,4)",
      sp.simplify(c_sym - (13 - 3 * a**2 - 12 / a**2)) == 0)
Npt, spt = curve_point(7)
check("curve_point(7) = (1/2, -7/2): the Ising central charge",
      Npt == sp.Rational(1, 2) and spt == sp.Rational(-7, 2))
Nm, sm = curve_point(-7)
check("curve_point(-7) = (1/2, +7/2)",
      Nm == sp.Rational(1, 2) and sm == sp.Rational(7, 2))
A_T7 = sp.sqrt(6) / 2                       # a^2 = 3/2
A_TM7 = 2 * sp.sqrt(6) / 3                  # a^2 = 8/3
for aval, sheet_s, lbl in ((A_T7, sp.Rational(-7, 2), "t=+7"),
                           (A_TM7, sp.Rational(7, 2), "t=-7")):
    check(f"{lbl}: c = 1/2 and sheet coordinate = s",
          sp.simplify(c_sym.subs(a, aval) - sp.Rational(1, 2)) == 0
          and sp.simplify((3 * a**2 - 12 / a**2).subs(a, aval)
                          - sheet_s) == 0)

# ---- composite basis -----------------------------------------------
TT = normal_prod(T_SYM, T_SYM)
C1 = normal_prod(T_SYM, TT)                 # :T :TT::
C2 = normal_prod(dsc(T_SYM), dsc(T_SYM))    # :dT dT:
C3 = normal_prod(dsc(dsc(T_SYM)), T_SYM)    # :d^2T T:
C4 = dsc(dsc(dsc(dsc(T_SYM))))              # d^4 T
COMPS = [C1, C2, C3, C4]

basis6 = gen_basis(6, 1, z2=False)
idx6 = {m: i for i, m in enumerate(basis6)}


def vec6(cur_b, aval=None):
    v = [sp.Integer(0)] * len(basis6)
    for m, c in cur_b.items():
        v[idx6[m]] += c
    return sp.Matrix([v])


def comp_vecs(aval):
    return [vec6(bos(at_a(c, aval))) for c in COMPS]


check("composite basis independent (a = 3 control point)",
      sp.Matrix.vstack(*comp_vecs(sp.Integer(3))).rank() == 4)


# ---- 2. the level-6 vacuum null field ------------------------------
def null_rays(aval):
    """Coefficient rays (x1..x4) making x.C primary: all poles 3..8
    of T(z) N(w) vanish."""
    Ta = at_a(T_SYM, aval)
    Cs = [at_a(c, aval) for c in COMPS]
    xs = sp.symbols('x0:4')
    conds = []
    for k in range(3, 9):
        acc = {}
        for x, comp in zip(xs, Cs):
            for m, cf in ope_pole(Ta, comp, k).items():
                acc[m] = sp.expand(acc.get(m, 0) + x * cf)
        conds += [v for v in acc.values() if v != 0]
    M, _ = sp.linear_eq_to_matrix(conds, xs)
    return M.nullspace()


ns_t7 = null_rays(A_T7)
ns_tm7 = null_rays(A_TM7)
ns_ctrl = null_rays(sp.Integer(3))
check("level-6 null: dim 1 at t=+7 point (a^2=3/2)", len(ns_t7) == 1)
check("level-6 null: dim 1 at t=-7 point (a^2=8/3)", len(ns_tm7) == 1)
check("level-6 null: dim 0 at control a=3 (c=-46/3)",
      len(ns_ctrl) == 0)

# ---- mod-d reduction -----------------------------------------------
drows = []
for bmono in gen_basis(5, 1, z2=False):
    v = [sp.Integer(0)] * len(basis6)
    for m, c in d_mono(bmono).items():
        v[idx6[m]] += c
    drows.append(v)
DMAT = sp.Matrix(drows)


def charge_class_rank(rows):
    """rank of the span of the given class vectors together with
    Im d, minus rank(Im d): the charge-class rank."""
    base = DMAT.rank()
    return sp.Matrix.vstack(DMAT, *rows).rank() - base


PHI12 = (12 * a**4 - 5 * a**2 + 3) / (8 * a**2)
PHI21 = (3 * a**4 - 20 * a**2 + 192) / (32 * a**2)


def p6_vec(phi, aval):
    cur = {}
    for m, cf in C1.items():
        cur[m] = sp.expand(8 * cf)
    for m, cf in C2.items():
        cur[m] = sp.expand(cur.get(m, 0) + 8 * phi * cf)
    return vec6(bos(at_a(cur, aval)))


# ---- 3. the verdict ------------------------------------------------
print()
verdicts = {}
PHYS = {("t=+7", "sheet21"): "physical t=+7 (s=-7/2)",
        ("t=+7", "sheet12"): "physical t=-7 (s=+7/2)",
        ("t=-7", "sheet12"): "physical t=+7 (s=-7/2)",
        ("t=-7", "sheet21"): "physical t=-7 (s=+7/2)"}
REALZ = {"t=+7": "realization a^2=3/2", "t=-7": "realization a^2=8/3"}
for aval, tlbl, ns in ((A_T7, "t=+7", ns_t7), (A_TM7, "t=-7", ns_tm7)):
    # build the null current from the ray
    ray = ns[0]
    ncur = {}
    for x, comp in zip(ray, [at_a(c, aval) for c in COMPS]):
        for m, cf in comp.items():
            ncur[m] = sp.expand(ncur.get(m, 0) + x * cf)
    nvec = vec6(bos({m: c for m, c in ncur.items() if c != 0}))
    check(f"{tlbl}: null charge class is itself nonzero mod d",
          charge_class_rank([nvec]) == 1)
    check(f"{tlbl}: discriminating control -- :T^3: alone is NOT "
          f"proportional to the null class",
          charge_class_rank([vec6(bos(at_a(C1, aval))), nvec]) == 2)
    for phi, slbl in ((PHI12, "sheet12"), (PHI21, "sheet21")):
        pv = p6_vec(sp.simplify(phi.subs(a, aval)), aval)
        pr = charge_class_rank([pv])
        joint = charge_class_rank([pv, nvec])
        verdicts[(tlbl, slbl)] = (pr, joint)
        word = ("DIES in the Ising quotient" if (pr, joint) == (1, 1)
                else "SURVIVES as a nonzero class"
                if (pr, joint) == (1, 2) else "UNEXPECTED / inspect")
        print(f"{REALZ[tlbl]}, {PHYS[(tlbl, slbl)]}: P6 class "
              f"rank {pr}, span with null rank {joint} -> {word}",
              flush=True)

# the discovered pattern, baked after computation: at each dual
# Coulomb-gas realization exactly ONE sheet charge is null, and via
# the certified phi<->s map (phi21 <-> s = +sheet_coord, phi12 <->
# s = -sheet_coord) BOTH dying charges are Vitchev's s = -7/2 sheet
# (t = +7), with K = 225/8; the s = +7/2 charge (K = 15) survives.
check("dying entries have exact ranks (1,1); surviving entries "
      "(1,2) -- nonzero class NOT in the null ray (review: pr=0 "
      "must not count as survival)",
      verdicts[("t=+7", "sheet21")] == (1, 1)
      and verdicts[("t=-7", "sheet12")] == (1, 1)
      and verdicts[("t=+7", "sheet12")] == (1, 2)
      and verdicts[("t=-7", "sheet21")] == (1, 2))
check("both dying charges are the SAME Vitchev sheet: K = 225/8 "
      "(s = -7/2, t = +7)",
      sp.simplify(8 * PHI21.subs(a, A_T7)
                  - sp.Rational(225, 8)) == 0
      and sp.simplify(8 * PHI12.subs(a, A_TM7)
                      - sp.Rational(225, 8)) == 0)
check("the surviving sheet has K = 15 (s = +7/2, t = -7)",
      sp.simplify(8 * PHI12.subs(a, A_T7) - 15) == 0
      and sp.simplify(8 * PHI21.subs(a, A_TM7) - 15) == 0)

print()
if failures:
    print(f"FAILURES: {failures}")
    print("FF ISING QUOTIENT EXPERIMENT FAILED")
    sys.exit(1)

print("VERDICT: ONE SHEET DIES.  At both dual Coulomb-gas "
      "realizations of the Ising point, Vitchev's s = -7/2 spin-5 "
      "charge (t = +7; P6 = 8:T^3: + (225/8):(dT)^2:) is "
      "proportional to the level-6 vacuum null field modulo total "
      "derivatives -- a null charge, vanishing in the Ising "
      "quotient -- while the deck-conjugate s = +7/2 charge "
      "(K = 15) remains a nonzero class.  This is exactly the "
      "Feigin--Frenkel degeneration the E8 spin list predicts: the "
      "generic A2(2) spins 6n+-1 contain 5, the Ising phi_12 (E8) "
      "list 1,7,11,13,... omits it, and the certificate exhibits "
      "the mechanism -- the spin-5 charge class equals the null-"
      "field class modulo total derivatives.  (FF Sec. 4.8.3 report "
      "the degree-5 dropout, citing [Zam, fz] for its discovery, "
      "with the E8 pattern conjectural; the explicit density/null-"
      "class exhibit in Vitchev normalization and the one-sheet "
      "resolution are new IN THIS PROJECT -- earlier literature "
      "has not been exhaustively checked for the same exhibit.)")
print("FF ISING QUOTIENT EXPERIMENT COMPLETE")
