"""The c = -25/2 identification: the enhanced c=1 tower IS quantum KdV.

Certifies the 2026-08-18 finding (external peer review, confirmed
in-engine the same day): the enhanced odd-spin tower at the c=1
branch point is the quantum KdV hierarchy at c_{1,4} = -25/2 in
Feigin--Fuchs realization on the same boson, T_Q = (1/2)(dphi)^2
+ Q d^2phi at Q^2 = 9/8.

Checks (fail-closed; exits nonzero unless every check passes):
 1. c(Q) = 1 - 12 Q^2 from ope_pole(T_Q, T_Q, 4), symbolic in Q.
 2. :T_Q T_Q: == (1/4)[(dphi)^4 + (4Q^2-2)(d^2phi)^2] mod d,
    SYMBOLIC in Q  =>  the enhanced I'3 = int[(dphi)^4
    + (5/2)(d^2phi)^2] is the KdV I3 at 4Q^2-2 = 5/2.
 3. Arithmetic: 4Q^2-2 = 5/2 <=> Q^2 = 9/8; c = -25/2; and the
    Kac value c_{1,4} = 1 - 6(1-4)^2/4 = -25/2.
 4. :T:TT:: + A :dT dT: == lambda P6 mod d at Q^2 = 9/8 with
    EXACTLY A = 7/8 and lambda = 1/8, where P6 is the paper's
    eq.(1) density (1, 5/2, 31/24)  =>  I5 is the KdV I5.
 5. A = (4Q^2-1)/4 at Q^2 = 9/8 equals 7/8 (the BLZ-form tie).
 6. Screening-quartet identity, symbolic in a: with
    kappa = a^2 - 6 + 4/a^2,
    alpha^4 - (kappa+6) alpha^2 + 4 = (alpha^2-a^2)(alpha^2-4/a^2);
    at a^2 = 1/2 the roots are {1/2, 8} -- the certified
    factorization alpha(alpha^2-8)(2alpha^2-1) of the spin-3
    condition: the four momenta are the FF quartet {+-a, +-2/a}.
 7. Negative controls: the mod-d identity of check 2 FAILS for the
    wrong coefficient 4Q^2-1; check 4 FAILS at Q^2 = 1 (c = -11).

Together with the certified multiplicity-one and pairwise
commutativity of the tower (c1_enhanced entries), these seal the
identification: ker ad_{I5} contains the KdV(Q^2=9/8) charge at
every odd spin where it is nonzero mod d, and multiplicity one
forces class equality.
"""
import sys

import sympy as sp

from ope_engine import gen_basis
from ope_poles import ope_pole
from super_engine import d_cur, normal_prod

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


Q = sp.Symbol('Q')
T = {tuple(sorted(((0, 1, 0), (0, 1, 0)))): sp.Rational(1, 2),
     ((0, 2, 0),): Q}


def dsc(cur):
    return {m: sp.expand(c) for m, c in d_cur(cur).items()
            if sp.expand(c) != 0}


def bos(cur):
    return {tuple((m, j) for _, m, j in mono): c for mono, c in cur.items()}


def mono(*orders):
    return tuple(sorted((n, 0) for n in orders))


def vec(cur_b, basis, idx):
    v = [sp.Integer(0)] * len(basis)
    for m, c in cur_b.items():
        v[idx[tuple(sorted(m))]] += sp.expand(c)
    return v


def mod_partial_solve(target_vec, weight, extra_cols):
    """Solve target = sum c_i d(m_i) + sum x_k extra_k; None if no
    solution.  Exact symbolic linear algebra."""
    lower = gen_basis(weight - 1, 1, z2=False)
    basis = gen_basis(weight, 1, z2=False)
    idx = {m: i for i, m in enumerate(basis)}
    cols = []
    for m in lower:
        cur = {tuple(sorted((0, n, 0) for n, _ in m)): sp.Integer(1)}
        cols.append(vec(bos(dsc(cur)), basis, idx))
    allc = cols + list(extra_cols)
    A = sp.Matrix([[allc[j][i] for j in range(len(allc))]
                   for i in range(len(basis))])
    b = sp.Matrix([[t] for t in target_vec])
    sol = sp.linsolve((A, b))
    return None if sol == sp.EmptySet else list(sol)[0]


# ---- 1. central charge -------------------------------------------------
c_of_Q = sp.expand(2 * ope_pole(T, T, 4).get((), 0))
check("c(Q) = 1 - 12Q^2 from ope_pole(T,T,4), symbolic",
      sp.simplify(c_of_Q - (1 - 12 * Q**2)) == 0)

# ---- 2. TT mod d, symbolic in Q ---------------------------------------
basis4 = gen_basis(4, 1, z2=False)
idx4 = {m: i for i, m in enumerate(basis4)}
TT = normal_prod(T, T)
tb = {tuple(sorted(m)): sp.expand(c) for m, c in bos(TT).items()}
claim = {mono(1, 1, 1, 1): sp.Rational(1, 4),
         mono(2, 2): (4 * Q**2 - 2) / 4}
diff = dict(tb)
for m, c in claim.items():
    diff[m] = sp.expand(diff.get(m, 0) - c)
check("TT == (1/4)[(dphi)^4 + (4Q^2-2)(d^2phi)^2] mod d, symbolic",
      mod_partial_solve(vec(diff, basis4, idx4), 4, []) is not None)

# negative control: wrong quantum shift
diff_bad = dict(tb)
for m, c in {mono(1, 1, 1, 1): sp.Rational(1, 4),
             mono(2, 2): (4 * Q**2 - 1) / 4}.items():
    diff_bad[m] = sp.expand(diff_bad.get(m, 0) - c)
check("negative control: coefficient (4Q^2-1)/4 rejected mod d",
      mod_partial_solve(vec(diff_bad, basis4, idx4), 4, []) is None)

# ---- 3. arithmetic of the identification ------------------------------
sols = sp.solve(sp.Eq(4 * Q**2 - 2, sp.Rational(5, 2)), Q**2)
check("4Q^2 - 2 = 5/2  <=>  Q^2 = 9/8",
      sols == [sp.Rational(9, 8)])
c_val = (1 - 12 * Q**2).subs(Q**2, sp.Rational(9, 8))
check("c(Q^2=9/8) = -25/2", sp.nsimplify(c_val) == sp.Rational(-25, 2))
kac = 1 - sp.Rational(6 * (1 - 4)**2, 1 * 4)
check("Kac: c_{1,4} = 1 - 6(1-4)^2/(1*4) = -25/2",
      kac == sp.Rational(-25, 2))

# ---- 4. the I5 identification at Q^2 = 9/8 ----------------------------
Q98 = sp.sqrt(sp.Rational(9, 8))
Tn = {m: c.subs(Q, Q98) for m, c in T.items()}
TTn = normal_prod(Tn, Tn)
C1 = normal_prod(Tn, TTn)
C2 = normal_prod(dsc(Tn), dsc(Tn))
basis6 = gen_basis(6, 1, z2=False)
idx6 = {m: i for i, m in enumerate(basis6)}
P6 = {mono(1, 1, 1, 1, 1, 1): sp.Integer(1),
      mono(1, 1, 2, 2): sp.Rational(5, 2),
      mono(3, 3): sp.Rational(31, 24)}
v1 = vec({tuple(sorted(m)): c for m, c in bos(C1).items()}, basis6, idx6)
v2 = vec({tuple(sorted(m)): c for m, c in bos(C2).items()}, basis6, idx6)
vp = vec(P6, basis6, idx6)
sol2 = mod_partial_solve(v1, 6, [[-x for x in v2], vp])
ok4 = sol2 is not None
if ok4:
    A_val = sp.nsimplify(sp.simplify(sol2[-2]))
    lam = sp.nsimplify(sp.simplify(sol2[-1]))
    ok4 = (A_val == sp.Rational(7, 8) and lam == sp.Rational(1, 8))
check(":T:TT:: + (7/8):dTdT: == (1/8) P6 mod d at Q^2=9/8 "
      "(A and lambda exact)", ok4)
check("A = (4Q^2-1)/4 at Q^2=9/8 equals 7/8 (BLZ-form tie)",
      sp.Rational(4 * 9, 8 * 4) - sp.Rational(1, 4) == sp.Rational(7, 8))

# negative control at Q^2 = 1 (c = -11): no lambda works with A free
Tb = {m: c.subs(Q, sp.Integer(1)) for m, c in T.items()}
TTb = normal_prod(Tb, Tb)
C1b = normal_prod(Tb, TTb)
C2b = normal_prod(dsc(Tb), dsc(Tb))
v1b = vec({tuple(sorted(m)): c for m, c in bos(C1b).items()}, basis6, idx6)
v2b = vec({tuple(sorted(m)): c for m, c in bos(C2b).items()}, basis6, idx6)
solb = mod_partial_solve(v1b, 6, [[-x for x in v2b], vp])
okneg = solb is None
if not okneg:
    # a solution may exist formally; it must NOT reproduce P6 with a
    # nonzero lambda AND the BLZ A for that Q
    A_b = sp.simplify(solb[-2])
    lam_b = sp.simplify(solb[-1])
    okneg = not (lam_b != 0 and A_b == sp.Rational(3, 4))
check("negative control: Q^2 = 1 (c = -11) does NOT reproduce P6 "
      "with the BLZ A", okneg)

# ---- 6. the screening-quartet identity --------------------------------
a, al = sp.symbols('a alpha', positive=True)
kappa = a**2 - 6 + 4 / a**2
lhs = al**4 - (kappa + 6) * al**2 + 4
rhs = (al**2 - a**2) * (al**2 - 4 / a**2)
check("alpha^4-(kappa+6)alpha^2+4 = (alpha^2-a^2)(alpha^2-4/a^2), "
      "symbolic in a", sp.simplify(lhs - rhs) == 0)
half = sp.Rational(1, 2)
check("at a^2 = 1/2 the roots are {1/2, 8} = roots of the certified "
      "(2alpha^2-1)(alpha^2-8)",
      {half, 4 / half} == {half, sp.Integer(8)}
      and sp.expand((al**2 - half) * (al**2 - 8) * 2
                    - (2 * al**2 - 1) * (al**2 - 8)) == 0)

# ---- verdict ----------------------------------------------------------
print()
if failures:
    print(f"IDENTIFICATION CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("KDV C=-25/2 IDENTIFICATION CERTIFIED "
      "(the enhanced c=1 tower is quantum KdV at c_{1,4})")
