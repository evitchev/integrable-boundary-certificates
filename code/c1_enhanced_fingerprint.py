"""Screening fingerprint of the enhanced c=1 tower (KL-menagerie
identification, stage 1).

The A2(2) family at the branch point N=1 (c=1, s=0) acquires an
enhanced all-odd centralizer; its spin-3 charge I3' has density
(dphi)^4 - (5/2) dphi d3phi == (dphi)^4 + (5/2)(d2phi)^2 mod d.
This script certifies, in exact arithmetic:

 1. I3' is the unique genuine spin-3 charge in ker ad_{I5} at N=1,
    with exactly that ray;
 2. I3' is NOT the c=1 quantum KdV ray ((dphi)^4 - 2(d2phi)^2), and
    [I3_KdV, I5(c=1)] is not a total derivative;
 3. generic screening fit: the alpha-conditions of I3' factor as
    a (a^2-8)(2a^2-1) -- and do NOT contain the KdV factor (a^2-2);
 4. branch-safe exact verification (numeric algebraic alpha, no
    parameter pivots) with a non-root control:
      +-1/sqrt2  annihilate I3' and I5   (the tower's screenings);
      +-2 sqrt2  annihilate I3' and I5   (their Liouville duals);
      +-sqrt2    annihilate I5 ONLY      (the KdV pair!);
      sqrt3      annihilates neither     (control).

All membership statements are finite-range (charges I3'..I9' as
computed); no all-spin claim is made or implied.  CANDIDATE
mechanism (asserted arithmetic, causality not proven): at c=1 the
Coulomb-gas parameter satisfies a^2 = 2, a double root of the A2(2)
dictionary, where gamma_12 = 1/a and gamma_21 = -a/2 merge into the
symmetric pair +-1/sqrt2.  The self-dual degeneracy is PROVEN here
(review round): the joint {+-sqrt2} kernel at density weight 6 has
charge dimension exactly 2, containing the independent classes of
I5_KdV and I5_A2.  Identification target (KL 1910.05947, their eq.
labeled Get): the Getmanov integrable structure, all-odd IM on the
cigar W-infinity algebra with I3^G ~ (k+6) W4 + 4(2k+1) W2^2; the
quantum-AKNS candidate is disfavored (all-integer spins, two-boson
hairpin algebra vs our all-odd one-boson tower)."""
import sys
from fractions import Fraction as F

import sympy as sp

from families import L
from ope_engine import Sector, gen_basis, residue_cur
from screening_fit import screening_conditions
from solve_charges import genuine_kernel

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def acc(pairs):
    out = {}
    for mono, c in pairs:
        out[mono] = out.get(mono, 0) + c
    return {m: c for m, c in out.items() if c}


# A2(2) I5 density at N=1 (eq. (34) with the (22)(11)==(12)(12)
# collision accumulated; cf. a2_branch_check.py)
A2_P6 = acc([(L((3, 0), (3, 0)), F(124, 96)),
             (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
             (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-68, 8)),
             (L((1, 0),) * 6, F(1))])

# 1. I3'
ker3 = genuine_kernel(A2_P6, gen_basis(4, 1, True),
                      Sector(9, 1, True), Sector(4, 1, True))
check("multiplicity one at sigma=3", len(ker3) == 1)
I3p = ker3[0]
norm = I3p.get(L((1, 0), (1, 0), (1, 0), (1, 0)), 0)
check("I3' ray = (1111) - 5/2 (13)",
      norm != 0 and {m: c / norm for m, c in I3p.items()}
      == {L((1, 0),) * 4: F(1), L((1, 0), (3, 0)): F(-5, 2)})

# 2. distinct from KdV
KDV = {L((1, 0),) * 4: F(1), L((2, 0), (2, 0)): F(-2)}
sec4 = Sector(4, 1, True)
diff = dict(I3p)
for m, c in KDV.items():
    diff[m] = diff.get(m, 0) - norm * c
diff = {m: c for m, c in diff.items() if c}
check("I3' is not the c=1 KdV ray",
      bool(diff) and not sec4.is_total_derivative(diff))
S = residue_cur(KDV, A2_P6)
check("[I3_KdV, I5(c=1)] is NOT a total derivative",
      bool(S) and not Sector(9, 1, True).is_total_derivative(S))

# 3. generic fit: the condition equals c * a(a^2-8)(2a^2-1) EXACTLY
a = sp.Symbol('a')
c3 = screening_conditions(I3p, (a,), 1)
target = a * (a**2 - 8) * (2 * a**2 - 1)
ok3 = (len(c3) == 1
       and sp.simplify(sp.cancel(c3[0] / target)).is_constant()
       and sp.simplify(sp.cancel(c3[0] / target)) != 0)
check("I3' condition = const * a(a^2-8)(2a^2-1) exactly", bool(ok3))
kdv_factor = sp.div(sp.Poly(sp.expand(c3[0]), a),
                    sp.Poly(a**2 - 2, a))[1].is_zero
check("I3' condition does NOT carry the KdV factor (a^2-2)",
      not kdv_factor)

# merger arithmetic (candidate mechanism): c(a) = 1 has the DOUBLE
# root a^2 = 2, where gamma_12 = 1/a and -gamma_21 = a/2 coincide
ca = 13 - 3 * a**2 - 12 / a**2
check("c(a) - 1 = -3 (a^2-2)^2 / a^2 (double root at a^2=2)",
      sp.simplify(ca - 1 + 3 * (a**2 - 2)**2 / a**2) == 0)
check("gamma_12 = -gamma_21 = 1/sqrt2 at a = sqrt2",
      sp.simplify(1 / sp.sqrt(2) - sp.sqrt(2) / 2) == 0)

# 4. branch-safe exact table, BOTH signs explicit, with control
def screens(P, v):
    return not [x for x in screening_conditions(P, (v,), 1)
                if sp.simplify(x) != 0]

expect = [(sp.sqrt(2) / 2, True, True), (2 * sp.sqrt(2), True, True),
          (sp.sqrt(2), False, True), (sp.sqrt(3), False, False)]
for v, want3, want5 in expect:
    for sgn in (1, -1):
        u = sgn * v
        check(f"alpha={u}: I3' {want3}, I5 {want5}",
              screens(I3p, u) == want3 and screens(A2_P6, u) == want5)

# 5. the self-dual degeneracy, proven (review round): the joint
# {+-sqrt2} kernel at density weight 6 has charge dimension exactly 2,
# containing the independent classes of I5_KdV and I5_A2
from families import KDV_P6
b6 = gen_basis(6, 1, True)
cs = [sp.Symbol(f'c{i}') for i in range(len(b6))]
Csym = {m: c for m, c in zip(b6, cs)}
conds = []
for sgn in (1, -1):
    conds += screening_conditions(Csym, (sgn * sp.sqrt(2),), 1)
sol = sp.solve(conds, cs, dict=True)
check("joint {+-sqrt2} conditions solvable", len(sol) == 1)
sub = sol[0]
free = [c for c in cs if c not in sub]
sec6 = Sector(6, 1, True)
vecs = []
for f in free:
    vals = {c: (1 if c == f else 0) for c in free}
    cur = {}
    for m, c in zip(b6, cs):
        v = sp.simplify(c.subs(sub).subs(vals)) if c in sub else vals[c]
        v = sp.nsimplify(v)
        if v != 0:
            cur[m] = F(int(sp.numer(v)), int(sp.denom(v)))
    if cur:
        vecs.append(sec6.reduce_mod_d(sec6.vec(cur)))
from ope_engine import rref
rows = [list(v) for v in vecs]
_, piv = rref(rows, len(sec6.basis))
check("joint kernel charge dimension = 2", len(piv) == 2)
check("KdV I5 is screened by +-sqrt2",
      screens(KDV_P6, sp.sqrt(2)) and screens(KDV_P6, -sp.sqrt(2)))
r1 = sec6.reduce_mod_d(sec6.vec(KDV_P6))
r2 = sec6.reduce_mod_d(sec6.vec(A2_P6))
rows2 = [list(r1), list(r2)]
_, piv2 = rref(rows2, len(sec6.basis))
check("I5_KdV and I5_A2 classes are linearly independent",
      len(piv2) == 2)

# 6 (finding 1): the four roots also screen I7' and I9' -- extends
# the finite evidence (review-verified, now retained).  Review round
# 2 additionally found the +-sqrt2 pattern is spin-dependent:
# it annihilates I7' but NOT I9'.
reps = {}
for sigma in (7, 9):
    kerx = genuine_kernel(A2_P6, gen_basis(sigma + 1, 1, True),
                          Sector(sigma + 6, 1, True),
                          Sector(sigma + 1, 1, True))
    check(f"multiplicity one at sigma={sigma}", len(kerx) == 1)
    reps[sigma] = kerx[0]
    ok = all(screens(reps[sigma], sgn * v)
             for v in (sp.sqrt(2) / 2, 2 * sp.sqrt(2))
             for sgn in (1, -1))
    check(f"I{sigma}' screened by all four momenta", ok)
check("+-sqrt2 annihilates I7' (both signs)",
      screens(reps[7], sp.sqrt(2)) and screens(reps[7], -sp.sqrt(2)))
check("+-sqrt2 does NOT annihilate I9' (either sign)",
      not screens(reps[9], sp.sqrt(2))
      and not screens(reps[9], -sp.sqrt(2)))
check("sqrt3 control: annihilates neither I7' nor I9'",
      not screens(reps[7], sp.sqrt(3))
      and not screens(reps[9], sp.sqrt(3)))

# 7 (review round 2): the weight-8 joint {+-sqrt2} kernel has charge
# dimension exactly 3, with the KdV I7 and enhanced I7' classes
# independent -- the self-dual degeneracy holds at spins 5 AND 7
# (and, by the check above, not at 9 for the enhanced ray)
def joint_dim(w):
    bw = gen_basis(w, 1, True)
    xs = [sp.Symbol(f'x{i}') for i in range(len(bw))]
    conds = []
    for sgn in (1, -1):
        conds += screening_conditions(dict(zip(bw, xs)),
                                      (sgn * sp.sqrt(2),), 1)
    M, _ = sp.linear_eq_to_matrix(conds, xs)
    sec = Sector(w, 1, True)
    return len(bw) - M.rank() - len(sec._img_pivots)

# the full certified table (review round 3): dims 1,2,3,4 at weights
# 4,6,8,10 -- the joint-kernel degeneracy GROWS through the computed
# range; what happens at spin 9 is that the enhanced I9' ray exits
# the joint kernel, not that the degeneracy ends
# complete table through weight 10 (review round 4: the full
# character includes a spin-1 class (oint T) and a spin-8 class at
# weight 9; the 1,2,3,4 progression is the odd-spin subsequence)
for w, want in ((2, 1), (3, 0), (4, 1), (5, 0), (6, 2), (7, 0),
                (8, 3), (9, 1), (10, 4)):
    check(f"joint {{+-sqrt2}} charge dimension at weight {w} = {want}",
          joint_dim(w) == want)
sec8 = Sector(8, 1, True)
KDV_I7 = genuine_kernel(KDV, gen_basis(8, 1, True),
                        Sector(11, 1, True), Sector(8, 1, True))[0]
check("KdV I7 screened by +-sqrt2",
      screens(KDV_I7, sp.sqrt(2)) and screens(KDV_I7, -sp.sqrt(2)))
ra = sec8.reduce_mod_d(sec8.vec(KDV_I7))
rb = sec8.reduce_mod_d(sec8.vec(reps[7]))
_, piv8 = rref([list(ra), list(rb)], len(sec8.basis))
check("KdV I7 and I7' classes are linearly independent",
      len(piv8) == 2)

print()
if failures:
    print(f"FAILURES: {failures}")
    print("C1 FINGERPRINT FAILED")
    sys.exit(1)
print("C1 ENHANCED-TOWER FINGERPRINT CERTIFIED: exactly four "
      "nonzero screening momenta {+-1/sqrt2, +-2 sqrt2} through the "
      "computed range I3..I9; the KdV pair additionally annihilates "
      "I5 and I7' but not I9'; complete joint-kernel table at charge "
      "spins 1..9: dims 1,0,1,0,2,0,3,1,4 (series "
      "q+q^3+2q^5+3q^7+q^8+4q^9, incl. the even spin-8 class; the "
      "enhanced ray exits at spin 9, the degeneracy does not); "
      "merger arithmetic asserted (candidate mechanism)")
