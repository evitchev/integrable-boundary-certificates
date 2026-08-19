"""The odd-parity sector of the centralizer: the spin-6 charge J.

Round-2 review finding (arm A, 2026-08-19), confirmed in-engine:
the FULL-space centralizer ker ad_{I5} is strictly larger than the
Z2-even tower -- it contains an odd-parity spin-6 charge

    J = oint[ (64/7)(dphi)^7 - 40 (dphi)^3 (d2phi)^2
              + dphi (d3phi)^2 ],

identified by the review as proportional to the singlet-algebra
commutator [Q_{2sqrt2}, Q_{-2sqrt2}] (the h_{3,1} = 7 generator at
c_{1,4}), whence by the Jacobi identity it commutes with every
charge annihilated by both +-2sqrt2 conditions -- the entire
quartet kernel, at all spins.  This certificate pins the certified
part of that picture (the residue identification itself is the
review's, consistent with everything below):

 1. J is a nonzero charge class (nonzero mod total derivatives).
 2. [J, I3'] = [J, I5] = [J, I7'] = [J, I9'] = 0 mod d, with the
    tower representatives recomputed from scratch.
 3. Conservation pattern: J is conserved to first order under
    e^{+-phi/sqrt2} and NOT under e^{+-2sqrt2 phi} -- the pattern
    forced by the [Q_+, Q_-] identification.
 4. Negative control: perturbing one coefficient of J breaks
    commutativity with I3'.
 5. The c_{1,p} kernel-jump fact (sharpens the multiplicity-one
    caution): at a = 1 (c = c_{1,2} = -2) the even-sector
    sine-Gordon pair kernel has K-dims 1,1,1,2 at spins 1,3,5,7
    (Q-dims 2,2,2,4) -- enlargements at degenerate couplings are
    real, so all-spin multiplicity-one at a = 1/sqrt2 legitimately
    remains open.

Fail-closed: exits nonzero unless every check passes.
"""
import sys
from fractions import Fraction as F

import sympy as sp

from families import L
from joint_kernel import joint
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
    for m, c in pairs:
        out[m] = out.get(m, 0) + c
    return {m: c for m, c in out.items() if c}


SQ2 = sp.sqrt(2)

J = acc([(L(*((1, 0),) * 7), F(64, 7)),
         (L((1, 0), (1, 0), (1, 0), (2, 0), (2, 0)), F(-40)),
         (L((1, 0), (3, 0), (3, 0)), F(1))])

A2_P6 = acc([(L((3, 0), (3, 0)), F(124, 96)),
             (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
             (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-68, 8)),
             (L((1, 0),) * 6, F(1))])

ker3 = genuine_kernel(A2_P6, gen_basis(4, 1, True),
                      Sector(9, 1, True), Sector(4, 1, True))
check("I'3 unique at sigma=3", len(ker3) == 1)
I3p = ker3[0]

# tower representatives at sigma = 7, 9
reps = {}
for sigma in (7, 9):
    w = sigma + 1
    kerx = genuine_kernel(A2_P6, gen_basis(w, 1, True),
                          Sector(w + 5, 1, True), Sector(w, 1, True))
    check(f"tower multiplicity one at sigma={sigma}", len(kerx) == 1)
    reps[sigma] = kerx[0]

# ---- 1. J is a nonzero class ------------------------------------------
sec7 = Sector(7, 1, False)
check("J nonzero mod d", any(sec7.reduce_mod_d(sec7.vec(J))))


def commutes(P, wP, Jc, wJ):
    r = residue_cur(P, Jc)
    sec = Sector(wP + wJ - 1, 1, False)
    return not any(sec.reduce_mod_d(sec.vec(r)))


# ---- 2. J commutes with the tower charges ----------------------------
check("[I'3, J] = 0 mod d", commutes(I3p, 4, J, 7))
check("[I5, J] = 0 mod d", commutes(A2_P6, 6, J, 7))
check("[I'7, J] = 0 mod d", commutes(reps[7], 8, J, 7))
check("[I'9, J] = 0 mod d", commutes(reps[9], 10, J, 7))


# ---- 3. conservation pattern ------------------------------------------
def screens(P, v):
    return not [x for x in screening_conditions(P, (v,), 1)
                if sp.simplify(x) != 0]


for v, want, lbl in ((SQ2 / 2, True, "+1/sqrt2"),
                     (-SQ2 / 2, True, "-1/sqrt2"),
                     (2 * SQ2, False, "+2sqrt2"),
                     (-2 * SQ2, False, "-2sqrt2")):
    check(f"J conserved under e^({lbl} phi) is {want}",
          screens(J, v) == want)

# ---- 4. negative control ----------------------------------------------
Jbad = dict(J)
Jbad[L(*((1, 0),) * 7)] = F(65, 7)
check("negative control: perturbed J fails [I'3, J'] = 0",
      not commutes(I3p, 4, Jbad, 7))

# ---- 5. the c_{1,p} kernel-jump fact ----------------------------------
one = sp.Integer(1)
dims = [joint(w, ((one,), (-one,)), cand_z2=True) for w in (2, 4, 6, 8)]
check("a=1 (c_{1,2}): even-sector pair kernel Q-dims 2,2,2,4 at "
      "spins 1,3,5,7 (K-dim 2 at spin 7)", dims == [2, 2, 2, 4])

print()
if failures:
    print(f"ODD-SECTOR CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("ODD-SECTOR J CERTIFIED (nonzero class commuting with the tower; "
      "conserved under +-1/sqrt2 only; c_{1,2} kernel jump real)")
