"""Class-by-class: the enhanced tower equals ker ad_{I'3} = KdV classes.

Campaign item 2 (+ the finite-range half of item 4) for the
c = -25/2 identification (see kdv_identification.py):

 1. I'3 (recomputed as the unique sigma=3 member of ker ad_{I5})
    has ker ad_{I'3} of charge multiplicity EXACTLY 1 at odd spins
    5, 7, 9 and 0 at even spins 4, 6, 8.
 2. At sigma = 5, 7, 9 the ker ad_{I'3} class EQUALS the tower
    class (ker ad_{I5} representative) mod total derivatives.
    Since kdv_identification.py certifies I'3 = KdV I3 at
    Q^2 = 9/8, the tower is the KdV commutant class by class in
    the certified range.
 3. Sine-Gordon-pair kernel (infinitude anchor, campaign item 4):
    on the Z2-EVEN sector (where the tower lives) the joint kernel
    of the first-order conditions under e^{+phi/sqrt2} AND
    e^{-phi/sqrt2} has K-dimension over Q(sqrt2) EXACTLY 1 at odd
    spins 1, 3, 5, 7, 9 and 0 at even spins (a K-line counts as
    Q-dim 2 in the split representation, so the certified Q-dims
    are 2/0).  Membership of the tower representatives in this
    kernel is certified in c1_enhanced_fingerprint.py; dimension
    one therefore forces the pair-kernel class to BE the tower
    class.  INFORMATIONAL: on the FULL space the pair kernel
    additionally contains one odd-parity class at sigma = 6
    (observed Q-dims 2,0,2,0,2,2,2,0,2 at sigma 1..9) -- outside
    the even sector, killed by the full quartet
    (c1_four_condition_char.py), recorded here as an observed
    phenomenon and flagged for study.
 4. Genericity control: the even-sector pair {e^{+-2phi}} (a = 2,
    c = -2, a generic sine-Gordon point) also has K-dim 1 at
    sigma = 3 -- the dimension the semicontinuity argument
    transports along the family (Feigin--Frenkel generic a + upper
    semicontinuity => at least one charge per odd spin at EVERY a,
    in particular a = 1/sqrt2; multiplicity one in range is this
    certificate).

Fail-closed: exits nonzero unless every check passes.
"""
import sys
from fractions import Fraction as F

import sympy as sp

from families import L
from ope_engine import Sector, d_mono, gen_basis
from screening_fit import residue_screen, twisted_d
from solve_charges import genuine_kernel
from sparse_linalg import SparsePivots, lcm

from joint_kernel import joint

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


SQ2 = sp.sqrt(2)

# the A2(2) I5 density at N=1 (as in c1_enhanced_fingerprint.py)
A2_P6 = acc([(L((3, 0), (3, 0)), F(124, 96)),
             (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
             (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-68, 8)),
             (L((1, 0),) * 6, F(1))])

ker3 = genuine_kernel(A2_P6, gen_basis(4, 1, True),
                      Sector(9, 1, True), Sector(4, 1, True))
check("I'3 unique at sigma=3 in ker ad_{I5}", len(ker3) == 1)
I3p = ker3[0]

# ---- 1+2: ker ad_{I'3} spins 4..9, with class equality at 5,7,9 ------
for sigma in (4, 5, 6, 7, 8, 9):
    w = sigma + 1
    kerK = genuine_kernel(I3p, gen_basis(w, 1, True),
                          Sector(w + 3, 1, True), Sector(w, 1, True))
    want = 1 if sigma % 2 == 1 else 0
    check(f"ker ad_{{I'3}} multiplicity at sigma={sigma} is {want}",
          len(kerK) == want)
    if sigma % 2 == 0 or len(kerK) != 1:
        continue
    kerA = genuine_kernel(A2_P6, gen_basis(w, 1, True),
                          Sector(w + 5, 1, True), Sector(w, 1, True))
    check(f"tower multiplicity at sigma={sigma} is 1", len(kerA) == 1)
    if len(kerA) != 1:
        continue
    sec = Sector(w, 1, True)
    va = sec.reduce_mod_d(sec.vec(kerA[0]))
    vk = sec.reduce_mod_d(sec.vec(kerK[0]))
    ok = any(va) and any(vk)
    if ok:
        # proportionality of the two reduced class vectors
        ia = next(i for i, x in enumerate(va) if x)
        ok = (vk[ia] != 0 and
              all(va[i] * vk[ia] == vk[i] * va[ia]
                  for i in range(len(va))))
    check(f"sigma={sigma}: ker ad_{{I'3}} class == tower class mod d",
          ok)

# ---- 3: sine-Gordon pair kernel at a = 1/sqrt2 -----------------------
PAIR = ((SQ2 / 2,), (-SQ2 / 2,))
for w in range(2, 11):
    sigma = w - 1
    want = 2 if sigma % 2 == 1 else 0
    qd = joint(w, PAIR, cand_z2=True)
    check(f"even-sector pair {{e^{{+-phi/sqrt2}}}} kernel Q-dim at "
          f"sigma={sigma} is {want} (K-dim {want // 2})", qd == want)

# informational full-space table, observed values baked
FULL_OBS = [2, 0, 2, 0, 2, 2, 2, 0, 2]
for w in range(2, 11):
    sigma = w - 1
    check(f"full-space pair kernel Q-dim at sigma={sigma} is "
          f"{FULL_OBS[w - 2]} (observed; odd-parity stray at sigma=6)",
          joint(w, PAIR) == FULL_OBS[w - 2])

# ---- 4: genericity control at a = 2 ----------------------------------
check("generic control: even-sector pair {e^{+-2phi}} K-dim at "
      "sigma=3 is 1",
      joint(4, ((sp.Integer(2),), (sp.Integer(-2),)), cand_z2=True) == 2)

print()
if failures:
    print(f"TOWER-CLASS CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("KDV TOWER CLASS-BY-CLASS CERTIFIED "
      "(ker ad_{I'3} = tower classes; sine-Gordon pair kernel "
      "multiplicity one, odd spins 1-9)")
