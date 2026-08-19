"""The Z2 mechanism behind the c=1 enhancement (campaign item 3).

Why the A2(2) joint kernel at the branch point falls inside the
KdV/sine-Gordon hierarchy at a = 1/sqrt2 (see kdv_identification.py):
at Q = 0 every tower charge is Z2-EVEN (phi -> -phi), and for an
even charge the first-order conservation condition under e^{-alpha
phi} is the parity image of the condition under e^{+alpha phi}, so
the two conditions have the SAME kernel on the even sector -- one
twisted A2(2) screening condition symmetrizes into the sine-Gordon
pair.

Certified, fail-closed:

 1. MECHANISM (even sector): for alpha in {1/sqrt2, sqrt2, 2} and
    density weights 2..9, the even-sector kernels of the +alpha
    condition, the -alpha condition, and the {+-alpha} pair all
    have EQUAL dimension -- which forces the three kernels to
    coincide (the pair kernel is the intersection; equality of its
    dimension with each single kernel's is containment both ways).
 2. SYMMETRIZATION CONSEQUENCE: on the even sector the c=1 A2(2)
    pair {e^{sqrt2 phi}, e^{-phi/sqrt2}} has the same joint-kernel
    dimension as the full symmetrized quartet
    {e^{+-sqrt2 phi}, e^{+-phi/sqrt2}} at every weight 2..9 --
    the containment driver A2(2)|_{c=1} into KdV at c = -25/2.
 3. NEGATIVE CONTROL (parity is load-bearing): on the FULL space
    at weight 6 the single-condition kernels under e^{+-phi/sqrt2}
    have Q-dim 4 (K-dim 2) each while the pair has Q-dim 2 --
    off the even sector the +-alpha conditions genuinely differ,
    so the mechanism is the Z2 parity of the charges, not an
    identity of the conditions.

Q-dim convention: a K-line over Q(sqrt2) counts as Q-dim 2 (split
representation of joint_kernel.py).
"""
import sys

import sympy as sp

from joint_kernel import joint

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


SQ2 = sp.sqrt(2)

# ---- 1. the mechanism on the even sector -----------------------------
for label, al in (("1/sqrt2", SQ2 / 2), ("sqrt2", SQ2),
                  ("2", sp.Integer(2))):
    for w in range(2, 10):
        dp = joint(w, ((al,),), cand_z2=True)
        dm = joint(w, ((-al,),), cand_z2=True)
        dj = joint(w, ((al,), (-al,)), cand_z2=True)
        check(f"even sector, alpha={label}, w={w}: dim(+)=dim(-)="
              f"dim(pair) [{dp},{dm},{dj}]", dp == dm == dj)

# ---- 2. A2(2) pair == symmetrized quartet on the even sector ---------
A2PAIR = ((SQ2,), (-SQ2 / 2,))
QUARTET = ((SQ2,), (-SQ2,), (SQ2 / 2,), (-SQ2 / 2,))
for w in range(2, 10):
    da = joint(w, A2PAIR, cand_z2=True)
    dq = joint(w, QUARTET, cand_z2=True)
    check(f"even sector, w={w}: A2(2) pair kernel == symmetrized "
          f"quartet kernel [{da},{dq}]", da == dq)

# ---- 3. negative control: parity is load-bearing ---------------------
dp6 = joint(6, ((SQ2 / 2,),))
dm6 = joint(6, ((-SQ2 / 2,),))
dj6 = joint(6, ((SQ2 / 2,), (-SQ2 / 2,)))
check("negative control: FULL space w=6 singles have Q-dim 4, pair 2 "
      f"[{dp6},{dm6},{dj6}]",
      dp6 == 4 and dm6 == 4 and dj6 == 2)

print()
if failures:
    print(f"Z2 MECHANISM CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("Z2 MECHANISM CERTIFIED (even-sector +-alpha kernel coincidence; "
      "A2(2) pair symmetrizes to the quartet at c=1)")
