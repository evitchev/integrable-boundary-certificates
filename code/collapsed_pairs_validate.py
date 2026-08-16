"""Validation of the Phase-2 pair-collapse (collapsed_pairs.py)."""
import random, time
from fractions import Fraction as F
from invariant_engine import residue_inv_mono, canon, canon_pair
import collapsed_pairs as cp

random.seed(20260803)
bad = 0
t0 = time.time()
ntest = 0
for _ in range(250):
    P = canon(tuple(canon_pair(random.randint(1, 3), random.randint(1, 3))
                    for _ in range(random.randint(1, 3))))
    # Q with deliberate repetitions to exercise the collapse
    base = [canon_pair(random.randint(1, 2), random.randint(1, 3))
            for _ in range(random.randint(1, 2))]
    Q = canon(tuple(random.choice(base)
                    for _ in range(random.randint(1, 5))))
    a = residue_inv_mono(P, Q)
    b = cp.residue_inv_mono_fast(P, Q)
    if a != b:
        bad += 1
        print("MISMATCH", P, Q)
        if bad > 2:
            break
    ntest += 1
print(f"random equality: {'PASS' if not bad else 'FAIL'} "
      f"({ntest} pairs, {time.time()-t0:.0f}s)")
assert bad == 0

# medium all-(1,1) case: both engines feasible, collapse regime
P = canon(((1, 1), (1, 1), (1, 1)))
Q = canon(((1, 1),) * 5)
t0 = time.time(); a = residue_inv_mono(P, Q); t1 = time.time()
b = cp.residue_inv_mono_fast(P, Q); t2 = time.time()
print(f"medium (11)^3 x (11)^5: {'PASS' if a == b else 'FAIL'} "
      f"slow={t1-t0:.1f}s fast={t2-t1:.2f}s")
assert a == b

# physics regression: KdV and A2 through the smart router
import invariant_engine
from invariant_engine import (InvSector, gen_inv_basis, genuine_kernel_inv,
                              kdv_P4, kdv_P6, a2_P6, curve_point,
                              residue_inv_cur)
N = F(3)
r_slow = residue_inv_cur(kdv_P4(), kdv_P6(3), N)
invariant_engine.residue_inv_mono = cp.residue_inv_mono_smart
r_fast = residue_inv_cur(kdv_P4(), kdv_P6(3), N)
assert r_slow == r_fast
Nv, sv = curve_point(F(2))
g = genuine_kernel_inv(a2_P6(Nv, sv), gen_inv_basis(8), Nv)
assert len(g) == 1
invariant_engine.residue_inv_mono = cp._ORIG
print("physics regression (KdV residue, A2 sigma-7 kernel): PASS")
print("ALL PHASE-2 VALIDATIONS PASSED")
