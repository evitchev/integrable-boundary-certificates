"""Validation of collapsed.py.

The router sends only multiplicity-heavy inputs to the collapsed path;
so we validate: (a) exact equality original-vs-collapsed on
all-ones-dominated MEDIUM monomials (collapsed regime, original still
feasible); (b) a physics regression through residue_mixed_cur with the
router monkeypatched in.  (Mixed-order small pairs, where the router
uses the original engine verbatim, were additionally spot-validated on
29 random pairs during development: all equal.)
"""
import time
from fractions import Fraction as F
from mixed_engine import (residue_mixed_mono, canon_mixed, gen_mixed_basis,
                          cyl_P4, genuine_kernel_mixed, residue_mixed_cur)
from invariant_engine import curve_point
import collapsed, mixed_engine

cases = [
    (canon_mixed((1,)*6, ((1,1),)), canon_mixed((1,)*8, ())),
    (canon_mixed((1,)*8, ()), canon_mixed((1,)*8, ((1,1),))),
    (canon_mixed((1,1,1,1,2,2), ()), canon_mixed((1,)*8, ())),
    (canon_mixed((1,)*6, ((1,2),)), canon_mixed((1,)*8, ((1,1),))),
]
for P, Q in cases:
    t0 = time.time(); a = residue_mixed_mono(P, Q); t1 = time.time()
    b = collapsed.residue_mixed_mono_fast(P, Q); t2 = time.time()
    ok = a == b
    print(f"medium all-ones: {'PASS' if ok else 'FAIL'} "
          f"slow={t1-t0:.1f}s fast={t2-t1:.2f}s", flush=True)
    assert ok

N, s = curve_point(F(2))
P4 = cyl_P4(1, N, s)
g = genuine_kernel_mixed(P4, gen_mixed_basis(6), N)[0]
r_slow = residue_mixed_cur(P4, g, N)
orig = mixed_engine.residue_mixed_mono
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
r_fast = residue_mixed_cur(P4, g, N)
mixed_engine.residue_mixed_mono = orig
print("physics regression:", "PASS" if r_slow == r_fast else "FAIL", flush=True)
assert r_slow == r_fast
print("ALL COLLAPSED VALIDATIONS PASSED", flush=True)

# retained randomized equality sweep (guarded to sizes where both
# engines are fast; seed fixed so the sweep is reproducible)
import random as _r
_r.seed(20260802)
_n = 0
for _ in range(40):
    xs1 = tuple(sorted(_r.choice([1, 1, 2]) for _ in range(_r.choice([2, 4, 6]))))
    xs2 = tuple(sorted(_r.choice([1, 1, 2]) for _ in range(_r.choice([2, 4]))))
    yp1 = tuple(sorted(tuple(sorted((1, _r.choice([1, 2])))) for _ in range(_r.randint(0, 1))))
    P, Q = canon_mixed(xs1, yp1), canon_mixed(xs2, ())
    if collapsed._raw_count(len(P[0]), len(Q[0])) > 100_000:
        continue
    assert residue_mixed_mono(P, Q) == collapsed.residue_mixed_mono_fast(P, Q), (P, Q)
    _n += 1
print(f"randomized retained sweep: PASS ({_n} pairs)")
