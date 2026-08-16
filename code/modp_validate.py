"""Validation of the GF(p) vacancy certifier (genuine_kernel_modp)
against known exact verdicts at lower spins, across all three engines.

Expected values come from exact computations already in the repo:
the cylindrical sigma<=11 diagnosis, the A2(2) generic scan (spin-5
charge exists, sigma=8,9 vacant), and the fresh exact paperclip
sigma=12 vacancy.  For vacant sectors the certifier must return
tags == drank (a rigorous certificate); for dim-1 sectors it must
return the upper bound tags - drank == 1 with vacant=False."""
import sys
from fractions import Fraction as F

import collapsed
import collapsed_pairs
import invariant_engine
import mixed_engine
invariant_engine.residue_inv_mono = collapsed_pairs.residue_inv_mono_smart
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from families import paperclip_P4
from invariant_engine import (a2_P6, curve_point, d_inv_mono,
                              gen_inv_basis, residue_inv_cur, spin_inv)
from mixed_engine import (cyl_P4, d_mixed_mono, gen_mixed_basis,
                          residue_mixed_cur, spin_mixed)
from ope_engine import d_mono as pc_d
from ope_engine import gen_basis as pc_gen
from ope_engine import residue_cur as pc_res
from ope_engine import spin as pc_spin
from sparse_linalg import genuine_kernel_modp

N, S = curve_point(F(2))


def cyl(sol, sigma):
    return genuine_kernel_modp(cyl_P4(sol, N, S),
                               gen_mixed_basis(sigma + 1), N,
                               residue_mixed_cur, gen_mixed_basis,
                               d_mixed_mono, spin_mixed)


def a2(sigma):
    return genuine_kernel_modp(a2_P6(N, S), gen_inv_basis(sigma + 1),
                               N, residue_inv_cur, gen_inv_basis,
                               d_inv_mono, spin_inv)


def pc(sigma, nval):
    return genuine_kernel_modp(paperclip_P4(nval),
                               pc_gen(sigma + 1, 2, True), None,
                               lambda P, Q, _: pc_res(P, Q),
                               lambda s: pc_gen(s, 2, True),
                               pc_d, pc_spin)


CASES = [
    ("cyl1 s4", 0, lambda: cyl(1, 4)),
    ("cyl1 s5", 1, lambda: cyl(1, 5)),
    ("cyl2 s5", 1, lambda: cyl(2, 5)),
    ("cyl3 s7", 1, lambda: cyl(3, 7)),
    ("cyl1 s6", 0, lambda: cyl(1, 6)),
    ("a2 s5", 1, lambda: a2(5)),
    ("a2 s7", 1, lambda: a2(7)),
    ("a2 s8", 0, lambda: a2(8)),
    ("a2 s9", 0, lambda: a2(9)),
    ("pc s8 n=7/3", 0, lambda: pc(8, F(7, 3))),
    ("pc s12 n=7/3", 0, lambda: pc(12, F(7, 3))),
]

failures = []
for name, expected, fn in CASES:
    tags, drank, vacant = fn()
    bound = tags - drank
    ok = (vacant == (expected == 0)) and bound == expected
    print(f"{'PASS' if ok else 'FAIL'}  {name}: tags {tags}, d-rank "
          f"{drank}, bound {bound} (expected dim {expected}, "
          f"certified vacant: {vacant})", flush=True)
    if not ok:
        failures.append(name)

if failures:
    print(f"{len(failures)} FAILURES: {failures}")
    sys.exit(1)
print("ALL MOD-P VALIDATIONS PASSED")
