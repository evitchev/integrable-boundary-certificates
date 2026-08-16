"""Validation: sparse kernel pipeline == dense engines on knowns."""
import time
from fractions import Fraction as F
import invariant_engine as ie
import mixed_engine as me
import collapsed, collapsed_pairs
ie.residue_inv_mono = collapsed_pairs.residue_inv_mono_smart
me.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from invariant_engine import (gen_inv_basis, genuine_kernel_inv, a2_P6,
                              kdv_P4, curve_point, d_inv_mono, spin_inv,
                              residue_inv_cur, InvSector)
from mixed_engine import (gen_mixed_basis, genuine_kernel_mixed, cyl_P4,
                          d_mixed_mono, spin_mixed, residue_mixed_cur,
                          MixedSector)
from sparse_linalg import genuine_kernel_sparse

Nv, sv = curve_point(F(2))

def canon_class(cur, Sec, s):
    sec = Sec(s)
    v = sec.reduce_mod_d(sec.vec(cur))
    nz = next((x for x in v if x), None)
    return None if nz is None else tuple(x / nz for x in v)

checks = [
    ("a2 s3 (dim0)", a2_P6(Nv, sv), 3, "inv"),
    ("a2 s5 (dim1)", a2_P6(Nv, sv), 5, "inv"),
    ("a2 s7 (dim1)", a2_P6(Nv, sv), 7, "inv"),
    ("kdv s6 (dim0)", kdv_P4(), 6, "inv"),
    ("cyl1 s4 (dim0)", cyl_P4(1, Nv, sv), 4, "mix"),
    ("cyl1 s5 (dim1)", cyl_P4(1, Nv, sv), 5, "mix"),
    ("cyl2 s5 (dim1)", cyl_P4(2, Nv, sv), 5, "mix"),
]
for label, P, sig, fam in checks:
    if fam == "inv":
        dense = genuine_kernel_inv(P, gen_inv_basis(sig + 1), Nv)
        sparse = genuine_kernel_sparse(P, gen_inv_basis(sig + 1), Nv,
                                       residue_inv_cur, gen_inv_basis,
                                       d_inv_mono, spin_inv)
        Sec = InvSector
    else:
        dense = genuine_kernel_mixed(P, gen_mixed_basis(sig + 1), Nv)
        sparse = genuine_kernel_sparse(P, gen_mixed_basis(sig + 1), Nv,
                                       residue_mixed_cur, gen_mixed_basis,
                                       d_mixed_mono, spin_mixed)
        Sec = MixedSector
    ok = len(dense) == len(sparse)
    if ok and dense:
        ok = all(canon_class(a, Sec, sig + 1) == canon_class(b, Sec, sig + 1)
                 for a, b in zip(dense, sparse))
    print(f"{label}: {'PASS' if ok else 'FAIL'} "
          f"(dense {len(dense)}, sparse {len(sparse)})")
    assert ok
print("ALL SPARSE VALIDATIONS PASSED")
