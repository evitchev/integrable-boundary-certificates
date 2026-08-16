"""Spectral data for ad_{I_3}: kernel dimensions per spin.

For the bootstrap lemma we need, in each symmetry sector:
  (V) even charge spins sigma: no genuine charge commutes with I_3;
  (E) odd charge spins sigma: exactly one (multiplicity one).
Verified here for the N=1 KdV branch and for the paperclip branch at
several generic rational n (full rank at one point implies the generic
statement; kernel dims can only jump UP at special points).

Also: direct check [I_5, I_7] = 0 for the paperclip at many sample n.
With denominators cleared, every entry of the reduced commutator vector
is a polynomial in n of degree <= 5 + 9 = 14, so vanishing at >= 15
distinct sample points proves the identity for ALL n.
"""

from fractions import Fraction as F

from ope_engine import Sector, gen_basis, residue_cur, cur_scale
from families import KDV_P4, paperclip_P4, paperclip_P6
from solve_charges import genuine_kernel


def kernel_dim(P_low, sigma, nfields, z2):
    """# genuine charges of charge-spin sigma commuting with oint P_low."""
    cand = gen_basis(sigma + 1, nfields, z2)
    if not cand:
        return 0, 0
    cand_sec = Sector(sigma + 1, nfields, z2)
    target = Sector(sigma + 4, nfields, z2)
    gen = genuine_kernel(P_low, cand, target, cand_sec)
    return len(gen), cand_sec.q_dim()


print("== KdV-Virasoro, N = 1 (Z2 sector) ==")
print("sigma | dim Q^sigma | dim ker ad_I3")
for sigma in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
    k, q = kernel_dim(KDV_P4, sigma, 1, True)
    print(f"  {sigma:2d}  |     {q:3d}     |      {k}")

print("\n== paperclip (Z2 x Z2 sector) ==")
NS = [F(2), F(3), F(7, 3), F(-5, 2), F(13, 5)]
print("sigma | dim Q^sigma | dim ker ad_I3 at n in "
      + ", ".join(str(x) for x in NS))
for sigma in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
    dims = []
    for nval in NS:
        k, q = kernel_dim(paperclip_P4(nval), sigma, 2, True)
        dims.append(k)
    print(f"  {sigma:2d}  |     {q:3d}     |      {dims}")

# ------------------------------------------------ direct [I_5, I_7] check
print("\n== direct [I_5, I_7] = 0, paperclip ==")
import json
from compute_p8 import parse_mono
import sympy as sp

with open("../results/paperclip_P8.json") as f:
    DATA = json.load(f)
nsym = sp.Symbol('n')
P8_EXPR = {parse_mono(s): sp.cancel(sp.sympify(e, locals={'n': nsym}))
           for s, e in zip(DATA["basis"], DATA["coefficients"])}


def paperclip_P8(nval):
    out = {}
    for mono, e in P8_EXPR.items():
        v = e.subs(nsym, sp.Rational(nval))
        v = sp.Rational(v)
        if v != 0:
            out[mono] = F(int(sp.numer(v)), int(sp.denom(v)))
    return out


sec13 = Sector(13, 2, z2=True)
count = 0
n_samples = [F(k) for k in range(1, 13) if F(k) != 0] + \
            [F(-1), F(-3), F(-4), F(-5), F(1, 2), F(3, 2), F(7, 3), F(13, 5)]
for nval in n_samples:
    res = residue_cur(paperclip_P6(nval), paperclip_P8(nval))
    ok = sec13.is_total_derivative(res)
    if not ok:
        print(f"  n = {nval}: FAIL")
        raise SystemExit(1)
    count += 1
print(f"  [I_5, I_7] = 0 at {count} distinct n values "
      f"(> 14 = degree bound after clearing denominators),")
print("  hence [I_5, I_7] = 0 identically in n.")
