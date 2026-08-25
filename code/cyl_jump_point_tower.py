"""The jump point (N, s) = (-25/2, -45/2) of the cylindrical curve, up the tower.

cyl_hidden_kdv_scan.py found, at this point, an I_5 kernel of dimension 2
(a multiplicity jump) with a one-dimensional e^{sqrt2 X}-screened
subspace.  This certificate follows the point up the tower (weights 6,
8, 10, 12 = spins 5, 7, 9, 11), computing the FULL ad_{I_3} kernel at
each weight (the promoted generic densities cannot be trusted to span
an enlarged kernel) and the screened subspace inside it, exactly over
Q(sqrt2) (doubled-field fraction-free sparse elimination).

Findings, all asserted:
 1. MERGER ON ONE SHEET: Solutions 2 and 3 have the SAME I_3 density at
    this point -- (dX)^4 - 2(d2X)^2 - 3(dX)^2(dY.dY) + (dY.dY)^2 - 2(d2Y.d2Y),
    i.e. 4[:T_X^2: - 3 T_X T_Y + :T_Y^2:] -- unlike the N=28 merger,
    which pairs deck-conjugate points.
 2. Kernel dimensions [2, 1, 1, 1] at spins [5, 7, 9, 11]: the
    multiplicity jump is confined to spin 5 in the known range.
 3. Screened-subspace dimensions [1, 1, 1, 1] at a^2 = 2 and [0, 0, 0, 0]
    at a^2 = 3: the self-dual screened line runs through the whole
    known tower; at spins 7-11 the unique class IS screened.  So this is
    a tower-level L(1,0)_X-type embedding with exactly one extra
    (unscreened) charge at spin 5.
 4. The extra spin-5 class is genuinely MIXED: no combination of the
    spin-5 kernel is pure-Y (this is not a decoupling point).
 5. The promoted generic I_7, I_9, I_11 of BOTH Solutions 2 and 3
    specialize pole-free to the unique kernel class at this point: the
    two generic towers meet here through I_11, and the generic line is
    the screened line.
Finite-range statements throughout (spins <= 11).  Fail-closed.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

from invariant_engine import _rref
from mixed_engine import (MixedSector, _M, canon_mixed, cyl_P4, d_mixed_mono,
                          gen_mixed_basis, genuine_kernel_mixed)
from sparse_linalg import SparsePivots, lcm

# reuse the Q(sqrt r) residue machinery and the promotion loader
_n28 = Path(__file__).resolve().parent / "cyl_n28_tower_embedding.py"
_src = _n28.read_text()
_ns = {"__file__": str(_n28)}
exec(_src[:_src.index("POINTS = [(1, F(28), F(-9)), (2, F(28), F(9))]")], _ns)
residue_x_quadratic = _ns["residue_x_quadratic"]
gen_mixed_full = _ns["gen_mixed_full"]
promoted_density = _ns["promoted_density"]

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def screened_dim(currents, w, r):
    """dim of the subspace of span(currents) conserved under e^{sqrt(r) X}:
    one twisted-derivative pivot structure (doubled field), exact
    reductions of every residue, nullity of the stacked reductions."""
    s_ = w - 1
    tgt = []
    for t in range(0, s_ + 1):
        tgt += gen_mixed_full(t)
    lower = gen_mixed_full(s_ - 1)
    index = {m: i for i, m in enumerate(tgt)}
    n = len(tgt)
    piv = SparsePivots()
    for u in lower:
        du = d_mixed_mono(u)
        xu = canon_mixed(u[0] + (1,), u[1])
        col = {}
        for m, c in du.items():
            col[index[m]] = col.get(index[m], 0) + int(c)
        col[n + index[xu]] = col.get(n + index[xu], 0) + 1
        piv.insert({k: c for k, c in col.items() if c})
        col2 = {index[xu]: r}
        for m, c in du.items():
            col2[n + index[m]] = col2.get(n + index[m], 0) + int(c)
        piv.insert({k: c for k, c in col2.items() if c})
    rows = []
    for P in currents:
        res = residue_x_quadratic(P, r)
        den = 1
        for u, v in res.values():
            den = lcm(den, u.denominator)
            den = lcm(den, v.denominator)
        vec = {}
        for m, (u, v) in res.items():
            if u:
                vec[index[m]] = int(u * den)
            if v:
                vec[n + index[m]] = int(v * den)
        rem, c = piv.reduce_scaled(vec)
        rows.append({k: F(x) / c / den for k, x in rem.items()})
    keys = sorted({k for row in rows for k in row})
    if not keys:
        return len(currents)
    M = [[row.get(k, F(0)) for k in keys] for row in rows]
    _, pivots = _rref([list(m) for m in M], len(keys))
    return len(currents) - len(pivots)


def mixed_projection_ranks(ker, w):
    """(rank of the X-sector projection, rank of the Y-sector projection) of
    the reduced kernel vectors.  Full rank of the X-projection means NO
    combination is pure-Y; full rank of the Y-projection means NO
    combination is pure-X -- basis-independent (Codex audit 2026-08-23)."""
    sec = MixedSector(w)
    vs = [sec.reduce_mod_d(sec.vec(K)) for K in ker]
    xidx = [i for i, m in enumerate(sec.basis) if len(m[0]) > 0]
    yidx = [i for i, m in enumerate(sec.basis) if len(m[1]) > 0]
    rx = len(_rref([[v[i] for i in xidx] for v in vs], len(xidx))[1])
    ry = len(_rref([[v[i] for i in yidx] for v in vs], len(yidx))[1])
    return rx, ry


N, s = F(-25, 2), F(-45, 2)
P2, P3 = cyl_P4(2, N, s), cyl_P4(3, N, s)
check("J1: Solutions 2 and 3 have the SAME I_3 at (-25/2,-45/2) (same-sheet "
      "merger)", P2 == P3)
check("J1: the common I_3 is (dX)^4 - 2(d2X)^2 - 3(dX)^2(dY.dY) + (dY.dY)^2 "
      "- 2(d2Y.d2Y)",
      P2 == {_M((1, 1, 1, 1)): F(1), _M((2, 2)): F(-2),
             _M((1, 1), [(1, 1)]): F(-3), _M((), [(1, 1), (1, 1)]): F(1),
             _M((), [(2, 2)]): F(-2)})

EXPECT = {6: (2, 1, 0), 8: (1, 1, 0), 10: (1, 1, 0), 12: (1, 1, 0)}
kernels = {}
for w, (kdim, s2, s3) in EXPECT.items():
    ker = genuine_kernel_mixed(P2, gen_mixed_basis(w), N)
    kernels[w] = ker
    d2 = screened_dim(ker, w, 2) if ker else 0
    d3 = screened_dim(ker, w, 3) if ker else 0
    check(f"J2: spin {w-1}: ad_{{I_3}} kernel dim {len(ker)}, e^{{sqrt2 X}}-screened "
          f"dim {d2}, e^{{sqrt3 X}}-screened dim {d3}  (expected {kdim},{s2},{s3})",
          (len(ker), d2, d3) == (kdim, s2, s3))

# the extra spin-5 class is genuinely mixed: NO combination of the kernel
# is pure-Y or pure-X (projection ranks, basis-independent)
rx, ry = mixed_projection_ranks(kernels[6], 6)
check(f"J3: spin-5 kernel: X-projection rank {rx}, Y-projection rank {ry} of 2 "
      "-- no pure-Y and no pure-X combination (not a decoupling point)",
      (rx, ry) == (2, 2))

# promoted generic towers of both solutions meet the unique kernel class
for sigma in (7, 9, 11):
    w = sigma + 1
    sec = MixedSector(w)
    vk = sec.reduce_mod_d(sec.vec(kernels[w][0]))
    jk = next(i for i, x in enumerate(vk) if x)
    vk = [x / vk[jk] for x in vk]
    for sol in (2, 3):
        try:
            D = promoted_density(sol, sigma, N, s)
        except (ZeroDivisionError, ValueError) as exc:
            check(f"J4: Solution {sol} promoted I_{sigma} pole-free at the point "
                  f"({exc})", False)
            continue
        v = sec.reduce_mod_d(sec.vec(D))
        ok = any(v)
        if ok:
            j = next(i for i, x in enumerate(v) if x)
            ok = [x / v[j] for x in v] == vk
        check(f"J4: Solution {sol} promoted generic I_{sigma} specializes to the "
              f"unique (screened) kernel class at (-25/2,-45/2)", ok)

print()
if failures:
    print(f"JUMP-POINT TOWER CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("JUMP-POINT TOWER CERTIFIED ((-25/2,-45/2): Solutions 2 and 3 merge on one "
      "sheet; kernel dims 2,1,1,1 at spins 5,7,9,11; the self-dual screened "
      "line runs through I_11 with one extra mixed charge at spin 5; "
      "finite-range)")
