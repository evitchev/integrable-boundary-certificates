"""The t_2 shadow test -- outcome form (preregistration + refutation
addendum: results/cyl_t2_shadow_preregistration.md).

Point t_2 = 7/5, (N, s) = (-24, -35) = c_{1,6}, Solution 3.  The
PREREGISTERED expectation (a screened c_{1,6} line through the tower,
extra spin-9 class as its unscreened shadow) was REFUTED by the
development smoke -- foreseeably: the certified record
(results/cyl_strata_lab.md) already classes this jump as having NO
screening structure.  This certificate pins the observed result:

  T1  ad_{I_3} kernel dims [1, 1, 2, 1] at weights 6, 8, 10, 12
      (spins 5, 7, 9, 11) -- reconfirming the record.
  T2  the no-screening-structure statement EXTENDS to the c_{1,6} Kac
      pair: screened dims under a^2 = 1/3 and a^2 = 12 are [0, 0, 0, 0]
      each (the pair solves c(a) - 1 = -3(a^2-2)^2/a^2 at c = -24,
      product a_+ a_- = 2), and the joint screened dim at weight 10
      is 0.
  T3  a^2 = 2, 3 likewise screen nothing (cross-validated against the
      t_1 template's screened_dim at the integer values).
  T4  the spin-9 kernel is genuinely mixed (projection ranks 2, 2).

Machinery loaded from the certified cyl_jump_point_tower.py prefix; a
generalized screened machinery handles fractional a^2 and joint
screening.  Fail-closed.  Marker: T2 SHADOW TEST CERTIFIED."""
import sys
from fractions import Fraction as F
from pathlib import Path

_tpl = Path(__file__).resolve().parent / "cyl_jump_point_tower.py"
_src = _tpl.read_text()
_ns = {"__file__": str(_tpl)}
exec(_src[:_src.index("N, s = F(-25, 2)")], _ns)
screened_dim_tpl = _ns["screened_dim"]
mixed_projection_ranks = _ns["mixed_projection_ranks"]
gen_mixed_full = _ns["gen_mixed_full"]
residue_x_quadratic = _ns["residue_x_quadratic"]
d_mixed_mono = _ns["d_mixed_mono"]
canon_mixed = _ns["canon_mixed"]
SparsePivots = _ns["SparsePivots"]
lcm = _ns["lcm"]
_rref = _ns["_rref"]
cyl_P4 = _ns["cyl_P4"]
gen_mixed_basis = _ns["gen_mixed_basis"]
genuine_kernel_mixed = _ns["genuine_kernel_mixed"]

fails = []


def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


def screened_rows(currents, w, r):
    """Per-current reduced residue rows against the e^{sqrt(r) X} twisted
    image; generalizes the template to Fraction r (rows scaled integral)."""
    s_ = w - 1
    tgt = []
    for t in range(0, s_ + 1):
        tgt += gen_mixed_full(t)
    lower = gen_mixed_full(s_ - 1)
    index = {m: i for i, m in enumerate(tgt)}
    n = len(tgt)
    rN, rD = F(r).numerator, F(r).denominator
    piv = SparsePivots()
    for u in lower:
        du = d_mixed_mono(u)
        xu = canon_mixed(u[0] + (1,), u[1])
        col = {}
        for m, c in du.items():
            col[index[m]] = col.get(index[m], 0) + int(c)
        col[n + index[xu]] = col.get(n + index[xu], 0) + 1
        piv.insert({k: c for k, c in col.items() if c})
        col2 = {index[xu]: rN}
        for m, c in du.items():
            col2[n + index[m]] = col2.get(n + index[m], 0) + int(c) * rD
        piv.insert({k: c for k, c in col2.items() if c})
    rows = []
    for P in currents:
        res = residue_x_quadratic(P, F(r))
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
    return rows


def screened_dim_gen(currents, w, rs):
    """dim of the subspace of span(currents) screened under e^{sqrt(r) X}
    for EVERY r in rs (joint kernel of the stacked reduced rows)."""
    stacked = [{} for _ in currents]
    off = 0
    for r in rs:
        rows = screened_rows(currents, w, r)
        mx = 0
        for i, row in enumerate(rows):
            for k, x in row.items():
                stacked[i][off + k] = x
                mx = max(mx, k)
        off += mx + 1
    keys = sorted({k for row in stacked for k in row})
    if not keys:
        return len(currents)
    M = [[row.get(k, F(0)) for k in keys] for row in stacked]
    _, pivots = _rref([list(m) for m in M], len(keys))
    return len(currents) - len(pivots)


if __name__ == "__main__":
    N, s = F(-24), F(-35)
    require(s * s == (N - 25) * (N - 1), "(-24, -35) on the curve")
    a2p, a2m = F(1, 3), F(12)
    for a2 in (a2p, a2m):
        require(-3 * (a2 - 2) ** 2 / a2 == N - 1,
                f"a^2 = {a2} solves c(a) - 1 = -3(a^2-2)^2/a^2 at c = -24 (c_{{1,6}} pair)")
    require(a2p * a2m == 4, "dual-pair product a_+^2 a_-^2 = 4 (a_+ a_- = 2)")
    P3 = cyl_P4(3, N, s)
    EXPECT = {6: 1, 8: 1, 10: 2, 12: 1}
    kernels = {}
    for w, kdim in EXPECT.items():
        ker = genuine_kernel_mixed(P3, gen_mixed_basis(w), N)
        kernels[w] = ker
        require(len(ker) == kdim,
                f"T1: Solution 3 at t_2, spin {w-1}: kernel dim {len(ker)} (expected {kdim})")
        dp = screened_dim_gen(ker, w, [a2p])
        dm = screened_dim_gen(ker, w, [a2m])
        require(dp == 0 and dm == 0,
                f"T2: spin {w-1}: screened dims a^2=1/3: {dp}, a^2=12: {dm} (expected 0, 0 -- "
                "no screening structure, extended to the Kac pair)")
        d2 = screened_dim_gen(ker, w, [F(2)])
        d3 = screened_dim_gen(ker, w, [F(3)])
        d2t = screened_dim_tpl(ker, w, 2)
        d3t = screened_dim_tpl(ker, w, 3)
        require(d2 == d3 == 0 and (d2, d3) == (d2t, d3t),
                f"T3: spin {w-1}: controls a^2=2: {d2} (tpl {d2t}), a^2=3: {d3} (tpl {d3t}) "
                "(expected 0, 0, generalized == template)")
    dj = screened_dim_gen(kernels[10], 10, [a2p, a2m])
    require(dj == 0, f"T2: spin 9: JOINT screened dim {dj} (expected 0)")
    rx, ry = mixed_projection_ranks(kernels[10], 10)
    require((rx, ry) == (2, 2),
            f"T4: spin-9 kernel projection ranks ({rx}, {ry}) -- genuinely mixed")
    print()
    if fails:
        print(f"T2 SHADOW: {len(fails)} PIN(S) NOT MET -- see FINDING lines")
        sys.exit(1)
    print("T2 SHADOW TEST CERTIFIED (Solution 3 at t_2 = c_{1,6}: kernel dims [1,1,2,1] at "
          "spins 5-11; NO screening structure, now including the c_{1,6} Kac pair a^2 = 1/3, 12 "
          "-- individually and jointly zero at every weight; kernel genuinely mixed; the "
          "preregistered screened-shadow reading is refuted and recorded)")
