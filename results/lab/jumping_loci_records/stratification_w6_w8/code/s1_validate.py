"""stratification_w6_w8 S1 validation: at 5 seeded random rational t per (sheet, W), plus t = 3 and t = 5/3, compare
  (a) nullity of my M_poly(t0)                      with len(genuine_kernel_mixed(cyl_P4(sol, N0, s0), gen_mixed_basis(W), N0));
  (b) every entry of M(t0) = M_poly(t0) / ((t0^2-1)^k c0) with the engine's numeric column pipeline at t0
      (reduce_mod_d of residue_mixed_cur(cyl_P4 numeric, {b: 1}, N0), class coordinates);
  (c) my Q(t) I_3 coefficients evaluated at t0     with the engine's numeric cyl_P4(sol, N0, s0).
Seeds: 68000 + 10*sheet + W (random.Random), t0 = p/q with p in [-99, 99], q in [1, 49], t0 not in {-1, 0, 1}."""
import os
import random
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402
import linalg as la  # noqa: E402

me = fam.me


def N_s(t0):
    return (t0 * t0 - 25) / (t0 * t0 - 1), F(-24) * t0 / (t0 * t0 - 1)


def evK(e, t0):
    """Evaluate an element of K = Q(t) at a Fraction t0 (denominator nonzero there)."""
    num = sum((F(int(c.numerator), int(c.denominator)) * t0**d for (d,), c in e.numer.terms()), F(0))
    den = sum((F(int(c.numerator), int(c.denominator)) * t0**d for (d,), c in e.denom.terms()), F(0))
    return num / den


def points(sheet, W):
    seed = 68000 + 10 * sheet + W
    rnd = random.Random(seed)
    pts = []
    while len(pts) < 5:
        t0 = F(rnd.randint(-99, 99), rnd.randint(1, 49))
        if t0 in (-1, 0, 1) or t0 in pts or t0 in (3, F(5, 3)):
            continue
        pts.append(t0)
    return seed, pts + [F(3), F(5, 3)]


allok = True
print("sheet W  seed   t0        N0           s0           nullity(M)  len(genuine)  entries-mismatch  I3-coeffs  verdict")
for sheet in (1, 2, 3):
    P4K = fam.p4_family(sheet)
    for W in (6, 8):
        d = fam.load(os.path.join(fam.BASE, "data", f"M_sheet{sheet}_W{W}.json"))
        k, c0, Mpoly = d["k"], d["c0"], d["Mpoly"]
        cs, tgt, cols, rows = fam.domain_target(W)
        assert [fam.mono_str(cs.basis[j]) for j in cols] == d["cols"]
        assert [fam.mono_str(tgt.basis[i]) for i in rows] == d["rows"]
        seed, pts = points(sheet, W)
        for t0 in pts:
            N0, s0 = N_s(t0)
            assert s0 * s0 == (N0 - 25) * (N0 - 1)
            Mt = fam.ev_matrix(Mpoly, t0)
            nullity = len(cols) - la.rank(Mt)
            P4num = me.cyl_P4(sheet, N0, s0)
            g = me.genuine_kernel_mixed(P4num, me.gen_mixed_basis(W), N0)
            # (b) entrywise
            scale = (t0 * t0 - 1) ** k * c0
            mism = 0
            for jj, j in enumerate(cols):
                col = tgt.reduce_mod_d(tgt.vec(me.residue_mixed_cur(P4num, {cs.basis[j]: F(1)}, N0)))
                for ii, i in enumerate(rows):
                    if Mt[ii][jj] / scale != col[i]:
                        mism += 1
                assert not any(col[p] for p in tgt._img_pivots)
            # (c) I_3 coefficients
            P4ev = {m: evK(e, t0) for m, e in P4K.items()}
            P4ev = {m: c for m, c in P4ev.items() if c}
            coeff_ok = P4ev == P4num
            ok = (nullity == len(g)) and mism == 0 and coeff_ok
            allok &= ok
            print(f"{sheet:5d} {W:2d} {seed:5d}  {str(t0):8s}  {str(N0):11s}  {str(s0):11s}  {nullity:10d}  {len(g):12d}  "
                  f"{mism:16d}  {'equal' if coeff_ok else 'DIFFER':9s}  {'OK' if ok else 'FAIL'}", flush=True)
print("ALL OK" if allok else "SOME FAILED")
