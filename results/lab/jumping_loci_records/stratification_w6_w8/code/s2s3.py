"""stratification_w6_w8 S2 + S3 for one matrix file data/M_*.json.

S2  r = rank of M over Q(t):
      lower bound: rank M_poly(t*) at seeded random rational points; the pivot block A = M_poly[R, C] has det A(t*) != 0,
                   so det A is a nonzero polynomial and r >= r0;
      upper bound: for every non-pivot column j, the Cramer vector w^(j) = det A * (A^{-1} M_poly[R, j] on C, -1 at j)
                   (entries = r0 x r0 minors of M_poly[R], interpolated exactly under a degree bound) satisfies
                   M_poly * w^(j) = 0 identically in Q[t] (checked by exact polynomial multiplication, all rows), and the
                   w^(j) are independent (w^(j)_j = -det A), so r <= n - #free = r0.
    Over-approximation: for several row sets R (seeded random row orders), all the r x r minors of M_poly[R] met above
    (det A and the entries of the w^(j)); G = gcd of all of them. Complete: if rank M(t0) < r at some t0 in Q-bar,
    t0 != +-1, every r x r minor of M_poly = (t^2-1)^k c0 M vanishes at t0, so t0 is a common root, i.e. G(t0) = 0.
S3  For every irreducible factor f of G over Q: exact rank of M_poly over Q[t]/(f) (for linear f, at the rational root);
    nullity = n - rank; a drop iff rank < r. Roots t = +-1 are excluded points (reported, not decided).
"""
import argparse
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as F

from sympy import Poly, Rational, Symbol

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402
import linalg as la  # noqa: E402
import numfield as nf  # noqa: E402

t = Symbol("t")


def to_poly(coeffs):
    cs = [Rational(F(c).numerator, F(c).denominator) for c in coeffs] or [Rational(0)]
    return Poly(list(reversed(cs)), t, domain="QQ")


def from_poly(p):
    return [F(int(c.p), int(c.q)) for c in reversed(p.all_coeffs())]


def pmul_int(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def deg(e):
    return len(e) - 1 if e and any(e) else -1


def int_points():
    k = 0
    while True:
        yield F(k)
        if k:
            yield F(-k)
        k += 1


def fmt_poly(p):
    return str(p.as_expr())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--nsets", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    T0 = time.time()
    d = fam.load(args.matrix)
    Mpoly = d["Mpoly"]
    m, n = len(Mpoly), len(Mpoly[0])
    print(f"matrix {args.matrix}: sheet {d['sheet']}, W = {d['W']}, {m} x {n}, {d['scalar']}; note: {d.get('note', '')}")
    rnd = random.Random(args.seed)
    print(f"seed {args.seed}")

    # ---- generic rank at random points
    pts = []
    while len(pts) < args.nsets:
        t0 = F(rnd.randint(-999, 999), rnd.randint(1, 999))
        if t0 not in (-1, 1) and t0 not in pts:
            pts.append(t0)
    Mts = [fam.ev_matrix(Mpoly, t0) for t0 in pts]
    ranks = [la.rank(Mt) for Mt in Mts]
    r0 = max(ranks)
    print(f"random points {[str(p) for p in pts]}: ranks {ranks}; r0 = {r0}")

    # ---- row sets, Cramer vectors, minors
    allminors = []
    G = None
    sets = []
    for si in range(args.nsets):
        Mt = Mts[si]
        if la.rank(Mt) != r0:
            print(f"  row set {si}: rank at {pts[si]} is below r0, skipped")
            continue
        perm = list(range(m))
        rnd.shuffle(perm)
        rk, prow, pcol = la.rank_and_pivots([Mt[i] for i in perm])
        assert rk == r0
        R = sorted(perm[i] for i in prow)
        Cset = list(pcol)
        J = [c for c in range(n) if c not in Cset]
        A_star = [[Mt[r][c] for c in Cset] for r in R]
        dA_star = la.det(A_star)
        assert dA_star != 0
        # degree bound for every r0 x r0 minor of M_poly[R] with columns in C u {j}
        S = Cset + J
        rowb = sum(max(deg(Mpoly[r][c]) for c in S) for r in R)
        colmax = sorted((max(deg(Mpoly[r][c]) for r in R) for c in S), reverse=True)
        colb = sum(colmax[:r0])
        B = min(rowb, colb)
        npts = B + 1 + 2           # two extra points check the interpolation
        xs, vals = [], []          # vals[i] = {j: w^(j)(x_i)} and det A(x_i)
        gen = int_points()
        while len(xs) < npts:
            x = next(gen)
            Mx = fam.ev_matrix([Mpoly[r] for r in R], x)
            A = [[Mx[k][c] for c in Cset] for k in range(len(R))]
            try:
                sol = {}
                dA = None
                for j in J:
                    xsol, dA = la.solve(A, [Mx[k][j] for k in range(len(R))])
                    sol[j] = (xsol, dA)
                if not J:
                    dA = la.det(A)
                    if dA == 0:
                        raise ZeroDivisionError
            except ZeroDivisionError:
                continue
            xs.append(x)
            vals.append((dA, sol))
        # interpolate det A and every w^(j)
        dApoly = la.interpolate(xs[:B + 1], [v[0] for v in vals[:B + 1]])
        for x, v in zip(xs[B + 1:], vals[B + 1:]):
            assert fam.ev([c for c in dApoly], x) == v[0], "interpolation check failed (det A)"
        ws = {}
        for j in J:
            w = []
            for c in range(n):
                if c in Cset:
                    k = Cset.index(c)
                    ys = [v[1][j][1] * v[1][j][0][k] for v in vals]
                elif c == j:
                    ys = [-v[0] for v in vals]
                else:
                    ys = [F(0)] * len(vals)
                p = la.interpolate(xs[:B + 1], ys[:B + 1])
                for x, y in zip(xs[B + 1:], ys[B + 1:]):
                    assert fam.ev(p, x) == y, "interpolation check failed (w)"
                w.append(p)
            # clear to integers and verify M_poly * w = 0 identically (all m rows)
            den = 1
            for p in w:
                for c in p:
                    den = math.lcm(den, c.denominator)
            wi = [[int(c * den) for c in p] if any(p) else [] for p in w]
            for r in range(m):
                acc = []
                for c in range(n):
                    pr = pmul_int(Mpoly[r][c], wi[c])
                    if len(pr) > len(acc):
                        acc += [0] * (len(pr) - len(acc))
                    for i, a in enumerate(pr):
                        acc[i] += a
                assert not any(acc), f"M_poly * w^({j}) != 0 at row {r}"
            ws[j] = wi
        minors = [to_poly(dApoly)] + [to_poly(p) for j in J for p in ws[j] if any(p)]
        GR = None
        for p in minors:
            GR = p if GR is None else GR.gcd(p)
        GR = GR.monic()
        print(f"  row set {si}: point {pts[si]}, R = {R}, C = {Cset}, free columns {J}; det A(t*) = {dA_star} != 0; "
              f"degree bound {B} (row {rowb}, column {colb}); {len(xs)} interpolation points; deg det A = {deg(dApoly)}; "
              f"M_poly * w^(j) = 0 identically for every free j: verified; gcd of this set's {len(minors)} minors: degree {GR.degree()}",
              flush=True)
        sets.append({"point": str(pts[si]), "R": R, "C": Cset, "J": J, "detA_tstar": str(dA_star), "B": B,
                     "detA": [str(c) for c in dApoly], "gcd_deg": GR.degree(),
                     "w": {str(j): ws[j] for j in J}})
        G = GR if G is None else G.gcd(GR).monic()
    r = r0
    print(f"r = rank over Q(t) = {r} (n = {n}; generic nullity n - r = {n - r})")
    D = to_poly([F(c) for c in [F(x) for x in sets[0]['detA']]])
    print(f"D(t) = det A of row set 0: degree {D.degree()}")
    print(f"G(t) = gcd over {len(sets)} row sets: degree {G.degree()}: {fmt_poly(G)}")
    cG, facs = G.factor_list()
    print(f"factorization of G over Q: {len(facs)} irreducible factor(s)")
    decisions = []
    fpolys = []
    for f, mult in facs:
        f = f.monic()
        fpolys.append(f)
        fc = from_poly(f)
        rec = {"factor": fmt_poly(f), "degree": f.degree(), "multiplicity_in_G": mult}
        if f.degree() == 1:
            a = -fc[0] / fc[1]
            rec["root"] = str(a)
            if a in (1, -1):
                rec["decision"] = "excluded point t = +-1"
                print(f"  factor {fmt_poly(f)} (mult {mult}): root t = {a}: EXCLUDED (pole)")
                decisions.append(rec)
                continue
            rk = la.rank(fam.ev_matrix(Mpoly, a))
            Nv = (a * a - 25) / (a * a - 1)
            sv = F(-24) * a / (a * a - 1)
            rec.update({"beta2": str((a - 1) / (a + 1)), "N": str(Nv), "s": str(sv)})
        else:
            fz = [c for c in fc]
            rk = nf.rank_mod(Mpoly, fz)
        rec.update({"rank": rk, "nullity": n - rk, "drop": rk < r})
        print(f"  factor {fmt_poly(f)} (deg {f.degree()}, mult {mult}): exact rank {rk} over Q[t]/(f) -> nullity {n - rk}; "
              f"{'DROP' if rk < r else 'no drop'}" + (f"; t = {rec['root']}, beta^2 = {rec['beta2']}, N = {rec['N']}, s = {rec['s']}"
                                                       if 'beta2' in rec else ""), flush=True)
        decisions.append(rec)
    # jump set: nullity >= 2
    jumps = [x for x in decisions if x.get("nullity", 0) >= 2]
    print(f"jump set (m >= 2): {[x['factor'] for x in jumps]}")
    for x in jumps:
        f = fpolys[decisions.index(x)]
        fneg = Poly(f.as_expr().subs(t, -t), t, domain="QQ").monic()
        hit = next((y for y, fp in zip(decisions, fpolys) if fp == fneg), None)
        if hit is None:
            # -t0 is not a root of G, so rank M(-t0) = r: not even a drop
            x["deck_image"] = f"{fmt_poly(fneg)}: not a factor of G -> rank r, nullity {n - r}: not a jump point"
        else:
            x["deck_image"] = (f"{fmt_poly(fneg)}: factor of G, nullity {hit.get('nullity')}"
                               + (": JUMP" if hit.get("nullity", 0) >= 2 else ": not a jump point"))
        print(f"  deck image of {x['factor']} = 0: {x['deck_image']}")
    out = {"matrix": args.matrix, "sheet": d["sheet"], "W": d["W"], "m": m, "n": n, "seed": args.seed,
           "random_points": [str(p) for p in pts], "ranks": ranks, "r": r, "generic_nullity": n - r,
           "D_degree": D.degree(), "G": fmt_poly(G), "G_degree": G.degree(), "decisions": decisions,
           "jumps": jumps, "sets": sets, "seconds": round(time.time() - T0, 1)}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(f"saved {args.out} ({time.time() - T0:.1f} s)")


if __name__ == "__main__":
    main()
