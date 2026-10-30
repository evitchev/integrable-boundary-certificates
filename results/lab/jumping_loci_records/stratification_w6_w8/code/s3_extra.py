"""stratification_w6_w8 S3 extra checks (exact, direct): nullity of M(t0) at the deck images t0 = -3, -5/3 and at t0 = 0, 9, -9,
5/3, 3 for every (sheet, W); and the factorization over Q of the single minor D(t) = det A of row set 0 (from results/),
to show which of its factors the gcd step removed."""
import json
import os
import sys
from fractions import Fraction as F

from sympy import Poly, Rational, Symbol

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402
import linalg as la  # noqa: E402

t = Symbol("t")
PTS = [F(-3), F(-5, 3), F(0), F(9), F(-9), F(5, 3), F(3)]
print("sheet W   " + "  ".join(f"t={str(p):5s}" for p in PTS) + "   (nullity of M(t0); generic nullity 1)")
for sheet in (1, 2, 3):
    for W in (6, 8):
        d = fam.load(os.path.join(fam.BASE, "data", f"M_sheet{sheet}_W{W}.json"))
        n = len(d["cols"])
        nul = [n - la.rank(fam.ev_matrix(d["Mpoly"], p)) for p in PTS]
        print(f"{sheet:5d} {W:2d}  " + "  ".join(f"{x:7d}" for x in nul), flush=True)
for sheet in (1, 2, 3):
    for W in (6, 8):
        res = json.load(open(os.path.join(fam.BASE, "results", f"S23_sheet{sheet}_W{W}.json")))
        dA = [F(c) for c in res["sets"][0]["detA"]]
        D = Poly([Rational(c.numerator, c.denominator) for c in reversed(dA)], t, domain="QQ")
        c, facs = D.factor_list()
        desc = " * ".join(f"({f.as_expr()})" + (f"^{m}" if m > 1 else "") for f, m in facs)
        print(f"sheet {sheet}, W = {W}: D(t) (degree {D.degree()}) = {c} * {desc}", flush=True)
