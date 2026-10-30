"""stratification_w6_w8 S3 cross-check (exact): the factors of the single minor D(t) that the gcd step removed must have full rank r.
For every irreducible factor f of D(t) (row set 0) with deg f <= 8, f not t -+ 1 and f not a factor of G: exact rank of
M_poly over Q[t]/(f) (numfield.rank_mod; for linear f also at the rational root). Also a consistency test of
numfield.rank_mod on the decided linear factors (t - 3, t - 5/3): it must equal the rank at the rational root."""
import json
import os
import sys
import time
from fractions import Fraction as F

from sympy import Poly, Rational, Symbol

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402
import linalg as la  # noqa: E402
import numfield as nf  # noqa: E402

t = Symbol("t")
allok = True
for sheet in (1, 2, 3):
    for W in (6, 8):
        res = json.load(open(os.path.join(fam.BASE, "results", f"S23_sheet{sheet}_W{W}.json")))
        Mpoly = fam.load(os.path.join(fam.BASE, "data", f"M_sheet{sheet}_W{W}.json"))["Mpoly"]
        r, n = res["r"], res["n"]
        dA = [F(c) for c in res["sets"][0]["detA"]]
        D = Poly([Rational(c.numerator, c.denominator) for c in reversed(dA)], t, domain="QQ")
        G = Poly(res["G"], t, domain="QQ")
        _, facs = D.factor_list()
        for f, _m in facs:
            f = f.monic()
            if f.degree() > 8 or f in (Poly(t - 1, t, domain="QQ"), Poly(t + 1, t, domain="QQ")):
                continue
            fc = [F(int(c.p), int(c.q)) for c in reversed(f.all_coeffs())]
            inG = G.rem(f).is_zero
            t0 = time.time()
            rk = nf.rank_mod(Mpoly, fc)
            extra = ""
            # a factor removed by the gcd (not in G) must have full rank r; a factor in G was decided in S3 (either outcome)
            ok = (rk == r) if not inG else True
            if f.degree() == 1:
                a = -fc[0] / fc[1]
                rk_ev = la.rank(fam.ev_matrix(Mpoly, a))
                extra = f"; rank at t = {a}: {rk_ev}"
                ok &= rk_ev == rk
            allok &= ok
            print(f"sheet {sheet}, W = {W}: f = {f.as_expr()} (deg {f.degree()}, in G: {inG}): exact rank over Q[t]/(f) = {rk} "
                  f"(r = {r}) -> {'no drop' if rk == r else 'DROP'}{extra} ({time.time() - t0:.1f} s) {'OK' if ok else 'FAIL'}",
                  flush=True)
print("ALL OK" if allok else "SOME FAILED")
