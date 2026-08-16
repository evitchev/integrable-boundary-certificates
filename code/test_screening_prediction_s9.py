"""Spin-9 extension of the decisive screening test: the four vectors
alpha = (±sqrt(-n/2), ±sqrt((n+2)/2)), fitted only on I_3 and I_5,
must also annihilate the independently computed spin-9 charge
oint P_10 (upgrades the LVZ four-screening verification from spin 7,
cf. test_screening_prediction.py).

Checked at n = 2 and n = 3 for all four sign choices."""
import json
import sys
from fractions import Fraction as F

import sympy as sp

from screening_fit import screening_conditions

with open("../results/paperclip_P10.json") as f:
    DATA = json.load(f)
nsym = sp.Symbol('n')
P10_EXPR = {tuple(tuple(l) for l in m): sp.cancel(
                sp.sympify(e, locals={'n': nsym}))
            for m, e in zip(DATA["basis_monos"], DATA["coefficients"])}


def paperclip_P10(nval):
    out = {}
    for mono, e in P10_EXPR.items():
        v = sp.Rational(e.subs(nsym, sp.Rational(nval)))
        if v != 0:
            out[mono] = F(int(sp.numer(v)), int(sp.denom(v)))
    return out


def alpha_vectors(nval):
    a = sp.sqrt(sp.Rational(-nval, 2))
    b = sp.sqrt(sp.Rational(nval + 2, 2))
    return [(sa * a, sb * b) for sa in (1, -1) for sb in (1, -1)]


ok = True
for nval in [F(2), F(3)]:
    P10 = paperclip_P10(nval)
    for alph in alpha_vectors(nval):
        conds = screening_conditions(P10, alph, 2)
        residual = [c for c in conds if sp.simplify(c) != 0]
        verdict = not residual
        print(f"n = {nval}, alpha = {alph}: annihilates P10: "
              f"{verdict}", flush=True)
        ok = ok and verdict

print("SPIN-9 SCREENING TEST " + ("PASSED" if ok else "FAILED"),
      flush=True)
sys.exit(0 if ok else 1)
