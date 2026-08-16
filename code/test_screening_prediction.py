"""Decisive test: the screening vectors alpha = (±sqrt(-n/2),
±sqrt((n+2)/2)), fitted ONLY on [Q_alpha, I_3] = [Q_alpha, I_5] = 0,
must also annihilate the independently computed spin-7 charge oint P_8.

Checked at n = 2 (alpha = (±i, ±sqrt(2))) and n = 3
(alpha = (±i sqrt(6)/2, ±sqrt(10)/2)) for all sign choices.
"""

import json
from fractions import Fraction as F

import sympy as sp

from families import paperclip_P4, paperclip_P6
from screening_fit import screening_conditions
from compute_p8 import parse_mono

with open("../results/paperclip_P8.json") as f:
    DATA = json.load(f)
nsym = sp.Symbol('n')
P8_EXPR = {parse_mono(s): sp.cancel(sp.sympify(e, locals={'n': nsym}))
           for s, e in zip(DATA["basis"], DATA["coefficients"])}


def paperclip_P8(nval):
    out = {}
    for mono, e in P8_EXPR.items():
        v = sp.Rational(e.subs(nsym, sp.Rational(nval)))
        if v != 0:
            out[mono] = F(int(sp.numer(v)), int(sp.denom(v)))
    return out


def alpha_vectors(nval):
    a = sp.sqrt(sp.Rational(-nval, 2))
    b = sp.sqrt(sp.Rational(nval + 2, 2))
    return [(sa * a, sb * b) for sa in (1, -1) for sb in (1, -1)]


for nval in [F(2), F(3)]:
    charges = [("I_3", paperclip_P4(nval)),
               ("I_5", paperclip_P6(nval)),
               ("I_7 (PREDICTION)", paperclip_P8(nval))]
    print(f"n = {nval}: alpha in {alpha_vectors(nval)}")
    for label, P in charges:
        for al in alpha_vectors(nval):
            conds = screening_conditions(P, al, 2)
            conds = [c for c in conds if sp.simplify(c) != 0]
            status = "annihilated" if not conds else f"FAILS: {conds}"
            print(f"  [Q_alpha, {label:18s}] alpha={str(al):24s} {status}")
        if any(screening_conditions(P, al, 2) for al in alpha_vectors(nval)):
            raise SystemExit(1)
print("\nAll four screenings annihilate I_3, I_5 AND the predicted I_7.")
