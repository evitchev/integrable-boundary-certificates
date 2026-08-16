"""Symbolic-in-n verification that the paperclip charges commute:
[I_3, I_5] = 0 and [I_3, I_7] = 0 identically in n (charges built from
eq. (50), corrected eq. (51), and the newly computed P_8).

Trick: the residue coefficients residue_mono are field-independent exact
Fractions and the total-derivative subspace has a NUMERIC basis, so with
denominators cleared all symbolic work is polynomial arithmetic in n.
"""

import json
from fractions import Fraction as F

import sympy as sp

from ope_engine import Sector, residue_mono
from families import L, P4B, P6B, paperclip_P4, paperclip_P6

n = sp.Symbol('n')


def paperclip_P4_sym():
    c = [
        -(4 * n**2 + 7 * n + 2) / (3 * (3 * n + 2) * (3 * n + 4)),
        -(4 * n**2 + 9 * n + 4) / (3 * (3 * n + 2) * (3 * n + 4)),
        (n + 2) / (6 * (3 * n + 4)),
        n * (n + 2) / ((3 * n + 2) * (3 * n + 4)),
        n / (6 * (3 * n + 2)),
    ]
    return {b: sp.cancel(ci) for b, ci in zip(P4B, c)}


def paperclip_P6_sym():
    c = [
        -(96 + 500 * n + 860 * n**2 + 575 * n**3 + 131 * n**4)
        / (120 * n * (2 + n)),
        (2 + 11 * n + 6 * n**2) / 4,
        (4 + 14 * n + 7 * n**2) * sp.Integer(1),
        -(32 + 232 * n + 554 * n**2 + 473 * n**3 + 131 * n**4)
        / (120 * n * (2 + n)),
        (4 + 13 * n + 6 * n**2) / 4,
        (2 + 5 * n) * (4 + 9 * n + 4 * n**2) / (4 * n),
        (8 + 5 * n) * (2 + 7 * n + 4 * n**2) / (4 * (2 + n)),
        -(2 + n) * (2 + 5 * n) * (4 + 5 * n) / (120 * n),
        -(2 + n) * (2 + 5 * n) / 8,
        -n * (8 + 5 * n) / 8,
        -n * (6 + 5 * n) * (8 + 5 * n) / (120 * (2 + n)),  # corrected
    ]
    return {b: sp.cancel(ci) for b, ci in zip(P6B, c)}


def paperclip_P8_sym():
    with open("../results/paperclip_P8.json") as f:
        data = json.load(f)
    from compute_p8 import parse_mono
    return {parse_mono(s): sp.cancel(sp.sympify(e, locals={'n': n}))
            for s, e in zip(data["basis"], data["coefficients"])}


def residue_cur_sym(P, Q):
    out = {}
    for mp, cp in P.items():
        for mq, cq in Q.items():
            for mono, c in residue_mono(mp, mq).items():
                out[mono] = out.get(mono, 0) + \
                    sp.Rational(c.numerator, c.denominator) * cp * cq
    return out


def clear_denominators(cur):
    den = sp.lcm([sp.denom(sp.cancel(c)) for c in cur.values()])
    return {m: sp.expand(sp.cancel(c * den)) for m, c in cur.items()}, den


def is_total_derivative_sym(sec, cur):
    v = [sp.Integer(0)] * len(sec.basis)
    for mono, c in cur.items():
        v[sec.index[mono]] = c
    for row, p in zip(sec._img_rows, sec._img_pivots):
        f = v[p]
        if f == 0:
            continue
        v = [sp.expand(a - f * sp.Rational(b.numerator, b.denominator))
             for a, b in zip(v, row)]
    return all(sp.expand(x) == 0 for x in v)


if __name__ == "__main__":
    P4s, _ = clear_denominators(paperclip_P4_sym())
    P6s, _ = clear_denominators(paperclip_P6_sym())
    P8s, _ = clear_denominators(paperclip_P8_sym())

    res = residue_cur_sym(P4s, P6s)
    res = {m: sp.expand(c) for m, c in res.items()}
    ok = is_total_derivative_sym(Sector(9, 2, z2=True), res)
    print(f"[I_3, I_5] = 0 identically in n: {ok}")
    assert ok

    res = residue_cur_sym(P4s, P8s)
    res = {m: sp.expand(c) for m, c in res.items()}
    ok = is_total_derivative_sym(Sector(11, 2, z2=True), res)
    print(f"[I_3, I_7] = 0 identically in n: {ok}")
    assert ok

    print("\nSymbolic verification complete: the paperclip charges "
          "I_3, I_5, I_7 commute with I_3 for ALL n (poles excluded).")
