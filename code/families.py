"""Vitchev's (hep-th/0404195) published families, transcribed.

Fields: X = 0, Y = 1.  Monomial notation (n1..na|m1..mb) as in the paper.
"""

from fractions import Fraction as F


def L(*pairs):
    """monomial from (m, j) pairs"""
    return tuple(sorted(pairs))


# ---------------------------------------------------- KdV-Virasoro, N = 1

KDV_P4 = {  # eq. (28): -2*(22) + (11)(11)
    L((2, 0), (2, 0)): F(-2),
    L((1, 0), (1, 0), (1, 0), (1, 0)): F(1),
}

KDV_P6 = {}  # eq. (29) at N=1; note (22)(11) == (12)(12) as N=1 monomials
for _mono, _c in [
    (L((3, 0), (3, 0)), F(56 + 1, 18)),
    (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
    (L((1, 0), (2, 0), (1, 0), (2, 0)), F(-2 * (20 + 1), 3)),
    (L((1, 0),) * 6, F(1)),
]:
    KDV_P6[_mono] = KDV_P6.get(_mono, 0) + _c


# ------------------------------------------------- paperclip, Z2 x Z2

P4B = [  # eq. (47)
    L((2, 1), (2, 1)),                       # (|22)
    L((2, 0), (2, 0)),                       # (22|)
    L((1, 1), (1, 1), (1, 1), (1, 1)),       # (|1111)
    L((1, 0), (1, 0), (1, 1), (1, 1)),       # (11|11)
    L((1, 0), (1, 0), (1, 0), (1, 0)),       # (1111|)
]

P6B = [  # eq. (48)
    L((3, 0), (3, 0)),                                 # (33|)
    L((1, 0), (1, 0), (2, 1), (2, 1)),                 # (11|22)
    L((1, 0), (2, 0), (1, 1), (2, 1)),                 # (12|12)
    L((3, 1), (3, 1)),                                 # (|33)
    L((2, 0), (2, 0), (1, 1), (1, 1)),                 # (22|11)
    L((1, 1), (1, 1), (2, 1), (2, 1)),                 # (|1122)
    L((1, 0), (1, 0), (2, 0), (2, 0)),                 # (1122|)
    L((1, 1),) * 6,                                    # (|111111)
    L((1, 0), (1, 0)) + L((1, 1),) * 4,                # (11|1111)
    L((1, 0),) * 4 + L((1, 1), (1, 1)),                # (1111|11)
    L((1, 0),) * 6,                                    # (111111|)
]


def paperclip_P4(n):
    """eq. (50); its right side is printed with P_6^{[i]} but must read the
    spin-4 basis P_4^{[i]} of eq. (47) -- confirmed by [I_3, I_5] = 0."""
    n = F(n)
    c = [
        -(4 * n**2 + 7 * n + 2) / (3 * (3 * n + 2) * (3 * n + 4)),
        -(4 * n**2 + 9 * n + 4) / (3 * (3 * n + 2) * (3 * n + 4)),
        (n + 2) / (6 * (3 * n + 4)),
        n * (n + 2) / ((3 * n + 2) * (3 * n + 4)),
        n / (6 * (3 * n + 2)),
    ]
    return {b: ci for b, ci in zip(P4B, c) if ci}


def paperclip_P6(n):
    """eq. (51), with one correction: the paper prints the P_6^{[10]}
    coefficient as -(6+5n)(8+5n)/(120(2+n)); solving [I_3, I_5] = 0 from
    scratch at n = 2, 3, 7/3 shows it must be -n(6+5n)(8+5n)/(120(2+n))
    (they agree only at n = 1).  Erratum found 2026-07-27, see
    notes/ and code/diagnose_p6.py."""
    n = F(n)
    c = [
        -(96 + 500 * n + 860 * n**2 + 575 * n**3 + 131 * n**4)
        / (120 * n * (2 + n)),
        (2 + 11 * n + 6 * n**2) / 4,
        F(4 + 14 * n + 7 * n**2),
        -(32 + 232 * n + 554 * n**2 + 473 * n**3 + 131 * n**4)
        / (120 * n * (2 + n)),
        (4 + 13 * n + 6 * n**2) / 4,
        (2 + 5 * n) * (4 + 9 * n + 4 * n**2) / (4 * n),
        (8 + 5 * n) * (2 + 7 * n + 4 * n**2) / (4 * (2 + n)),
        -(2 + n) * (2 + 5 * n) * (4 + 5 * n) / (120 * n),
        -(2 + n) * (2 + 5 * n) / 8,
        -n * (8 + 5 * n) / 8,
        -n * (6 + 5 * n) * (8 + 5 * n) / (120 * (2 + n)),  # erratum: paper omits n
    ]
    return {b: ci for b, ci in zip(P6B, c) if ci}
