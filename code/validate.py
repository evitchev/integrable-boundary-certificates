"""Validation of the OPE engine against known results.

V1: Res T(z)T(w) = dT  (so [I_1, I_1] = 0), and [I_1, J] = 0 for arbitrary
    currents J (T generates translations).
V2: Vitchev eqs. (28)-(29): the N=1 KdV-Virasoro charges I_3, I_5 commute.
V3: Vitchev eqs. (47)-(51): the paperclip Z2xZ2 charges I_3 = oint P_4 and
    I_5 = oint P_6 commute, for several rational values of the parameter n.
    (Eq. (50) has an obvious typo: its right side is written with P_6^{[i]}
    but must be the spin-4 basis P_4^{[i]} of eq. (47).)
"""

from fractions import Fraction as F

from ope_engine import Sector, residue_cur, d_cur
from families import L, KDV_P4, KDV_P6, paperclip_P4, paperclip_P6


def check(label, ok):
    print(("PASS  " if ok else "FAIL  ") + label)
    if not ok:
        raise SystemExit(1)


def main():
    # ------------------------------------------------------------ V1
    T = {L((1, 0), (1, 0)): F(1, 2)}
    check("V1a: Res T(z)T(w) == dT", residue_cur(T, T) == d_cur(T))

    J = {L((2, 0), (2, 0), (1, 0)): F(3, 7),
         L((1, 0),) * 5: F(-2),
         L((3, 0), (1, 0), (1, 0)): F(5, 3)}
    check("V1b: [oint T, oint J] = 0 for a random J",
          Sector(6, 1, z2=False).is_total_derivative(residue_cur(T, J)))

    # ------------------------------------------------------------ V2
    check("V2: KdV N=1 [I_3, I_5] = 0",
          Sector(9, 1, z2=False).is_total_derivative(
              residue_cur(KDV_P4, KDV_P6)))

    # ------------------------------------------------------------ V3
    sec9 = Sector(9, 2, z2=True)
    for n in [F(1), F(2), F(3), F(7, 3), F(-1), F(13, 5)]:
        res = residue_cur(paperclip_P4(n), paperclip_P6(n))
        check(f"V3: paperclip [I_3, I_5] = 0 at n = {n}",
              sec9.is_total_derivative(res))

    # ------------------------------------------------------------ V4
    # negative control: a perturbed coefficient must break commutativity
    bad = dict(KDV_P6)
    bad[L((3, 0), (3, 0))] = bad[L((3, 0), (3, 0))] + F(1, 18)
    check("V4: perturbed KdV I_5 fails to commute (negative control)",
          not Sector(9, 1, z2=False).is_total_derivative(
              residue_cur(KDV_P4, bad)))

    # ------------------------------------------------------------ V5
    # antisymmetry mod d: Res(P,Q) + Res(Q,P) is always a total derivative
    s = residue_cur(KDV_P4, KDV_P6)
    for m, c in residue_cur(KDV_P6, KDV_P4).items():
        s[m] = s.get(m, 0) + c
    check("V5: Res(P,Q) + Res(Q,P) is a total derivative",
          Sector(9, 1, z2=False).is_total_derivative(
              {m: c for m, c in s.items() if c}))

    print("\nAll validations passed.")


if __name__ == "__main__":
    main()
