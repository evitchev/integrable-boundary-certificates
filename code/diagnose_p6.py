"""Solve for the genuine spin-5 charge commuting with the paperclip I_3
at sample n, express it in Vitchev's Q^5 representative basis (eq. 48),
and compare against the transcription of eq. (51)."""

from fractions import Fraction as F

from ope_engine import (Sector, gen_basis, residue_cur, nullspace, rref,
                        fmt_mono)
from families import P4B, P6B, paperclip_P4, paperclip_P6, L


def genuine_kernel(P4, basis6, sec9, sec6):
    cols = []
    for b in basis6:
        r = residue_cur(P4, {b: F(1)})
        cols.append(sec9.reduce_mod_d(sec9.vec(r)))
    ker = nullspace(cols, len(sec9.basis))
    # remove total-derivative directions: reduce each kernel current mod d,
    # keep independent nonzero images
    reduced = []
    genuine = []
    for v in ker:
        cur = {b: c for b, c in zip(basis6, v) if c}
        w = sec6.reduce_mod_d(sec6.vec(cur))
        if not any(w):
            continue
        # independence test against previously kept
        test = reduced + [w]
        rows = [list(x) for x in test]
        _, piv = rref(rows, len(w))
        if len(piv) == len(test):
            reduced.append(w)
            genuine.append(cur)
    return genuine


def in_vitchev_basis(cur, sec6):
    """Coefficients of a spin-6 current in Vitchev's representatives P6B,
    modulo total derivatives."""
    cols = [sec6.reduce_mod_d(sec6.vec({b: F(1)})) for b in P6B]
    target = sec6.reduce_mod_d(sec6.vec(cur))
    nrows = len(target)
    aug = [[cols[c][r] for c in range(len(cols))] + [target[r]]
           for r in range(nrows)]
    aug, piv = rref(aug, len(cols) + 1)
    if len(cols) in piv:
        raise ValueError("current not in span of Vitchev representatives")
    sol = [F(0)] * len(cols)
    for i, p in enumerate(piv):
        sol[p] = aug[i][len(cols)]
    return sol


basis6 = gen_basis(6, 2, z2=True)
sec9 = Sector(9, 2, z2=True)
sec6 = Sector(6, 2, z2=True)
print(f"full B^6 Z2xZ2 basis: {len(basis6)} monomials, "
      f"dim Q^5 = {sec6.q_dim()}")

for n in [F(2), F(3), F(7, 3)]:
    gen = genuine_kernel(paperclip_P4(n), basis6, sec9, sec6)
    print(f"\nn = {n}: {len(gen)} genuine commuting spin-5 charge(s)")
    for cur in gen:
        sol = in_vitchev_basis(cur, sec6)
        # normalize like eq. (51): coefficient of P6^[2] = 4+14n+7n^2
        scale = (4 + 14 * n + 7 * n**2) / sol[2] if sol[2] else F(1)
        sol = [s * scale for s in sol]
        mine = paperclip_P6(n)
        my_sol = in_vitchev_basis(mine, sec6)
        for i, (a, b) in enumerate(zip(sol, my_sol)):
            flag = "   " if a == b else "<<<"
            print(f"  P6[{i:2d}] solved: {str(a):>16}   "
                  f"transcribed: {str(b):>16} {flag}")
