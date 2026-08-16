"""Compute the (previously unpublished) spin-7 paperclip charge density P_8.

Method: at ~40 sample rational values of n, solve the exact linear system
[oint P_4(n), oint P_8] = 0 for P_8 in the 25-dimensional space Q^7 of
Z2xZ2-symmetric spin-8 densities modulo total derivatives; check the
solution space is one-dimensional (multiplicity one); normalize by the
coefficient of (11111111|) = (dX)^8; reconstruct each coefficient as a
rational function of n from the samples with held-out verification.

Output: results/paperclip_P8.md (human-readable) and .json (exact).
"""

import json
from fractions import Fraction as F

import sympy as sp

from ope_engine import Sector, gen_basis, fmt_mono, rref
from families import L, paperclip_P4
from solve_charges import genuine_kernel, express_in, rat_reconstruct, rat_eval

# Vitchev's Q^7 representatives, eq. (49)
P8B_NOTATION = [
    "(33|11)", "(44|)", "(13|13)", "(1133|)", "(1223|)", "(|44)",
    "(14|12)", "(|111122)", "(|1133)", "(11|1122)", "(12|1112)",
    "(11|33)", "(22|1111)", "(22|22)", "(1111|22)", "(1112|12)",
    "(23|12)", "(1122|11)", "(|1223)", "(111122|)", "(|11111111)",
    "(11|111111)", "(1111|1111)", "(111111|11)", "(11111111|)",
]


def parse_mono(s):
    xs, ys = s.strip("()").split("|")
    return L(*[(int(ch), 0) for ch in xs], *[(int(ch), 1) for ch in ys])


P8B = [parse_mono(s) for s in P8B_NOTATION]


if __name__ == "__main__":
    sec8 = Sector(8, 2, z2=True)
    sec11 = Sector(11, 2, z2=True)
    basis8 = gen_basis(8, 2, z2=True)
    print(f"dim B^8 = {len(basis8)}, dim Q^7 = {sec8.q_dim()}, "
          f"dim B^11 = {len(sec11.basis)}")

    # check the transcribed eq.-(49) monomials form a basis of Q^7
    cols = [sec8.reduce_mod_d(sec8.vec({b: F(1)})) for b in P8B]
    rows = [[cols[c][r] for c in range(len(cols))] for r in range(len(cols[0]))]
    _, piv = rref(rows, len(cols))
    assert len(piv) == len(P8B) == sec8.q_dim(), "eq. (49) transcription bad"
    print("eq. (49) representatives verified as a basis of Q^7")

    NORM = P8B.index(parse_mono("(11111111|)"))

    samples = []
    bad = {F(0), F(-2), F(-2, 3), F(-4, 3)}
    n_values = [F(k) for k in range(1, 16)] + [F(-k) for k in range(1, 10)] + \
               [F(1, 2), F(3, 2), F(5, 2), F(7, 2), F(-1, 2), F(-3, 2),
                F(7, 3), F(-7, 3), F(11, 4), F(2, 5), F(13, 5), F(1, 3),
                F(2, 3) + F(1, 100), F(17, 7), F(-11, 7), F(23, 9)]
    n_values = [n for n in n_values if n not in bad]

    for n in n_values:
        gen = genuine_kernel(paperclip_P4(n), basis8, sec11, sec8)
        if len(gen) != 1:
            print(f"  n = {n}: kernel dimension {len(gen)} != 1, skipping")
            continue
        coeffs = express_in(gen[0], P8B, sec8)
        if coeffs[NORM] == 0:
            print(f"  n = {n}: normalizing coefficient vanishes, skipping")
            continue
        coeffs = [c / coeffs[NORM] for c in coeffs]
        samples.append((n, coeffs))

    print(f"{len(samples)} usable samples; kernel is 1-dimensional at each")

    # ---------------------------------------------------- reconstruction
    nsym = sp.Symbol('n')
    exprs = []
    ok = True
    for i in range(len(P8B)):
        if i == NORM:
            exprs.append(sp.Integer(1))
            continue
        pts = [(n, c[i]) for n, c in samples]
        got = None
        for d in [(4, 4), (6, 6), (8, 8), (10, 10)]:
            need = d[0] + d[1] + 2
            train, held = pts[:need + 2], pts[need + 2:]
            pq = rat_reconstruct(train, *d)
            if pq is None:
                continue
            if all(rat_eval(pq, n) == f for n, f in held):
                got = pq
                break
        if got is None:
            print(f"  coefficient {i} NOT reconstructed")
            ok = False
            exprs.append(None)
            continue
        p, q = got
        num = sum(sp.Rational(c) * nsym**k for k, c in enumerate(p))
        den = sum(sp.Rational(c) * nsym**k for k, c in enumerate(q))
        exprs.append(sp.cancel(sp.factor(num) / sp.factor(den)))

    if ok:
        print("\nAll coefficients reconstructed and verified on held-out samples.")
        print("\nP_8 (normalized: coefficient of (11111111|) = 1):\n")
        lines = []
        for s, e in zip(P8B_NOTATION, exprs):
            lines.append(f"  [{s:>12}]  {sp.factor(sp.cancel(e))}")
            print(lines[-1])
        # common-denominator presentation
        common = sp.lcm([sp.denom(sp.cancel(e)) for e in exprs])
        print(f"\ncommon denominator: {sp.factor(common)}")
        with open("../results/paperclip_P8.json", "w") as f:
            json.dump({
                "normalization": "coefficient of (11111111|) = 1",
                "basis": P8B_NOTATION,
                "coefficients": [str(e) for e in exprs],
                "samples_used": len(samples),
            }, f, indent=1)
        print("\nwritten to results/paperclip_P8.json")
