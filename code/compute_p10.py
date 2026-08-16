"""Compute the spin-9 paperclip charge density P_10.

Same method as compute_p8.py one spin up: at ~45 sample rational n,
solve [oint P_4(n), oint P_10] = 0 in the 55-dimensional Q^9 of
Z2xZ2-symmetric weight-10 densities mod d (sparse kernels per
sample), check multiplicity one, normalize by (dX_1)^10, and
reconstruct each coefficient as a rational function of n with
held-out verification.  Vitchev prints no Q^9 basis, so
representatives are auto-selected greedily (normalization monomial
first) and recorded in the artifact as raw monomial tuples.

Output: results/paperclip_P10.{md,json}; fail-closed."""
import json
import sys
from fractions import Fraction as F

import sympy as sp

from compute_p8 import parse_mono
from families import paperclip_P4
from ope_engine import Sector, d_mono, fmt_mono, gen_basis, residue_cur, spin
from solve_charges import express_in, rat_reconstruct, rat_eval
from sparse_linalg import SparsePivots, genuine_kernel_sparse, lcm

if __name__ == "__main__":
    sec10 = Sector(10, 2, z2=True)
    basis10 = gen_basis(10, 2, z2=True)
    qdim = sec10.q_dim()
    print(f"dim B^10 = {len(basis10)}, dim Q^9 = {qdim}", flush=True)

    NORM_MONO = parse_mono("(1111111111|)")

    def int_vec(v):
        L = 1
        for c in v.values():
            L = lcm(L, c.denominator)
        return {k: int(c * L) for k, c in v.items()}

    piv = SparsePivots()
    reps = []
    for b in [NORM_MONO] + [b for b in basis10 if b != NORM_MONO]:
        if len(reps) == qdim:
            break
        v = sec10.reduce_mod_d(sec10.vec({b: F(1)}))
        sv = {i: c for i, c in enumerate(v) if c}
        if sv and piv.insert(int_vec(sv)):
            reps.append(b)
    assert len(reps) == qdim and reps[0] == NORM_MONO
    print(f"auto-selected {len(reps)} Q^9 representatives", flush=True)

    bad = {F(0), F(-2), F(-2, 3), F(-4, 3)}
    n_values = [F(k) for k in range(1, 16)] + \
               [F(-k) for k in range(1, 10)] + \
               [F(1, 2), F(3, 2), F(5, 2), F(7, 2), F(-1, 2), F(-3, 2),
                F(7, 3), F(-7, 3), F(11, 4), F(2, 5), F(13, 5), F(1, 3),
                F(2, 3) + F(1, 100), F(17, 7), F(-11, 7), F(23, 9),
                F(19, 3), F(-13, 4), F(29, 8), F(3, 7), F(31, 9),
                F(-17, 5), F(37, 10), F(5, 7), F(41, 11), F(-2, 7)]
    n_values = [n for n in n_values if n not in bad]

    samples = []
    for n in n_values:
        gen = genuine_kernel_sparse(
            paperclip_P4(n), basis10, None,
            lambda P, Q, N: residue_cur(P, Q),
            lambda s: gen_basis(s, 2, True),
            d_mono, spin)
        if len(gen) != 1:
            print(f"  n = {n}: kernel dimension {len(gen)} != 1, "
                  f"skipping", flush=True)
            continue
        coeffs = express_in(gen[0], reps, sec10)
        if coeffs[0] == 0:
            print(f"  n = {n}: normalizing coefficient vanishes, "
                  f"skipping", flush=True)
            continue
        coeffs = [c / coeffs[0] for c in coeffs]
        samples.append((n, coeffs))
        if len(samples) % 10 == 0:
            print(f"  {len(samples)} samples done", flush=True)

    print(f"{len(samples)} usable samples; kernel 1-dimensional at "
          f"each", flush=True)

    nsym = sp.Symbol('n')
    exprs = []
    ok = True
    for i in range(qdim):
        if i == 0:
            exprs.append(sp.Integer(1))
            continue
        pts = [(n, c[i]) for n, c in samples]
        got = None
        for d in [(4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]:
            need = d[0] + d[1] + 2
            train, held = pts[:need + 2], pts[need + 2:]
            if len(held) < 3:
                continue
            pq = rat_reconstruct(train, *d)
            if pq is None:
                continue
            if all(rat_eval(pq, n) == f for n, f in held):
                got = pq
                break
        if got is None:
            print(f"  coefficient {i} NOT reconstructed", flush=True)
            ok = False
            exprs.append(None)
            continue
        p, q = got
        num = sum(sp.Rational(c) * nsym**k for k, c in enumerate(p))
        den = sum(sp.Rational(c) * nsym**k for k, c in enumerate(q))
        exprs.append(sp.cancel(sp.factor(num) / sp.factor(den)))

    if not ok:
        print("P10 RECONSTRUCTION FAILED", flush=True)
        sys.exit(1)

    print("\nAll coefficients reconstructed and verified on held-out "
          "samples.", flush=True)
    common = sp.lcm([sp.denom(sp.cancel(e)) for e in exprs])
    print(f"common denominator: {sp.factor(common)}", flush=True)

    with open("../results/paperclip_P10.json", "w") as f:
        json.dump({
            "normalization": "coefficient of (dX_1)^10 = 1",
            "basis_monos": [[list(l) for l in m] for m in reps],
            "basis_pretty": [fmt_mono(m, 2) for m in reps],
            "coefficients": [str(e) for e in exprs],
            "n_samples": len(samples),
        }, f, indent=1)
    with open("../results/paperclip_P10.md", "w") as f:
        f.write("# Paperclip P_10 (spin-9 charge density)\n\n"
                "Computed by `code/compute_p10.py` (see docstring); "
                "normalized so (dX_1)^10 has coefficient 1.\n\n"
                f"{len(samples)} samples, multiplicity one at each; "
                "all 55 coefficients rational-reconstructed with "
                "held-out verification.\n\n"
                f"Common denominator: `{sp.factor(common)}`\n\n")
        for m, e in zip(reps, exprs):
            f.write(f"- `{fmt_mono(m, 2)}`: "
                    f"`{sp.factor(sp.cancel(e))}`\n")
    print("P10 COMPLETE", flush=True)
