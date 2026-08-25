"""Lab opener, Python side: the commutator matrix of ad_{I_3} for the
cylindrical families as an EXACT polynomial matrix in the curve
parameter t, exported for the Mathematica stratification script.

ad_{I_3}: Q^w (charge classes at weight w) -> Q^{w+3} (residue classes).
Both sides are taken in class coordinates (non-pivot monomials of the
mixed sectors modulo total derivatives), so the nullity of the matrix
IS the charge-class kernel dimension (multiplicity) at that spin.

Linearity: I_3 = sum_k c_k(N,s) mono_k, so M(N,s) = sum_k c_k(N,s) M_k(N),
where M_k is the residue map of the single monomial mono_k and depends
on N only through O(N) traces -- polynomially.  Each M_k(N) is built by
exact evaluation at NSAMP rational N and Lagrange interpolation, then
VERIFIED at held-out N (fail-closed).  On the curve N = (t^2-25)/(t^2-1),
s = -24t/(t^2-1); the common denominator is cleared so every entry is a
polynomial in t with integer coefficients.

Output: wolfram/lab/cyl_strata_w{w}_sol{sol}.m  (Mathematica list of
rows of polynomials in t) plus a .json sidecar with sizes, the
generic-rank check points, and the script hash.
Usage: python3 cyl_strata_build.py W [SOL ...]
"""
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "code"))
from mixed_engine import (MixedSector, _M, gen_mixed_basis,           # noqa: E402
                          residue_mixed_cur)

MONOS = {                       # I_3 monomials and their coefficients c_k(N, s)
    "X4":    _M((1, 1, 1, 1)),
    "X22":   _M((2, 2)),
    "cross": _M((1, 1), [(1, 1)]),
    "Y11":   _M((), [(1, 1), (1, 1)]),
    "Y22":   _M((), [(2, 2)]),
}
Nsym, ssym, t = sp.symbols('N s t')
COEFF = {
    1: {"X4": (2 + Nsym + ssym) / 3, "X22": (-9 - 3 * Nsym - ssym) / 6,
        "cross": sp.Integer(6), "Y11": sp.Integer(1), "Y22": sp.Integer(-2)},
    2: {"X4": (117 + 2 * Nsym + Nsym**2 + (15 + Nsym) * ssym) / 192,
        "X22": -(409 + 46 * Nsym + Nsym**2 + (Nsym - 5) * ssym) / 192,
        "cross": (11 + Nsym + ssym) / 8, "Y11": sp.Integer(1), "Y22": sp.Integer(-2)},
    3: {"X4": sp.Integer(1),
        "X22": (13 - 27 * Nsym + 2 * Nsym**2 + (5 - 2 * Nsym) * ssym) / 6,
        "cross": 7 - Nsym + ssym, "Y11": sp.Integer(1), "Y22": sp.Integer(-2)},
}
NSAMP = [F(k) for k in (2, 3, 5, 7, 11, 13, 17)]      # interpolation nodes
NCHECK = [F(19), F(-4), F(1, 2)]                       # held-out verification
# RIGOROUS DEGREE BOUND (audit 2026-08-23): the only N-dependence of a residue
# entry comes from full traces delta_ii = N, each produced by contracting both
# letters of one Y-pair of the I_3 monomial with both letters of one Y-pair of
# the candidate; a Y-pair of the monomial enters at most one trace, and there
# are no self-contractions.  Hence deg_N <= (# Y-pairs of the monomial) <= 2,
# and len(NSAMP) = 7 nodes determine a degree-<=6 polynomial EXACTLY.  The
# held-out checks are then a consistency control, not the proof.
YPAIRS = {"X4": 0, "X22": 0, "cross": 1, "Y11": 2, "Y22": 1}
if len(NSAMP) < max(YPAIRS.values()) + 1:
    sys.exit(f"FAIL  interpolation nodes {len(NSAMP)} < degree bound + 1")


def class_coords(sector):
    piv = set(sector._img_pivots)
    return [i for i in range(len(sector.basis)) if i not in piv]


def lagrange(nodes, values, x):
    """Exact Lagrange interpolation polynomial (sympy) through (nodes, values)."""
    P = sp.Integer(0)
    for i, (xi, yi) in enumerate(zip(nodes, values)):
        if yi == 0:
            continue
        term = sp.Rational(yi.numerator, yi.denominator)
        for j, xj in enumerate(nodes):
            if j != i:
                term *= (x - sp.Rational(xj.numerator, xj.denominator)) / \
                        sp.Rational((xi - xj).numerator, (xi - xj).denominator)
        P += term
    return sp.expand(P)


def build(w, sol):
    S, T = MixedSector(w), MixedSector(w + 3)
    cols_idx, rows_idx = class_coords(S), class_coords(T)
    cand = [S.basis[i] for i in cols_idx]
    nrow, ncol = len(rows_idx), len(cols_idx)

    def column(mono, b, Nval):
        v = T.reduce_mod_d(T.vec(residue_mixed_cur({mono: F(1)}, {b: F(1)}, Nval)))
        return [v[i] for i in rows_idx]

    Mk = {}
    for name, mono in MONOS.items():
        samples = {Nv: [column(mono, b, Nv) for b in cand] for Nv in NSAMP + NCHECK}
        mat = [[None] * ncol for _ in range(nrow)]
        for j in range(ncol):
            for i in range(nrow):
                vals = [samples[Nv][j][i] for Nv in NSAMP]
                poly = lagrange(NSAMP, vals, Nsym)
                for Nv in NCHECK:                      # fail-closed verification
                    got = poly.subs(Nsym, sp.Rational(Nv.numerator, Nv.denominator))
                    want = samples[Nv][j][i]
                    if sp.Rational(got) != sp.Rational(want.numerator, want.denominator):
                        raise SystemExit(f"interpolation degree too low: {name} "
                                         f"entry ({i},{j}) at N={Nv}")
                mat[i][j] = poly
        Mk[name] = mat
        degN = max((sp.Poly(e, Nsym).degree() for r in mat for e in r if e != 0),
                   default=0)
        if degN > YPAIRS[name]:
            raise SystemExit(f"degree bound violated: M_{name} has degree {degN} "
                             f"> {YPAIRS[name]} Y-pairs")
        print(f"  M_{name}: max degree in N = {degN} (bound {YPAIRS[name]}, "
              f"nodes {len(NSAMP)})", flush=True)
    # assemble on the curve
    Nt = (t**2 - 25) / (t**2 - 1)
    st = -24 * t / (t**2 - 1)
    M = [[sp.Integer(0)] * ncol for _ in range(nrow)]
    for name, mat in Mk.items():
        c = COEFF[sol][name].subs({Nsym: Nt, ssym: st})
        for i in range(nrow):
            for j in range(ncol):
                if mat[i][j] != 0:
                    M[i][j] += c * mat[i][j].subs(Nsym, Nt)
    # clear the common denominator -> integer polynomial matrix
    den = sp.Integer(1)
    for r in M:
        for e in r:
            if e != 0:
                den = sp.lcm(den, sp.denom(sp.together(e)))
    Mpoly = [[sp.expand(sp.cancel(e * den)) for e in r] for r in M]
    for r in Mpoly:
        for e in r:
            if sp.denom(sp.together(e)) != 1:
                sys.exit(f"FAIL  non-polynomial interpolant {e}")
    return Mpoly, cand, rows_idx, S.q_dim(), T.q_dim()


def main():
    w = int(sys.argv[1])
    sols = [int(x) for x in sys.argv[2:]] or [1, 2, 3]
    for sol in sols:
        print(f"weight {w}, Solution {sol}:", flush=True)
        M, cand, rows_idx, qw, qt = build(w, sol)
        out = HERE / f"cyl_strata_w{w}_sol{sol}.m"
        rows = ["{" + ", ".join(sp.mathematica_code(e) for e in r) + "}" for r in M]
        out.write_text("{\n" + ",\n".join(rows) + "\n}\n")
        meta = {"weight": w, "solution": sol, "rows": len(M), "cols": len(M[0]),
                "class_dim_w": qw, "class_dim_target": qt,
                "curve": "N=(t^2-25)/(t^2-1), s=-24t/(t^2-1)",
                "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "matrix_sha256": hashlib.sha256(out.read_bytes()).hexdigest()}
        (HERE / f"cyl_strata_w{w}_sol{sol}.json").write_text(json.dumps(meta, indent=1))
        print(f"  -> {out.name}: {len(M)}x{len(M[0])}", flush=True)


if __name__ == "__main__":
    main()
