"""Spin-21 null cells, gate G1 (classical, own code).  At weight W = 22 (spin 21) the classical multiplicity -- the nullity of the
conservation matrix jl_lib.matrix(22, H), built from the PDE H_xr (F_xx - F_rr) - (H_xx - H_rr) F_xr = 0 -- of the committed
I_3's weight-4 symbol H is required to be
  2 at S3: sheet 3 (solution 3) at t = 11/5, (N, s) = (-21/4, -55/4);
  2 at S2: sheet 2 (solution 2) at t = 2,    (N, s) = (-7, -16);
  1 at one control fiber per sheet, the other cell's fiber: sheet 3 at t = 2 and sheet 2 at t = 11/5 (sheet 2's normalization
    pole: its symbol is on the line at infinity).
Printed too (not gate conditions): G2's fiber (sheet 3 at t_5 = 13/11, (N, s) = (-119/2, -143/2)); the lattice indices (j, j')
of S3 and S2 from the closed forms (jl_lib.sheet3_points / sheet2_points, eq. (3) with 3j != s) checked against the symbol.
The symbol is read from cyl_P4 of the byte copy of the repository's code/mixed_engine.py (code/vendor/code, sha256-checked
against the repository at this run): the coefficients of (dX)^4, (dX)^2 (dY.dY), (dY.dY)^2 -> x^4, x^2 r^2, r^4.
Exact arithmetic (Fractions, sympy Rationals).  Exit 0 iff every gate condition holds."""
import os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "vendor", "code"))
from fractions import Fraction as Fr
import sympy as sp

REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))   # the repository root (this file: results/lab/jumping_loci_records/spin21_null_cells/code/)
JL_LIB_SHA256 = "840360c5180447fff81a0e39ee68520d98ce7f98f101a864c2325f62ee3d15a6"   # code/jl_lib.py as shipped (code/SOURCE_SHA256.txt)
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)
def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

print("== byte copies")
for rel in ("code/mixed_engine.py", "code/invariant_engine.py"):
    hv, hr = sha(os.path.join(HERE, "vendor", rel)), sha(os.path.join(REPO, rel))
    require(hv == hr, f"code/vendor/{rel} is a byte copy of the repository's {rel} (sha256 {hr[:16]})")
require(sha(os.path.join(HERE, "jl_lib.py")) == JL_LIB_SHA256, f"code/jl_lib.py is the pinned copy of the shared library (sha256 {JL_LIB_SHA256[:16]})")
if fails:
    sys.exit(1)
import jl_lib as L
import mixed_engine as me
require(os.path.realpath(me.__file__) == os.path.realpath(os.path.join(HERE, "vendor", "code", "mixed_engine.py"))
        and os.path.realpath(L.__file__) == os.path.realpath(os.path.join(HERE, "jl_lib.py")), "the imported modules are these copies")

W = 22; SPIN = W - 1
KEYS = (me._M((1, 1, 1, 1)), me._M((1, 1), [(1, 1)]), me._M((), [(1, 1), (1, 1)]))

def rat(v):
    return sp.Rational(v.numerator, v.denominator)
def symbol_at(sheet, tv):
    N, s = L.N_of(tv), L.s_of(tv)
    P = me.cyl_P4(sheet, N, s)
    return N, s, tuple(P.get(k, Fr(0)) for k in KEYS)
def chart(A, B, C):
    return f"(b, c) = ({B / A}, {C / A})" if A != 0 else f"at infinity, C' = c/b = {C / B}"

CASES = [  # label, sheet, t, required nullity, gate condition?
    ("S3", 3, Fr(11, 5), 2, True),
    ("S2", 2, Fr(2), 2, True),
    ("control: sheet 3 at S2's fiber", 3, Fr(2), 1, True),
    ("control: sheet 2 at S3's fiber", 2, Fr(11, 5), 1, True),
    ("G2's fiber t_5 (not a gate condition)", 3, Fr(13, 11), 2, False),
]
print(f"== the classical multiplicity at W = {W} (spin {SPIN}); the conservation matrix is {W // 2} x {W // 2 + 1}")
notes = []
symb = {}
for label, sheet, tv, want, gate in CASES:
    N, s = L.N_of(tv), L.s_of(tv)
    require(s * s == (N - 25) * (N - 1), f"{label}: t = {tv} gives (N, s) = ({N}, {s}), on the curve s^2 = (N - 25)(N - 1)")
    N, s, (A, B, C) = symbol_at(sheet, tv)
    symb[label] = (A, B, C)
    nul = L.nullity_at(W, rat(A), rat(B), rat(C))
    text = (f"{label}: sheet {sheet} at t = {tv}, tool arguments '{W} {sheet} {N} {s}'; symbol [x^4 : x^2 r^2 : r^4] = "
            f"[{A} : {B} : {C}], {chart(A, B, C)}; nullity {nul} (required {want})")
    if gate:
        require(nul == want, text)
    else:
        print(("NOTE  " if nul == want else "NOTE (differs)  ") + text, flush=True)
        notes.append((label, nul, want))

print(f"== the lattice indices from the closed forms (spin {SPIN})")
aff = L.lattice_affine(SPIN); inf = L.lattice_infinity(SPIN)
s3 = {j: tv for j, tv in L.sheet3_points(SPIN).items() if tv == Fr(11, 5)}
s2 = {jj: tv for jj, tv in L.sheet2_points(SPIN).items() if tv == Fr(2)}
require(list(s3) == [3], f"sheet 3: t = 11/5 is the point j = j' = 3 of sheet3_points({SPIN}) (found {sorted(s3)})")
require(list(s2) == [(6, 3)], f"sheet 2: t = 2 is the point (j, j') = (6, 3) of sheet2_points({SPIN}) (found {sorted(s2)})")
A, B, C = symb["S3"]
require(A != 0 and aff.get((3, 3)) == (B / A, C / A), f"S3's symbol is the eq. (3) lattice point (j, j') = (3, 3): (b, c) = {aff.get((3, 3))}")
A, B, C = symb["S2"]
require(A != 0 and aff.get((6, 3)) == (B / A, C / A), f"S2's symbol is the eq. (3) lattice point (j, j') = (6, 3): (b, c) = {aff.get((6, 3))}")
A, B, C = symb["control: sheet 3 at S2's fiber"]
require(A != 0 and B != 0 and (B / A, C / A) not in set(aff.values()),
        f"control sheet 3 at t = 2: (b, c) = ({B / A}, {C / A}) is off the line b = 0 and not an eq. (3) point at spin {SPIN}")
A, B, C = symb["control: sheet 2 at S3's fiber"]
require(A == 0 and B != 0 and C / B not in set(inf.values()),
        f"control sheet 2 at t = 11/5: at infinity with C' = {C / B}, not among the spin-{SPIN} points at infinity {sorted(set(inf.values()))}")
A, B, C = symb["G2's fiber t_5 (not a gate condition)"]
print(f"NOTE  G2's fiber: (b, c) = ({B / A}, {C / A}); eq. (3) point (j, j') = (5, 5): {aff.get((5, 5))}", flush=True)

print()
if fails:
    print(f"G1: {len(fails)} FAILURE(S): stop and report"); sys.exit(1)
print(f"G1 PASSED: classical multiplicity 2 at S3 and S2, 1 at both control fibers (W = {W}); notes: {notes}")
