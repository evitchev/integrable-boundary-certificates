"""The OTHER half of the four-screening characterization at spins 7, 9.

R4 verifies that the four paperclip screenings annihilate oint P_8 and
oint P_10 (containment in the joint kernel).  The characterization of
the tower (eq. (84) of LVZ) also asserts the joint kernel is no BIGGER:
one dimension per spin.  This certificate computes that joint kernel
exactly at the sampled parameters n = 2, 3 and charge spins 7, 9
(density weights 8, 10): all four screening conditions are imposed
simultaneously (no symmetry shortcut), each as "residue is a twisted
total derivative" via the generic-stratum elimination of
screening_fit.screening_conditions, on the Z2 x Z2-invariant two-boson
density space; solutions are then counted modulo ORDINARY total
derivatives.  Expected: dimension exactly 1, spanned by the known
charge -- membership plus one-dimensionality plus the charge being a
NONZERO class (each checked explicitly; membership alone would not
rule out a total-derivative representative).

Fail-closed: any dimension other than 1, at any of the four (n, sigma)
points, fails the run.  Sampled-n scope: this is exact at n = 2, 3,
not an identity in n -- the same scope as R4.
"""
import json
import sys
from fractions import Fraction as F

import sympy as sp

from compute_p8 import parse_mono
from ope_engine import d_mono, gen_basis
from screening_fit import screening_conditions

failures = []


def _dm(M):
    from sympy.polys.matrices import DomainMatrix
    return DomainMatrix.from_Matrix(M, extension=True)


def dm_rank(M):
    return _dm(M).rank() if M.cols and M.rows else 0


def dm_nullspace(M):
    """Nullspace column vectors of M via DomainMatrix over the algebraic
    extension field of the entries (orders of magnitude faster than
    Matrix.nullspace on algebraic entries)."""
    if not M.cols:
        return []
    ns = _dm(M).nullspace().to_Matrix()   # rows span the nullspace
    return [ns[i, :].T for i in range(ns.rows)]


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def alpha_vectors(nval):
    a = sp.sqrt(sp.Rational(-nval, 2))
    b = sp.sqrt(sp.Rational(nval + 2, 2))
    return [(sa * a, sb * b) for sa in (1, -1) for sb in (1, -1)]


def load_charge(path, nval):
    with open(path) as f:
        data = json.load(f)
    basis = data.get("basis", data.get("basis_monos"))
    if basis is None:                       # fail closed on schema drift
        raise SystemExit(f"no basis/basis_monos key in {path}")
    if len(basis) != len(data["coefficients"]):   # zip would truncate
        raise SystemExit(f"basis/coefficients length mismatch in {path}")
    nsym = sp.Symbol('n')
    out = {}
    for s, e in zip(basis, data["coefficients"]):
        mono = parse_mono(s) if isinstance(s, str) else tuple(
            tuple(l) for l in s)
        v = sp.Rational(sp.sympify(e, locals={'n': nsym})
                        .subs(nsym, sp.Rational(nval)))
        if v != 0:
            out[mono] = F(int(sp.numer(v)), int(sp.denom(v)))
    return out


def genuine_joint_kernel_dim(nval, w):
    """dim of (joint kernel of all four screenings) / (total derivatives)
    on the Z2 x Z2-invariant two-boson density space of weight w.

    Per screening: build the residue map Res_alpha and the twisted
    derivative D_alpha on the full residue target space; the screening
    condition is "Res_alpha(P) lies in Im D_alpha", i.e. the projection
    of Res_alpha(P) to coker D_alpha vanishes.  The cokernel projection
    is computed by row reduction; the four small projected condition
    blocks are stacked and the joint nullspace taken -- no witness
    unknowns, so the matrices stay small."""
    from screening_fit import residue_screen, twisted_d
    basis = gen_basis(w, 2, z2=True)
    nP = len(basis)
    s = w - 1
    tgt = [()]
    for tt in range(1, s + 1):
        tgt += gen_basis(tt, 2, z2=False)
    tix = {m: i for i, m in enumerate(tgt)}
    lower = [()] if s - 1 == 0 else gen_basis(s - 1, 2, z2=False)
    blocks = []
    for al in alpha_vectors(nval):
        Dm = sp.zeros(len(tgt), len(lower))
        for k, u in enumerate(lower):
            for mm, c in twisted_d(u, al, 2).items():
                Dm[tix[mm], k] = c
        R = sp.zeros(len(tgt), nP)
        for j, m in enumerate(basis):
            for mm, c in residue_screen({m: F(1)}, al).items():
                R[tix[mm], j] += c
        # cokernel projection: the left-null space of D_alpha (computed
        # as the nullspace of D^T) gives the functionals vanishing on
        # Im D_alpha; applying them to Res_alpha is the condition block
        ln = dm_nullspace(Dm.T)
        if ln:
            Pc = sp.Matrix.vstack(*[v.T for v in ln])
            blocks.append(sp.expand(Pc * R))
    A = sp.Matrix.vstack(*blocks) if blocks else sp.zeros(1, nP)
    ker = dm_nullspace(A)
    # ordinary total derivatives inside the weight-w invariant space
    idx = {m: i for i, m in enumerate(basis)}
    dcols = []
    for u in gen_basis(w - 1, 2, z2=True):
        col = [sp.Integer(0)] * len(basis)
        img = d_mono(u)
        if any(m not in idx for m in img):
            continue                      # leaves the invariant sector
        for m, c in img.items():
            col[idx[m]] = sp.Integer(c)
        dcols.append(col)
    D = sp.Matrix.zeros(len(basis), 0) if not dcols else sp.Matrix(
        [list(r) for r in zip(*dcols)])
    K = sp.Matrix.hstack(*ker) if ker else sp.zeros(len(basis), 0)
    rD = dm_rank(D)
    rKD = dm_rank(sp.Matrix.hstack(K, D)) if ker else rD
    return rKD - rD, basis, ker


def joint_kernel_dim_witness_block(nval, w):
    """CONTROL implementation: the original explicit-witness construction
    (unknowns = P coordinates + one twisted witness per screening; one
    joint nullspace).  Used only at low weight to validate the cokernel-
    projection method against an independent construction."""
    from screening_fit import residue_screen, twisted_d
    basis = gen_basis(w, 2, z2=True)
    nP = len(basis)
    alphas = alpha_vectors(nval)
    s = w - 1
    tgt = [()]
    for tt in range(1, s + 1):
        tgt += gen_basis(tt, 2, z2=False)
    tix = {m: i for i, m in enumerate(tgt)}
    lower = [()] if s - 1 == 0 else gen_basis(s - 1, 2, z2=False)
    nW = len(lower)
    rows_per = len(tgt)
    A = sp.zeros(4 * rows_per, nP + 4 * nW)
    for ai, al in enumerate(alphas):
        r0 = ai * rows_per
        for j, m in enumerate(basis):
            for mm, c in residue_screen({m: F(1)}, al).items():
                A[r0 + tix[mm], j] += c
        for k, u in enumerate(lower):
            for mm, c in twisted_d(u, al, 2).items():
                A[r0 + tix[mm], nP + ai * nW + k] -= c
    kf = dm_nullspace(A)
    Pcols = [v[:nP, 0] for v in kf]
    K = sp.Matrix.hstack(*Pcols) if Pcols else sp.zeros(nP, 0)
    idx = {m: i for i, m in enumerate(basis)}
    dcols = []
    for u in gen_basis(w - 1, 2, z2=True):
        col = [sp.Integer(0)] * len(basis)
        img = d_mono(u)
        if any(m not in idx for m in img):
            continue
        for m, c in img.items():
            col[idx[m]] = sp.Integer(c)
        dcols.append(col)
    D = sp.Matrix([list(r) for r in zip(*dcols)]) if dcols         else sp.zeros(nP, 0)
    rD = dm_rank(D)
    rKD = dm_rank(sp.Matrix.hstack(K, D)) if Pcols else rD
    return rKD - rD


from pathlib import Path as _P
from launch_provenance import capture_launch
_ROOT = _P(__file__).resolve().parent.parent
_LAUNCH = capture_launch(_ROOT, _P(__file__).resolve())

# METHOD CONTROL (review requirement): at the small weight 6 (charge
# spin 5) the cokernel-projection method must agree with the original
# witness-block construction, at both sampled n.
RECORDS = {"controls": [], "points": []}
for nval in (F(2), F(3)):
    d_fast, _, _ = genuine_joint_kernel_dim(nval, 6)
    d_ctrl = joint_kernel_dim_witness_block(nval, 6)
    check(f"method control at n={nval}, weight 6: cokernel-projection "
          f"dim {d_fast} == witness-block dim {d_ctrl}", d_fast == d_ctrl)
    RECORDS["controls"].append({"n": str(nval), "weight": 6,
                                "cokernel_dim": d_fast,
                                "witness_block_dim": d_ctrl})

for nval, w, path in [(F(2), 8, "../results/paperclip_P8.json"),
                      (F(3), 8, "../results/paperclip_P8.json"),
                      (F(2), 10, "../results/paperclip_P10.json"),
                      (F(3), 10, "../results/paperclip_P10.json")]:
    dim, basis, ker = genuine_joint_kernel_dim(nval, w)
    check(f"joint kernel of ALL FOUR screenings at n={nval}, weight {w} "
          f"(charge spin {w-1}) has class dimension {dim} == 1 -- the "
          "characterization's uniqueness half at this sample", dim == 1)
    # membership sanity: the known charge solves every condition
    P = load_charge(path, nval)
    bad = [al for al in alpha_vectors(nval)
           if any(sp.simplify(c) != 0
                  for c in screening_conditions(P, al, 2))]
    check(f"the known spin-{w-1} charge at n={nval} lies in the joint "
          "kernel", not bad)
    # membership alone would also pass for a total-derivative class;
    # spanning needs the charge to be a NONZERO class (review, 2026-08-26)
    basis_w = gen_basis(w, 2, z2=True)
    vec = sp.Matrix([sp.Rational(P.get(m, 0)) for m in basis_w])
    idx = {m: i for i, m in enumerate(basis_w)}
    dcols = []
    for u in gen_basis(w - 1, 2, z2=True):
        col = [sp.Integer(0)] * len(basis_w)
        img = d_mono(u)
        if any(m not in idx for m in img):
            continue
        for m, c in img.items():
            col[idx[m]] = sp.Integer(c)
        dcols.append(col)
    Dw = sp.Matrix([list(r) for r in zip(*dcols)])
    nz = dm_rank(sp.Matrix.hstack(Dw, vec)) == dm_rank(Dw) + 1
    check(f"the known spin-{w-1} charge at n={nval} is a NONZERO class "
          "(not a total derivative) -- membership + dim 1 => it spans",
          nz)
    RECORDS["points"].append({"n": str(nval), "weight": w,
                              "charge_spin": w - 1, "class_dim": dim,
                              "known_charge_in_kernel": not bad,
                              "known_charge_nonzero_class": bool(nz)})

print()
if failures:
    print(f"PAPERCLIP JOINT-KERNEL TEST FAILED ({len(failures)})")
    sys.exit(1)
from pathlib import Path
from certificate import stamp
RECORDS["launch"] = _LAUNCH
RECORDS["script_sha256"] = _LAUNCH["script_sha256_at_start"]
RECORDS["object"] = ("joint kernel of the four paperclip screenings on "
                     "the Z2xZ2-invariant sector, modulo total "
                     "derivatives, at sampled n = 2, 3, spins 7 and 9")
RECORDS["acceptance"] = ("class dimension exactly 1 at each point, known "
                         "charge a member AND a nonzero class (hence "
                         "spans); method control: "
                         "cokernel-projection == witness-block at weight 6")
art = Path(__file__).resolve().parent.parent / "results" / \
    "paperclip_joint_kernel_79.json"


def _canon(d):
    d = {k: v for k, v in d.items()
         if k not in ("provenance", "launch", "script_sha256",
                      "invocation")}   # provenance-class keys: a
    # provenance-only refactor must not fake a scientific change
    # (review, 2026-08-27)
    return json.dumps(json.loads(json.dumps(d)), sort_keys=True)


if art.exists() and _canon(json.loads(art.read_text())) == _canon(RECORDS):
    print("artifact of record unchanged (payload agrees); "
          "provenance preserved")
else:
    art.write_text(json.dumps(stamp(RECORDS), indent=1) + "\n")
    print("artifact written")
print("PAPERCLIP JOINT KERNEL CERTIFIED (all four screening conditions "
      "imposed simultaneously: class dimension exactly 1 at n=2,3, "
      "spins 7 and 9, spanned by the known charges -- the uniqueness "
      "half of the four-screening characterization at the sampled n)")
