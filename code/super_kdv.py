"""Phase 2 calibration: quantum super-KdV (one N=1 superfield, c = 3/2).

AUDIT-CORRECTED (2026-08-04, external audit finding 1): the composite
pencil is built with the POINT-SPLIT normal product
:AB:(w) = oint dz (z-w)^-1 A(z)B(w), not naive concatenation.  The two
differ by non-d-exact correction terms, and only the point-split
pencil contains the standard quantum super-KdV ray
    I3 = :TT: - 1/4 :G dG:
of Kulish-Zeitlin (SUSY N=1 KdV hierarchy I & II: hep-th/0407154,
hep-th/0501019) and Chistyakova-Litvinov-Orlov (arXiv:2110.05870,
where I3 and the fermionic-screening commutant characterization
appear explicitly; they cite the two KZ papers for the hierarchy).

Selection of the ray inside the pencil :TT: + kappa :G dG::
  - ad_G on charges is a differential (ad_G^2 = ad_{oint T} = 0), so
    its kernel is huge (G-cohomology) and cannot select;
  - at c = 3/2 the dual screening pair psi e^{+-X} (CLO conventions:
    c-hat = 1 forces Q = b + 1/b = 0, b = -1/b = +-i, alpha = -ib =
    +-1) annihilates EVERY member of the point-split pencil, so
    screenability does not select either (super_screening_fit.py);
  - G-invariance of the density ray DOES select: Res(G, :TT: +
    kappa :G dG:) is a total derivative iff kappa = -1/4.

All spectral statements are scoped to the computed range sigma <= MAX
(default 18, i.e. sigma <= 8; pass 24 for sigma <= 11).  No all-spin
theorem is claimed.

--naive reproduces the legacy concatenation-pencil family: a DIFFERENT
exactly-commuting tower (its brackets vanish exactly; it does not
commute with the point-split I3), conjecturally a background-charged
member of the super-KdV family — see results/super_kdv_calibration.md
section 4.  Its artifacts carry a _naive suffix.

Outputs results/super_kdv_calibration[_naive][_ext].{json}; exits
nonzero (and writes .partial) on any anomaly.  The success JSON is
deterministic (no volatile fields)."""
import json
import resource
import sys
import time
from fractions import Fraction as F

from ope_engine import nullspace
from sparse_linalg import lcm
from super_engine import (SuperSector, canonicalize, cur_add, d_cur,
                          d_mono, gen_basis_super, genuine_kernel_super,
                          mono_weight2, normal_prod, q_dim_super,
                          residue_cur)

lim = 12 * 1024**3
resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

NF = 1
DX, PSI = (0, 1, 0), (1, 0, 0)
T = {(DX, DX): F(1, 2), (PSI, (1, 1, 0)): F(-1, 2)}
G = {(DX, PSI): F(1)}


def sig(w2):
    """Charge-spin label for density doubled-weight w2."""
    return str((w2 - 2) // 2) if w2 % 2 == 0 else f"{w2 - 2}/2"


def fmt_letter(l):
    p, m, j = l
    core = ("psi" if p else "X") + (f"_{j}" if NF > 1 else "")
    return core if (p and m == 0) else (f"d{core}" if m == 1
                                        else f"d{m}{core}")


def fmt_cur(cur):
    terms = []
    for mono in sorted(cur):
        c = cur[mono]
        terms.append(f"{'+' if c > 0 else '-'} {abs(c)} "
                     + " ".join(fmt_letter(l) for l in mono))
    return " ".join(terms)


def primitive(cur):
    """Scale to coprime integers, first (sorted) coefficient positive."""
    from math import gcd
    L = 1
    for c in cur.values():
        L = lcm(L, c.denominator)
    g = 0
    for c in cur.values():
        g = gcd(g, abs(int(c * L)))
    s = 1 if cur[sorted(cur)[0]] > 0 else -1
    return {m: F(int(c * L) * s, g) for m, c in cur.items()}


def naive_mult(A, B):
    """Concatenation (contraction-free) product — NOT the quantum
    normal-ordered product; kept for the --naive legacy family."""
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            cm, sg = canonicalize(ma + mb)
            if cm is not None:
                v = out.get(cm, 0) + sg * ca * cb
                if v:
                    out[cm] = v
                else:
                    out.pop(cm, None)
    return out


def select_I3(prod):
    """The unique G-invariant ray TTc + kappa GdGc in the pencil built
    with `prod`; returns (kappa, I3 current).  Raises on anomaly."""
    TTc, GdGc = prod(T, T), prod(G, d_cur(G))
    w9 = gen_basis_super(9, NF)
    idx9 = {m: i for i, m in enumerate(w9)}
    cols = []
    for R in (residue_cur(G, TTc), residue_cur(G, GdGc)):
        v = [F(0)] * len(w9)
        for m, c in R.items():
            v[idx9[m]] = c
        cols.append(v)
    for b in gen_basis_super(7, NF):
        v = [F(0)] * len(w9)
        for m, c in d_mono(b).items():
            v[idx9[m]] = F(c)
        cols.append(v)
    sols = [v for v in nullspace(cols, len(w9)) if v[0] or v[1]]
    if len(sols) != 1 or not sols[0][0]:
        raise RuntimeError(f"G-invariant ray not unique/regular: "
                           f"{len(sols)} solutions")
    kappa = sols[0][1] / sols[0][0]
    I3 = {}
    cur_add(I3, TTc)
    cur_add(I3, GdGc, kappa)
    if SuperSector(8, NF).is_total_derivative(I3):
        raise RuntimeError("I3 candidate is d-exact")
    return kappa, I3


def derive_tower(prod=normal_prod, weights=(12, 16)):
    """{8: I3, w2: I_sigma ...} via the joint (ad_G, ad_I3) kernel.
    Importable (used by super_screening_fit.py)."""
    _, I3 = select_I3(prod)
    tower = {8: I3}
    for w2 in weights:
        greps = genuine_kernel_super(
            G, [{b: F(1)} for b in gen_basis_super(w2, NF)], NF)
        joint = genuine_kernel_super(I3, greps, NF)
        if len(joint) != 1:
            raise RuntimeError(f"joint kernel at w2={w2}: dim "
                               f"{len(joint)} != 1")
        tower[w2] = joint[0]
    return tower


def scan(P, w2s, cand_map=None):
    """dims and genuine representatives of ker ad_P per doubled weight.
    cand_map: optional {w2: [currents]} restricting candidates (joint
    kernels); default = full monomial basis."""
    dims, reps = {}, {}
    for w2 in w2s:
        cands = (cand_map.get(w2, []) if cand_map is not None
                 else [{b: F(1)} for b in gen_basis_super(w2, NF)])
        g = genuine_kernel_super(P, cands, NF) if cands else []
        dims[w2], reps[w2] = len(g), g
        print(f"  w2={w2:2d} (sigma={sig(w2):>4}): dim {len(g)}",
              flush=True)
    return dims, reps


def check_bracket(name, A, B, results):
    S = residue_cur(A, B)
    ok = (not S) or SuperSector(
        mono_weight2(next(iter(S))), NF).is_total_derivative(S)
    print(f"  [{name}] total derivative: {ok}", flush=True)
    results[name] = ok
    return ok


def main():
    t0 = time.time()
    args = [a for a in sys.argv[1:] if a != "--naive"]
    naive = "--naive" in sys.argv[1:]
    MAX = int(args[0]) if args else 18
    prod = naive_mult if naive else normal_prod
    w2s = list(range(3, MAX + 1))
    out = {"max_w2": MAX,
           "product": "naive-concatenation" if naive else "point-split",
           "qdim": {}, "adG": {}, "adI3": {}, "joint": {},
           "checks": {}, "anomalies": []}

    for w2 in w2s:
        out["qdim"][w2] = q_dim_super(w2, NF)

    print("ad_G spectrum (G-cohomological: kernel contains G-exact "
          "charges):", flush=True)
    gdims, greps = scan(G, w2s)
    out["adG"] = gdims

    try:
        kappa, I3 = select_I3(prod)
    except RuntimeError as e:
        out["anomalies"].append(str(e))
        finish(out, t0, naive, failed=True)
    out["kappa"] = str(kappa)
    I3 = primitive(I3)
    out["I3"] = fmt_cur(I3)
    print(f"\nI3 = (TT) + {kappa}*(G dG) [{out['product']} pencil, "
          f"selected by G-invariance]\nI3 density = {out['I3']}\n",
          flush=True)

    print("ad_I3 spectrum (full free-field centralizer, computed "
          "range only):", flush=True)
    i3dims, i3reps = scan(I3, w2s)
    out["adI3"] = i3dims

    print("joint (ad_G, ad_I3) kernel:", flush=True)
    jdims, jreps = scan(I3, w2s, cand_map=greps)
    out["joint"] = jdims

    for w2 in (12, 16):
        if jdims.get(w2) != 1:
            out["anomalies"].append(f"joint kernel at w2={w2} has dim "
                                    f"{jdims.get(w2)} (expected 1)")
    if out["anomalies"]:
        finish(out, t0, naive, failed=True)
    I5, I7 = primitive(jreps[12][0]), primitive(jreps[16][0])
    out["I5"], out["I7"] = fmt_cur(I5), fmt_cur(I7)

    print("\ndirect bracket verifications:", flush=True)
    ok = True
    for name, A, B in (("I5,I7", I5, I7), ("G,I5", G, I5),
                       ("G,I7", G, I7), ("I3,I5", I3, I5),
                       ("I3,I7", I3, I7)):
        ok = check_bracket(name, A, B, out["checks"]) and ok
    if jdims.get(20) == 1:
        I9 = primitive(jreps[20][0])
        out["I9"] = fmt_cur(I9)
        for name, A, B in (("I5,I9", I5, I9), ("I7,I9", I7, I9)):
            ok = check_bracket(name, A, B, out["checks"]) and ok
    if not ok:
        out["anomalies"].append("a direct bracket check failed")

    finish(out, t0, naive, failed=bool(out["anomalies"]))


def finish(out, t0, naive, failed):
    # elapsed goes to stdout only: the success JSON is tracked and must
    # stay deterministic (no volatile fields; cf. the manifest rule)
    elapsed = round(time.time() - t0, 1)
    tag = ("_naive" if naive else "") + ("_ext" if out["max_w2"] > 18
                                         else "")
    path = (f"../results/super_kdv_calibration{tag}.json" if not failed
            else f"../results/super_kdv_calibration{tag}.partial.json")
    json.dump(out, open(path, "w"), indent=1)
    print(f"\n{'ANOMALIES: ' + '; '.join(out['anomalies']) if failed else 'SUPER-KDV CALIBRATION COMPLETE'}"
          f"  ({elapsed}s)", flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
