"""Exact stratification over Q-bar at weights 10, 12 and 14 (spins 9, 11, 13) -- fail-closed
certificate with open tools only (Python, sympy, numpy; FLINT through a small C driver).

Tests the preregistered expectations of results/cyl_discriminant_preregistration.md (H_F1-H_F3)
and results/cyl_discriminant_w14_preregistration.md (H_W14), which the certificates
cyl_discriminant_w1012 and cyl_discriminant_w14 test with wolframscript; this one needs no
Mathematica.  For each weight W and each sheet 1-3:

  1. The family matrix M_poly(t) = (t^2-1)^k c0 M(t) of ad I_3 in class coordinates is REBUILT
     here from the engine (strat1014_family: mixed_engine's cyl_P4 evaluated over Q(t),
     residue_mixed_mono, MixedSector).  The reduction shortcut is compared with the engine's
     reduce_mod_d (every domain column at W <= 12, the first 40 at W = 14), and the sha256 of
     each matrix file must equal the pinned value (strat1014_verify.EXPECTED).
  2. The archived claims in results/cyl_stratification_w1014/claims/ (every file listed in its
     SHA256SUMS, whose own sha256 is pinned below) are re-verified by strat1014_verify.main:
     the generic rank r, proved both ways (rank modulo a prime at a point; exact null vectors);
     two maximal minors P1, P2 (row and column content and powers of t -+ 1 removed),
     recomputed modulo every listed 62-bit prime by a FLINT driver compiled here
     (strat1014_detmodp.c) and matched with the archived exact polynomials under proven degree
     and Hadamard bounds; G = gcd(P1o, P2o), checked exactly (divisibility, cofactors coprime
     modulo a certificate prime, irreducible factors); every rank drop away from t = +-1 is a
     root of G, and every factor of G is decided -- no drop (rank at a root modulo a prime), or
     a drop at a rational root with exact null vectors.
  3. The outcome is pinned (t != +-1):
       W = 10: sheet 1 no jump, sheet 2 {3}, sheet 3 {3, 7/5};
       W = 12: sheet 1 no jump, sheet 2 {3, 9/7}, sheet 3 {3};
       W = 14: sheet 1 no jump, sheet 2 {3}, sheet 3 {3, 9/7};
     generic rank r = n - 1 and nullity 2 at every jump; beta^2 = (t-1)/(t+1) images exactly
     1/2, 1/6, 1/8; no jump at t = 4/3, 3/2, 2 (beta^2 = 1/7, 1/5, 1/3); the deck image -t0 of
     every jump is not a jump.

Needs, beyond the fast suite: numpy; bash and coreutils timeout; a C compiler (CC, default gcc)
and FLINT with GMP and MPFR (tested with FLINT 3.5.0; set FLINT_PREFIX=<prefix> when its headers
and library are not on the default paths).  About 50 minutes; at most three determinant
processes at once, each under a 16 GB address-space limit.  Refuses to run under python -O.

Writes results/cyl_stratification_w1014.json (deterministic) only when every check passes at
all three weights.  Drill options, which never write the artifact or print the marker:
--weights W [W ...] (a subset of 10 12 14) and --claims DIR (an unpinned claims directory).
Marker: STRATIFICATION W10-W14 CERTIFIED."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction as F

import strat1014_family as fam
import strat1014_s2 as s2
import strat1014_verify as s5

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CLAIMS = os.path.join(ROOT, "results", "cyl_stratification_w1014", "claims")
OUT = os.path.join(ROOT, "results", "cyl_stratification_w1014.json")
CSRC = os.path.join(HERE, "strat1014_detmodp.c")
CLAIMS_SUMS_SHA256 = "265a6cbef0e63aa6d626d17e65c717e4cad3722fb8ec23fd1348156b6434a43b"
WEIGHTS = (10, 12, 14)
SIZES = {10: (185, 68), 12: (461, 160), 14: (1082, 365)}          # (rows m, columns n)
EXPECT = {10: {1: {}, 2: {"3": 2}, 3: {"3": 2, "7/5": 2}},
          12: {1: {}, 2: {"3": 2, "9/7": 2}, 3: {"3": 2}},
          14: {1: {}, 2: {"3": 2}, 3: {"3": 2, "9/7": 2}}}          # jump t -> nullity
BETA2 = {"3": F(1, 2), "7/5": F(1, 6), "9/7": F(1, 8)}
ODD_Q = (F(4, 3), F(3, 2), F(2))                                    # beta^2 = 1/7, 1/5, 1/3

fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


def sha256_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def stop(msg):
    print(f"STRATIFICATION: {len(fails)} CHECK(S) NOT MET -- {msg}", flush=True)
    sys.exit(1)


def verify_claims_dir(path):
    """Every file of the claims directory is listed in its SHA256SUMS with the right hash, and nothing else is there."""
    sums = os.path.join(path, "SHA256SUMS")
    check(os.path.isfile(sums) and sha256_file(sums) == CLAIMS_SUMS_SHA256,
          f"claims SHA256SUMS has the pinned sha256 ({CLAIMS_SUMS_SHA256[:12]})")
    if fails:
        return
    listed = {}
    for line in open(sums, encoding="utf-8"):
        h, name = line.rstrip("\n").split("  ", 1)
        listed[name] = h
    present = set(os.listdir(path)) - {"SHA256SUMS"}
    check(present == set(listed), f"claims directory holds exactly the {len(listed)} listed files")
    bad = [n for n, h in sorted(listed.items()) if not os.path.isfile(os.path.join(path, n))
           or sha256_file(os.path.join(path, n)) != h]
    check(not bad, f"every listed claim file has its listed sha256 (differing: {bad})")


def compile_driver(work):
    cc = os.environ.get("CC", "gcc")
    prefix = os.environ.get("FLINT_PREFIX")
    flags = ([f"-I{prefix}/include", f"-L{prefix}/lib", f"-Wl,-rpath,{prefix}/lib"] if prefix else [])
    libs = ["-lflint", "-lgmp", "-lmpfr"]
    exe = os.path.join(work, "detmodp")
    cp = subprocess.run([cc, "-O2", "-Wall", "-o", exe, CSRC] + flags + libs, capture_output=True, text=True)
    check(cp.returncode == 0, f"the FLINT driver compiles ({cc}; exit {cp.returncode}) {cp.stderr.strip()[-300:]}")
    ver = os.path.join(work, "flintver.c")
    with open(ver, "w") as fh:
        fh.write('#include <stdio.h>\n#include <flint/flint.h>\nint main(void) { printf("%s\\n", flint_version); return 0; }\n')
    cv = subprocess.run([cc, "-o", ver[:-2], ver] + flags + libs, capture_output=True, text=True)
    fv = subprocess.run([ver[:-2]], capture_output=True, text=True).stdout.strip() if cv.returncode == 0 else "unknown"
    print(f"      driver {os.path.basename(CSRC)} {sha256_file(CSRC)[:12]} compiled; binary sha256 "
          f"{sha256_file(exe)[:12] if os.path.isfile(exe) else '-'} (platform-dependent, not pinned); FLINT {fv}", flush=True)
    return exe


def build_matrices(weights, data):
    """The matrices M_sheet{S}_W{W}.json, exactly as the original build step (all three sheets per weight)."""
    sheets = (1, 2, 3)
    P4s = {sh: fam.p4_family(sh) for sh in sheets}
    keys = sorted({m for sh in sheets for m in P4s[sh]})
    for W in weights:
        cs, tgt, cols, rows = fam.domain_target(W)
        rowpos, red = fam.target_reduction(tgt, rows)
        ncheck = len(cols) if W <= 12 else 40
        nchk = fam.check_reduction(keys, W, cs, tgt, cols, rows, rowpos, red, ncheck)
        print(f"      W = {W}: classes n = {len(cols)}, target m = {len(rows)}; reduction shortcut = engine reduce_mod_d "
              f"on {nchk} residue vectors", flush=True)
        R = fam.residue_columns(keys, W, cs, tgt, cols, rowpos, red)
        for sh in sheets:
            Rsh = {mq: v for mq, v in R.items() if mq[0] in P4s[sh]}
            k, c0, ent, _ = fam.assemble(P4s[sh], Rsh)
            path = os.path.join(data, f"M_sheet{sh}_W{W}.json")
            fam.save(path, sh, W, k, c0, len(rows), len(cols), ent,
                     [tgt.basis[i] for i in rows], [cs.basis[j] for j in cols])
            h = sha256_file(path)
            check(h == s5.EXPECTED[W][sh], f"W = {W} sheet {sh}: rebuilt M_poly {len(rows)} x {len(cols)} has the pinned "
                                           f"sha256 ({h[:12]})")


def pins(W, claims):
    """The preregistered classification, read from the claims that strat1014_verify has just verified."""
    out = {}
    for sh in (1, 2, 3):
        tag = f"sheet{sh}_W{W}"
        S2 = json.load(open(os.path.join(claims, f"S2_{tag}.json")))
        S3 = json.load(open(os.path.join(claims, f"S3_{tag}.json")))
        m, n = SIZES[W]
        check((S2["m"], S2["n"], S2["r"]) == (m, n, n - 1) and S3["r"] == n - 1,
              f"W = {W} sheet {sh}: {m} x {n}, generic rank r = n - 1 = {n - 1}")
        jumps = {str(F(j["t"])): j["m"] for j in S3["jump_set"]}
        check(jumps == EXPECT[W][sh], f"W = {W} sheet {sh}: jump set over Q-bar (t != +-1) exactly "
                                      f"{EXPECT[W][sh] or 'empty'} (t: nullity)")
        check(all(t in BETA2 and (F(t) - 1) / (F(t) + 1) == BETA2[t] for t in jumps),
              f"W = {W} sheet {sh}: beta^2 images {sorted(str(BETA2[t]) for t in jumps if t in BETA2)}")
        if fails:
            return out
        check(not any(q in {F(t) for t in jumps} for q in ODD_Q), f"W = {W} sheet {sh}: no jump at t = 4/3, 3/2, 2")
        check(not any(str(-F(t)) in jumps for t in jumps), f"W = {W} sheet {sh}: no deck image -t0 of a jump is a jump")
        decs = {}
        for d in S3["decisions"]:
            decs[d["decision"]] = decs.get(d["decision"], 0) + 1
        out[f"sheet{sh}"] = {
            "m": S2["m"], "n": S2["n"], "r": S2["r"], "matrix_sha256": S2["matrix_sha256"], "scalar": S2["scalar"],
            "P1_degree": S2["P1_degree"], "P1o_degree": S2["P1o_degree"], "P1o_sha256": S2["P1o_sha256"],
            "P2_degree": S2["P2_degree"], "P2o_degree": S2["P2o_degree"], "P2o_sha256": S2["P2o_sha256"],
            "G_degree": S2["G_degree"],
            "G_factors": [[f["factor"], f["multiplicity"]] for f in S2["G_factors"]],
            "decisions": dict(sorted(decs.items())),
            "jump_set": [{"t": str(F(t)), "nullity": jumps[t], "beta2": str(BETA2[t])}
                         for t in sorted(jumps, key=F)],
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", type=int, nargs="+", default=list(WEIGHTS))
    ap.add_argument("--claims", default=None)
    a = ap.parse_args()
    if not __debug__:
        print("FAIL  run without -O: the build step's sanity checks are assert statements")
        sys.exit(1)
    weights = list(a.weights)
    drill = sorted(weights) != list(WEIGHTS) or a.claims is not None
    check(set(weights) <= set(WEIGHTS) and len(set(weights)) == len(weights), f"weights {weights} within {list(WEIGHTS)}")
    for t, x in BETA2.items():
        check((F(t) - 1) / (F(t) + 1) == x, f"beta^2 table: t = {t} -> {x}")
    if a.claims is None:
        claims = CLAIMS
        verify_claims_dir(claims)
    else:
        claims = os.path.abspath(a.claims)
        print(f"NOTE  drill: unpinned claims directory {claims}", flush=True)
    if drill:
        print("NOTE  drill run: the artifact is not written and the marker is not printed", flush=True)
    if fails:
        stop("claims")
    with tempfile.TemporaryDirectory(prefix="strat1014_") as work:
        data = os.path.join(work, "data")
        os.makedirs(data)
        exe = compile_driver(work)
        if fails:
            stop("driver")
        build_matrices(weights, data)
        if fails:
            stop("matrices")
        s5.BASE = work
        s2.DETMODP = exe
        summary = {}
        for W in weights:
            try:
                s5.main([str(W), "--data", data, "--results", claims, "--jobs", "3"])
                check(True, f"W = {W}: strat1014_verify re-verified every claim on sheets 1, 2, 3")
            except s5.CheckFailed as e:
                check(False, f"W = {W}: strat1014_verify: {e}")
            if fails:
                stop(f"verification at W = {W}")
            summary[str(W)] = pins(W, claims)
            if fails:
                stop(f"pins at W = {W}")
    print()
    if drill:
        print(f"STRATIFICATION DRILL RUN: every check passed at W = {weights} (no artifact, no marker)")
        return
    record = {
        "certificate": "cyl_stratification_w1014",
        "claims_sha256sums": CLAIMS_SUMS_SHA256,
        "expectations": ["results/cyl_discriminant_preregistration.md (H_F1-H_F3)",
                         "results/cyl_discriminant_w14_preregistration.md (H_W14)"],
        "parametrization": "N = (t^2-25)/(t^2-1), s = -24t/(t^2-1), beta^2 = (t-1)/(t+1); t = +-1 excluded",
        "weights": summary,
    }
    with open(OUT, "w") as fh:
        json.dump(record, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("STRATIFICATION W10-W14 CERTIFIED (exact classification over Q-bar away from t = +-1, open tools: "
          "W = 10 sheet 2 {3}, sheet 3 {3, 7/5}; W = 12 sheet 2 {3, 9/7}, sheet 3 {3}; W = 14 sheet 2 {3}, "
          "sheet 3 {3, 9/7}; sheet 1 never; nullity 2 at every jump; beta^2 images 1/2, 1/6, 1/8; "
          "no jump at beta^2 = 1/7, 1/5, 1/3; no deck image jumps)")


if __name__ == "__main__":
    main()
