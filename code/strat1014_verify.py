"""S5: fail-closed re-verification of every exact S2-S3 claim at one weight, from the saved matrix files (called by the
certificate cyl_stratification_w1014, which sets BASE to its work directory and strat1014_s2.DETMODP to the driver it compiles).
Usage: strat1014_verify.py W [--data DIR] [--results DIR] [--jobs K]
Exit status 0 iff every check passes for sheets 1, 2, 3; any failed check or exception exits nonzero.

Per sheet (claims read from results/S2_*.json, S3_*.json and the files they name):
  M  the matrix file's sha256 equals the value fixed below (EXPECTED);
  R1 rank lower bound: rank of M_poly(t*) mod p recomputed = r;
  R2 rank upper bound: the n - r saved null vectors satisfy M_poly w = 0 exactly (integer polynomial products), each nonzero at its
     own free column and zero at the other free columns;
  P1 the submatrix A1 = M_poly[R1, pivot columns] is rebuilt and content-reduced; its degree bound and Hadamard bound B are
     recomputed; det A1' mod every listed prime is recomputed by bin/detmodp in a fresh directory; the saved exact P1 must agree
     with every image, have degree <= the bound and |coefficients| <= B, and prod p > 2 B 2^64 must hold; then P1 = det A1';
  P1o P1 = (t-1)^a (t+1)^b P1o exactly, P1o(+-1) != 0;
  P2 the same for A2' = content-reduced M_poly[R2, C2] and the saved exact P2 (all its primes recomputed), P2 = (t-1)^b (t+1)^b' P2o,
     and the saved pure-Python spot value of det A2'(t) mod one 62-bit prime equals the fresh FLINT image there;
  G  G divides P1o and P2o exactly; the cofactors are coprime (gcd mod the certificate prime q' = 1, q' dividing neither leading
     coefficient); the saved factor list multiplies to G and each factor is irreducible (sympy);
  S3 every irreducible factor of G has a decision: no drop -> f(alpha) = 0 mod p, p does not divide lc(f), rank M_poly(alpha) mod p
     = r; drop (linear f, root t0) -> the saved exact null vectors satisfy M_poly(t0) v = 0 exactly, are in RREF normalization on
     their free columns, and rank M_poly(t0) mod each listed prime = n - nullity; the jump set and multiplicities equal the saved
     ones; each jump's deck image -t0 is either not a root of G or is decided;
  then prints the verdict line."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from fractions import Fraction as F

from sympy import Poly, factor_list

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat1014_family as fam  # noqa: E402
import strat1014_modla as modla  # noqa: E402
import strat1014_s2 as s2  # noqa: E402

BASE = fam.BASE
EXPECTED = {
    10: {1: "86791da0754b204a9dac8fc950d70ebaff1d44a7dcf33be9c37e68175c9149ec",
         2: "90b68ffa717b0f498c658fb55a4d487aa038819ce2c82979483c1d0d5b09b9d6",
         3: "4f7bbba93047e3bb1175c3471c5fdeac0165419df4d932c1105e8424b08b7bd8"},
    12: {1: "d0a897d9ac6d746ae77d1b7669805fa97906acc3a482e8cce186c521775d928a",
         2: "1e333f94753809bc1a9158c522a1cce6a37968ab9e07a921f49bbc70f58fc9b1",
         3: "9c8e57d80a7c6036f7e843eaaa665d73f4022bf26e039660a1b4e333f8efcef9"},
    14: {1: "7b2e8ae47226156965dbc27af722df21da7e8061067b3a87f1190e6a0a019811",
         2: "1cac02ab7a3d0c5c8c67ec1089c2d401c29ebccb128b46277cfb5395dcbae2d8",
         3: "615780bbd8723a6fb25437a19148b15296634dfae6626bd85a9450a61abfe41c"},
}


class CheckFailed(Exception):
    pass


def check(cond, msg):
    if not cond:
        raise CheckFailed(msg)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("W", type=int)
    ap.add_argument("--data", default=os.path.join(BASE, "data"))
    ap.add_argument("--results", default=os.path.join(BASE, "results"))
    ap.add_argument("--jobs", type=int, default=3)
    a = ap.parse_args(argv)
    W = a.W
    T0 = time.time()
    stamp = time.strftime("%Y%m%dT%H%M%S")
    wd = os.path.join(BASE, "s5_runs", f"W{W}_{stamp}")
    os.makedirs(wd, exist_ok=True)
    print(f"S5 verify W = {W}: data {a.data}, results {a.results}, fresh detmodp directory {os.path.relpath(wd, BASE)}", flush=True)
    sheets = {}
    procs = []
    # phase 1: matrices, claims, and the detmodp recomputations (run concurrently, at most --jobs at once)
    for sheet in (1, 2, 3):
        tag = f"sheet{sheet}_W{W}"
        mpath = os.path.join(a.data, f"M_{tag}.json")
        raw = open(mpath, "rb").read()
        h = hashlib.sha256(raw).hexdigest()
        check(h == EXPECTED[W][sheet], f"{tag}: matrix sha256 {h} != expected {EXPECTED[W][sheet]}")
        Md = json.loads(raw)
        S2 = json.load(open(os.path.join(a.results, f"S2_{tag}.json")))
        S3 = json.load(open(os.path.join(a.results, f"S3_{tag}.json")))
        check(S2["matrix_sha256"] == h and S3["matrix_sha256"] == h, f"{tag}: results name another matrix")
        A1r, _ = s2.content_reduce(s2.submatrix(Md, S2["R1"], S2["pivot_cols"]))
        D1, B2, _, _ = s2.bounds(A1r)
        check(D1 == S2["P1_degree_bound"] and B2 == int(S2["P1_B2"]), f"{tag}: P1 bounds differ")
        check(s2.primes_for(B2) == S2["P1_primes"], f"{tag}: P1 prime list is not the one the bound requires")
        in1, out1 = os.path.join(wd, f"A1_{tag}.in"), os.path.join(wd, f"A1_{tag}.out")
        s2.write_detmodp_input(in1, A1r, S2["P1_primes"])
        A2r, _ = s2.content_reduce(s2.submatrix(Md, S2["R2"], S2["C2"]))
        check(len(S2["R2"]) == S2["r"] and len(S2["C2"]) == S2["r"], f"{tag}: A2 is not r x r")
        D2, B22, _, _ = s2.bounds(A2r)
        check(D2 == S2["P2_degree_bound"] and B22 == int(S2["P2_B2"]), f"{tag}: P2 bounds differ")
        check(s2.primes_for(B22) == S2["P2_primes"], f"{tag}: P2 prime list is not the one the bound requires")
        in2, out2 = os.path.join(wd, f"A2_{tag}.in"), os.path.join(wd, f"A2_{tag}.out")
        s2.write_detmodp_input(in2, A2r, S2["P2_primes"])
        sheets[sheet] = dict(tag=tag, Md=Md, S2=S2, S3=S3, D1=D1, B2=B2, D2=D2, B22=B22, A2r=A2r, out1=out1, out2=out2)
        for inp, outp in ((in1, out1), (in2, out2)):
            while sum(1 for pr in procs if pr.poll() is None) >= a.jobs:
                time.sleep(1)
            procs.append(subprocess.Popen(["bash", "-c", f"ulimit -v 16000000; exec timeout 86400 {s2.DETMODP} {inp} {outp}"]))
    for pr in procs:
        check(pr.wait() == 0, "a detmodp run failed")
    print(f"  detmodp recomputations done ({time.time() - T0:.1f} s)", flush=True)
    # phase 2: checks
    for sheet in (1, 2, 3):
        T1 = time.time()
        z = sheets[sheet]
        tag, Md, S2, S3 = z["tag"], z["Md"], z["S2"], z["S3"]
        m, n, r = Md["m"], Md["n"], S2["r"]
        p, ts = S2["rank_prime"], S2["tstar"]
        C31 = modla.coeff_arrays(Md, p)
        check(modla.rank_mod(modla.eval_mod(C31, ts, p), p) == r, f"{tag}: rank at t* differs from r")
        nv = json.load(open(os.path.join(a.results, os.path.basename(S2["nullvectors_file"]))))
        J = S2["free_cols"]
        check(len(J) == n - r and sorted(int(j) for j in nv["w"]) == sorted(J), f"{tag}: wrong number of null vectors")
        for j in J:
            w = nv["w"][str(j)]
            check(s2.verify_null(Md, w), f"{tag}: null vector {j} fails M_poly w = 0")
            check(all(bool(w[jj]) == (jj == j) for jj in J), f"{tag}: null vectors not independent")
        P = json.load(open(os.path.join(a.results, os.path.basename(S2["P1_file"]))))
        P1, P1o = P["P1"], P["P1o"]
        imgs = s2.read_detmodp(z["out1"], S2["P1_primes"])
        for q in S2["P1_primes"]:
            check(len(imgs[q]) - 1 <= z["D1"], f"{tag}: image mod {q} exceeds the degree bound")
            red = s2.ip_trim([c % q for c in P1])
            check(red == imgs[q], f"{tag}: P1 mod {q} differs from det A1' mod {q}")
        check(len(P1) - 1 <= z["D1"] and all(c * c <= z["B2"] for c in P1) and any(P1), f"{tag}: P1 bound check")
        prod = 1
        for q in S2["P1_primes"]:
            prod *= q
        check(prod * prod > 4 * z["B2"] * 2**128, f"{tag}: primes do not cover 2 B 2^64")
        rec = list(P1o)
        for _ in range(P["v_tm1"]):
            rec = s2.ip_mul(rec, [-1, 1])
        for _ in range(P["v_tp1"]):
            rec = s2.ip_mul(rec, [1, 1])
        check(s2.ip_trim(rec) == P1 and s2.ip_eval(P1o, 1) != 0 and s2.ip_eval(P1o, -1) != 0, f"{tag}: P1o")
        P2d = json.load(open(os.path.join(a.results, os.path.basename(S2["P2_file"]))))
        P2, P2o = P2d["P2"], P2d["P2o"]
        img2 = s2.read_detmodp(z["out2"], S2["P2_primes"])
        for q in S2["P2_primes"]:
            check(len(img2[q]) - 1 <= z["D2"], f"{tag}: P2 image mod {q} exceeds the degree bound")
            check(s2.ip_trim([c % q for c in P2]) == img2[q], f"{tag}: P2 mod {q} differs from det A2' mod {q}")
        check(len(P2) - 1 <= z["D2"] and all(c * c <= z["B22"] for c in P2) and any(P2), f"{tag}: P2 bound check")
        prod = 1
        for q in S2["P2_primes"]:
            prod *= q
        check(prod * prod > 4 * z["B22"] * 2**128, f"{tag}: P2 primes do not cover 2 B 2^64")
        rec = list(P2o)
        for _ in range(P2d["v_tm1"]):
            rec = s2.ip_mul(rec, [-1, 1])
        for _ in range(P2d["v_tp1"]):
            rec = s2.ip_mul(rec, [1, 1])
        check(s2.ip_trim(rec) == P2 and s2.ip_eval(P2o, 1) != 0 and s2.ip_eval(P2o, -1) != 0, f"{tag}: P2o")
        sp3 = S2["spot3"]
        check(s2.det_mod_py(z["A2r"], sp3["t"], sp3["prime"]) == modla.pmod_eval(img2[sp3["prime"]], sp3["t"], sp3["prime"])
              == sp3["value"], f"{tag}: pure-Python det check of a FLINT image fails")
        G = S2["G"]
        q1, r1 = s2.q_divmod(P1o, G)
        q2, r2 = s2.q_divmod(P2o, G)
        check(not r1 and not r2 and all(F(x).denominator == 1 for x in q1 + q2), f"{tag}: G does not divide P1o and P2o")
        H1, H2 = [int(x) for x in q1], [int(x) for x in q2]
        qc = S2["certificate_prime"]
        check(H1[-1] % qc != 0 and H2[-1] % qc != 0, f"{tag}: certificate prime divides a leading coefficient")
        check(modla.pgcdp([c % qc for c in H1], [c % qc for c in H2], qc) == [1], f"{tag}: cofactors not coprime mod q'")
        prodf = [int(S2["G_content"])]
        for f in S2["G_factors"]:
            fp = Poly(list(reversed(f["coeffs"])), s2.TSYM)
            fl = factor_list(fp.as_expr(), s2.TSYM)
            check(len(fl[1]) == 1 and fl[1][0][1] == 1 and abs(int(fl[0])) == 1, f"{tag}: factor {f['factor']} not irreducible")
            for _ in range(f["multiplicity"]):
                prodf = s2.ip_mul(prodf, f["coeffs"])
        check(s2.ip_trim(prodf) == G, f"{tag}: factor list does not multiply to G")
        # S3
        decs = {d["factor"]: d for d in S3["decisions"]}
        check(sorted(decs) == sorted(f["factor"] for f in S2["G_factors"]), f"{tag}: decisions do not match the factors")
        jumps = []
        roots_G = [F(-f["coeffs"][0], f["coeffs"][1]) for f in S2["G_factors"] if len(f["coeffs"]) == 2]
        Cache = {}
        for f in S2["G_factors"]:
            d = decs[f["factor"]]
            fc = f["coeffs"]
            if d["decision"] == "no drop" and "alpha" in d:
                pp, al = d["p"], d["alpha"]
                check(fc[-1] % pp != 0 and sum(c * pow(al, k, pp) for k, c in enumerate(fc)) % pp == 0, f"{tag}: alpha")
                if pp not in Cache:
                    Cache[pp] = modla.coeff_arrays(Md, pp)
                check(modla.rank_mod(modla.eval_mod(Cache[pp], al, pp), pp) == r, f"{tag}: no-drop rank at {f['factor']}")
            else:
                check(len(fc) == 2, f"{tag}: a nonlinear factor is undecided")
                t0 = F(-fc[0], fc[1])
                check(str(t0) == d["t0"], f"{tag}: root mismatch")
                nvd = json.load(open(os.path.join(a.results, os.path.basename(d["null_vectors_file"]))))
                Mint = modla.exact_int_matrix(Md, t0)
                vecs, fr = nvd["vectors"], nvd["free_cols"]
                nul = len(vecs)
                for v in vecs:
                    check(modla.sparse_matvec_is_zero(Mint, m, v), f"{tag}: null vector at {t0} fails")
                check(all((v[c] != 0) == (k == kk) for k, v in enumerate(vecs) for kk, c in enumerate(fr)),
                      f"{tag}: null vectors at {t0} not in RREF normalization")
                for pp in d["rank_primes"]:
                    Ap = modla.eval_mod(modla.coeff_arrays(Md, pp), modla.tmod(t0, pp), pp)
                    check(modla.rank_mod(Ap, pp) == n - nul, f"{tag}: rank at {t0} mod {pp} differs")
                check(nul == d["nullity"], f"{tag}: nullity at {t0} differs")
                if nul > n - r:
                    jumps.append((t0, nul))
        saved = sorted((F(x["t"]), x["m"]) for x in S3["jump_set"])
        check(sorted(jumps) == saved, f"{tag}: jump set differs from the saved one")
        deck = {str(-t0): ("not a root of G, so no drop" if -t0 not in roots_G else
                           ("a jump" if any(-t0 == jt for jt, _ in jumps) else "a root of G, decided no drop"))
                for t0, _ in jumps}
        js = ", ".join(f"t = {t0} (m = {mm})" for t0, mm in sorted(jumps)) or "empty"
        print(f"  {tag}: r = {r} (n = {n}, m = {m}); G of degree {len(G) - 1} = "
              + " ".join(f"({f['factor']})^{f['multiplicity']}" for f in S2["G_factors"])
              + f"; jump set {js}; deck images of the jumps {deck}; ALL CHECKS PASS ({time.time() - T1:.1f} s)", flush=True)
    print(f"S5 W = {W}: ALL CHECKS PASS for sheets 1, 2, 3 ({time.time() - T0:.1f} s)", flush=True)


if __name__ == "__main__":
    try:
        main()
    except CheckFailed as e:
        print(f"S5 CHECK FAILED: {e}", flush=True)
        sys.exit(1)
