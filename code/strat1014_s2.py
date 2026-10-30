"""S2 for one (sheet, W) (the step that produced the claims verified by strat1014_verify; the certificate
cyl_stratification_w1014 imports its helpers only -- run directly, it reads data/ and writes work/ and results/ under the
parent of the directory of strat1014_family.py, so run it from a copy): the generic rank r of M_poly over Q(t), proved both ways, and a complete over-approximation of
the locus {t : rank M(t) < r}, exact over Q.

Steps (all choices seeded; every intermediate saved under work/ or results/, so a rerun resumes):
 1. Rank lower bound: M_poly(t*) mod p = 2^31 - 1 at a seeded t*, reduced row echelon form -> rank r, pivot columns Cc, free
    columns J, pivot rows.  An r x r minor of M_poly is nonzero mod p at t*, hence nonzero in Z[t]: rank over Q(t) >= r.
 2. Rank upper bound: for each free column j, the null vector x with x_j = 1, x_(J\\j) = 0 over Q(t), found by modular solving at
    points, rational function reconstruction (common denominator), CRT and rational number reconstruction over 31-bit primes, then
    cleared to a primitive vector w^(j) in Z[t]^n and VERIFIED EXACTLY: M_poly * w^(j) = 0 identically (integer polynomial products).
    The w^(j) are independent (w^(j) is nonzero at column j and zero at the other free columns), so rank <= n - |J| = r.
 3. Minor P1 = det A1', A1' = M_poly[R1, Cc] with the integer content and the (t-1), (t+1) powers of each row and then each column
    divided out (rows R1: the lowest-norm rows independent at t*).  P1 differs from that r x r minor of M_poly by a nonzero integer
    and powers of t -+ 1 only, so every t != +-1 with rank M(t) < r is a root of P1.  P1 is computed EXACTLY by CRT from det A1' mod
    62-bit primes (bin/detmodp, FLINT nmod_poly_mat_det_interpolate) with the proven bounds
       deg P1 <= Dbound = min(sum of row degrees, sum of column degrees),
       |coefficient| <= B, B^2 = min over rows/columns of prod_i sum_j ||A_ij||_1^2   (Hadamard's inequality on |t| = 1 and
       Cauchy's coefficient estimate: |c_k| <= max_{|t|=1} |P1(t)| <= prod_i ||row_i(t)||_2 <= prod_i (sum_j ||A_ij||_1^2)^(1/2)),
    with prod p > 2 B 2^64 (a 64-bit margin, so a corrupted image fails the |c_k| <= B check with probability >= 1 - 2^-64).
    Checks: every image has degree <= Dbound; every reconstructed |c_k| <= B; a spot value mod a fresh 31-bit prime at a random
    point equals numpy's det of A1'(t) there.  P1o = P1 / ((t-1)^a (t+1)^b), exact.
 4. Second minor P2 = det A2', A2' from M_poly[R2, C2] (C2 = Cc with one column j1 exchanged for the free column j0, chosen
    with gcd(w_j0, w_j1) = 1 where possible; R2 greedy in a seeded random order), content-reduced the same way and computed EXACTLY
    by CRT with its own degree and Hadamard bounds, exactly as P1 (spot checks: numpy det mod a fresh 31-bit prime, and a pure-Python
    det mod one 62-bit prime against the FLINT image).  P2o = P2 / ((t-1)^b (t+1)^b').
 5. G = gcd(P1o, P2o), exact: monic modular gcds at the primes of least gcd degree, CRT + rational reconstruction to a primitive G
    in Z[t]; then EXACTLY G | P1o and G | P2o, and the cofactors H1 = P1o/G, H2 = P2o/G are coprime: a certificate prime q' dividing
    neither lc(H1) nor lc(H2) with gcd(H1 mod q', H2 mod q') = 1 (a common factor h over Q, primitive in Z[t], keeps its degree mod
    q' because lc(h) | lc(H1), and would divide both reductions).  Hence G = gcd(P1o, P2o) over Q (up to a unit); where the
    degrees are small, sympy's exact gcd over Z is compared as well.  Every t != +-1 with rank M(t) < r is a common root of P1o and
    P2o, hence a root of G: G is the complete over-approximation away from +-1.
 6. G factored over Q (sympy factor_list).
Usage: strat1014_s2.py SHEET W"""
import hashlib
import json
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as F

import numpy as np
from sympy import Poly, Symbol, factor_list

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat1014_family as fam  # noqa: E402
import strat1014_modla as modla  # noqa: E402

BASE = fam.BASE
DETMODP = os.path.join(BASE, "bin", "detmodp")
TSYM = Symbol("t")
P62_BOUND = 2**62


def log(msg):
    print(msg, flush=True)


# ------------------------------------------------------------------ exact integer polynomials (lists low -> high)

def ip_trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def ip_mul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def ip_add_into(acc, a):
    if len(acc) < len(a):
        acc.extend([0] * (len(a) - len(acc)))
    for i, x in enumerate(a):
        acc[i] += x


def ip_eval(a, x):
    acc = 0
    for c in reversed(a):
        acc = acc * x + c
    return acc


def ip_div_linear(a, s):
    """a / (t - s) for an integer polynomial with a(s) = 0 and integer s (synthetic division, exact)."""
    n = len(a) - 1
    q = [0] * n
    acc = 0
    for i in range(n, 0, -1):
        acc = acc * s + a[i]
        q[i - 1] = acc
    assert acc * s + a[0] == 0
    return q


def ip_val(a, s):
    """(t - s)-adic valuation of a nonzero integer polynomial, s in {1, -1}; returns (v, a / (t-s)^v)."""
    v = 0
    while a and ip_eval(a, s) == 0:
        a = ip_div_linear(a, s)
        v += 1
    return v, a


def ip_content(a):
    g = 0
    for x in a:
        g = math.gcd(g, x)
    return g


def q_divmod(a, b):
    """Division of polynomials with rational coefficients (Fractions)."""
    a = [F(x) for x in a]
    b = ip_trim([F(x) for x in b])
    q = [F(0)] * max(len(a) - len(b) + 1, 0)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] / b[-1]
        q[i] = c
        if c:
            for j, y in enumerate(b):
                a[i + j] -= c * y
    return ip_trim(q), ip_trim(a)


def to_sympy(a):
    return Poly(list(reversed(a)), TSYM, domain="ZZ") if a else Poly(0, TSYM, domain="ZZ")


def canon(a):
    """Canonical string of an integer polynomial (coefficients low -> high) for hashing."""
    return ",".join(str(x) for x in a)


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


# ------------------------------------------------------------------ submatrices

def submatrix(Md, rows, cols):
    """Dense list-of-lists of integer coefficient lists for M_poly[rows, cols]."""
    ci = {c: k for k, c in enumerate(cols)}
    ri = {r: k for k, r in enumerate(rows)}
    A = [[[] for _ in cols] for _ in rows]
    for ii, jj, e in Md["entries"]:
        if ii in ri and jj in ci:
            A[ri[ii]][ci[jj]] = list(e)
    return A


def content_reduce(A):
    """Divide each row, then each column, by its integer content and its powers of (t-1), (t+1) (two passes).
    Returns (A', record) with det A = (prod of integer contents) (t-1)^(sum) (t+1)^(sum) det A'."""
    A = [[list(e) for e in row] for row in A]
    rec = {"int": 1, "v_tm1": 0, "v_tp1": 0}
    r = len(A)
    for _pass in range(2):
        for axis in ("row", "col"):
            for k in range(r):
                ents = [A[k][j] for j in range(r)] if axis == "row" else [A[i][k] for i in range(r)]
                nz = [e for e in ents if e]
                if not nz:
                    continue
                g = 0
                for e in nz:
                    g = math.gcd(g, ip_content(e))
                vm = min(ip_val(e, 1)[0] for e in nz)
                vp = min(ip_val(e, -1)[0] for e in nz)
                if g == 1 and vm == 0 and vp == 0:
                    continue
                new = []
                for e in ents:
                    if e:
                        e = [x // g for x in e]
                        for _ in range(vm):
                            e = ip_div_linear(e, 1)
                        for _ in range(vp):
                            e = ip_div_linear(e, -1)
                    new.append(e)
                for idx, e in enumerate(new):
                    if axis == "row":
                        A[k][idx] = e
                    else:
                        A[idx][k] = e
                rec["int"] *= g
                rec["v_tm1"] += vm
                rec["v_tp1"] += vp
    return A, rec


def bounds(A):
    r = len(A)
    rowdeg = [max((len(e) - 1 for e in A[i] if e), default=-1) for i in range(r)]
    coldeg = [max((len(A[i][j]) - 1 for i in range(r) if A[i][j]), default=-1) for j in range(r)]
    D = min(sum(rowdeg), sum(coldeg))
    l1 = [[sum(abs(x) for x in e) for e in row] for row in A]
    B2row = 1
    for i in range(r):
        B2row *= sum(x * x for x in l1[i])
    B2col = 1
    for j in range(r):
        B2col *= sum(l1[i][j] ** 2 for i in range(r))
    return D, min(B2row, B2col), sum(rowdeg), sum(coldeg)


def primes_for(B2, extra_bits=64, skip=()):
    """62-bit primes (descending from 2^62) until prod p > 2 B 2^extra_bits, i.e. prod^2 > 4 B^2 2^(2 extra_bits)."""
    target = 4 * B2 * 2 ** (2 * extra_bits)
    ps, prod2 = [], 1
    q = P62_BOUND - 1
    while prod2 <= target:
        while not modla.isprime(q) or q in skip:
            q -= 1
        ps.append(q)
        prod2 *= q * q
        q -= 1
    return ps


def write_detmodp_input(path, A, primes):
    with open(path, "w") as fh:
        fh.write(f"{len(A)}\n")
        for row in A:
            fh.write(" ".join(f"{len(e)} " + " ".join(str(x) for x in e) if e else "0" for e in row) + "\n")
        fh.write(f"{len(primes)}\n" + "\n".join(str(p) for p in primes) + "\n")


def run_detmodp(inp, out, timeout):
    t0 = time.time()
    cp = subprocess.run(["bash", "-c", f"ulimit -v 16000000; exec timeout {timeout} {DETMODP} {inp} {out}"],
                        capture_output=True, text=True)
    if cp.returncode != 0:
        raise RuntimeError(f"detmodp exit {cp.returncode}: {cp.stderr}")
    return time.time() - t0


def read_detmodp(out, primes):
    res = {}
    for line in open(out):
        tok = line.split()
        if len(tok) < 2:
            continue
        p, D = int(tok[0]), int(tok[1])
        coeffs = [int(x) for x in tok[2:]]
        assert len(coeffs) == D + 1, f"incomplete line for prime {p}"
        assert all(0 <= c < p for c in coeffs)
        res[p] = coeffs
    missing = [p for p in primes if p not in res]
    assert not missing, f"{len(missing)} primes missing in {out}"
    return res


def eval_dense_mod(A, t0, p):
    r = len(A)
    X = np.zeros((r, r), dtype=np.int64)
    for i in range(r):
        for j in range(r):
            if A[i][j]:
                X[i, j] = modla.pmod_eval([c % p for c in A[i][j]], t0, p)
    return X


def det_mod_py(A, t0, q):
    """det A(t0) mod q in pure Python integers (Gaussian elimination); A dense list of integer coefficient lists."""
    M = [[modla.pmod_eval([c % q for c in e], t0, q) if e else 0 for e in row] for row in A]
    n = len(M)
    det = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c]), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det % q
        pv = M[c][c]
        det = det * pv % q
        inv = pow(pv, q - 2, q)
        rowc = M[c]
        for i in range(c + 1, n):
            f = M[i][c]
            if f:
                f = f * inv % q
                ri = M[i]
                M[i] = [(x - f * y) % q for x, y in zip(ri, rowc)]
    return det


# ------------------------------------------------------------------ null vectors

def solve_mod(A, b, p):
    """Solve A x = b mod p (A square int64, nonsingular); None if singular."""
    n = A.shape[0]
    Ab = np.concatenate([A % p, (b % p)[:, None]], axis=1)
    for c in range(n):
        nz = np.nonzero(Ab[c:, c])[0]
        if nz.size == 0:
            return None
        i = c + int(nz[0])
        if i != c:
            Ab[[c, i]] = Ab[[i, c]]
        inv = pow(int(Ab[c, c]), p - 2, p)
        Ab[c] = (Ab[c] * inv) % p
        col = Ab[:, c].copy()
        col[c] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            Ab[nzr] = (Ab[nzr] - (col[nzr, None] * Ab[c][None, :]) % p) % p
    return Ab[:, n].copy()


def nullvec_mod_prime(Cp, p, Rrows, Ccols, j, n, rng, L0=24):
    """Null vector x (x_j = 1, x_(J\\j) = 0, x_Cc = -A^-1 M[R, j]) as rational functions mod p:
    returns (q monic, {col: numerator}) or None.  Doubles the number of points until reconstruction checks out."""
    L = L0
    while L <= 4096:
        xs, vals, seen = [], [], set()
        while len(xs) < L + 6:
            t0 = rng.randrange(2, p - 2)
            if t0 in seen:
                continue
            seen.add(t0)
            Mt = modla.eval_mod(Cp, t0, p)
            A = Mt[np.ix_(Rrows, Ccols)]
            x = solve_mod(A, (-Mt[Rrows, j]) % p, p)
            if x is None:
                continue
            xs.append(t0)
            vals.append(x)
        V = np.array(vals, dtype=np.int64)          # (L+6) x r
        lam = [rng.randrange(1, p) for _ in range(V.shape[1])]
        z = [int(sum(int(V[k, i]) * lam[i] for i in range(V.shape[1])) % p) for k in range(V.shape[0])]
        D = (L - 1) // 2
        rr = modla.ratrecon_poly(xs[:L], z[:L], p, D)
        if rr is not None:
            a, q = rr
            if all(modla.pmod_eval(a, xs[k], p) == z[k] * modla.pmod_eval(q, xs[k], p) % p for k in range(L, L + 6)):
                qv = [modla.pmod_eval(q, x0, p) for x0 in xs]
                nums = {}
                ok = True
                for i, c in enumerate(Ccols):
                    ys = [int(V[k, i]) * qv[k] % p for k in range(L + 6)]
                    num = modla.interp_mod(xs[:L], ys[:L], p)
                    if len(num) - 1 > D or any(modla.pmod_eval(num, xs[k], p) != ys[k] for k in range(L, L + 6)):
                        ok = False
                        break
                    nums[c] = num
                if ok:
                    nums[j] = q
                    return q, nums, L
        L *= 2
    return None


def exact_nullvector(Md, Rrows, Ccols, j, seed, max_primes=400):
    """Exact primitive w in Z[t]^n with M_poly w = 0 identically (verified) and w_(J\\j) = 0; None on failure."""
    n, m = Md["n"], Md["m"]
    rng = random.Random(seed)
    primes = modla.primes_below(modla.P31 - 1000, max_primes)
    acc = None
    shape = None
    used = 0
    info = {}
    for p in primes:
        Cp = modla.coeff_arrays(Md, p)
        res = nullvec_mod_prime(Cp, p, Rrows, Ccols, j, n, rng)
        if res is None:
            continue
        q, nums, L = res
        sh = (len(q), tuple(len(nums[c]) for c in sorted(nums)))
        if shape is None or sh[0] > shape[0]:
            shape, acc, used = sh, {c: [(x, p) for x in nums[c]] + [(0, p)] * 0 for c in nums}, 1
            info["points_per_prime"] = L
        elif sh != shape:
            continue                                    # degenerate image at this prime
        else:
            acc = {c: [modla.crt_pair(r, mm, x, p) for (r, mm), x in zip(acc[c], nums[c])] for c in nums}
            used += 1
        if used >= 2:
            rec = {}
            ok = True
            for c, lst in acc.items():
                vals = [modla.ratrecon(r, mm) for r, mm in lst]
                if any(v is None for v in vals):
                    ok = False
                    break
                rec[c] = vals
            if not ok:
                continue
            den = 1
            for vals in rec.values():
                for v in vals:
                    den = math.lcm(den, v.denominator)
            w = [[] for _ in range(n)]
            for c, vals in rec.items():
                w[c] = ip_trim([int(v * den) for v in vals])
            g = 0
            for e in w:
                g = math.gcd(g, ip_content(e))
            w = [[x // g for x in e] for e in w]
            if verify_null(Md, w):
                info.update({"primes_used": used, "last_prime": p, "deg": max(len(e) - 1 for e in w if e),
                             "max_bits": max(abs(x).bit_length() for e in w for x in e)})
                return w, info
    return None, info


def verify_null(Md, w):
    """Exact: every row of M_poly * w is the zero polynomial."""
    acc = [[] for _ in range(Md["m"])]
    for ii, jj, e in Md["entries"]:
        if w[jj]:
            ip_add_into(acc[ii], ip_mul(e, w[jj]))
    return all(not any(a) for a in acc)


# ------------------------------------------------------------------ main

def main():
    sheet, W = int(sys.argv[1]), int(sys.argv[2])
    suffix = sys.argv[3] if len(sys.argv) > 3 else ""
    T0 = time.time()
    tag = f"sheet{sheet}_W{W}" + (f"_{suffix}" if suffix else "")
    mpath = os.path.join(BASE, "data", f"M_{tag}.json")
    msha = hashlib.sha256(open(mpath, "rb").read()).hexdigest()
    Md = fam.load(mpath)
    m, n = Md["m"], Md["n"]
    os.makedirs(os.path.join(BASE, "work"), exist_ok=True)
    out = {"sheet": sheet, "W": W, "matrix": os.path.relpath(mpath, BASE), "matrix_sha256": msha, "m": m, "n": n,
           "scalar": Md["scalar"]}
    seed = 1014100 + 10 * sheet + W
    rng = random.Random(seed)
    p = modla.P31
    C31 = modla.coeff_arrays(Md, p)
    tstar = rng.randrange(2, p - 2)
    R, piv, prow = modla.rref_mod(modla.eval_mod(C31, tstar, p), p)
    r = len(piv)
    Ccols = list(piv)
    J = [c for c in range(n) if c not in set(piv)]
    log(f"{tag}: {m} x {n}; seed {seed}; t* = {tstar} mod p = {p}: rank {r} (pivot rows {len(prow)}), free columns {J} "
        f"({Md['cols'][J[0]] if J else ''}) ({time.time() - T0:.1f} s)")
    out.update({"seed": seed, "rank_prime": p, "tstar": tstar, "r": r, "pivot_cols": Ccols, "pivot_rows": prow, "free_cols": J})
    # ---- 2. null vectors
    nv_path = os.path.join(BASE, "results", f"S2_nullvectors_{tag}.json")
    ws = {}
    if os.path.exists(nv_path):
        saved = json.load(open(nv_path))
        for j, w in saved["w"].items():
            assert verify_null(Md, w), "saved null vector fails verification"
            ws[int(j)] = w
        log(f"  null vectors loaded from {os.path.relpath(nv_path, BASE)} and re-verified exactly")
    else:
        for j in J:
            t1 = time.time()
            w, info = exact_nullvector(Md, prow, Ccols, j, seed + 7 * (j + 1))
            assert w is not None, f"null vector reconstruction failed for column {j}: {info}"
            ws[j] = w
            log(f"  null vector for free column {j} ({Md['cols'][j]}): degree {info['deg']}, max coefficient bits "
                f"{info['max_bits']}, {info['primes_used']} primes x {info['points_per_prime']}+6 points; "
                f"M_poly * w = 0 VERIFIED EXACTLY ({time.time() - t1:.1f} s)")
        json.dump({"sheet": sheet, "W": W, "matrix_sha256": msha, "w": {str(j): ws[j] for j in J}}, open(nv_path, "w"))
    indep = all(bool(ws[j][jj]) == (j == jj) for j in J for jj in J)
    assert indep
    out["nullvectors"] = {str(j): {"deg": max(len(e) - 1 for e in ws[j] if e),
                                   "max_bits": max(abs(x).bit_length() for e in ws[j] for x in e),
                                   "nonzero_entries": sum(1 for e in ws[j] if e), "sha256": sha(json.dumps(ws[j]))} for j in J}
    out["nullvectors_file"] = os.path.relpath(nv_path, BASE)
    log(f"  rank over Q(t): r = {r} (lower bound: nonzero minor mod p at t*; upper bound: {len(J)} independent exact null "
        f"vectors); generic nullity n - r = {n - r}")
    # ---- 3. P1
    j0 = J[0] if J else None
    # row order by content-reduced log-norm on the pivot columns
    Cset = set(Ccols)
    rowl2 = [0] * m
    rowg = [0] * m
    for ii, jj, e in Md["entries"]:
        if jj in Cset:
            rowg[ii] = math.gcd(rowg[ii], ip_content(e))
    for ii, jj, e in Md["entries"]:
        if jj in Cset and rowg[ii]:
            l1 = sum(abs(x) for x in e) // rowg[ii]
            rowl2[ii] += l1 * l1
    order = sorted((i for i in range(m) if rowl2[i]), key=lambda i: (rowl2[i], i))
    Mt = modla.eval_mod(C31, tstar, p)
    _, pv, _ = modla.rref_mod(Mt[np.ix_(order, Ccols)].T.copy(), p)
    R1 = sorted(order[k] for k in pv)
    assert len(R1) == r
    A1 = submatrix(Md, R1, Ccols)
    A1r, rec1 = content_reduce(A1)
    D1, B2, sr, sc = bounds(A1r)
    primes1 = primes_for(B2)
    log(f"  P1: rows R1 = the {r} lowest-norm rows independent at t*; content removed: integer {rec1['int'].bit_length()} bits, "
        f"(t-1)^{rec1['v_tm1']}, (t+1)^{rec1['v_tp1']}; degree bound {D1} (rows {sr}, columns {sc}); log2 B = "
        f"{B2.bit_length() / 2:.0f}; primes needed {len(primes1)} (62-bit, prod > 2 B 2^64)")
    in1 = os.path.join(BASE, "work", f"A1_{tag}.in")
    out1 = os.path.join(BASE, "work", f"A1_{tag}.out")
    write_detmodp_input(in1, A1r, primes1)
    dt = run_detmodp(in1, out1, 86400)
    imgs = read_detmodp(out1, primes1)
    assert all(len(imgs[q]) - 1 <= D1 for q in primes1), "an image exceeds the degree bound"
    log(f"  detmodp over {len(primes1)} primes: {dt:.1f} s (input sha256 {hashlib.sha256(open(in1, 'rb').read()).hexdigest()})")
    # CRT
    deg = max(len(imgs[q]) for q in primes1) - 1
    P1 = []
    for k in range(deg + 1):
        x, mm = 0, 1
        for q in primes1:
            c = imgs[q][k] if k < len(imgs[q]) else 0
            x, mm = modla.crt_pair(x, mm, c, q)
        P1.append(modla.symmetric(x, mm))
    P1 = ip_trim(P1)
    assert all(c * c <= B2 for c in P1), "a reconstructed coefficient exceeds the Hadamard bound"
    # spot check vs numpy
    ps_spot = modla.primes_below(modla.P31 - 5000, 1)[0]
    ts = rng.randrange(2, ps_spot - 2)
    spot_np = modla.det_mod(eval_dense_mod(A1r, ts, ps_spot), ps_spot)
    spot_P1 = modla.pmod_eval([c % ps_spot for c in P1], ts, ps_spot)
    assert spot_np == spot_P1, "spot check failed"
    a_m, P1o = ip_val(P1, 1)
    a_p, P1o = ip_val(P1o, -1)
    p1path = os.path.join(BASE, "results", f"S2_P1_{tag}.json")
    json.dump({"P1": P1, "P1o": P1o, "v_tm1": a_m, "v_tp1": a_p}, open(p1path, "w"))
    log(f"  P1 exact: degree {len(P1) - 1}, max coefficient bits {max(abs(c).bit_length() for c in P1)} (bound "
        f"{B2.bit_length() // 2}); spot check mod {ps_spot} at t = {ts}: numpy det = P1(t) = {spot_np}; "
        f"P1 = (t-1)^{a_m} (t+1)^{a_p} P1o, deg P1o = {len(P1o) - 1}; sha256(canonical P1o) {sha(canon(P1o))}")
    out.update({"R1": R1, "A1_content": {"int_bits": rec1["int"].bit_length(), "v_tm1": rec1["v_tm1"], "v_tp1": rec1["v_tp1"]},
                "P1_degree_bound": D1, "P1_log2_B": B2.bit_length() / 2, "P1_B2": str(B2), "P1_primes": primes1,
                "P1_degree": len(P1) - 1, "P1_v_tm1": a_m, "P1_v_tp1": a_p, "P1o_degree": len(P1o) - 1,
                "P1o_sha256": sha(canon(P1o)), "P1_file": os.path.relpath(p1path, BASE), "detmodp_P1_seconds": round(dt, 1),
                "spot": {"prime": ps_spot, "t": ts, "value": spot_np}})
    # ---- 4. P2
    if J:
        w0 = to_sympy(ws[j0][j0])
        j1, best = None, None
        for c in Ccols:
            if ws[j0][c]:
                gdeg = w0.gcd(to_sympy(ws[j0][c])).degree()
                if best is None or gdeg < best:
                    j1, best = c, gdeg
                if gdeg == 0:
                    break
        C2 = sorted([c for c in Ccols if c != j1] + [j0])
    else:
        j1, best = None, None                       # full column rank: P2 uses all columns and other rows
        C2 = list(Ccols)
    perm = list(range(m))
    random.Random(seed + 2).shuffle(perm)
    _, pv2, _ = modla.rref_mod(Mt[np.ix_(perm, C2)].T.copy(), p)
    R2 = sorted(perm[k] for k in pv2)
    assert len(R2) == r
    A2r, rec2 = content_reduce(submatrix(Md, R2, C2))
    D2, B22, sr2, sc2 = bounds(A2r)
    primes2 = primes_for(B22)
    in2 = os.path.join(BASE, "work", f"A2x_{tag}.in")
    out2 = os.path.join(BASE, "work", f"A2x_{tag}.out")
    write_detmodp_input(in2, A2r, primes2)
    dt2 = run_detmodp(in2, out2, 86400)
    img2 = read_detmodp(out2, primes2)
    assert all(len(img2[q]) - 1 <= D2 for q in primes2), "a P2 image exceeds the degree bound"
    deg2 = max(len(img2[q]) for q in primes2) - 1
    P2 = []
    for k in range(deg2 + 1):
        x, mm = 0, 1
        for q in primes2:
            x, mm = modla.crt_pair(x, mm, img2[q][k] if k < len(img2[q]) else 0, q)
        P2.append(modla.symmetric(x, mm))
    P2 = ip_trim(P2)
    assert all(c * c <= B22 for c in P2), "a reconstructed P2 coefficient exceeds the Hadamard bound"
    ts2 = rng.randrange(2, ps_spot - 2)
    spot2_np = modla.det_mod(eval_dense_mod(A2r, ts2, ps_spot), ps_spot)
    assert spot2_np == modla.pmod_eval([c % ps_spot for c in P2], ts2, ps_spot), "P2 spot check failed"
    q62 = primes2[0]
    ts3 = random.Random(seed + 3).randrange(2, q62 - 2)
    spot3_py = det_mod_py(A2r, ts3, q62)
    assert spot3_py == modla.pmod_eval(img2[q62], ts3, q62), "pure-Python check of a FLINT image failed"
    b_m, P2o = ip_val(P2, 1)
    b_p, P2o = ip_val(P2o, -1)
    p2path = os.path.join(BASE, "results", f"S2_P2_{tag}.json")
    json.dump({"P2": P2, "P2o": P2o, "v_tm1": b_m, "v_tp1": b_p}, open(p2path, "w"))
    log(f"  P2: columns C2 = pivot columns with {j1} ({Md['cols'][j1] if j1 is not None else '-'}) exchanged for {j0} "
        f"(deg gcd(w_j0, w_j1) = {best}); rows R2 greedy in a seeded random order; content removed: integer "
        f"{rec2['int'].bit_length()} bits, (t-1)^{rec2['v_tm1']}, (t+1)^{rec2['v_tp1']}; degree bound {D2}, log2 B = "
        f"{B22.bit_length() / 2:.0f}, {len(primes2)} primes ({dt2:.1f} s); P2 exact: degree {len(P2) - 1}, max coefficient "
        f"bits {max(abs(c).bit_length() for c in P2)}; spot check mod {ps_spot} at t = {ts2}: numpy det = P2(t) = {spot2_np}; "
        f"pure-Python det A2'(t) mod {q62} at t = {ts3} = FLINT image = {spot3_py}; P2 = (t-1)^{b_m} (t+1)^{b_p} P2o, "
        f"deg P2o = {len(P2o) - 1}; sha256(canonical P2o) {sha(canon(P2o))}")
    # ---- 5. G = gcd(P1o, P2o), exact
    gq = {}
    for q in primes2:
        if P1o[-1] % q == 0 or P2o[-1] % q == 0:
            continue
        gq[q] = modla.pgcdp([c % q for c in P1o], [c % q for c in P2o], q)
    dmin = min(len(g) for g in gq.values())
    lucky = [q for q in primes2 if q in gq and len(gq[q]) == dmin]
    G = None
    used = []
    acc = None
    cert = None
    for q in lucky:
        used.append(q)
        acc = [(c, q) for c in gq[q]] if acc is None else [modla.crt_pair(r_, mm, c, q) for (r_, mm), c in zip(acc, gq[q])]
        vals = [modla.ratrecon(r_, mm) for r_, mm in acc]
        if any(v is None for v in vals):
            continue
        den = 1
        for v in vals:
            den = math.lcm(den, v.denominator)
        cand = [int(v * den) for v in vals]
        g_ = ip_content(cand)
        cand = [x // g_ for x in cand]
        if cand[-1] < 0:
            cand = [-x for x in cand]
        q1, r1 = q_divmod(P1o, cand)
        q2_, r2 = q_divmod(P2o, cand)
        if r1 or r2 or any(F(x).denominator != 1 for x in q1 + q2_):
            continue
        H1, H2 = [int(x) for x in q1], [int(x) for x in q2_]
        assert ip_mul(cand, H1) == P1o and ip_mul(cand, H2) == P2o
        for qc in [x for x in primes2 if x not in used]:
            if H1[-1] % qc == 0 or H2[-1] % qc == 0:
                continue
            if modla.pgcdp([c % qc for c in H1], [c % qc for c in H2], qc) == [1]:
                cert = qc
                break
        if cert is not None:
            G = cand
            break
    assert G is not None and cert is not None, "no certified gcd from the available primes"
    cross = None
    if len(P1o) + len(P2o) <= 800:
        gs = to_sympy(P1o).gcd(to_sympy(P2o)).primitive()[1]       # sympy keeps the common integer content; compare primitive parts
        cross = [int(c) for c in reversed(gs.all_coeffs())]
        cross = cross if cross[-1] > 0 else [-c for c in cross]
        assert cross == G, "sympy's exact gcd differs"
    log(f"  G = gcd(P1o, P2o) EXACT: degree {len(G) - 1} (modular gcd degree {dmin - 1} at {len(lucky)} of {len(gq)} primes; "
        f"reconstructed from {len(used)} primes); G | P1o and G | P2o exactly; cofactors coprime: gcd(P1o/G mod q', P2o/G mod q')"
        f" = 1 at q' = {cert} (q' divides neither leading coefficient)"
        + (f"; sympy's exact gcd over Z agrees" if cross is not None else ""))
    fl = factor_list(to_sympy(G).as_expr(), TSYM)
    facs = []
    for f, e in fl[1]:
        fp = Poly(f, TSYM)
        facs.append({"factor": str(f), "degree": fp.degree(), "multiplicity": e,
                     "coeffs": [int(c) for c in reversed(fp.all_coeffs())]})
    log(f"  G over Q: content {fl[0]}; factors: " + "; ".join(f"({x['factor']})^{x['multiplicity']}" for x in facs))
    out.update({"j0": j0, "j1": j1, "gcd_deg_w_j0_w_j1": best, "C2": C2, "R2": R2,
                "A2_content": {"int_bits": rec2["int"].bit_length(), "v_tm1": rec2["v_tm1"], "v_tp1": rec2["v_tp1"]},
                "P2_degree_bound": D2, "P2_log2_B": B22.bit_length() / 2, "P2_B2": str(B22), "P2_primes": primes2,
                "P2_degree": len(P2) - 1, "P2_v_tm1": b_m, "P2_v_tp1": b_p, "P2o_degree": len(P2o) - 1,
                "P2o_sha256": sha(canon(P2o)), "P2_file": os.path.relpath(p2path, BASE), "detmodp_P2_seconds": round(dt2, 1),
                "spot2": {"prime": ps_spot, "t": ts2, "value": spot2_np}, "spot3": {"prime": q62, "t": ts3, "value": spot3_py},
                "gcd_mod_degree": dmin - 1, "gcd_primes_used": used, "G": G, "G_degree": len(G) - 1,
                "certificate_prime": cert, "sympy_gcd_crosscheck": cross is not None,
                "G_content": str(fl[0]), "G_factors": facs, "seconds": round(time.time() - T0, 1)})
    rpath = os.path.join(BASE, "results", f"S2_{tag}.json")
    json.dump(out, open(rpath, "w"), indent=1)
    log(f"  saved {os.path.relpath(rpath, BASE)} ({time.time() - T0:.1f} s)")


if __name__ == "__main__":
    main()
