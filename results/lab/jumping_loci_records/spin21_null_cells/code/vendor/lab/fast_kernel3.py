"""Lab: genuine kernel of ad_{I_3} by the QUOTIENT formulation, modular linear algebra, CRT lifting and exact
verification (2026-09-06, v2: one multi-prime sweep in the worker pool; per-prime eliminations in parallel; v4 2026-09-07: a prime dividing a denominator is skipped instead of aborting the build -- review finding).
Same output record format as lab/fast_kernel2.py; lab record, not a certificate.

Why a new path.  fast_kernel2 eliminates the full [A | D_t] system (target x (candidates + d-source)) with a dense
mod-p rref, keeps every exact column in memory across primes, and recovers each class vector's d-part by a second
full-size elimination; at weight 22 (target index 25: 78,771 x (25,258 + 54,370)) that is two 50 GB matrices and
about 30 h per prime on top of an 8 h single-core exact build.  Here:

  1. TRIANGULAR d.  Order monomials by the descending multiset of their letter indices (tie-break: the canonical
     tuple).  For a source monomial b the max-key term of d(b) raises a maximal letter; it contains the unique letter
     M+1, so b -> LT(d(b)) is injective, and every other term of d(b) has a smaller key.  So Im d is triangular with
     pivot set {LT(d(b))}: the quotient T/Im d has basis = the NON-PIVOT monomials, and reducing a vector modulo Im d
     is one descending sweep.  Both the target sector and the candidate sector get this structure; the pivot
     coefficients are multiplicities (small positive integers), so the structure is the same over Q and mod p.
  2. THE INDUCED MAP.  ad(d b') = d(ad b'), so A(Im D_c) is inside Im D_t and A descends to
     A-hat : V/Im D_c -> T/Im D_t, whose matrix in the non-pivot bases has as columns the normal forms of A_b for the
     NON-PIVOT candidates b only.  The genuine kernel (classes of charges modulo total derivatives) is exactly
     ker A-hat: a kernel vector supported on non-pivot candidates cannot lie in Im D_c (every nonzero element of
     Im D_c has a nonzero pivot entry), so independence is plain independence, and the representatives are canonical.
  3. The exact residues of the non-pivot candidates are built in a fork pool, reduced modulo the whole prime set at
     once and the rationals dropped (a denominator divisible by a prime marks it bad).  The normal forms are computed
     in the pool by ONE sweep carrying all primes' residues (the elimination order is the pivot structure, prime-
     independent; a zero residue gives a zero factor).  The dense rrefs of A-hat, one per prime, run in parallel
     processes; the kernel comes from the free columns, the d-part S of each kernel vector from the sweep factors of
     A v (A v = -D_t S).  Good-prime checks: pivot columns and nullity must agree across primes.
  4. CRT + rational reconstruction over growing prime subsets until stable; then EXACT verification with the
     certified engine: sum_b v_b ad(b) + sum_j S_j d(b'_j) = 0 identically (exact residues recomputed for the
     support), and the upper bound nullity_p >= nullity_Q closes the count: genuine dimension = U, fail-closed.

Usage: fast_kernel3.py W SOL N S [--primes=K] [--workers=n] [--solvers=m]   (m parallel eliminations; ~9 GB each at weight 22)   (record path from FASTK_RECORD, as fast_kernel2)

Distributed phases (v3; one prime per node for the big eliminations):
  --phase=build --dir=DIR        sectors, exact build, normal forms; writes DIR/meta.json, DIR/cols.pkl, DIR/nf_p<k>.pkl
  --phase=solve --dir=DIR --prime-index=k   rref of A-hat mod primes[k] from DIR/nf_p<k>.pkl; writes DIR/sol_p<k>.pkl
  --phase=merge --dir=DIR        consistency, d-parts, CRT + rational reconstruction, exact verification, record
The phases reproduce the single-process path exactly (same pivot structure, same canonical representatives)."""
import hashlib
import heapq
import json
import os
import sys
import time
from fractions import Fraction as F
from math import gcd, isqrt
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("IB_CODE", str(ROOT / "code")))      # IB_CODE: run from outside the repository
import flint                                                   # noqa: E402
from mixed_engine import (cyl_P4, d_mixed_mono, gen_mixed_basis,   # noqa: E402
                          residue_mixed_cur, spin_mixed)


def _more_primes(n):
    ps, c = [], (1 << 61)
    while len(ps) < n:
        c -= 1
        if flint.fmpz(c).is_prime():
            ps.append(c)
    return ps


def _ratrec(a, m):
    a %= m
    bound = isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gcd(r1, s1) != 1:
        return None
    return F(r1, s1) if s1 > 0 else F(-r1, -s1)


def mkey(m):
    """descending multiset of letter indices, then the canonical tuple: a total order under which d is triangular"""
    xs, yp = m
    letters = list(xs) + [a for pr in yp for a in pr]
    return (tuple(sorted(letters, reverse=True)), m)


def triangular(src_basis, tgt_index):
    """pivot structure of d: src_basis -> target.  Returns (pivot: {target row -> (src idx, coef)},
    images: list of {target row: coef}) with the fail-closed distinctness check."""
    pivot = {}
    images = []
    for j, b in enumerate(src_basis):
        img = d_mixed_mono(b)
        rows = {tgt_index[m]: c for m, c in img.items()}
        lead = max(img, key=mkey)
        r = tgt_index[lead]
        if r in pivot:
            raise SystemExit(f"FAIL  triangular structure broken: two d-images share the leading term {lead}")
        pivot[r] = (j, img[lead])
        images.append(rows)
    return pivot, images


def sweep(col, pivot, images, keyrow, primes):
    """reduce a sparse column {row: residue tuple (one per prime)} modulo Im d for all primes at once; returns
    (normal form, factors {src idx: residue tuple}) with col - sum_j f_j d(b_j) = normal form."""
    col = dict(col)
    np_ = len(primes)
    heap = [(keyrow[r], r) for r in col if r in pivot]          # keyrow = rank in DESCENDING key order (0 = largest)
    heapq.heapify(heap)
    factors = {}
    while heap:
        _, r = heapq.heappop(heap)
        v = col.get(r)
        if v is None:
            continue
        j, c = pivot[r]
        f = tuple((v[k] * pow(c, -1, primes[k])) % primes[k] for k in range(np_))
        if j in factors:
            g = factors[j]
            factors[j] = tuple((g[k] + f[k]) % primes[k] for k in range(np_))
        else:
            factors[j] = f
        for rr, cc in images[j].items():
            old = col.get(rr)
            if old is None:
                nv = tuple((-f[k] * cc) % primes[k] for k in range(np_))
                if any(nv):
                    if rr in pivot:
                        heapq.heappush(heap, (keyrow[rr], rr))
                    col[rr] = nv
            else:
                nv = tuple((old[k] - f[k] * cc) % primes[k] for k in range(np_))
                if any(nv):
                    col[rr] = nv
                else:
                    del col[rr]
        assert r not in col
    return col, factors


# ---- pool tasks (fork: module globals are inherited) ----
_G = {}


def _build_columns(idx_list):
    P_low, Nval, cb, tidx, primes = _G["P_low"], _G["Nval"], _G["cb"], _G["tidx"], _G["primes"]
    out = {}
    denom_lcm = 1
    for i in idx_list:
        col = residue_mixed_cur(P_low, {cb[i]: F(1)}, Nval)
        entries = []
        for m, c in col.items():
            if c:
                denom_lcm = denom_lcm * c.denominator // gcd(denom_lcm, c.denominator)
                # a prime dividing the denominator is a BAD prime (excluded by the parent through denom_lcm); store 0 there
                entries.append((tidx[m], tuple(((c.numerator * pow(c.denominator, -1, p)) % p) if c.denominator % p else 0 for p in primes)))
        out[i] = entries
    return out, denom_lcm


def _exact_columns(idx_list):
    P_low, Nval, cb = _G["P_low"], _G["Nval"], _G["cb"]
    return {i: residue_mixed_cur(P_low, {cb[i]: F(1)}, Nval) for i in idx_list}


def _sweep_columns(idx_list):
    cols, piv_t, img_t, keyrow_t, primes = _G["cols"], _G["piv_t"], _G["img_t"], _G["keyrow_t"], _G["primes"]
    out = {}
    for i in idx_list:
        nf, _ = sweep({r: vals for r, vals in cols[i]}, piv_t, img_t, keyrow_t, primes)
        out[i] = nf
    return out


def _solve_prime(pi):
    """rref of A-hat mod primes[pi]: returns (pi, pivot columns, kernel vectors {cand idx: residue})"""
    p = _G["primes"][pi]
    nf_cols, nonpiv_t, nonpiv_c, qrow = _G["nf_cols"], _G["nonpiv_t"], _G["nonpiv_c"], _G["qrow"]
    nq, mq = len(nonpiv_t), len(nonpiv_c)
    qcol = {i: j for j, i in enumerate(nonpiv_c)}
    Ahat = flint.nmod_mat(nq, mq, p)
    for i, nf in nf_cols.items():
        j = qcol[i]
        for r, vals in nf.items():
            if vals[pi]:
                Ahat[qrow[r], j] = vals[pi]
    R, rank = Ahat.rref()
    piv, row = [], 0
    for c in range(mq):
        if row < rank and int(R[row, c]) != 0:
            piv.append(c)
            row += 1
    pset = set(piv)
    kers = []
    for f in range(mq):
        if f in pset:
            continue
        v = {nonpiv_c[f]: 1}
        for r, pc in enumerate(piv):
            e = int(R[r, f])
            if e:
                v[nonpiv_c[pc]] = (-e) % p
        kers.append(v)
    return pi, piv, kers


def fast_genuine_kernel3(P_low, cand_basis, Nval, n_primes=4, workers=14, solvers=4, log=print):
    t0 = time.time()
    s_low = spin_mixed(next(iter(P_low)))
    s_c = spin_mixed(cand_basis[0])
    s_t = s_low + s_c - 1
    tb = gen_mixed_basis(s_t)
    tidx = {m: i for i, m in enumerate(tb)}
    nt = len(tb)
    cb = list(cand_basis)
    nc = len(cb)
    src_t = gen_mixed_basis(s_t - 1)                 # D_t source
    src_c = gen_mixed_basis(s_c - 1)                 # D_c source
    cidx = {m: i for i, m in enumerate(cb)}
    piv_t, img_t = triangular(src_t, tidx)
    piv_c, _ = triangular(src_c, cidx)
    kt, kc = len(src_t), len(src_c)
    nonpiv_t = [r for r in range(nt) if r not in piv_t]
    qrow = {r: i for i, r in enumerate(nonpiv_t)}
    nonpiv_c = [i for i in range(nc) if i not in piv_c]
    order_t = sorted(range(nt), key=lambda i: mkey(tb[i]), reverse=True)
    keyrow_t = {r: k for k, r in enumerate(order_t)}            # rank in descending key order
    log(f"  sectors: cand {nc} (d-image {kc}, non-pivot {len(nonpiv_c)}), target {nt} (d-image {kt}, non-pivot {len(nonpiv_t)}); "
        f"structure {time.time() - t0:.1f}s")
    assert len(nonpiv_t) == nt - kt and len(nonpiv_c) == nc - kc
    primes = _more_primes(n_primes)
    _G.update(P_low=P_low, Nval=Nval, cb=cb, tidx=tidx, primes=primes)
    t1 = time.time()
    chunks = [nonpiv_c[k::workers * 4] for k in range(workers * 4)]
    with Pool(workers) as pool:
        parts = pool.map(_build_columns, chunks)
    cols, dens = {}, 1
    for part, dl in parts:
        cols.update(part)
        dens = dens * dl // gcd(dens, dl)
    log(f"  exact build of {len(nonpiv_c)} non-pivot columns in {time.time() - t1:.0f}s "
        f"({sum(len(v) for v in cols.values())} nonzeros)")
    good = [pi for pi, p in enumerate(primes) if dens % p != 0]
    for pi, p in enumerate(primes):
        if pi not in good:
            log(f"  prime {p}: divides a denominator, skipped")
    t1 = time.time()
    _G.update(cols=cols, piv_t=piv_t, img_t=img_t, keyrow_t=keyrow_t)
    with Pool(workers) as pool:
        parts = pool.map(_sweep_columns, chunks)
    nf_cols = {}
    for part in parts:
        nf_cols.update(part)
    log(f"  normal forms of {len(nf_cols)} columns for {len(primes)} primes at once in {time.time() - t1:.0f}s")
    _G.update(nf_cols=nf_cols, nonpiv_t=nonpiv_t, nonpiv_c=nonpiv_c, qrow=qrow)
    t1 = time.time()
    with Pool(min(solvers, len(good))) as pool:                 # each solver holds its own copy of the tables (~9 GB at weight 22)
        solved = pool.map(_solve_prime, good, chunksize=1)
    log(f"  rref of A-hat ({len(nonpiv_t)} x {len(nonpiv_c)}) for {len(good)} primes in parallel in {time.time() - t1:.0f}s")
    residues = {}
    pivots_ref, U = None, None
    primes_used = []
    for pi, piv, kers in solved:
        p = primes[pi]
        if pivots_ref is None:
            pivots_ref, U = piv, len(kers)
            log(f"  prime {p}: rank {len(piv)}, nullity U = {U}")
        elif piv != pivots_ref or len(kers) != U:
            log(f"  prime {p}: pivot structure differs, skipped (bad prime)")
            continue
        primes_used.append(p)
        for r_i, v in enumerate(kers):
            Av = {}
            for i, vi in v.items():
                for row_, vals in cols[i]:
                    if vals[pi]:
                        nv = (Av.get(row_, 0) + vi * vals[pi]) % p
                        if nv:
                            Av[row_] = nv
                        else:
                            Av.pop(row_, None)
            nf, factors = sweep({r: (val,) for r, val in Av.items()}, piv_t, img_t, keyrow_t, (p,))
            if nf:
                raise SystemExit("FAIL  A v not in Im D_t mod p for a kernel vector")
            vec = dict(v)
            vec.update({nc + j: (-fj[0]) % p for j, fj in factors.items() if fj[0]})
            residues.setdefault(r_i, {})[p] = vec
    for n_use in range(1, len(primes_used) + 1):
        primes_sub = primes_used[:n_use]
        Mprod = 1
        for q in primes_sub:
            Mprod *= q
        lifted, ok = [], True
        for r_i in range(U):
            vec, entries = {}, set()
            for q in primes_sub:
                entries |= set(residues[r_i][q])
            for idx in entries:
                x, mod = 0, 1
                for q in primes_sub:
                    a = residues[r_i][q].get(idx, 0)
                    h = ((a - x) * pow(mod, -1, q)) % q
                    x, mod = x + mod * h, mod * q
                rr = _ratrec(x, Mprod)
                if rr is None:
                    ok = False
                    break
                if rr:
                    vec[idx] = rr
            if not ok:
                break
            lifted.append(vec)
        if not ok:
            log(f"  {n_use} prime(s): reconstruction not yet stable")
            continue
        support = sorted({idx for vec in lifted for idx in vec if idx < nc})
        t3 = time.time()
        chunks_s = [support[k::workers * 2] for k in range(workers * 2)]
        with Pool(workers) as pool:
            parts = pool.map(_exact_columns, chunks_s)
        Aex = {}
        for part in parts:
            Aex.update(part)
        allgood = True
        for vec in lifted:
            acc = {}
            for idx, coef in vec.items():
                src = Aex[idx] if idx < nc else d_mixed_mono(src_t[idx - nc])
                for m, c in src.items():
                    val = acc.get(m, 0) + coef * c
                    if val:
                        acc[m] = val
                    else:
                        acc.pop(m, None)
            if acc:
                allgood = False
                break
        if not allgood:
            log(f"  {n_use} prime(s): lifted vector fails the exact identity ({time.time() - t3:.0f}s), adding primes")
            continue
        Vex = [{cb[i]: c for i, c in vec.items() if i < nc} for vec in lifted]
        # independence: the U vectors are distinct free-column vectors (unit in their own free coordinate), supported on
        # non-pivot candidates, hence independent modulo Im D_c (no nonzero element of Im D_c is non-pivot-supported)
        log(f"  CERTIFIED genuine dim {U} with {n_use} prime(s): {U} vectors verified by exact identities "
            f"({time.time() - t3:.0f}s), canonical non-pivot representatives; total {time.time() - t0:.0f}s")
        return U, Vex, {"kc": kc, "nullity_p": U + kc, "primes": n_use, "full_witnesses": lifted, "nc": nc}
    raise SystemExit("FAIL  could not certify within the prime budget")


# ---- distributed phases (v3) ----
import pickle


def _structure(P_low, cand_basis):
    s_low = spin_mixed(next(iter(P_low)))
    s_c = spin_mixed(cand_basis[0])
    s_t = s_low + s_c - 1
    tb = gen_mixed_basis(s_t)
    tidx = {m: i for i, m in enumerate(tb)}
    cb = list(cand_basis)
    cidx = {m: i for i, m in enumerate(cb)}
    src_t = gen_mixed_basis(s_t - 1)
    src_c = gen_mixed_basis(s_c - 1)
    piv_t, img_t = triangular(src_t, tidx)
    piv_c, _ = triangular(src_c, cidx)
    nt, nc = len(tb), len(cb)
    nonpiv_t = [r for r in range(nt) if r not in piv_t]
    nonpiv_c = [i for i in range(nc) if i not in piv_c]
    order_t = sorted(range(nt), key=lambda i: mkey(tb[i]), reverse=True)
    keyrow_t = {r: k for k, r in enumerate(order_t)}
    return dict(tb=tb, tidx=tidx, cb=cb, src_t=src_t, piv_t=piv_t, img_t=img_t, nt=nt, nc=nc, kt=len(src_t), kc=len(src_c),
                nonpiv_t=nonpiv_t, nonpiv_c=nonpiv_c, keyrow_t=keyrow_t)


def phase_build(W, sol, N, s, n_primes, workers, outdir, log=print):
    t0 = time.time()
    P_low, cb = cyl_P4(sol, N, s), gen_mixed_basis(W)
    S = _structure(P_low, cb)
    log(f"  sectors: cand {S['nc']} (d-image {S['kc']}, non-pivot {len(S['nonpiv_c'])}), target {S['nt']} (d-image {S['kt']}, non-pivot {len(S['nonpiv_t'])}); structure {time.time() - t0:.1f}s")
    primes = _more_primes(n_primes)
    _G.update(P_low=P_low, Nval=N, cb=S['cb'], tidx=S['tidx'], primes=primes)
    t1 = time.time()
    chunks = [S['nonpiv_c'][k::workers * 4] for k in range(workers * 4)]
    with Pool(workers) as pool:
        parts = pool.map(_build_columns, chunks)
    cols, dens = {}, 1
    for part, dl in parts:
        cols.update(part)
        dens = dens * dl // gcd(dens, dl)
    log(f"  exact build of {len(S['nonpiv_c'])} non-pivot columns in {time.time() - t1:.0f}s ({sum(len(v) for v in cols.values())} nonzeros)")
    good = [pi for pi, p in enumerate(primes) if dens % p != 0]
    t1 = time.time()
    _G.update(cols=cols, piv_t=S['piv_t'], img_t=S['img_t'], keyrow_t=S['keyrow_t'])
    with Pool(workers) as pool:
        parts = pool.map(_sweep_columns, chunks)
    nf_cols = {}
    for part in parts:
        nf_cols.update(part)
    log(f"  normal forms of {len(nf_cols)} columns for {len(primes)} primes at once in {time.time() - t1:.0f}s")
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    qrow = {r: i for i, r in enumerate(S['nonpiv_t'])}
    for pi in good:
        rows, colsq, vals = [], [], []
        # compact per-prime triplets (quotient row, non-pivot column position, value)
        qcol = {i: j for j, i in enumerate(S['nonpiv_c'])}
        for i, nf in nf_cols.items():
            j = qcol[i]
            for r, v in nf.items():
                if v[pi]:
                    rows.append(qrow[r]); colsq.append(j); vals.append(v[pi])
        with open(out / f"nf_p{pi}.pkl", "wb") as fh:
            pickle.dump((primes[pi], len(S['nonpiv_t']), len(S['nonpiv_c']), rows, colsq, vals), fh, protocol=pickle.HIGHEST_PROTOCOL)
        log(f"  wrote nf_p{pi}.pkl ({len(vals)} nonzeros)")
    with open(out / "cols.pkl", "wb") as fh:
        pickle.dump(cols, fh, protocol=pickle.HIGHEST_PROTOCOL)
    meta = {"W": W, "sol": sol, "N": str(N), "s": str(s), "primes": primes, "good": good, "nonpiv_t": S['nonpiv_t'], "nonpiv_c": S['nonpiv_c'],
            "nc": S['nc'], "kc": S['kc'], "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "build_seconds": round(time.time() - t0)}
    (out / "meta.json").write_text(json.dumps(meta) + "\n")
    log(f"  build phase done in {time.time() - t0:.0f}s -> {out}")


def phase_solve(outdir, pi, log=print):
    t0 = time.time()
    with open(Path(outdir) / f"nf_p{pi}.pkl", "rb") as fh:
        p, nq, mq, rows, colsq, vals = pickle.load(fh)
    Ahat = flint.nmod_mat(nq, mq, p)
    for r, j, v in zip(rows, colsq, vals):
        Ahat[r, j] = v
    log(f"  prime index {pi} ({p}): matrix {nq} x {mq} filled ({len(vals)} nonzeros) in {time.time() - t0:.0f}s; rref ...")
    t1 = time.time()
    R, rank = Ahat.rref()
    piv, row = [], 0
    for c in range(mq):
        if row < rank and int(R[row, c]) != 0:
            piv.append(c)
            row += 1
    pset = set(piv)
    kers = []
    for f in range(mq):
        if f in pset:
            continue
        v = {f: 1}                                   # in non-pivot column POSITIONS; mapped to candidate indices at merge
        for r, pc in enumerate(piv):
            e = int(R[r, f])
            if e:
                v[pc] = (-e) % p
        kers.append(v)
    meta = json.loads((Path(outdir) / "meta.json").read_text())
    binding = {"solver": "fast_kernel3 phase_solve (dense rref)", "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "build_tool_sha256": meta.get("tool_sha256"), "nq": nq, "mq": mq, "nnz": len(vals),
               "input_files": [f"nf_p{pi}.pkl"], "input_sha256": hashlib.sha256((Path(outdir) / f"nf_p{pi}.pkl").read_bytes()).hexdigest()}
    with open(Path(outdir) / f"sol_p{pi}.pkl", "wb") as fh:
        pickle.dump((pi, p, piv, kers, binding), fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"  prime index {pi}: rank {rank}, nullity {len(kers)}; rref {time.time() - t1:.0f}s; total {time.time() - t0:.0f}s -> sol_p{pi}.pkl")


def _solution_record_defect(rec, pi, primes, mq, meta):
    """Round-29 (Codex) High: the merge used to take the nullity from the first record's vector count and certify a dimension from
    exact residue identities alone -- a duplicated vector certified dim 2, an empty list dim 0.  A record is now accepted only if it
    is (pi, p, piv, kers[, binding]) with: the prime index and value those of the build; pivots a sorted list of distinct column
    positions in range; len(kers) == mq - len(piv) (the modular nullity, an upper bound on the rational one: the certification is
    'U independent vectors verified exactly' == 'U = modular nullity'); the vectors the CANONICAL free-coordinate basis (the i-th
    vector has residue 1 at the i-th non-pivot column, no other non-pivot entry, pivot entries in [0, p)) -- which makes them
    independent by construction; and, when the optional binding is present, the build's tool hash and the matrix sizes.  Returns
    None or the defect (a string).  Explicit checks, not asserts: they survive python -O."""
    if not isinstance(rec, tuple) or len(rec) not in (4, 5):
        return "not a (pi, p, piv, kers[, binding]) tuple"
    rpi, p, piv, kers = rec[:4]
    if not (isinstance(rpi, int) and rpi == pi):
        return f"prime index {rpi!r} is not the record's index {pi}"
    if not (isinstance(p, int) and 0 <= pi < len(primes) and primes[pi] == p):
        return f"prime {p} is not primes[{pi}] = {primes[pi] if 0 <= pi < len(primes) else None} of the build"
    if not (isinstance(piv, list) and all(isinstance(c, int) and 0 <= c < mq for c in piv) and piv == sorted(set(piv))):
        return "pivot list is not a sorted list of distinct column positions in range"
    if not isinstance(kers, list) or len(kers) != mq - len(piv):
        return f"vector count {len(kers) if isinstance(kers, list) else '?'} != columns {mq} - rank {len(piv)} (modular nullity)"
    pset = set(piv)
    free = [c for c in range(mq) if c not in pset]
    for i, (v, fcol) in enumerate(zip(kers, free)):
        if not isinstance(v, dict) or v.get(fcol) != 1:
            return f"vector {i} is not the canonical representative of free column {fcol} (entry 1 there)"
        for c, val in v.items():
            if c != fcol and c not in pset:
                return f"vector {i} has an entry at another free column {c}: not the canonical basis"
            if not (isinstance(val, int) and 0 <= val < p):
                return f"vector {i}: residue {val!r} at column {c} is not in [0, {p})"
    if len(rec) == 5:
        b = rec[4]
        if not isinstance(b, dict):
            return "binding is not a dict"
        if b.get("mq") not in (None, mq):
            return f"binding column count {b.get('mq')} != {mq}"
        if b.get("nq") not in (None, len(meta["nonpiv_t"])):
            return f"binding row count {b.get('nq')} != {len(meta['nonpiv_t'])}"
        if b.get("build_tool_sha256") not in (None, meta.get("tool_sha256")):
            return "binding names a different build tool"
    return None


def _load_triplets(out, pi, meta):
    """(p, nq, mq, rows, cols, vals) of A-hat mod primes[pi] from the build directory: the v3 pickle or the v4 chunk files; None if absent."""
    f = out / f"nf_p{pi}.pkl"
    if f.exists():
        with open(f, "rb") as fh:
            p, nq, mq, rows, colsq, vals = pickle.load(fh)
        return p, nq, mq, list(rows), list(colsq), list(vals)
    if meta.get("chunked"):
        import glob
        files = sorted(glob.glob(str(out / f"nf_p{pi}_chunk*.npz")), key=lambda x: int(x.rsplit("chunk", 1)[1].split(".")[0]))
        if len(files) != meta.get("chunks"):
            return None
        import numpy as np
        rows, colsq, vals = [], [], []
        for x in files:
            q = np.load(x); rows += q["rows"].tolist(); colsq += q["cols"].tolist(); vals += q["vals"].tolist()
        return meta["primes"][pi], len(meta["nonpiv_t"]), len(meta["nonpiv_c"]), rows, colsq, vals
    return None


def _provenance_check(rec, pi, out, log):
    """Re-hash the record's input files when they are present.  True = verified, False = MISMATCH (fatal), None = not verifiable."""
    if len(rec) < 5 or not isinstance(rec[4], dict) or "input_sha256" not in rec[4]:
        log(f"  record {pi}: legacy record without input provenance (accepted; upper bound asserted)")
        return None
    b = rec[4]; files = b.get("input_files") or []
    # the claimed files must be THIS prime's matrix files (Kimi ce0885e, C9 side finding: a record could otherwise name another
    # prime's file and pass the hash check with it -- provenance alone did not bind a record to its own prime's matrix)
    import re as _re
    own = _re.compile(rf"^nf_p{pi}(\.pkl|_chunk\d+\.npz)$")
    if files and not all(own.match(x) for x in files):
        raise SystemExit(f"FAIL  solve record {pi}: input provenance names files {files} that are not prime index {pi}'s matrix files")
    if not files or not all((out / x).exists() for x in files):
        log(f"  record {pi}: input files not present in the build directory; provenance not re-verified")
        return None
    h = hashlib.sha256()
    for x in files:
        h.update((out / x).read_bytes())
    if h.hexdigest() != b["input_sha256"]:
        raise SystemExit(f"FAIL  solve record {pi}: input provenance MISMATCH (the record was not computed from these matrix files)")
    return True


def _rank_certificate(rec, pi, out, meta, log):
    """Independent upper bound on the rational nullity: recompute the rank of the PIVOT columns of A-hat mod p (they must be
    independent, so rank_p >= len(piv) and nullity_p <= mq - len(piv) = U; nullity_Q <= nullity_p).  Done when the matrix files
    are present and the pivot submatrix is within the replay budget (FK3_VERIFY_RANK_LIMIT entries, default 3e8; FK3_VERIFY_RANK=0
    disables, =1 forces).  True = verified, False = the claim is FALSE, None = not replayed (asserted)."""
    mode = os.environ.get("FK3_VERIFY_RANK", "auto")
    if mode == "0":
        return None
    data = _load_triplets(out, pi, meta)
    if data is None:
        log(f"  record {pi}: matrix files absent; modular rank not re-verified (asserted by the solve record)")
        return None
    p, nq, mq, rows, colsq, vals = data
    piv = rec[2]; r = len(piv)
    limit = float(os.environ.get("FK3_VERIFY_RANK_LIMIT", "3e8"))
    if mode != "1" and nq * r > limit:
        log(f"  record {pi}: pivot submatrix {nq} x {r} above the replay budget; modular rank not re-verified (asserted)")
        return None
    t1 = time.time()
    pos = {c: j for j, c in enumerate(piv)}
    M = flint.nmod_mat(nq, r, p)
    for a, b, v in zip(rows, colsq, vals):
        j = pos.get(b)
        if j is not None:
            M[a, j] = v
    rk = M.rank()
    log(f"  record {pi}: rank of the {r} claimed pivot columns mod {p} recomputed = {rk} ({time.time() - t1:.0f}s): "
        f"{'INDEPENDENT -- the modular nullity is at most ' + str(mq - r) if rk == r else 'NOT independent: the claimed rank is false'}")
    return rk == r


def phase_merge(outdir, workers, log=print):
    t0 = time.time()
    out = Path(outdir)
    meta = json.loads((out / "meta.json").read_text())
    W, sol, N, s = meta["W"], meta["sol"], F(meta["N"]), F(meta["s"])
    primes, nonpiv_c = meta["primes"], meta["nonpiv_c"]
    P_low, cb = cyl_P4(sol, N, s), gen_mixed_basis(W)
    S = _structure(P_low, cb)
    if not (S["nonpiv_c"] == nonpiv_c and S["nonpiv_t"] == meta["nonpiv_t"]):
        raise SystemExit("FAIL  structure differs from the build (meta.json does not describe this point)")
    nc, kc, src_t = S["nc"], S["kc"], S["src_t"]
    with open(out / "cols.pkl", "rb") as fh:
        cols = pickle.load(fh)
    mq = len(nonpiv_c)
    sols = []; attest = []
    for pi in meta["good"]:
        f = out / f"sol_p{pi}.pkl"
        if f.exists():
            with open(f, "rb") as fh:
                rec = pickle.load(fh)
            why = _solution_record_defect(rec, pi, primes, mq, meta)
            if why:
                raise SystemExit(f"FAIL  solve record {f.name} rejected: {why}")
            prov = _provenance_check(rec, pi, out, log)
            rankv = _rank_certificate(rec, pi, out, meta, log)
            if rankv is False:
                raise SystemExit(f"FAIL  solve record {f.name} rejected: its pivot columns are NOT independent modulo its prime "
                                 f"(the claimed rank is false)")
            sols.append(rec[:4]); attest.append((pi, prov, rankv))
    log(f"  merge: {len(sols)} solved prime(s) of {len(meta['good'])} good, every record validated (schema, prime binding, "
        f"nullity = columns - rank, canonical free-coordinate basis); structure + load {time.time() - t0:.0f}s")
    n_rank_verified = sum(1 for _, _, rv in attest if rv is True); n_prov = sum(1 for _, pv, _ in attest if pv is True)
    log(f"  modular rank independently re-verified for {n_rank_verified} of {len(sols)} record(s); input provenance re-verified for "
        f"{n_prov} of {len(sols)} (a record without an independently verified rank contributes an ASSERTED upper bound only)")
    # The reference pivot structure is taken from a record whose rank was independently verified (the largest verified rank
    # first), otherwise from the record with the largest claimed rank (a bad prime has a SMALLER rank; a forged record may claim a
    # larger one).  A verified rank r of some record bounds the rational nullity by mq - r whether or not that record is used
    # below; the certified dimension counts as independently verified only if the REFERENCE record's rank was verified, i.e.
    # the verified upper bound equals the dimension (Codex bf7c317 follow-up, finding 2: previously any verified record, even
    # one skipped as a "bad prime", labelled the dimension of the selected structure as verified).
    verified = {pi: (rv is True) for pi, _, rv in attest}
    sols.sort(key=lambda rec: (not verified[rec[0]], -len(rec[2]), rec[0]))
    verified_upper = min([mq - len(rec[2]) for rec in sols if verified[rec[0]]], default=None)
    residues, pivots_ref, U, primes_used, ref_verified, ref_p, skipped = {}, None, None, [], False, None, []
    for pi, p, piv, kers in sols:
        if pivots_ref is None:
            pivots_ref, U, ref_verified, ref_p = piv, len(kers), verified[pi], p
            log(f"  prime {p}: rank {len(piv)} ({'independently verified' if ref_verified else 'ASSERTED'}), nullity U = {U} -- reference structure")
        elif piv != pivots_ref or len(kers) != U:
            log(f"  prime {p}: pivot structure differs (claimed rank {len(piv)}, nullity {len(kers)}, rank "
                f"{'independently verified' if verified[pi] else 'asserted'}), skipped (bad prime)")
            skipped.append({"prime": p, "prime_index": pi, "claimed_rank": len(piv), "claimed_nullity": len(kers), "rank_verified": verified[pi]})
            continue
        primes_used.append(p)
        for r_i, vpos in enumerate(kers):
            v = {nonpiv_c[j]: val for j, val in vpos.items()}
            Av = {}
            for i, vi in v.items():
                for row_, vals in cols[i]:
                    if vals[pi]:
                        nv = (Av.get(row_, 0) + vi * vals[pi]) % p
                        if nv:
                            Av[row_] = nv
                        else:
                            Av.pop(row_, None)
            nf, factors = sweep({r: (val,) for r, val in Av.items()}, S["piv_t"], S["img_t"], S["keyrow_t"], (p,))
            if nf:
                raise SystemExit("FAIL  A v not in Im D_t mod p for a kernel vector")
            vec = dict(v)
            vec.update({nc + j: (-fj[0]) % p for j, fj in factors.items() if fj[0]})
            residues.setdefault(r_i, {})[p] = vec
    _G.update(P_low=P_low, Nval=N, cb=S["cb"])
    for n_use in range(1, len(primes_used) + 1):
        primes_sub = primes_used[:n_use]
        Mprod = 1
        for q in primes_sub:
            Mprod *= q
        lifted, ok = [], True
        for r_i in range(U):
            vec, entries = {}, set()
            for q in primes_sub:
                entries |= set(residues[r_i][q])
            for idx in entries:
                x, mod = 0, 1
                for q in primes_sub:
                    a = residues[r_i][q].get(idx, 0)
                    h = ((a - x) * pow(mod, -1, q)) % q
                    x, mod = x + mod * h, mod * q
                rr = _ratrec(x, Mprod)
                if rr is None:
                    ok = False
                    break
                if rr:
                    vec[idx] = rr
            if not ok:
                break
            lifted.append(vec)
        if not ok:
            log(f"  {n_use} prime(s): reconstruction not yet stable")
            continue
        support = sorted({idx for vec in lifted for idx in vec if idx < nc})
        t3 = time.time()
        chunks_s = [support[k::workers * 2] for k in range(workers * 2)]
        with Pool(workers) as pool:
            parts = pool.map(_exact_columns, chunks_s)
        Aex = {}
        for part in parts:
            Aex.update(part)
        allgood = True
        for vec in lifted:
            acc = {}
            for idx, coef in vec.items():
                src = Aex[idx] if idx < nc else d_mixed_mono(src_t[idx - nc])
                for m, c in src.items():
                    val = acc.get(m, 0) + coef * c
                    if val:
                        acc[m] = val
                    else:
                        acc.pop(m, None)
            if acc:
                allgood = False
                break
        if not allgood:
            log(f"  {n_use} prime(s): lifted vector fails the exact identity ({time.time() - t3:.0f}s), adding primes")
            continue
        if verified_upper is not None and verified_upper < U:
            raise SystemExit(f"FAIL  {U} exact independent witnesses, but an independently verified modular rank bounds the rational "
                             f"nullity by {verified_upper}: the merge is inconsistent (report this)")
        # a skipped record claiming a rank ABOVE the certified rank mq - U is refuted: rank mod p never exceeds the rational rank,
        # and U exact independent witnesses prove the rational rank is at most mq - U
        refuted = [sk for sk in skipped if sk["claimed_rank"] > mq - U]
        conflict = [sk for sk in skipped if sk["claimed_rank"] != mq - U]
        if conflict and not ref_verified and os.environ.get("FK3_ALLOW_CONFLICT") != "1":
            raise SystemExit(f"FAIL  conflicting solve records: the reference prime {ref_p} claims rank {mq - U} (not independently verified) and "
                             f"prime(s) {[sk['prime'] for sk in conflict]} claim rank(s) {[sk['claimed_rank'] for sk in conflict]}; verify the "
                             f"reference rank (matrix files present, FK3_VERIFY_RANK=1) or solve more primes; FK3_ALLOW_CONFLICT=1 "
                             f"certifies {U} with an ASSERTED upper bound naming the conflict")
        rank_status = ("modular rank independently re-verified (reference record)" if ref_verified else
                       "modular rank ASSERTED by the solve records, not independently verified"
                       + (f"; best independently verified upper bound {verified_upper}" if verified_upper is not None else "")
                       + (f"; CONFLICTING skipped record(s) {[sk['prime'] for sk in conflict]} with claimed rank(s) "
                          f"{[sk['claimed_rank'] for sk in conflict]}" if conflict else ""))
        log(f"  dimension {U}: lower bound {U} by exact independent witnesses (independent); upper bound {U} by the modular nullity -- {rank_status}")
        if refuted:
            log(f"  WARNING  skipped record(s) {[sk['prime'] for sk in refuted]} claim rank(s) {[sk['claimed_rank'] for sk in refuted]} above the "
                f"certified rank {mq - U}: REFUTED by the exact witnesses (a corrupt or forged solve record); the record is flagged and the exit code is 3")
        # Kimi ce0885e (C9, F4): with the rank replay skipped and no conflicting second record, a lone solve record that over-claims
        # its rank (nullity too small) with genuine kernel vectors is certified at the under-stated dimension; the exact witnesses
        # are a true LOWER bound, the upper bound is the record's claim.  The headline says which of the two the dimension is.
        dim_status = (f"dimension {U} VERIFIED (lower bound by exact witnesses, upper bound by an independently re-verified modular rank)"
                      if ref_verified else
                      f"dimension >= {U} VERIFIED by exact witnesses; dimension = {U} ASSERTED (the upper bound is the solve records' claim, "
                      f"not independently re-verified)")
        log(f"  CERTIFIED genuine dim {U} with {n_use} prime(s): {U} vectors verified by exact identities ({time.time() - t3:.0f}s), canonical non-pivot representatives; merge {time.time() - t0:.0f}s")
        log(f"  {dim_status}")
        record_summaries = [{"prime_index": pi_, "prime": p_, "claimed_rank": len(piv_), "claimed_nullity": len(kers_), "rank_verified": verified[pi_],
                             "provenance_verified": next((pv is True) for q, pv, _ in attest if q == pi_), "used": p_ in primes_used}
                            for pi_, p_, piv_, kers_ in sols]
        return U, lifted, {"kc": kc, "nullity_p": U + kc, "primes": n_use, "full_witnesses": lifted, "nc": nc, "meta": meta,
                           "rank_verified_records": n_rank_verified, "provenance_verified_records": n_prov, "records": len(sols),
                           "upper_bound": rank_status, "reference_rank_verified": ref_verified, "verified_upper_bound": verified_upper,
                           "skipped_records": skipped, "refuted_records": refuted, "dimension_status": dim_status,
                           "primes_used": primes_used[:n_use], "solve_records": record_summaries}
    raise SystemExit(f"FAIL  could not certify with the {len(primes_used)} solved prime(s); solve more primes and merge again")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    phase = opts.get("--phase")
    if phase == "solve":
        phase_solve(opts["--dir"], int(opts["--prime-index"]))
        sys.exit(0)
    if phase == "merge":
        t0 = time.time()
        g, lifted, info = phase_merge(opts["--dir"], workers=int(opts.get("--workers", 14)))
        meta = info["meta"]; W, sol, N, s = meta["W"], meta["sol"], F(meta["N"]), F(meta["s"])
        print(f"sol={sol} N={N} s={s} dim={g} kc={info['kc']} primes={info['primes']} time={time.time() - t0:.0f}s (merge; build {meta['build_seconds']}s)", flush=True)
        out = os.environ.get("FASTK_RECORD")
        if out:
            tool_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            rec = {"tool": "lab/fast_kernel3.py (quotient formulation, v3 distributed phases)", "tool_sha256": tool_sha, "build_tool_sha256": meta["tool_sha256"],
                   "weight": W, "solution": sol, "N": str(N), "s": str(s), "genuine_dim": g, "kc": info["kc"],
                   "primes_used": info["primes"], "seconds": round(time.time() - t0), "build_seconds": meta["build_seconds"],
                   "witnesses": [{"charge_v": {str(k): str(c) for k, c in w.items() if k < info["nc"]},
                                  "d_coefficients_S": {str(k - info["nc"]): str(c) for k, c in w.items() if k >= info["nc"]}}
                                 for w in lifted],
                   "attestation": {"lower_bound": f"{g} exact independent witnesses (verified by exact identities in the merge)",
                                   "dimension_status": info.get("dimension_status"),
                                   "primes": info.get("primes_used"), "solve_records": info.get("solve_records"),
                                   "upper_bound": info.get("upper_bound"), "reference_rank_verified": info.get("reference_rank_verified"),
                                   "verified_upper_bound": info.get("verified_upper_bound"), "rank_verified_records": info.get("rank_verified_records"),
                                   "provenance_verified_records": info.get("provenance_verified_records"), "records": info.get("records"),
                                   "skipped_records": info.get("skipped_records"), "refuted_records": info.get("refuted_records")},
                   "status": "lab record: exact fast-path computation (quotient formulation, distributed phases), pending "
                             "independent certification by the registered pure-Python certificate"}
            Path(out).write_text(json.dumps(rec, indent=1) + "\n")
            print(f"record written: {out} (tool sha256 {tool_sha[:16]})")
        if info.get("refuted_records"):
            print(f"WARNING  refuted solve record(s) {[sk['prime'] for sk in info['refuted_records']]}: the dimension {g} is certified, but the record set is corrupt", flush=True)
            sys.exit(3)
        sys.exit(0)
    W, sol, N, s = int(args[0]), int(args[1]), F(args[2]), F(args[3])
    if s * s != (N - 25) * (N - 1):
        sys.exit("point not on the curve")
    if phase == "build":
        phase_build(W, sol, N, s, int(opts.get("--primes", 8)), int(opts.get("--workers", 14)), opts["--dir"])
        sys.exit(0)
    t0 = time.time()
    g, V, info = fast_genuine_kernel3(cyl_P4(sol, N, s), gen_mixed_basis(W), N,
                                      n_primes=int(opts.get("--primes", 4)), workers=int(opts.get("--workers", 14)),
                                      solvers=int(opts.get("--solvers", 4)))
    print(f"sol={sol} N={N} s={s} dim={g} kc={info['kc']} primes={info['primes']} time={time.time() - t0:.0f}s", flush=True)
    out = os.environ.get("FASTK_RECORD")
    if out:
        tool_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        rec = {"tool": "lab/fast_kernel3.py (quotient formulation, v2)", "tool_sha256": tool_sha, "weight": W, "solution": sol,
               "N": str(N), "s": str(s), "genuine_dim": g, "kc": info["kc"],
               "primes_used": info["primes"], "seconds": round(time.time() - t0),
               "witnesses": [{"charge_v": {str(k): str(c) for k, c in w.items() if k < info["nc"]},
                              "d_coefficients_S": {str(k - info["nc"]): str(c) for k, c in w.items() if k >= info["nc"]}}
                             for w in info["full_witnesses"]],
               "status": "lab record: exact fast-path computation (quotient formulation), pending "
                         "independent certification by the registered pure-Python certificate"}
        Path(out).write_text(json.dumps(rec, indent=1) + "\n")
        print(f"record written: {out} (tool sha256 {tool_sha[:16]})")
