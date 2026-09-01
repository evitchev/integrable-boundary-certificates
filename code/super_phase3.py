"""SUSY Phase 3 remainder (preregistered: results/super_phase3_preregistration.md).

S1  vacancy in the enlarged spaces at general lambda: the ad_{I3(lam)}
    centralizer spectrum (classes mod d) per doubled weight w2 = 1..MAX
    for lam in LAMBDAS, pinned to the free-point pattern (dim 1 at
    sigma in {0, 1/2, 1, 3, 5, 7, 9, 11}, 0 elsewhere in range);
S2  supermultiplet pairing: I5, I7, I9 unique as joint (ad_G, ad_I3)
    kernels, all direct brackets total derivatives, and NO fermionic
    partners (every odd w2 >= 5 vacant): a tower of G-singlets;
S3  super corank-character instance at lam = 3/4: joint kernels of the
    fermionic screenings psi e^{alpha X} for the dual pair F+ = {1/2, -2}
    and for all four roots {2, -1/2, 1/2, -2}, per w2 <= 16, against the
    centralizer spectrum -- the four roots pinned EQUAL at every
    3 <= w2 <= 16 (momentum exception oint dX pinned at sigma = 0); the
    dual pair alone pinned = NS super-Virasoro vacuum class count (the
    preregistered 'equal at the odd-spin slots' clause was refuted in
    the smoke run and is pinned as a record).

Everything exact (Fractions / integer sparse elimination); fail-closed:
a failed pin prints FINDING with the observed values and exits 1.
Usage: super_phase3.py [MAX_W2]   (default 24).  Only MAX_W2 == 24 is a
certification: the artifact results/super_phase3.json is written and
the CERTIFIED marker printed only then; smaller values are smoke runs
(no artifact, marker 'SMOKE RUN').  Finite-range statements only: five
sampled lambda, doubled weight <= 24."""
import json
import os
import sys
import time
from fractions import Fraction as F

import sympy as sp

from sparse_linalg import SparsePivots, _to_int_vec
from super_engine import (SuperSector, cur_add, d_cur, d_mono,
                          gen_basis_super, genuine_kernel_super,
                          mono_weight2, normal_prod, residue_cur)
from super_kdv import G, NF, primitive, sig
from super_lambda import G_l, T_l, to_fraction_cur
from super_screening_fit import residue_screen_super, twisted_d_super

MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
if MAX < 12 or MAX > 24:
    # Codex round 12: tiny ranges made every range-limited pin vacuous
    # ("every 3 <= w2 <= 2 ... all pins met").  Smoke runs need the
    # spin-5 tower slot (w2 = 12) at least; 24 is the certified range.
    print(f"SUPER PHASE 3: MAX_W2 must be between 12 and 24 (got {MAX}); nothing checked")
    sys.exit(2)
CERTIFY = (MAX == 24)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, 'results', 'super_phase3.json')
LAMBDAS = [F(0), F(3, 4), F(1, 2), F(1), F(2)]
def sigma(w2):
    return F(w2 - 2, 2)


EXPECT_SPIN_DIMS = {w2: (1 if sigma(w2) in (0, F(1, 2), 1, 3, 5, 7, 9, 11) else 0)
                    for w2 in range(1, MAX + 1)}


def require(cond, msg, findings):
    print(f"{'PASS' if cond else 'FAIL'}  {msg}", flush=True)
    if not cond:
        findings.append(msg)
    return cond


def I3_of(lam):
    Tl = to_fraction_cur(T_l(sp.Rational(lam.numerator, lam.denominator)))
    Gl = to_fraction_cur(G_l(sp.Rational(lam.numerator, lam.denominator)))
    cur = dict(normal_prod(Tl, Tl))
    for m, c in normal_prod(Gl, d_cur(Gl)).items():
        cur[m] = cur.get(m, 0) + F(-1, 4) * c
    return primitive({m: c for m, c in cur.items() if c})


def total_derivative(S):
    return (not S) or SuperSector(mono_weight2(next(iter(S))), NF).is_total_derivative(S)


def spectrum(P, w2s):
    return {w2: len(genuine_kernel_super(
        P, [{b: F(1)} for b in gen_basis_super(w2, NF)], NF)) for w2 in w2s}


def twisted_kernel(cands, alpha):
    """Genuine classes in span(cands) whose residue against
    V_alpha = psi e^{alpha X} is D_alpha-exact (D_alpha J = dJ + alpha dX J).
    Mirrors genuine_kernel_super with the twisted image as target span;
    every emitted word has doubled weight w2(P) - 1."""
    w2c = mono_weight2(next(iter(cands[0])))
    nc = len(cands)
    w2t = w2c - 1
    tgt_index = {m: i for i, m in enumerate(gen_basis_super(w2t, NF))}
    apiv = SparsePivots()
    if w2t >= 2:
        # gen_basis_super(0) is empty: the constant's twisted image
        # D_alpha 1 = alpha dX must be inserted by hand at w2t == 2.
        lower = gen_basis_super(w2t - 2, NF) if w2t > 2 else [()]
        for b in lower:
            iv, _ = _to_int_vec(twisted_d_super(b, alpha), tgt_index)
            col = {k + nc: c for k, c in iv.items()}
            if col:
                apiv.insert(col)
    tags, scales = [], []
    for j, cand in enumerate(cands):
        iv, L = _to_int_vec(residue_screen_super(cand, alpha), tgt_index)
        scales.append(L)
        aug = {k + nc: c for k, c in iv.items()}
        aug[j] = 1
        r = apiv.reduce(aug)
        if r and max(r) < nc:
            tags.append(r)
        elif r:
            apiv.rows[max(r)] = r
    cs_index = {m: i for i, m in enumerate(gen_basis_super(w2c, NF))}
    span = SparsePivots()
    for b in gen_basis_super(w2c - 2, NF):
        iv, _ = _to_int_vec({m: F(c) for m, c in d_mono(b).items()}, cs_index)
        if iv:
            span.insert(iv)
    genuine = []
    for t in tags:
        cur = {}
        for j, c in t.items():
            cur_add(cur, cands[j], F(c) * scales[j])
        cur = {m: c for m, c in cur.items() if c}
        iv, _ = _to_int_vec(cur, cs_index)
        if span.insert(iv):
            genuine.append(cur)
    return genuine


def joint_screening_dim(w2, alphas):
    kernel = [{b: F(1)} for b in gen_basis_super(w2, NF)]
    for a in alphas:
        if not kernel:
            break
        kernel = twisted_kernel(kernel, a)
    return len(kernel)


if __name__ == "__main__":
    t0 = time.time()
    findings, out = [], {"max_w2": MAX, "lambdas": [str(l) for l in LAMBDAS],
                         "spectra": {}, "tower": {}, "screening": {}}
    w2s = list(range(1, MAX + 1))
    towers = {}
    for lam in LAMBDAS:
        I3 = I3_of(lam)
        require(total_derivative(residue_cur(G, I3)),
                f"lam={lam}: I3(lam) is oint G-invariant", findings)
        spec = spectrum(I3, w2s)
        out["spectra"][str(lam)] = spec
        bad = {w2: d for w2, d in spec.items() if d != EXPECT_SPIN_DIMS[w2]}
        require(not bad, f"lam={lam}: S1 spectrum equals the free-point pattern "
                f"(FINDING if not: {bad})", findings)
        # S2: tower by joint (ad_G, ad_I3) kernel; brackets; no fermionic partners
        tower = {8: I3}
        for w2 in (12, 16, 20):
            if w2 > MAX:
                continue
            greps = genuine_kernel_super(G, [{b: F(1)} for b in gen_basis_super(w2, NF)], NF)
            joint = genuine_kernel_super(I3, greps, NF)
            require(len(joint) == 1, f"lam={lam}: joint (ad_G, ad_I3) kernel at w2={w2} "
                    f"(sigma={sig(w2)}) has dim 1 (got {len(joint)})", findings)
            if joint:
                tower[w2] = primitive(joint[0])
        for (a, b) in ((12, 16), (8, 12), (8, 16), (12, 20), (16, 20)):
            if a in tower and b in tower:
                require(total_derivative(residue_cur(tower[a], tower[b])),
                        f"lam={lam}: [I{sig(a)}, I{sig(b)}] is a total derivative", findings)
        for w2 in tower:
            if w2 != 8:
                require(total_derivative(residue_cur(G, tower[w2])),
                        f"lam={lam}: ad_(oint G) I{sig(w2)}(lam) = 0 (G-singlet)", findings)
        partners = {w2: spec[w2] for w2 in w2s if w2 % 2 == 1 and w2 >= 5 and spec[w2]}
        require(not partners, f"lam={lam}: no fermionic partners -- every half-integer "
                f"spin >= 3/2 vacant (FINDING if not: {partners})", findings)
        towers[str(lam)] = {str(w2): len(v) for w2, v in tower.items()}
    out["tower"] = towers
    # S3 at lam = 3/4
    spec34 = out["spectra"]["3/4"]
    pair = (F(1, 2), F(-2))
    four = (F(2), F(-1, 2), F(1, 2), F(-2))
    w2s3 = [w2 for w2 in range(2, min(MAX, 16) + 1)]
    dpair = {w2: joint_screening_dim(w2, pair) for w2 in w2s3}
    dfour = {w2: joint_screening_dim(w2, four) for w2 in w2s3}
    out["screening"] = {"lambda": "3/4", "dual_pair": {str(k): v for k, v in dpair.items()},
                        "all_four": {str(k): v for k, v in dfour.items()},
                        "centralizer": {str(k): spec34[k] for k in w2s3}}
    print("S3 at lam=3/4 (w2: pair / four / centralizer):", flush=True)
    for w2 in w2s3:
        print(f"  w2={w2:2d} sigma={sig(w2):>4}: {dpair[w2]} / {dfour[w2]} / {spec34[w2]}", flush=True)
    require(dfour[2] == 0 and dpair[2] == 0 and spec34[2] == 1,
            "S3: momentum exception at sigma=0 -- oint dX commutes with I3 (dim 1) but is "
            f"screened by no root (dims {dpair[2]}, {dfour[2]}): it carries screening charge",
            findings)
    bad4 = {w2: (dfour[w2], spec34[w2]) for w2 in w2s3 if w2 >= 3 and dfour[w2] != spec34[w2]}
    require(not bad4, f"S3: joint kernel of all four screenings equals the centralizer at every "
            f"3 <= w2 <= {w2s3[-1]} (FINDING if not: {bad4})", findings)
    # Preregistered H_S3(a) ("dual pair equal to the centralizer at the
    # odd-spin slots") was REFUTED in the w2 <= 12 smoke run: the pair's
    # commutant is larger (it contains the super-Virasoro classes and, at
    # this rational point, more).  Pinned as the observed record.
    SMOKE_PAIR_34 = {3: 1, 4: 1, 5: 0, 6: 0, 7: 1, 8: 2, 9: 1, 10: 0, 11: 2, 12: 4}
    badp = {w2: (dpair[w2], spec34[w2]) for w2 in w2s3 if w2 >= 3 and dpair[w2] < spec34[w2]}
    require(not badp, "S3: dual pair F+ kernel >= centralizer for w2 >= 3 "
            f"(FINDING if not: {badp})", findings)
    obs = {w2: dpair[w2] for w2 in SMOKE_PAIR_34 if w2 <= MAX}
    require(obs == {w2: v for w2, v in SMOKE_PAIR_34.items() if w2 <= MAX},
            "S3 (refutation record): the dual-pair kernel at lam=3/4 has dims "
            f"{SMOKE_PAIR_34} at w2=3..12 -- STRICTLY LARGER than the centralizer at sigma=3,5; "
            "the preregistered 'equal at odd-spin slots' clause is refuted", findings)
    # Exploratory (recorded, unpinned): dual-pair kernels at three rational
    # dual pairs vs the NS super-Virasoro vacuum class count p(h) - p(h-1)
    # (character prod_{n>=2} 1/(1-q^n) prod_{r>=3/2} (1+q^r); classes mod d).
    def ns_vacuum_dims(max_h2):
        p = {0: 1}
        for n in range(2, max_h2 // 2 + 1):
            for h2 in range(n * 2, max_h2 + 1):
                p[h2] = p.get(h2, 0) + p.get(h2 - 2 * n, 0)
        for r2 in range(3, max_h2 + 1, 2):
            for h2 in range(max_h2, r2 - 1, -1):
                p[h2] = p.get(h2, 0) + p.get(h2 - r2, 0)
        return {h2: p.get(h2, 0) - p.get(h2 - 2, 0) for h2 in range(3, max_h2 + 1)}
    nsv = ns_vacuum_dims(w2s3[-1])
    out["screening"]["ns_vacuum_classes"] = {str(k): v for k, v in nsv.items()}
    out["screening"]["dual_pairs"] = {}
    print("Dual-pair kernels vs NS vacuum classes (w2: NS / lam=3/4 / lam=4/3 / lam=5/12):", flush=True)
    pairs = {"3/4": (F(1, 2), F(-2)), "4/3": (F(1, 3), F(-3)), "5/12": (F(2, 3), F(-3, 2))}
    dp = {"3/4": dpair}
    for key in ("4/3", "5/12"):
        lam = F(key)
        a1, a2 = pairs[key]
        require(a1 * a1 + 2 * lam * a1 - 1 == 0 and a2 * a2 + 2 * lam * a2 - 1 == 0,
                f"S3c: {a1}, {a2} are the F+ roots a^2 + 2 lam a - 1 = 0 at lam={key}", findings)
        dp[key] = {w2: joint_screening_dim(w2, pairs[key]) for w2 in w2s3}
    for key in dp:
        out["screening"]["dual_pairs"][key] = {str(k): v for k, v in dp[key].items() if k >= 3}
    for w2 in w2s3:
        if w2 >= 3:
            print(f"  w2={w2:2d} sigma={sig(w2):>4}: {nsv[w2]} / {dp['3/4'][w2]} / {dp['4/3'][w2]} / {dp['5/12'][w2]}", flush=True)
    for key in ("3/4", "4/3", "5/12"):
        badn = {w2: (dp[key][w2], nsv[w2]) for w2 in w2s3 if w2 >= 3 and dp[key][w2] != nsv[w2]}
        require(not badn, f"S3c (prediction): dual-pair kernel at lam={key} equals the NS "
                f"super-Virasoro vacuum class count at every 3 <= w2 <= {w2s3[-1]} "
                f"(FINDING if not: {badn})", findings)
    seconds = round(time.time() - t0)
    out["findings"] = findings
    if CERTIFY:
        with open(OUT, "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
            fh.write("\n")
    print()
    if findings:
        print(f"SUPER PHASE 3: {len(findings)} PIN(S) NOT MET -- see FINDING lines")
        sys.exit(1)
    if not CERTIFY:
        print(f"SUPER PHASE 3 SMOKE RUN (w2 <= {MAX}): all pins met; NOT a certification, "
              f"no artifact written; {seconds} s")
        sys.exit(0)
    print(f"SUPER PHASE 3 CERTIFIED (S1 spectrum identical at the five sampled lambda "
          f"{out['lambdas']} through w2=24; S2 tower of G-singlets, no fermionic partners through "
          f"w2=24; S3 four-screening joint kernel = centralizer at lam=3/4 for 3 <= w2 <= 16, "
          f"dual pair = NS vacuum classes at three rational pairs; {seconds} s)")
