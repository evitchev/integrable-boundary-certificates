"""Check A5 (with the record parts of A1, A6, A9, A10): the lifting-rule ledger, odd spins 5-29.

Classical side (own code): every isolated non-decoupling classical jumping point per sheet and spin, from the closed forms
(jl_lib.sheet3_points / sheet2_points, checked in run_a4.py), each VERIFIED here by the exact nullity (= 2) of the conservation
matrix of the committed I_3's weight-4 symbol at that fiber (cyl_P4 at (N(t), s(t)), byte copy, Fractions), plus generic
controls (nullity 1).  Sheet 1 has none (run_a4.py, run_a2.py).
Quantum side: the records only, each parsed or transcribed with its file and line; nothing is recomputed.  Evidence classes
(the paper's labels; registry: code/certifications.json): certified, certified scan evidence, pinned lab record, lab result,
untested.  Output: LEDGER.json, LEDGER.md (written by this script) and the A1/A6/A9/A10 tables in the log."""
import os, sys, re, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "vendor"))
from fractions import Fraction as Fr
import sympy as sp
import jl_lib as L

REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))   # repository root (code/ -> classical_locus_checks/ -> jumping_loci_records/ -> lab/ -> results/ -> root)
if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]):
    sys.exit("usage: run_a5.py OUTPUT_DIR (an existing directory; LEDGER.json and LEDGER.md are written there)")
WS = os.path.abspath(sys.argv[1])   # output directory for LEDGER.json, LEDGER.md (run.sh passes its $1)
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)
def rp(path):
    return os.path.join(REPO, path)
for f in ("mixed_engine.py", "invariant_engine.py"):
    h1 = hashlib.sha256(open(os.path.join(HERE, "vendor", f), "rb").read()).hexdigest()
    h2 = hashlib.sha256(open(rp("code/" + f), "rb").read()).hexdigest()
    require(h1 == h2, f"vendor/{f} is a byte copy of the repository's code/{f}")
if fails:
    sys.exit(1)
import mixed_engine as me
KEYS = (me._M((1, 1, 1, 1)), me._M((1, 1), [(1, 1)]), me._M((), [(1, 1), (1, 1)]))

def rat(v):
    return sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else sp.Rational(v)
def symbol_at(sheet, tv):
    P = me.cyl_P4(sheet, L.N_of(tv), L.s_of(tv))
    return tuple(P.get(k, Fr(0)) for k in KEYS)
def nullity(sheet, tv, W):
    A, B, C = symbol_at(sheet, tv)
    return L.nullity_at(W, rat(A), rat(B), rat(C))
def F(*v):
    return Fr(*v)
def fs(v):
    return str(v)

CLASS_ORDER = ["certified", "certified scan evidence", "pinned lab record", "lab result", "untested"]

# ------------------------------------------------------------------ records: certified exact kernels (transcribed, cited)
CERT = []   # (sheet, t, W, dim, record)
def add(sheets, tv, Ws, dims, rec):
    for sh in sheets:
        for W, d in zip(Ws, dims):
            CERT.append((sh, F(tv), W, d, rec))
SP = "code/cyl_strata_points.py:63-71 (TABLE; certified)"
add([1], 3, (6, 8, 10), (1, 1, 1), SP); add([2, 3], 3, (6, 8, 10), (2, 2, 2), SP)
add([1], "5/3", (6, 8, 10), (1, 1, 1), SP); add([2, 3], "5/3", (6, 8, 10), (2, 1, 1), SP)
add([1, 2], "7/5", (6, 8, 10), (1, 1, 1), SP); add([3], "7/5", (6, 8, 10), (1, 1, 2), SP)
add([1, 2, 3], "1/3", (6, 8, 10), (1, 1, 1), SP)                       # (N, s) = (28, 9): t = 1/3
S11 = "code/cyl_strata_points_s11.py:50-68 (certified)"
add([2], "9/7", (6, 8, 10, 12), (1, 1, 1, 2), S11); add([1, 3], "9/7", (12,), (1,), S11)
add([3], "7/5", (12,), (1,), S11); add([1], 3, (12,), (1,), S11); add([2, 3], 3, (12,), (2,), S11)
JT = "code/cyl_jump_point_tower.py:127 (EXPECT; the merged I_3 of sheets 2 and 3 at t = 5/3; certified)"
add([2, 3], "5/3", (6, 8, 10, 12), (2, 1, 1, 1), JT)
S13 = "code/cyl_spin13_prediction.py:51-61 (certified; results/cyl_spin13_prediction_run.json pass)"
for sh in (1, 2, 3):
    add([sh], "11/9", (14,), (1,), S13); add([sh], 2, (14,), (1,), S13); add([sh], "5/4", (14,), (1,), S13)
    add([sh], "7/6", (14,), (1,), S13); add([sh], 3, (14,), (2 if sh > 1 else 1,), S13)
add([2], "9/7", (14,), (1,), S13); add([3], "7/5", (14,), (1,), S13); add([2, 3], "5/3", (14,), (1,), S13)
S13b = "code/cyl_spin13_sol3_sporadic.py:59-62 (certified; results/cyl_spin13_sol3_sporadic_run.json pass)"
add([1, 2], "9/7", (14,), (1,), S13b); add([3], "9/7", (14,), (2,), S13b)
S15 = "code/cyl_spin15_stratification.py:65-68 (certified; results/cyl_spin15_stratification_run.json pass)"
add([2, 3], "5/3", (16,), (2,), S15); add([2, 3], 3, (16,), (2,), S15); add([1], "5/3", (16,), (1,), S15); add([1], 3, (16,), (1,), S15)
S15p = "code/cyl_spin15_ladder_probe.py:57-61 (certified; results/cyl_spin15_ladder_probe_run.json pass)"
add([1, 2, 3], 2, (16,), (1,), S15p); add([2, 3], "9/7", (16,), (1,), S15p); add([2, 3], "11/9", (16,), (1,), S15p); add([3], "7/5", (16,), (1,), S15p)
S17 = "code/cyl_spin17_double_prediction.py:58-63 (certified; results/cyl_spin17_double_prediction_run.json pass)"
add([3], "11/9", (18,), (2,), S17); add([2], "13/11", (18,), (2,), S17); add([1, 2, 3], 2, (18,), (1,), S17)
add([2], "11/9", (18,), (1,), S17); add([3], "13/11", (18,), (1,), S17); add([2, 3], "9/7", (18,), (1,), S17)
add([1], "11/9", (18,), (1,), S17); add([1], "13/11", (18,), (1,), S17)
OQ = "code/cyl_oddq_controls.py:29-34 (certified; results/cyl_oddq_controls_run.json pass)"
for tv in (2, "3/2", "4/3"):
    add([2, 3], tv, (10, 12, 14), (1, 1, 1), OQ)
# the run records named above must say pass
for f in ("cyl_spin13_prediction", "cyl_spin13_sol3_sporadic", "cyl_spin15_stratification", "cyl_spin15_ladder_probe",
          "cyl_spin17_double_prediction", "cyl_oddq_controls", "cyl_discriminant_w1012", "cyl_discriminant_w14"):
    rj = json.load(open(rp(f"results/{f}_run.json")))
    require(all(x_["status"] == "pass" and x_["exit_code"] == 0 for x_ in rj["results"]), f"results/{f}_run.json: status pass, exit 0")
# consistency of the transcription: no (sheet, t, W) with two different dimensions
seen = {}
for sh, tv, W, d, rec in CERT:
    seen.setdefault((sh, tv, W), set()).add(d)
require(all(len(v) == 1 for v in seen.values()), f"the {len(CERT)} transcribed certified pins are mutually consistent ({len(seen)} cells)")

# ------------------------------------------------------------------ records: certified Q-bar classification at W = 10, 12, 14
DISC = {}
for fn, Ws in (("results/cyl_discriminant_w1012.json", (10, 12)), ("results/cyl_discriminant_w14.json", (14,))):
    dj = json.load(open(rp(fn)))
    for W in Ws:
        for sh in (1, 2, 3):
            e = dj[f"w{W}_sol{sh}"]
            DISC[(sh, W)] = ({F(j["t0"]) for j in e["jumps"]}, {F(n_["t0"]) for n_ in e["nonjump_factors"] if n_.get("degree") == 1}, fn)
# ------------------------------------------------------------------ records: lab strata (Mathematica) at W = 6..14
STRATA = {}
for W in (6, 8, 10, 12, 14):
    sj = json.load(open(rp(f"wolfram/lab/cyl_strata_w{W}.json")))
    for sh in (1, 2, 3):
        e = sj[f"sol{sh}"]
        jumps = {F(j["t0"]): j["nullity"] for j in e["jumps"]}
        rej = {F(c_["t0"]): c_["nullity"] for c_ in e["rejected_candidates"] if c_.get("degree") == 1}
        STRATA[(sh, W)] = (jumps, rej, f"wolfram/lab/cyl_strata_w{W}.json")
# ------------------------------------------------------------------ records: residue scans at W = 16, 18 (pinned) and 20 (certified scan evidence)
SCANS = {}   # (sheet, W, p) -> (evaluated, {residue: lift}, record, class)
for W in (16, 18):
    for p in (10007, 10009):
        fn = f"results/lab/modp_scan/w{W}_p{p}_summary.json"
        sj = json.load(open(rp(fn)))
        for sh in (1, 2, 3):
            e = sj["sheets"][str(sh)]
            SCANS[(sh, W, p)] = (e["residues_evaluated"], {d_["residue"]: d_["lift"] for d_ in e["rank_drops"]}, fn, "pinned lab record")
sj = json.load(open(rp("results/lab/spin19/summary.json")))
for e in sj["scans"]:
    SCANS[(e["sheet"], 20, e["prime"])] = (e["evaluated"], {d_["residue"]: d_["lift"] for d_ in e["drop_details"]},
                                           "results/lab/spin19/summary.json", "certified scan evidence")
require(all((sh, W, p) in SCANS for sh in (1, 2, 3) for W in (16, 18, 20) for p in (10007, 10009)),
        "scan summaries present for every sheet, weight 16/18/20 and both primes")
require(all(SCANS[k][0] == k[2] - 2 for k in SCANS), "every scan evaluated p - 2 residues (all but the poles 1 and p - 1)")
# ------------------------------------------------------------------ records: exact kernels at W = 22..28 (fast_kernel3)
FK3 = {}
for fn, cls in (("w22_t1311_sheet3.json", "pinned lab record"), ("w24_t1715_sheet2.json", "pinned lab record"),
                ("w22_sol1_t117.json", "lab result"), ("w26_t1513_sheet3.json", "lab result"), ("w26_t53_sheet2.json", "lab result"),
                ("w26_t53_sheet3.json", "lab result"), ("w28_t75_sheet3.json", "lab result")):
    kj = json.load(open(rp("results/lab/fast_kernel3/" + fn)))
    Nv, sv = F(kj["N"]), F(kj["s"])
    cands = [tv for tv in (F(5, 3), F(7, 5), F(13, 11), F(15, 13), F(17, 15), F(11, 7)) if (L.N_of(tv), L.s_of(tv)) == (Nv, sv)]
    require(len(cands) == 1, f"fast_kernel3/{fn}: (N, s) = ({Nv}, {sv}) is the fiber t = {cands}")
    FK3[(kj["solution"], cands[0], kj["weight"])] = (kj["genuine_dim"], "results/lab/fast_kernel3/" + fn, cls)
reg = open(rp("code/certifications.json")).read()
for (sh, tv, W), (d, fn, cls) in FK3.items():
    pinned = fn in reg
    require(pinned == (cls == "pinned lab record"), f"{fn}: {'listed' if pinned else 'not listed'} in code/certifications.json -> {cls}")
# ------------------------------------------------------------------ records: lifting_out.txt (pinned lab record)
LIFT = {}
for line in open(rp("results/lab/classical/lifting_out.txt")):
    m = re.match(r"Solution (\d), t = (\S+), weight (\d+) \[[^\]]*\]: classical conserved dim (\d+), quantum kernel dim (\d+)", line)
    if m:
        LIFT[(int(m.group(1)), F(m.group(2)), int(m.group(3)))] = (int(m.group(4)), int(m.group(5)))
require("results/lab/classical/lifting_out.txt" in reg and len(LIFT) == 6, f"lifting_out.txt: {len(LIFT)} entries, listed in the registry")

def residue(tv, p):
    if tv.denominator % p == 0:
        return None
    return tv.numerator * pow(tv.denominator, -1, p) % p

def evidence(sh, tv, W):
    """all record facts at (sheet, t, W): list of (outcome, dim or None, class, record)"""
    ev = []
    for csh, ctv, cW, d, rec in CERT:
        if (csh, ctv, cW) == (sh, tv, W):
            ev.append(("jump" if d >= 2 else "no jump", d, "certified", rec))
    if (sh, W) in DISC:
        jumps, nonj, fn = DISC[(sh, W)]
        ev.append(("jump" if tv in jumps else "no jump", None, "certified",
                   f"{fn} (Q-bar classification minus t = +-1; t {'a jump' if tv in jumps else ('a non-jump factor' if tv in nonj else 'not a root of the minor')})"))
    if (sh, W) in STRATA:
        jumps, rej, fn = STRATA[(sh, W)]
        if tv in jumps:
            ev.append(("jump", jumps[tv], "lab result", fn + " (jump root)"))
        else:
            ev.append(("no jump", rej.get(tv, 1), "lab result", fn + (" (rejected root, exact rank)" if tv in rej else " (not a root of the maximal minor)")))
    if (sh, tv, W) in LIFT:
        cd, qd = LIFT[(sh, tv, W)]
        ev.append(("jump" if qd >= 2 else "no jump", qd, "pinned lab record", "results/lab/classical/lifting_out.txt"))
    if W in (16, 18, 20):
        outs = []
        for p in (10007, 10009):
            evn, drops, fn, cls = SCANS[(sh, W, p)]
            rr = residue(tv, p)
            if rr is None or rr in (1, p - 1):
                outs.append((p, "not covered", fn, cls))
            else:
                outs.append((p, "drop" if rr in drops else "no drop", fn, cls))
        cls = outs[0][3]; fn = outs[0][2]
        if all(o[1] == "no drop" for o in outs):
            ev.append(("no jump", None, cls, f"{fn} (no drop at either prime)"))
        elif all(o[1] == "drop" for o in outs):
            ev.append(("jump", None, cls, f"{fn} (drop at both primes)"))
        else:
            ev.append(("mixed", None, cls, f"{fn} ({[(o[0], o[1]) for o in outs]})"))
    if (sh, tv, W) in FK3:
        d, fn, cls = FK3[(sh, tv, W)]
        ev.append(("jump" if d >= 2 else "no jump", d, cls, fn))
    return ev

def primary(ev):
    if not ev:
        return ("untested", None, "untested", "")
    return min(ev, key=lambda e: CLASS_ORDER.index(e[2]))

# ------------------------------------------------------------------ the ledger
rows = []
ctrl_bad = []
for s in range(5, 31, 2):
    W = s + 1
    pts = [(3, (j, j), tv, "affine") for j, tv in L.sheet3_points(s).items()]
    pts += [(2, jj, tv, "infinity" if jj[0] == 3 * jj[1] else "affine") for jj, tv in L.sheet2_points(s).items()]
    for sh, (j, jp), tv, chart in pts:
        xv = L.x_of(tv)
        nl = nullity(sh, tv, W)
        q = xv.denominator if xv.numerator == 1 else None
        rule_pos = q is not None and q % 2 == 0
        kind = "rule null"
        if rule_pos:
            k = (q - 2) // 2
            low = 4 * k + 1 if sh == 3 else 3 * k + 2
            kind = "ladder" if s == low else "recurrence"
        ev = evidence(sh, tv, W)
        pr = primary(ev)
        outcome = pr[0]
        consistent = None if outcome == "untested" else ((outcome == "jump") == rule_pos)
        conflict = len({e[0] for e in ev}) > 1
        rows.append(dict(spin=s, weight=W, sheet=sh, t=fs(tv), beta2=fs(xv), N=fs(L.N_of(tv)), s_curve=fs(L.s_of(tv)), j=j, jp=jp,
                         chart=chart, classical_nullity=nl, q=q, odd_q=(q is not None and q % 2 == 1), kind=kind,
                         rule_prediction="jump" if rule_pos else "no jump", outcome=outcome, evidence_class=pr[2],
                         dim=pr[1], record=pr[3], all_records=[dict(outcome=e[0], dim=e[1], cls=e[2], record=e[3]) for e in ev],
                         records_agree=not conflict, consistent_with_rule=consistent))
    # controls: two generic fibers per sheet and spin must have classical nullity 1
    for sh in (1, 2, 3):
        for tv in (F(37, 11), F(-29, 13)):
            if nullity(sh, tv, W) != 1:
                ctrl_bad.append((sh, tv, W))
require(all(r_["classical_nullity"] == 2 for r_ in rows), f"all {len(rows)} ledger points have classical nullity exactly 2 with the committed I_3's symbol at the fiber")
require(not ctrl_bad, f"generic controls t = 37/11, -29/13 on the three sheets at every W = 6..30 have classical nullity 1 ({ctrl_bad})")
require(all(r_["records_agree"] for r_ in rows), "where several records cover a point they agree")
incons = [r_ for r_ in rows if r_["consistent_with_rule"] is False]
require(not incons, f"no recorded outcome contradicts the lifting rule ({len(incons)}: {[(r_['spin'], r_['sheet'], r_['t']) for r_ in incons]})")

print("== the ledger (spin, sheet, t, beta^2, (j, j'), chart, kind, rule, outcome, class, record)")
for r_ in rows:
    print(f"  s={r_['spin']:2d} sheet {r_['sheet']} t={r_['t']:>6} x={r_['beta2']:>6} (j,j')=({r_['j']},{r_['jp']}) {r_['chart']:8} "
          f"{r_['kind']:10} rule:{r_['rule_prediction']:7} -> {r_['outcome']:8} [{r_['evidence_class']}] {r_['record'][:110]}", flush=True)
cnt = {}
for r_ in rows:
    cnt.setdefault((r_["kind"], r_["evidence_class"], r_["outcome"]), 0)
    cnt[(r_["kind"], r_["evidence_class"], r_["outcome"])] += 1
print("== counts (kind, class, outcome): " + "; ".join(f"{k}: {v}" for k, v in sorted(cnt.items())), flush=True)

print("== R2's three statements")
n21 = sorted(Fr(r_["beta2"]) for r_ in rows if r_["spin"] == 21 and r_["sheet"] == 3 and r_["kind"] == "rule null")
require(n21 == sorted([Fr(3, 8), Fr(5, 18), Fr(13, 30), Fr(17, 36)]), f"(i) the sheet-3 non-lifting classical points at spin 21: beta^2 = {[str(v) for v in n21]}")
tested_null_hi = [r_ for r_ in rows if r_["spin"] >= 21 and r_["kind"] == "rule null" and r_["outcome"] != "untested"]
require(not tested_null_hi, f"(ii) no rule-null classical point at spin >= 21 has a record ({len(tested_null_hi)})")
nonclass_hi = [(sh, str(tv), W, d, fn) for (sh, tv, W), (d, fn, cls) in FK3.items()
               if not any(r_["sheet"] == sh and Fr(r_["t"]) == tv and r_["weight"] == W for r_ in rows)]
print(f"      exact kernels at weight >= 22 at points that are not classical jumping points: {nonclass_hi}", flush=True)
pos_rule = {(r_["spin"], r_["sheet"], r_["t"]) for r_ in rows if r_["spin"] >= 21 and r_["rule_prediction"] == "jump"}
lad_rec = set()
for k in range(1, 40):
    tk = L.ladder(k)
    for i in range(0, 20):
        s3 = (4 * k + 1) * (2 * i + 1)
        if 21 <= s3 <= 29: lad_rec.add((s3, 3, fs(tk)))
        if k % 2 == 1:
            s2 = (3 * k + 2) * (2 * i + 1)
            if 21 <= s2 <= 29: lad_rec.add((s2, 2, fs(tk)))
require(pos_rule == lad_rec, f"(iii) the rule's positive cells at spins 21-29 = the two ladders' slots plus the classical recurrence spins "
        f"(4k+1)(2i+1) on sheet 3, (3k+2)(2i+1), k odd, on sheet 2: {sorted(pos_rule)}")
pos_out = {(r_["spin"], r_["sheet"], r_["t"]): (r_["outcome"], r_["evidence_class"]) for r_ in rows if r_["spin"] >= 21 and r_["rule_prediction"] == "jump"}
print(f"      their recorded outcomes: {sorted(pos_out.items())}", flush=True)
# the same equivalence for all s <= 61 (the rule's positives are exactly the t_k at their classical spins)
eq_all = True
for s in range(5, 63, 2):
    for sh, pts in ((3, L.sheet3_points(s).values()), (2, L.sheet2_points(s).values())):
        for tv in pts:
            xv = L.x_of(tv); rule_pos = xv.numerator == 1 and xv.denominator % 2 == 0
            is_tk = tv in {L.ladder(k) for k in range(1, 40)}
            eq_all &= rule_pos == is_tk
require(eq_all, "for every classical point at odd s <= 61: beta^2 = 1/(2k) iff t is a ladder point t_k (the rule's positive content = ladders + recurrences)")

print("== A1: the deck images -t0 of the jump points")
deck = []
for tv in [F(3)] + [L.ladder(k) for k in range(1, 8)]:
    for W in range(6, 30, 2):
        for sh in (1, 2, 3):
            ev = evidence(sh, -tv, W)
            pr = primary(ev)
            deck.append((str(-tv), W - 1, sh, pr[0], pr[2], pr[3]))
deck_jumps = [d_ for d_ in deck if d_[3] == "jump"]
require(not deck_jumps, f"no record shows a jump at any deck image -t0 ({len(deck_jumps)})")
by = {}
for d_ in deck:
    by.setdefault((d_[0], d_[4]), []).append(d_[1])
for (tt, cls), spins in sorted(by.items(), key=lambda z: (Fr(z[0][0]), CLASS_ORDER.index(z[0][1]))):
    print(f"      -t0 = {tt:>7}: {cls:24} at spins {sorted(set(spins))} (all three sheets)" if len(spins) == 3 * len(set(spins))
          else f"      -t0 = {tt:>7}: {cls:24} at (spin) {sorted(spins)}", flush=True)

print("== A6: t = 3, per sheet and spin")
t3 = {}
for W in range(6, 30, 2):
    for sh in (1, 2, 3):
        ev = [e for e in evidence(sh, F(3), W)]
        pr = primary(ev)
        t3[(sh, W - 1)] = (pr[0], pr[1], pr[2], pr[3], len(ev))
for sh in (1, 2, 3):
    print(f"      sheet {sh}: " + "; ".join(f"s={s}: {t3[(sh, s)][0]}{'' if t3[(sh, s)][1] is None else ' (dim ' + str(t3[(sh, s)][1]) + ')'} [{t3[(sh, s)][2]}]"
                                         for s in range(5, 29, 2)), flush=True)
require(all(t3[(1, s)][0] in ("no jump", "untested") for s in range(5, 29, 2)) and all(t3[(1, s)][0] == "no jump" for s in range(5, 21, 2)),
        "sheet 1 at t = 3: no jump at every recorded spin (5-19); nothing recorded at 21-27")
require(all(t3[(sh, s)][0] == "jump" for sh in (2, 3) for s in range(5, 21, 2)), "sheets 2 and 3 at t = 3: a jump at every spin 5-19 in the records")
require(all(t3[(sh, s)][0] == "untested" for sh in (1, 2, 3) for s in range(21, 29, 2)), "t = 3 at spins 21-27: no record on any sheet")
cl3 = {sh: [nullity(sh, F(3), W) for W in range(4, 32, 2)] for sh in (1, 2, 3)}
require(cl3[2] == [2] * 14 and cl3[3] == [2] * 14 and cl3[1] == [1] * 14,
        f"classically at t = 3 (W = 4..30): sheets 2, 3 nullity 2 (sheet 3 on b = 0, sheet 2 at r^4), sheet 1 nullity 1 ((b, c) = {L.MAPS[1][0].subs(L.t, 3)}, {L.MAPS[1][1].subs(L.t, 3)})")
sp7 = [e for e in evidence(2, F(3), 8) + evidence(3, F(3), 8) if e[2] == "certified"]
require(any(e[1] == 2 for e in sp7), f"Table 1's spin-7 row '1 everywhere': the certified record has dimension 2 at t = 3 on sheets 2, 3 at spin 7 ({[e[3][:40] for e in sp7]})")
sp7_cls = [(r_["sheet"], r_["t"], r_["outcome"], r_["evidence_class"], r_["record"]) for r_ in rows if r_["spin"] == 7]
print(f"      the spin-7 classical points and their records: {sp7_cls}", flush=True)

print("== A9: sheet 2 at t = 7/5 (q = 6)")
ev = evidence(2, F(7, 5), 8)
print(f"      records at spin 7: {[(e[0], e[1], e[2], e[3][:60]) for e in ev]}", flush=True)
jumps8, rej8, _ = STRATA[(2, 8)]
w8 = json.load(open(rp("wolfram/lab/cyl_strata_w8.json")))["sol2"]
facs = [c_.get("factor") for c_ in w8["rejected_candidates"]] + [j_["factor"] for j_ in w8["jumps"]]
require(not any(f_ and f_.replace(" ", "") in ("-7+5*t", "5*t-7") for f_ in facs) and all(e[0] == "no jump" for e in ev),
        f"sheet 2, spin 7, t = 7/5: no jump (records above); the weight-8 maximal minor of sheet 2 has no factor 5t - 7 (factors {facs})")
s2_q6 = [s for s in range(3, 200, 2) if any(tv == F(7, 5) for tv in L.sheet2_points(s).values())]
require(s2_q6 == [] and (3 * 6 - 2) % 2 == 0 and (3 * 6 - 2) // 2 == 8,
        "q = 6: the sheet-2 slot (3q-2)/2 = 8 is even, and t = 7/5 is a sheet-2 classical point at no odd spin < 200 (x = 1/6 needs s = 8j)")

print("== A10: t_2 = 7/5 on sheet 3")
sp_t2 = [s for s in range(3, 62, 2) if F(7, 5) in L.sheet3_points(s).values()]
require(sp_t2 == [9, 27, 45], f"sheet-3 classical spins of t_2 (s <= 61): {sp_t2} = 9(2i+1)")
d28 = FK3.get((3, F(7, 5), 28))
require(d28 is not None and d28[0] == 2, f"the weight-28 record at t_2 on sheet 3: genuine dimension {d28[0] if d28 else None} [{d28[2] if d28 else ''}] {d28[1] if d28 else ''}")

# ------------------------------------------------------------------ outputs
out = dict(generated_by="code/run_a5.py",
           classes=CLASS_ORDER, rows=rows,
           deck_images=[dict(t=d_[0], spin=d_[1], sheet=d_[2], outcome=d_[3], evidence_class=d_[4], record=d_[5]) for d_ in deck],
           t3=[dict(sheet=sh, spin=s, outcome=v[0], dim=v[1], evidence_class=v[2], record=v[3]) for (sh, s), v in sorted(t3.items())],
           nonclassical_kernels_w22plus=[dict(sheet=a, t=b_, weight=c_, dim=d_, record=e_) for a, b_, c_, d_, e_ in nonclass_hi])
with open(os.path.join(WS, "LEDGER.json"), "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
md = ["# LEDGER: the classical jumping points of the cylindrical sheets at odd spins 5-29 and their quantum outcomes", "",
      "Generated by code/run_a5.py (log logs/run_a5.log); machine-readable copy LEDGER.json.", "",
      "Points: every isolated non-decoupling classical jumping point (eq. (3) with 3j != s, and the point at infinity of sheet 2 at",
      "t = 11/5), each verified by the exact nullity 2 of the conservation matrix of the committed I_3's weight-4 symbol at the fiber.",
      "Sheet 1 has none at any odd spin. The decoupling point t = 3 (sheets 2 and 3) is not isolated and is listed in REPORT.md (A6).",
      "kind: ladder = t_k at its lowest classical spin on the sheet (4k+1 on sheet 3, 3k+2 for odd k on sheet 2); recurrence = t_k at a",
      "higher classical spin; rule null = beta^2 not of the form 1/(2k) (the lifting rule predicts no quantum jump). odd q: beta^2 = 1/q, q odd.",
      "Outcome and class: the strongest record (certified > certified scan evidence > pinned lab record > lab result > untested); every",
      "record covering the point is in LEDGER.json (all_records). Scan rows: no drop at either prime (10007, 10009) at the point's residue.", ""]
for s in range(5, 31, 2):
    rs = [r_ for r_ in rows if r_["spin"] == s]
    md.append(f"## spin {s} (weight {s + 1})")
    md.append("")
    md.append("| sheet | t | beta^2 | (N, s) | (j, j') | kind | rule | outcome | class | record |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    for r_ in sorted(rs, key=lambda z: (-z["sheet"], Fr(z["t"]))):
        kind = r_["kind"] + (" (odd q = %d)" % r_["q"] if r_["odd_q"] else "") + (" (at infinity)" if r_["chart"] == "infinity" else "")
        dim = "" if r_["dim"] is None else f", dim {r_['dim']}"
        md.append(f"| {r_['sheet']} | {r_['t']} | {r_['beta2']} | ({r_['N']}, {r_['s_curve']}) | ({r_['j']}, {r_['jp']}) | {kind} | "
                  f"{r_['rule_prediction']} | {r_['outcome']}{dim} | {r_['evidence_class']} | {r_['record']} |")
    md.append("")
md.append("## Summary counts (kind, class, outcome)")
md.append("")
for k, v in sorted(cnt.items()):
    md.append(f"- {k[0]}, {k[1]}, {k[2]}: {v}")
md.append("")
with open(os.path.join(WS, "LEDGER.md"), "w") as fh:
    fh.write("\n".join(md) + "\n")
print(f"wrote LEDGER.json ({len(rows)} rows) and LEDGER.md", flush=True)
print()
if fails:
    print(f"A5: {len(fails)} FAILURE(S)"); sys.exit(1)
print("A5 DONE: all checks pass")
