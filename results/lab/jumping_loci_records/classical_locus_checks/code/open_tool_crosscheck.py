"""Open-tool cross-check of the ledger.

Some ledger entries (LEDGER.json, written by run_a5.py) name as their record an earlier computation of the multiplicity
stratification at weights 6-14 made with Mathematica: results/cyl_discriminant_w1012.json, results/cyl_discriminant_w14.json
and wolfram/lab/cyl_strata_w*.json.  The same classification has since been computed with open-source tools only:
  - weights 10, 12, 14: the registered certificate cyl_stratification_w1014 (artifact results/cyl_stratification_w1014.json);
  - weights 6, 8: results/lab/jumping_loci_records/stratification_w6_w8/results/S23_sheet{1,2,3}_W{6,8}.json.
Both are complete classifications over Q-bar away from t = +-1: every point where the multiplicity jumps is listed, with its
multiplicity (nullity).  This script checks that every ledger entry resting on a Mathematica record states the outcome that the
open-source record gives: "jump" exactly when t is in the open record's jump set, with the same multiplicity where the ledger
gives one; "no jump" otherwise.  It fails (exit 1) on any disagreement and on any such entry it cannot place.

Usage: open_tool_crosscheck.py LEDGER.json
"""
import hashlib
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
W1014 = os.path.join(REPO, "results", "cyl_stratification_w1014.json")
W68 = os.path.join(REPO, "results", "lab", "jumping_loci_records", "stratification_w6_w8", "results")
MATHEMATICA_MARKS = ("results/cyl_discriminant_w1012.json", "results/cyl_discriminant_w14.json", "wolfram/lab/cyl_strata_w")


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def open_jumps():
    """{(W, sheet): {t: nullity}} from the two open-source records."""
    out = {}
    with open(W1014) as fh:
        a = json.load(fh)
    print(f"read {os.path.relpath(W1014, REPO)} sha256 {sha(W1014)}")
    for W in (10, 12, 14):
        for s in (1, 2, 3):
            out[(W, s)] = {Fraction(j["t"]): int(j["nullity"]) for j in a["weights"][str(W)][f"sheet{s}"]["jump_set"]}
    for W in (6, 8):
        for s in (1, 2, 3):
            p = os.path.join(W68, f"S23_sheet{s}_W{W}.json")
            with open(p) as fh:
                b = json.load(fh)
            assert b["W"] == W and b["sheet"] == s, p
            print(f"read {os.path.relpath(p, REPO)} sha256 {sha(p)}")
            out[(W, s)] = {Fraction(j["root"]): int(j["nullity"]) for j in b["jumps"]}
    return out


def entries(ledger):
    """Every (place, weight, sheet, t, outcome, dim, record) whose record is a Mathematica one."""
    for i, row in enumerate(ledger["rows"]):
        for k, rec in enumerate(row.get("all_records", [])):
            yield (f"rows[{i}].all_records[{k}]", row["weight"], row["sheet"], row["t"], rec["outcome"], rec.get("dim"),
                   rec["record"])
        yield (f"rows[{i}]", row["weight"], row["sheet"], row["t"], row["outcome"], row.get("dim"), row["record"])
    for key in ("deck_images", "t3"):
        for i, e in enumerate(ledger[key]):
            yield (f"{key}[{i}]", e["spin"] + 1, e["sheet"], e.get("t", "3"), e["outcome"], e.get("dim"), e["record"])
    for i, e in enumerate(ledger["nonclassical_kernels_w22plus"]):
        yield (f"nonclassical_kernels_w22plus[{i}]", e["weight"], e["sheet"], e["t"], "jump", e.get("dim"), e["record"])


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: open_tool_crosscheck.py LEDGER.json")
    with open(sys.argv[1]) as fh:
        ledger = json.load(fh)
    jumps = open_jumps()
    checked, bad, counts = 0, 0, {}
    for place, W, s, t, outcome, dim, record in entries(ledger):
        if not any(m in record for m in MATHEMATICA_MARKS):
            continue
        t = Fraction(str(t))
        key = (int(W), int(s))
        fam = next(m for m in MATHEMATICA_MARKS if m in record)
        if key not in jumps or t in (1, -1):
            print(f"CANNOT PLACE {place}: weight {W}, sheet {s}, t = {t}, record {record}")
            bad += 1
            continue
        js = jumps[key]
        expect = "jump" if t in js else "no jump"
        ok = (outcome == expect) and (expect == "no jump" or dim in (None, js[t]))
        checked += 1
        c = counts.setdefault((fam, int(W), expect), [0, 0])
        c[0] += 1
        if not ok:
            c[1] += 1
            bad += 1
            print(f"DISAGREES {place}: weight {W}, sheet {s}, t = {t}: ledger {outcome} (dim {dim}), "
                  f"open record {expect}" + (f" (nullity {js[t]})" if t in js else ""))
    for (fam, W, expect), (n, nbad) in sorted(counts.items()):
        print(f"{fam:<40} weight {W:>2}  {expect:<8} {n:>3} entries, {n - nbad:>3} agree")
    if bad or checked == 0:
        print(f"OPEN-TOOL CROSS-CHECK: FAIL ({bad} of {checked} entries disagree or cannot be placed)")
        sys.exit(1)
    print(f"OPEN-TOOL CROSS-CHECK: ALL {checked} ENTRIES AGREE")


if __name__ == "__main__":
    main()
