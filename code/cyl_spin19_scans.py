"""Spin-19 (weight-20) residue scans: machine-readable summary and consistency
certificate (Codex round 19, 2026-09-03).

results/lab/spin19/summary.json records, per (sheet, prime): the prime, the
16 residue chunks with their [T0, T1) ranges, evaluated counts and drops,
the drop details, the archived log path, and the scanner / matrix sidecar
hashes.  Provenance caveat (Codex round 20): the archived chunk logs carry
NO execution provenance (no scanner hash, matrix hash, invocation or launch
state), so the scanner and matrix hashes in the summary are an
operator-asserted, retrospective association with the files now in the
repository -- pinning prevents later alteration but does not establish
historically which scanner and matrix produced the logs.  Future chunks run
through lab/modp_scan_attested.sh, which prints exactly that provenance
before computing.

What this certificate ENFORCES against the archived logs: every log header
names the expected sheet and prime; the parsed chunk table (T0, T1,
evaluated, drops) equals the summary's chunk table entry by entry; ranges
are contiguous and cover [0, P) exactly; evaluated residues total P - 2
(poles t = 1, P-1 skipped); the generic rank is 3927 of 3928 at t = 2, 5, 7
in every chunk (and matches the summary's generic_rank field); every drop
line's modulus equals P; the drop list and drop details equal the recorded
ones; the recorded check points hold in every chunk; the summary's log path
is the file checked; the scanner and sidecar hashes equal the repository's.

What this certifies: finite-field EVIDENCE -- a residue scan detecting no
rank drop at a prime is not a characteristic-zero proof of the absence of
a jump; two primes strengthen the evidence.  Preregistration e6dc226.
"""
import sys, json, re, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)

def parse_log(path):
    """returns (headers, chunks): headers = list of (sheet, prime, T0) from '== scan_w20_...' lines"""
    chunks, headers = [], []
    cur = None
    for line in path.read_text().splitlines():
        m = re.match(r"== scan_w20_sol(\d)_p(\d+)_chunk_(\d+)\.out", line)
        if m:
            headers.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
            cur = {"T0": int(m.group(3)), "generic": None, "range": None, "evaluated": None, "drops": [], "checks": {}, "moduli": []}
            chunks.append(cur); continue
        if cur is None: continue
        m = re.match(r"generic ranks at t=2,5,7: \((\d+), (\d+), (\d+)\)", line)
        if m: cur["generic"] = tuple(int(x) for x in m.groups())
        m = re.match(r"scanned residues \[(\d+),(\d+)\): (\d+) evaluated", line)
        if m: cur["range"] = (int(m.group(1)), int(m.group(2))); cur["evaluated"] = int(m.group(3))
        m = re.match(r"  t = (\d+) \(mod (\d+)\)  rank (\d+) \(nullity (\d+)\)  lifts to t = (\S+)", line)
        if m:
            cur["drops"].append({"residue": int(m.group(1)), "rank": int(m.group(3)), "nullity": int(m.group(4)), "lift": m.group(5)})
            cur["moduli"].append(int(m.group(2)))
        m = re.match(r"  check t=(\S+): residue (\d+), rank (\d+)", line)
        if m: cur["checks"][m.group(1)] = (int(m.group(2)), int(m.group(3)))
        m = re.match(r"ATTEST (\S+)=(.*)", line)
        if m: cur.setdefault("attest", []).append(line)
    return headers, chunks

def residue_of(tstr, P):
    """the residue of the rational t = a/b modulo P"""
    a, b = (tstr.split("/") + ["1"])[:2]
    return (int(a) * pow(int(b), -1, P)) % P

def _public_edition():
    """A genuine public edition, by the runner's own criteria (run_certifications._public_edition, Codex rounds
    16 and 25): EXPORT_RECORD.json at the root parsing as JSON with kind == "public_export_record" and a 40-hex
    private_freeze_commit, the registry stamped with a _public_edition record, and the record bound to the
    edition -- tracked by git, or listed with a matching sha256 in RELEASE_MANIFEST.json of a git-free archive.
    A stray or empty EXPORT_RECORD.json in the development repository does not qualify."""
    import json as _json, subprocess as _sp
    rec = ROOT / "EXPORT_RECORD.json"
    if not rec.is_file():
        return False
    try:
        d = _json.loads(rec.read_text(encoding="utf-8"))
        reg = _json.loads((ROOT / "code" / "certifications.json").read_text(encoding="utf-8"))
    except Exception:
        return False
    fc = str(d.get("private_freeze_commit", ""))
    if (not isinstance(d, dict) or d.get("kind") != "public_export_record" or len(fc) != 40
            or any(ch not in "0123456789abcdef" for ch in fc) or not isinstance(reg.get("_public_edition"), dict)):
        return False
    try:
        tracked = _sp.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", "EXPORT_RECORD.json"], capture_output=True)
        if tracked.returncode == 0:
            return True
    except OSError:
        pass
    man = ROOT / "RELEASE_MANIFEST.json"
    if not man.is_file():
        return False
    try:
        m = _json.loads(man.read_text(encoding="utf-8"))
        files = m.get("files", m) if isinstance(m, dict) else {}
        pin = files.get("EXPORT_RECORD.json")
        want = pin.get("sha256") if isinstance(pin, dict) else pin
        have = "sha256:" + hashlib.sha256(rec.read_bytes()).hexdigest()
        return isinstance(want, str) and (want == have or want == have.split(":", 1)[1])
    except Exception:
        return False


def check_attest(lines, P, S, scanner_sha, side_sha, T0, T1, W=20):
    """exactly one attestation block (6 lines) whose fields are bound to THIS chunk: scanner hash, matrix
    filename and hash == sidecar hash, structural invocation 'modp_scan.py W S P T0 T1' matching the
    chunk's header and range, a UTC timestamp, python/flint versions, a 40-hex commit with a clean tree;
    when the commit is resolvable in this history, the scanner blob at that commit must have the attested hash."""
    if len(lines) != 6:
        return False, f"{len(lines)} ATTEST lines (expected 6)"
    kv = {}
    for l in lines:
        body = l[len("ATTEST "):]
        m = re.match(r"invocation=modp_scan\.py (\d+) (\d+) (\d+) (\d+) (\d+)\s+started=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)$", body)
        if m:
            kv["invocation"] = tuple(int(x) for x in m.groups()[:5]); kv["started"] = m.group(6); continue
        for tok in body.split():
            if "=" in tok:
                k, val = tok.split("=", 1); kv.setdefault(k, []).append(val)
    if kv.get("invocation") != (W, S, P, T0, T1):
        return False, f"invocation {kv.get('invocation')} != (20, {S}, {P}, {T0}, {T1})"
    if kv.get("scanner_sha256") != [scanner_sha[7:]]:
        return False, "scanner hash"
    if kv.get("matrix") != [f"cyl_strata_w{W}_sol{S}.m"] or kv.get("matrix_sha256") != [side_sha] or kv.get("sidecar_matrix_sha256") != [side_sha]:
        return False, "matrix filename / hash"
    if not kv.get("host") or not kv.get("python") or not kv.get("flint"):
        return False, "host/python/flint fields"
    commit = (kv.get("git") or [""])[0]
    if re.fullmatch(r"[0-9a-f]{40}", commit) is None or kv.get("dirty") != ["0"]:
        return False, "commit / dirty"
    # commit resolution (Codex round 23): in the DEVELOPMENT repository the attested commit must resolve
    # and the scanner blob at it must carry the attested hash; only a public edition (EXPORT_RECORD.json at
    # the root: exported clone or archive) may treat a private commit as an opaque archival reference.
    import subprocess
    public_edition = _public_edition()
    try:
        blob = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:lab/modp_scan.py"], capture_output=True)
    except OSError:
        return (public_edition, "git unavailable: commit not resolved" + ("" if public_edition else " -- fatal outside a public edition"))
    if blob.returncode == 0:
        if "sha256:" + hashlib.sha256(blob.stdout).hexdigest() != scanner_sha:
            return False, "scanner blob at the attested commit differs from the attested hash"
        return True, "commit resolved, scanner blob verified"
    if public_edition:
        return True, "commit not in this history (opaque archival reference, public edition)"
    return False, f"attested commit {commit[:12]} does not resolve in the development repository"

SCHEMA_SCAN = {"sheet", "prime", "matrix_sha256", "chunks", "evaluated", "drops", "drop_details", "checks", "log", "attested"}
INVENTORY = {(1, 10007), (1, 10009), (2, 10007), (2, 10009), (3, 10007), (3, 10009)}   # the exact set of (sheet, prime) scans this certificate covers; extend when new scans land
SCHEMA_CHUNK = {"T0", "T1", "evaluated", "drops"}

if __name__ == "__main__":
    summary = json.loads((ROOT / "results/lab/spin19/summary.json").read_text())
    require(set(summary) >= {"weight", "spin", "preregistration", "scanner", "scanner_sha256", "generic_rank", "evidence_class", "provenance", "scans"}
            and summary["weight"] == 20 and summary["spin"] == 19, "summary schema: top-level fields present, weight 20 / spin 19")
    scanner_sha = "sha256:" + hashlib.sha256((ROOT / summary["scanner"]).read_bytes()).hexdigest()
    require(summary["scanner_sha256"] == scanner_sha, f"scanner {summary['scanner']} hash in the repository equals the summary's (retrospective association, see docstring)")
    gr = summary["generic_rank"]
    require(gr == {"columns": 3928, "rank": 3927, "test_points": [2, 5, 7]}, "summary generic_rank field: 3927 of 3928 at t = 2, 5, 7")
    require({(sc.get("sheet"), sc.get("prime")) for sc in summary["scans"]} == INVENTORY and len(summary["scans"]) == len(INVENTORY),
            f"summary inventory is exactly {sorted(INVENTORY)}")
    for scan in summary["scans"]:
        require(set(scan) >= SCHEMA_SCAN, f"scan entry schema complete (sheet {scan.get('sheet')}, prime {scan.get('prime')})")
        S, P = scan["sheet"], scan["prime"]
        require(scan["log"] == f"results/lab/spin19/scan_w20_sol{S}_p{P}.txt", f"sheet {S} p={P}: summary log path names the archived log for this sheet and prime")
        side = json.loads((ROOT / f"wolfram/lab/cyl_strata_w20_sol{S}.json").read_text())
        require(side["matrix_sha256"] == scan["matrix_sha256"] and side["rows"] == 11618 and side["cols"] == 3928,
                f"sheet {S}: sidecar matrix hash {scan['matrix_sha256'][:12]} and shape 11618 x 3928 match the summary (retrospective association)")
        headers, chunks = parse_log(ROOT / scan["log"])
        require(len(headers) == 16 and all(h[0] == S and h[1] == P for h in headers),
                f"sheet {S} p={P}: all 16 log headers name sheet {S} and prime {P}")
        chunks = sorted(chunks, key=lambda c: c["T0"])
        require(len(chunks) == 16 and all(c["range"] for c in chunks), f"sheet {S} p={P}: 16 finished chunks in the archived log")
        require(all(set(c) >= SCHEMA_CHUNK for c in scan["chunks"]) and len(scan["chunks"]) == 16, f"sheet {S} p={P}: summary chunk table has 16 complete entries")
        table = sorted(scan["chunks"], key=lambda c: c["T0"])
        same = all(c["range"] == (e["T0"], e["T1"]) and c["evaluated"] == e["evaluated"] and c["drops"] == e["drops"]
                   for c, e in zip(chunks, table))
        require(same, f"sheet {S} p={P}: parsed chunk table (T0, T1, evaluated, drops) equals the summary's entry by entry")
        ranges = [c["range"] for c in chunks]
        contiguous = ranges[0][0] == 0 and ranges[-1][1] == P and all(ranges[i][1] == ranges[i+1][0] for i in range(15))
        require(contiguous, f"sheet {S} p={P}: chunk ranges contiguous and covering [0, {P}) exactly")
        ev = sum(c["evaluated"] for c in chunks)
        require(ev == P - 2 and ev == scan["evaluated"], f"sheet {S} p={P}: {ev} residues evaluated = P - 2 (poles t = 1, P-1 skipped), as recorded")
        require(all(c["generic"] == (gr["rank"],) * 3 for c in chunks), f"sheet {S} p={P}: generic rank {gr['rank']} of {gr['columns']} at t = 2, 5, 7 in every chunk")
        require(all(mod == P for c in chunks for mod in c["moduli"]), f"sheet {S} p={P}: every drop line's modulus equals {P}")
        drops = sorted(d["residue"] for c in chunks for d in c["drops"])
        details = [d for c in chunks for d in c["drops"]]
        detail = ", ".join(f"{d['residue']}: nullity {d['nullity']}, lifts to t = {d['lift']}" for d in details) or "none"
        require(drops == sorted(scan["drops"]) and details == scan["drop_details"], f"sheet {S} p={P}: rank drops and drop details as recorded ({detail})")
        for tpt, rank in scan["checks"].items():
            res = residue_of(tpt, P)
            require(all(c["checks"].get(tpt) == (res, rank) for c in chunks),
                    f"sheet {S} p={P}: check point t = {tpt} logged at residue {res} (= num * den^-1 mod P) with rank {rank} in every chunk")
        if scan["attested"]:
            verdicts = [check_attest(c.get("attest", []), P, S, scanner_sha, side["matrix_sha256"], c["range"][0], c["range"][1]) for c in chunks]
            require(all(ok for ok, _ in verdicts),
                    f"sheet {S} p={P}: every chunk carries exactly one attestation block bound to its own invocation (20 {S} {P} T0 T1), matrix, timestamp and a 40-hex clean commit"
                    + ("" if all(ok for ok, _ in verdicts) else f" -- {[why for ok, why in verdicts if not ok][:2]}")
                    + (f" [{verdicts[0][1]}]" if verdicts else ""))
        else:
            require(all("attest" not in c for c in chunks), f"sheet {S} p={P}: recorded as unattested (no attestation blocks; retrospective association)")
    print()
    if fails:
        print(f"SPIN-19 SCANS: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("SPIN-19 SCANS CONSISTENT (" + "; ".join(f"sheet {s['sheet']} p={s['prime']}: {s['evaluated']} residues, drops {s['drops'] or 'none'}" for s in summary["scans"]) + " -- finite-field evidence per preregistration e6dc226; scanner/matrix association retrospective)")
