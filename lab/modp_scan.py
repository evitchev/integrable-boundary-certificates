"""Lab: mod-p residue scan of a committed stratification matrix.
Parses wolfram/lab/cyl_strata_w{W}_sol{SOL}.m (rows of polynomials in t),
evaluates the matrix at every residue t in GF(p) except the parametrization
poles t = +-1 (t = 0, the ramification point (N,s) = (25,0), IS scanned --
Codex round 2, 2026-08-28: it was omitted before), computes the rank with
python-flint, and reports every residue where the rank drops below the
generic value, lifting each by rational reconstruction (small height).

What a single-prime scan proves and does not prove (Codex round 3
wording): at a calibrated good prime -- one at which the baseline ranks
equal the generic rank -- a rational jump point t0 = a/b whose reduction
is defined (p does not divide b) and non-polar (t0 is not +-1 mod p)
shows up at its residue; two rational points can share a residue (a
second prime separates them); an ALGEBRAIC jump point is seen only when
its minimal polynomial has a root mod p (use several primes).  It is a
discovery instrument: every candidate gets exact certification in the
open stack, and absence of drops is evidence, not proof.

Usage: modp_scan.py W SOL P [T0 T1]   (default residue range: all of
GF(P); the poles t = 1 and t = P-1 are skipped inside the loop)."""
import re, sys, time
import numpy as np
import flint
from fractions import Fraction as F

W, SOL, P = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
# optional residue range [T0, T1) for parallel chunks (default: all of GF(P); the poles 1 and P-1 are skipped in the loop)
T0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
T1 = int(sys.argv[5]) if len(sys.argv) > 5 else P
ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
path = f"{ROOT}/wolfram/lab/cyl_strata_w{W}_sol{SOL}.m"
t0 = time.time()
txt = open(path).read().strip()
if not (txt.startswith("{\n{") and txt.endswith("}\n}")):
    raise SystemExit("unexpected .m framing")           # explicit, survives -O
rows = txt[3:-3].split("},\n{")            # "{\n{" ... "}\n}" framing; the old [2:-2] left a stray "{"
# Strict grammar (Codex round 3: the old regex never proved it consumed
# the whole entry): entry = term ((+|-) term)*, term = int[*t[^k]] | t[^k]
TERM = r"(?:\d+\*t\^\d+|\d+\*t|\d+|t\^\d+|t)"
ENTRY = re.compile(rf"^[+-]?{TERM}(?:[+-]{TERM})*$")
term = re.compile(r"([+-]?)(?:(\d+)\*t\^(\d+)|(\d+)\*t|(\d+)|t\^(\d+)|t)")
nrows, ncols = len(rows), None
coef = {}                                   # degree -> dense int64 array (mod P)
for i, r in enumerate(rows):
    ents = r.split(", ")
    if ncols is None: ncols = len(ents)
    if len(ents) != ncols:
        raise SystemExit(f"row {i} has {len(ents)} entries, expected {ncols}")
    for j, e in enumerate(ents):
        if e == "0": continue
        s = e.replace(" ", "")
        if not ENTRY.fullmatch(s):
            raise SystemExit(f"entry ({i},{j}) does not parse as a polynomial in t: {e!r}")
        pos = 0
        for m in term.finditer(s):
            if m.start() != pos:
                raise SystemExit(f"entry ({i},{j}): unconsumed text at {pos}: {s!r}")
            pos = m.end()
            sgn = m.group(1)
            if m.group(2) is not None:   c, deg = int(m.group(2)), int(m.group(3))
            elif m.group(4) is not None: c, deg = int(m.group(4)), 1
            elif m.group(5) is not None: c, deg = int(m.group(5)), 0
            elif m.group(6) is not None: c, deg = 1, int(m.group(6))
            else:                        c, deg = 1, 1
            if sgn == "-": c = -c
            coef.setdefault(deg, np.zeros((nrows, ncols), dtype=np.int64))[i, j] += c % P
        if pos != len(s):
            raise SystemExit(f"entry ({i},{j}): trailing text after {pos}: {s!r}")
maxdeg = max(coef)
for k in coef: coef[k] %= P
print(f"parsed {nrows}x{ncols}, max degree {maxdeg}, parse {time.time()-t0:.1f}s", flush=True)

def rank_at(t):
    M = np.zeros((nrows, ncols), dtype=np.int64); tp = 1
    for k in range(maxdeg + 1):
        if k in coef: M = (M + coef[k] * tp) % P
        tp = (tp * t) % P
    return flint.nmod_mat(nrows, ncols, M.ravel().tolist(), P).rank()

def lift(r):
    """rational reconstruction of residue r with small height"""
    for b in range(1, 400):
        a = (r * b) % P
        if a > P // 2: a -= P
        if abs(a) < 400: return F(a, b)
    return None

t1 = time.time()
generic = rank_at(2), rank_at(5), rank_at(7)
print("generic ranks at t=2,5,7:", generic, f"({(time.time()-t1)/3:.3f}s each)")
g = max(generic)
drops = []
t1 = time.time()
evaluated = 0
for t in range(T0, T1):
    if t in (1, P - 1):                      # parametrization poles t = +-1
        continue
    evaluated += 1
    r = rank_at(t)
    if r < g: drops.append((t, r, lift(t)))
print(f"scanned residues [{T0},{T1}): {evaluated} evaluated (poles skipped) in {time.time()-t1:.0f}s; rank drops: {len(drops)}")
for t, r, q in drops:
    print(f"  t = {t} (mod {P})  rank {r} (nullity {ncols - r})  lifts to t = {q}")
for name, q in (("t=3", F(3)), ("t=9/7", F(9, 7)), ("t=11/9", F(11, 9)), ("t=2", F(2))):
    res = (q.numerator * pow(q.denominator, -1, P)) % P
    print(f"  check {name}: residue {res}, rank {rank_at(res)}")
