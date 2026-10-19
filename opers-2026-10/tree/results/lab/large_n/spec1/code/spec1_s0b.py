"""SPEC1 S0b (lead amendment): (H) at the added generic point Sol 1 a = 3 (t = 7) for momenta A, B, plus the parameter row.
Imports spec1_s0 unchanged.  NEW FILE.  Usage: spec1_s0b.py S0b"""
import hashlib
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_s0 as S  # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "S0b":
        raise SystemExit("usage: spec1_s0b.py S0b")
    for f in (__file__, S.__file__):
        print(f"imported {os.path.abspath(f)} sha256 {hashlib.sha256(open(f, 'rb').read()).hexdigest()}", flush=True)
    sol, a = "Sol1", 3
    t = S.t_of(sol, a)
    pr = S.params(sol, t)
    tt = F(t)
    N = (tt * tt - 25) / (tt * tt - 1)
    print(f"{sol} a={a}: t={t} N={N} n={pr['n']} s={pr['s']} sigma={pr['sigma']} c={pr['c']}", flush=True)
    for nm, (l0, l1) in (("A", (F(1, 3), F(3, 5))), ("B", (F(9, 10), F(1, 5)))):
        th, base, bad = S.check_H(sol, pr, l0, l1)
        print(f"  point {nm}: theta_j={[str(x) for x in th]}, frozen base={base} (lattice base + {pr['sigma']} Z); "
              f"(H) {'HOLDS' if not bad else 'FAILS ' + str(bad)}", flush=True)
    for sol2, a2 in (("Sol3", 1), ("Sol3", 3), ("Sol1", 1)):
        t2 = F(S.t_of(sol2, a2))
        print(f"  info: {sol2} a={a2} t={t2} N={(t2 * t2 - 25) / (t2 * t2 - 1)}", flush=True)
    print("S0b DONE", flush=True)
