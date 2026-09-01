"""Spin-17 (weight-18) ad_{I_3} kernel dimension of one cylindrical Solution
at one rational point of the curve, exact, open stack.  Usage:
  cyl_spin17_kdim.py SOL N S     (N, S as fractions a/b)
Prints one line: sol N s dim seconds."""
import sys, time
from fractions import Fraction as F
from mixed_engine import cyl_P4, gen_mixed_basis, genuine_kernel_mixed

sol = int(sys.argv[1]); N = F(sys.argv[2]); s = F(sys.argv[3])
if s * s != (N - 25) * (N - 1):
    sys.exit("point not on the curve")
t0 = time.time()
ker = genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(18), N)
print(f"sol={sol} N={N} s={s} dim={len(ker)} time={time.time()-t0:.0f}s", flush=True)
