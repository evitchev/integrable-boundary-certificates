"""control of the numerical route: one point, random complex starts -> the four exact solutions (z^2 = 1/324, 49/5832), r21 = -b."""
import sys, random
from mpmath import mp, mpf, mpc
sys.dont_write_bytecode = True
import num_two as NT
a0, a1, b = mpf(1) / 2, mpf(1) / 3, mpf(2) / 3
random.seed(7)
found = {}
for trial in range(120):
    x0 = [(mpc(random.uniform(-1, 1), random.uniform(-1, 1)),) for _ in range(7)]
    pt = tuple(v[0] for v in x0)
    sol = NT.newton([pt], a0, a1, b)
    if sol is None:
        continue
    z = sol[0][6]
    if abs(z) < mpf(10) ** (-12):
        continue
    key = (round(float(z.real), 8), round(float(z.imag), 8))
    found.setdefault(key, sol[0])
for key, s in sorted(found.items()):
    print('z = %s   z^2 = %s   r21 = %s' % (mp.nstr(s[6], 15), mp.nstr(s[6] ** 2, 15), mp.nstr(s[0], 12)))
print('distinct solutions found: %d (expected 4); 1/324 = %s, 49/5832 = %s' % (len(found), mp.nstr(mpf(1) / 324, 15), mp.nstr(mpf(49) / 5832, 15)))
