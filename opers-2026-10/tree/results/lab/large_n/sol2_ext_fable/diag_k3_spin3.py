"""SOL2F post-hoc diagnostic (labelled): k = 3 (t = 11/5) spin-3 single-coefficient difference.  Which coefficient; is the certified cyl(2,4) regular at
t = 11/5 (compare value vs limit t -> 11/5); the family's spin 3 as t varies along the family is not available (one member per k), so compare also the
family's k = 3 spin 3 with the certified spin 3 at t = 11/5 +- 1/1000 trend."""
import sys
from fractions import Fraction as F
import sympy as sp
from common import *
LOG = open(os.path.join(HERE, 'run_diag_k3_spin3.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
tv = sp.Rational(11, 5); cert = certified(2, tv)
g = classG(2, tv); M = F(1, 3); Cv = 2 * g['D'] * (1 + M)
sv = gwkb.wkb(blocks_of(g, M), g['s'], M, 4); r = thwkb.integrate(sv[4], g['D'] * M)
ch = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in r[2].items()))
m = monic(ch, 2)
say(f'family k = 3 spin 3 (monic): {m}')
say(f'certified Sol 2 t = 11/5 spin 3: {cert[3]}')
say(f'difference: {sp.expand(m - cert[3])}')
raw = sp.expand(cyl(2, 4))
say(f'cyl(2, 4) symbolic: {raw}')
for c, mon in zip(sp.Poly(raw, X, Y).coeffs(), sp.Poly(raw, X, Y).monoms()):
    say(f'   {mon}: {sp.factor(c)}  | at 11/5: {sp.factor(c).subs(tc, tv)}')
