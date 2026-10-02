"""SOL2F post-hoc (labelled): k = 3, t = 11/5 spin 3 compared PROJECTIVELY.  Certified = lim_(t -> 11/5) (5t - 11) cyl(2, 4), converted by the VIR3-R dictionary
(regular); WKB = raw family spin-3 charge.  Equal iff proportional (all 2x2 cross products vanish).  Also: does the WKB's P_X^4 coefficient vanish?"""
import sys
from fractions import Fraction as F
import sympy as sp
from common import *
LOG = open(os.path.join(HERE, 'run_diag_k3_projective.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
tv = sp.Rational(11, 5); p2 = 2 * (tv - 1) / (tv + 1); rho2 = (p2 - 2)**2 / (4 * p2)
reg = sp.expand(sum(sp.limit(sp.factor(c) * (5 * tc - 11), tc, tv) * Pa**(2 * i) * Qa**(2 * j) for (i, j), c in zip(sp.Poly(sp.expand(cyl(2, 4)), X, Y).monoms(), sp.Poly(sp.expand(cyl(2, 4)), X, Y).coeffs())))
cert = sp.expand(reg.subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(rho2 - PI**2)}))
g = classG(2, tv); M = F(1, 3); Cv = 2 * g['D'] * (1 + M)
sv = gwkb.wkb(blocks_of(g, M), g['s'], M, 4); r = thwkb.integrate(sv[4], g['D'] * M)
w = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in r[2].items()))
say(f'certified regularised spin 3: {cert}')
say(f'WKB family k = 3 spin 3 (raw): {w}')
pc = sp.Poly(cert, PX, PI).as_dict(); pw = sp.Poly(w, PX, PI).as_dict()
mons = sorted(set(pc) | set(pw)); ref = next(m for m in mons if pc.get(m, 0) != 0)
prop = all(sp.expand(pc.get(m, 0) * pw.get(ref, 0) - pw.get(m, 0) * pc[ref]) == 0 for m in mons) and pw.get(ref, 0) != 0
say(f'P_X^4 coefficient: certified {pc.get((4, 0), 0)}, WKB {pw.get((4, 0), 0)}')
say(f'PROPORTIONAL (projectively equal, {len(mons)} monomials): {prop}')
for M2 in (F(-1, 4),):
    Cv2 = 2 * g['D'] * (1 + M2); sv2 = gwkb.wkb(blocks_of(g, M2), g['s'], M2, 4); r2 = thwkb.integrate(sv2[4], g['D'] * M2)
    w2 = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * (sp.Rational(Cv2.numerator, Cv2.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv2.numerator, Cv2.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in r2[2].items()))
    pw2 = sp.Poly(w2, PX, PI).as_dict()
    say(f"dual M' = {M2}: PROPORTIONAL: {all(sp.expand(pc.get(m, 0) * pw2.get(ref, 0) - pw2.get(m, 0) * pc[ref]) == 0 for m in mons)}")
tam = cert + sp.Rational(1, 1000) * PI**2
say(f'scoring tamper (+pi^2/1000) proportional: {all(sp.expand(sp.Poly(tam, PX, PI).as_dict().get(m, 0) * pw.get(ref, 0) - pw.get(m, 0) * sp.Poly(tam, PX, PI).as_dict()[ref]) == 0 for m in mons)} (must be False)')
