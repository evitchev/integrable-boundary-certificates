"""NUM1 post-hoc (labelled): at M = 5, the extracted e^-3 coefficient vs  lim_(M->5) C_6^WKB(M)  +  alpha * Mellin(6),
Mellin(s) = Gamma(s/2) Gamma(s/2 + nu) Gamma(1/2 - s/2) / (4 sqrt(pi) Gamma(1 + nu - s/2)), nu = l + 1/2.  Limit by the symmetric average at M = 5 +- 1e-30."""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from rwkb import wkb, integral_coeff
mp.mp.dps = 90
LOG = open('run_diag_resonance.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
J = json.load(open('num1_dps60.json'))
PTS = {'A': (mp.mpf(1), mp.mpf(1) / 3), 'B': (mp.mpf(-3) / 2, mp.mpf(7) / 10), 'C': (mp.mpf(0), mp.mpf(1) / 3)}
eps = Fr(1, 10 ** 30); Sp = wkb(Fr(5) + eps, 6); Sm = wkb(Fr(5) - eps, 6); S5 = wkb(5, 6)
def mellin(s, nu): return mp.gamma(s / 2) * mp.gamma(s / 2 + nu) * mp.gamma(mp.mpf(1) / 2 - s / 2) / (4 * mp.sqrt(mp.pi) * mp.gamma(1 + nu - s / 2))
for nm, (al, l) in PTS.items():
    lam = l * (l + 1); nu = l + mp.mpf(1) / 2
    lim = -(integral_coeff(Sp[6], Fr(5) + eps, al, lam)[0] + integral_coeff(Sm[6], Fr(5) - eps, al, lam)[0]) / 2
    at5 = -integral_coeff(S5[6], 5, al, lam)[0]
    T = al * mellin(mp.mpf(6), nu); pred = lim + T; ext = mp.mpf(J['E2'][nm]['6'][0])
    dg = -mp.log10(abs(pred - ext) / abs(ext))
    say(f'{nm}: WKB at M = 5 exactly {mp.nstr(at5, 22)}; limit M -> 5 {mp.nstr(lim, 22)}; alpha Mellin(6) {mp.nstr(T, 22)}')
    say(f'    limit + Mellin = {mp.nstr(pred, 30)}  vs extracted {mp.nstr(ext, 30)}: {mp.nstr(dg, 4)} digits (E2 fit stability at grade 6 ~ 27-28 digits)')
    say(f'    closed form check: alpha Mellin(6) == -(4/15) alpha nu (nu^2-1)(nu^2-4): {mp.nstr(T + mp.mpf(4) / 15 * al * nu * (nu**2 - 1) * (nu**2 - 4), 3)}')
