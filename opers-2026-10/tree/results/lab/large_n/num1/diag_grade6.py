"""NUM1 post-hoc diagnostic (labelled): the grade-6 (e^-3, integer power) mismatch at alpha != 0.  Hypothesis: a 0 x infinity term at M = 5 exactly
(zero coefficient times a Gamma pole); the correct WKB value is the limit M -> 5 (as VIR5b's residue rule).  Computes C_6(5 + eps) for eps -> 0 and
compares with the M = 5 value and the numerically extracted coefficient (num1_dps60.json)."""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from rwkb import wkb, integral_coeff
mp.mp.dps = 50
LOG = open('run_diag_grade6.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
J = json.load(open('num1_dps60.json'))
PTS = {'A': (mp.mpf(1), mp.mpf(1) / 3), 'B': (mp.mpf(-3) / 2, mp.mpf(7) / 10), 'C': (mp.mpf(0), mp.mpf(1) / 3)}
S5 = wkb(5, 6)
for nm, (al, l) in PTS.items():
    lam = l * (l + 1); ext = mp.mpf(J['E2'][nm]['6'][0])
    c5 = -integral_coeff(S5[6], 5, al, lam)[0]
    say(f'point {nm}: extracted grade 6 = {mp.nstr(ext, 25)};  WKB at M = 5 exactly = {mp.nstr(c5, 25)}')
    for k in (6, 10, 14, 20):
        Me = Fr(5) + Fr(1, 10 ** k); Se = wkb(Me, 6)
        ce = -integral_coeff(Se[6], Me, al, lam)[0]
        say(f'   M = 5 + 1e-{k}: WKB grade 6 = {mp.nstr(ce, 25)};  |minus extracted| = {mp.nstr(abs(ce - ext), 3)}')
# which terms carry the pole near M = 5
Se = wkb(Fr(5) + Fr(1, 10**10), 6)
from fractions import Fraction
for (a, b), p in sorted(S5[6].items()):
    A_ = a / 10; B_ = -b - a / 10
    if (A_.denominator == 1 and A_ <= 0) or (B_.denominator == 1 and B_ <= 0): say(f'   M = 5 term at a Gamma pole: x^{a} W^{b}, coeff {p}')
say('   (terms whose coefficient vanishes identically at M = 5 do not appear in the M = 5 expression; their limit is the residue contribution)')
