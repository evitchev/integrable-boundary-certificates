"""EXC1 control: the Schrodinger-form engine with the alpha term (potential-term feature of mellin_pot.py) reproduces
the VIR7 Gamma-form S1 charges (odd spins 1..7) at t = 9/4 and t = 7/3 under lam^2 = (M+1) pi^2, alpha^2 = 4M(M+1) PX^2.
ODE-vs-ODE; no excited-state data.  Exit 0 iff all agree."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
import mellin_pot as MP
from fractions import Fraction as F
PX, PI = sp.symbols('PX PI')
ok = True
for tstr, fn in (('9/4', '../vir7/s1_pred_t9_4.json'), ('7/3', '../vir7/s1_pred_t7_3.json')):
    t = sp.Rational(tstr)
    M = (t + 3) / (t - 1)
    K = 8
    h = [MP.ONE, MP.ZERO, MP.P_(-MP.l1**2)] + [MP.ZERO] * (K - 1)
    eng = MP.MellinWKB(2, M, h)
    Mf = F(int(M.p), int(M.q))
    V = {1: {(Mf - 1, F(-1)): MP.P_(MP.l0)}}
    eng.solve(K, V)
    pred = json.load(open(fn))['odd_spins']
    for i in range(2, K + 1, 2):
        R, good = eng.R(i)
        assert good and eng.last_half == 0, (good, eng.last_half)     # alpha-odd part at an even order: a total derivative
        Pl = sp.Poly(R, MP.l0, MP.l1)
        ex = sum(c * (4 * M * (M + 1)) ** sp.Rational(a, 2) * (M + 1) ** sp.Rational(b, 2) * PX**a * PI**b for (a, b), c in Pl.terms())
        Q = sp.Poly(sp.expand(ex), PX, PI)
        Q = sp.Poly(sp.expand(Q.as_expr() / Q.coeff_monomial(PX**i)), PX, PI)
        ref = pred[str(i - 1)]['coefficients']
        same = all(Q.coeff_monomial(PX**int(k.split(',')[0]) * PI**int(k.split(',')[1])) == sp.Rational(v) for k, v in ref.items()) and len(Q.terms()) == len(ref)
        ok &= same
        print('t = %s (M = %s) spin %d: Schrodinger engine with alpha term == Gamma-form S1 prediction: %s' % (tstr, M, i - 1, same), flush=True)
print('vacuum control', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
