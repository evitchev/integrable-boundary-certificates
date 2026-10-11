"""Sol 1's certified spin-3 vacuum eigenvalue (the record's vev_w4_points.json, Solution 1, 29 fibres) vs
Fateev-Lukyanov hep-th/0510271v1 eq. (95), I_3^vac of quantum AKNS, at n = -(t-1)/(t+1), with ONE fixed
dictionary Q_FL^2 = 2 P_X^2 (= -2X), P_FL^2 = 2 pi^2 (= 2(rho^2 - Y)).  (95) is read with (5n+4)/16: the printed
(5n+4) lacks the 1/16 and then contradicts FL's own general top-term formula (checked below at three n).
Controls: the literal reading and a tamper n + 1/10 must fail."""
import json, sys, sympy as sp
D = json.load(open(sys.argv[1]))['1']
P2, Q2, X, Y, n = sp.symbols('P2 Q2 X Y n')
def fl(nv, c5=16):
    return (nv/16*(P2**2 - P2 + sp.Rational(1, 12)) + 3*nv/8*(P2 - sp.Rational(1, 6))*(Q2 - sp.Rational(1, 6))
            + (5*nv + 4)/c5*(Q2**2 - Q2 + sp.Rational(1, 12)) + (nv + 3)*(2*nv + 1)/240)
def run(nf, c5, label):
    ex = tot = 0
    for ts, coef in D.items():
        if not isinstance(coef, dict): continue
        t = sp.Rational(ts); nv = nf(t); rho2 = 2/(t**2 - 1)
        rec = sum(sp.Rational(v)*X**int(k.split(',')[0])*Y**int(k.split(',')[1]) for k, v in coef.items())
        f = sp.expand(fl(nv, c5).subs({Q2: -2*X, P2: 2*(rho2 - Y)})); lead = f.coeff(X, 2).subs(Y, 0)
        tot += 1; ex += sp.expand(f/lead - rec) == 0
    print(f'{label}: {ex}/{tot} fibres exact in all 6 coefficients')
run(lambda t: -(t - 1)/(t + 1), 16, '(95) with (5n+4)/16')
run(lambda t: -(t - 1)/(t + 1), 1, 'CONTROL (95) literal')
run(lambda t: -(t - 1)/(t + 1) + sp.Rational(1, 10), 16, 'CONTROL tamper n + 1/10')
T = {l: 2**(4 - 2*l)*(n + 2)**sp.Rational(4 - 2*l, 2)*sp.gamma(sp.Rational(1, 2) - l + 3*(n + 1)/(n + 2))/(n**(3 - l)*sp.factorial(4 - 2*l)*sp.factorial(l)) for l in range(3)}
for nv in (sp.Rational(-1, 3), sp.Rational(1, 2), sp.Rational(3)):
    print(f'FL general top formula at n={nv}: Q^4/P^4 = {sp.N((T[0]/T[2]).subs(n, nv), 15)}; (5n+4)/n = {sp.N((5*nv + 4)/nv, 15)}; literal 16(5n+4)/n = {sp.N(16*(5*nv + 4)/nv, 15)}')
