"""VIR5 pre-seal, ODE side only: the Mellin-WKB engine against the verified chain engine at
integer k (k = 2 and k = 4, spins 1..9; k = 3 for the missing spin).  Exit 0 iff identical."""
import sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
import wkb_lib as WL
PX, PI = sp.symbols('PX PI')
t0 = time.time()
fails = 0
print('Gamma-ratio coefficients c_j(k): k = 2:', MW.gamma_ratio_coeffs(2, 3), ' [exact product: (u - 1/2)(u + 1/2) = u^2 - 1/4]')
print('                                 k = 4:', MW.gamma_ratio_coeffs(4, 3), ' [(u^2 - 1/4)(u^2 - 9/4) = u^4 - (5/2) u^2 + 9/16]')
for k in (2, 4, 3):
    n = k + 4
    new, odd = MW.family_monic(k, 10, PX, PI)
    print('k = %d: Mellin-WKB solved to order 10 (%.0fs); odd orders vanish: %s' % (k, time.time() - t0, odd), flush=True)
    s_ = sp.Rational(n - 1, 2)
    frozen = [s_ + sp.Rational(n, k) * sp.Rational(2 * j - (k - 1), 2) for j in range(k)]
    ex = [s_ + WL.l0, s_ - WL.l0, s_ + WL.l1, s_ - WL.l1] + frozen
    pr = WL.Chain([('D', ex[i] - i) for i in range(n)] + [('p', -1)], '1', -1 if n % 2 == 0 else +1, n, (WL.F(1, n), -1))
    pr.solve(10)
    c2 = n**2 * sp.Rational(k + 1, k)
    for s in (1, 3, 5, 7, 9):
        Rk = sp.cancel(pr.R(s + 1)[0].subs(WL.M, sp.Rational(1, k)))
        Pl = sp.Poly(sp.expand(Rk), WL.l0, WL.l1)
        Rm = sp.Poly(sp.expand(sum(cf * c2 ** ((a + b) // 2) * PX**a * PI**b for (a, b), cf in Pl.terms())), PX, PI)
        lead = Rm.coeff_monomial(PX**(s + 1))
        if lead == 0:
            same = new[s][0] is None
            print('   spin %d: chain engine has no charge; Mellin engine: %s' % (s, 'none either' if same else 'HAS one'))
        else:
            old = sp.Poly(sp.expand(Rm.as_expr() / lead), PX, PI)
            same = new[s][0] is not None and sp.expand(old.as_expr() - new[s][0].as_expr()) == 0
            print('   spin %d: monic polynomials identical: %s (%d monomials)' % (s, same, len(old.terms())))
        fails += not same
    fails += not all(odd.values())
print('VALIDATION', 'FAILED (%d)' % fails if fails else 'OK: the Mellin-WKB engine reproduces the chain engine at integer k')
sys.exit(1 if fails else 0)
