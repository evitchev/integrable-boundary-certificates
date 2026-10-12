"""RED2 (B2, exact operator-level check): Sol 3's k = -2 member S and the gauged tensor product T of two rotated Coulomb problems are related by a
Mellin multiplier.  Convention: phi = int phihat(nu) x^nu dnu; an operator sum_i x^i p_i(theta) gives the recurrence sum_i p_i(nu - i) phihat(nu - i) = 0.
Claim: with phihat_T(nu) = m(nu) psihat_S(nu), m(nu)/m(nu-1) = rho (nu - 1/2)/(nu - 2) (m = rho^nu Gamma(nu + 1/2)/Gamma(nu - 1)),
the T-recurrence divided by m(nu) IS the S-recurrence, identically in nu, for a suitable x-scale (rho), with lam_pm = ((l0 -+ l1)^2 - 1)/4, g = 1/2."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
x, D, nu = sp.symbols('x D nu')
lp, lm, E, l0, l1, C, rho, mu = sp.symbols('lam_p lam_m E l0 l1 C rho mu')
T = (lm**2 - 2*lm*lp - 2*lm + lp**2 - 2*lp + 2*(2*lm + 2*lp - 1)*D + (-2*lm - 2*lp + 5)*D**2 - 4*D**3 + D**4
     + x*(2*D - 4*D**2) + x**2*(4*E*D + 4*E*D**2))          # copied from b2_coulomb.log
sub = {lp: ((l0 - l1)**2 - 1) / 4, lm: ((l0 + l1)**2 - 1) / 4}
Tg = sp.expand(sum(c * mu**i * x**i * (D + sp.Rational(1, 2))**k for (i, k), c in sp.Poly(sp.expand(T.subs(sub)), x, D).terms()))
S = sp.expand((((D - sp.Rational(1, 2))**2 - l0**2) * ((D - sp.Rational(1, 2))**2 - l1**2)) - C * x * D * (D - 1) + C * E * x**2 * D * (D - 1))
def coeffs(Op): P = sp.Poly(Op, x); return {i: sp.expand(P.coeff_monomial(x**i)) for i in range(3)}
pT, pS = coeffs(Tg), coeffs(S)
print('x^0 parts equal (P4):', sp.expand(pT[0] - pS[0]) == 0)
r = lambda n: rho * (n - sp.Rational(1, 2)) / (n - 2)       # m(n)/m(n-1)
# T-recurrence / m(nu):  pT0(nu) psihat(nu) + pT1(nu-1) [m(nu-1)/m(nu)] psihat(nu-1) + pT2(nu-2) [m(nu-2)/m(nu)] psihat(nu-2)
t1 = sp.cancel(pT[1].subs(D, nu - 1) / r(nu)); t2 = sp.cancel(pT[2].subs(D, nu - 2) / (r(nu) * r(nu - 1)))
s1 = sp.expand(pS[1].subs(D, nu - 1)); s2 = sp.expand(pS[2].subs(D, nu - 2))
# C is not independent: x -> nu x rescales (C, E) -> (nu C, E nu); solve with C free (C = 4 expected, mu = rho)
sol = [s_ for s_ in sp.solve([sp.numer(sp.together(t1 - s1)), sp.numer(sp.together(t2 - s2))], [rho, mu, C], dict=True) if s_.get(C, 1) != 0 and not any(v.has(nu) for v in s_.values())]
print('solve for (rho, mu) making the recurrences identical:', sol)
ok = False
for s_ in sol:
    a = sp.simplify(t1.subs(s_) - s1.subs(s_)); b = sp.simplify(t2.subs(s_) - s2.subs(s_))
    print('  identity in nu:', a == 0 and b == 0, '(residuals', a, b, ')'); ok |= (a == 0 and b == 0)
# negative control: multiplier with the wrong Gamma shift (nu + 1/2 -> nu + 1)
r_bad = lambda n: rho * (n - 1) / (n - 2)
t1b = sp.cancel(pT[1].subs(D, nu - 1) / r_bad(nu)); t2b = sp.cancel(pT[2].subs(D, nu - 2) / (r_bad(nu) * r_bad(nu - 1)))
solb = [s_ for s_ in sp.solve([sp.numer(sp.together(t1b - s1)), sp.numer(sp.together(t2b - s2))], [rho, mu, C], dict=True) if s_.get(rho, 1) != 0]
print('NEGATIVE CONTROL (wrong multiplier Gamma(nu+1)/Gamma(nu-1)) constant solutions (must be none):', solb)
# negative control 2: unrotated momenta lam_pm = (l0^2 - 1)/4, (l1^2 - 1)/4 -> P4 mismatch
sub2 = {lp: (l0**2 - 1) / 4, lm: (l1**2 - 1) / 4}
T2 = sp.expand(sum(c * x**i * (D + sp.Rational(1, 2))**k for (i, k), c in sp.Poly(sp.expand(T.subs(sub2)), x, D).terms()))
print('NEGATIVE CONTROL (unrotated momenta) x^0 parts equal:', sp.expand(coeffs(T2)[0] - pS[0]) == 0)
sys.exit(0 if ok else 1)
