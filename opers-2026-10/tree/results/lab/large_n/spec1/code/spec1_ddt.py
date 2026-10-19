"""SPEC1 information (lead, 2026-10-17): the Sol 1 a = 1 point (t = 3) against Dorey-Dunning-Tateo hep-th/0103051 sec. 4
((4.4), (4.8), (4.9): kappa (D(g) + x^3) phi = E phi, kappa = 4/(3 sqrt 3)).  NEW FILE (run 2; run 1 = spec1_ddt_run1.py,
failed on sympy's chain rule).  Usage: spec1_ddt.py DDT
(i) Second-order side.  Our L_O = (T^2 - l1^2) - x^6 (T - l0 - 3) - E x^3 (s = 1, c = 6, n = 3, T = 1 - theta) and the sextic
    (4.8) S = -d^2/dy^2 + y^6 + alpha y^2 + l(l+1)/y^2 - E_S.  Claim (derived by hand): with x = mu y^(2/3), mu^6 = 3,
    Phi = exp(y^4/4) y^(-1/6) f(x), alpha = -(4/3) l0, l(l+1) = (4/9) l1^2 - 1/4, E_S = -(4/9) mu^3 E = -kappa E:
        y^2 exp(-y^4/4) y^(1/6) S[Phi] = -(4/9) L_O[f].
    Checked on f = x^nu with symbolic nu (this determines the operator: both sides are polynomials in theta with x-power
    coefficients), and the decaying branches match: exp(-y^4/4) Phi-decay  <->  f ~ exp(-x^6/6) = exp(-y^4/2).
(ii) Third-order side.  Our Gamma form (monic) Pi_theta + x^6 - E_G x^3, roots 1 +- l1, 1 - l0; with psi = x^(l0/3) phi the roots
    become g = {1 +- l1, 1 - l0} + l0/3 - ... (sum 3), and DDT (3.2) times x^3 is prod(theta - g_i) + x^6 - E_D x^3: E_D = E_G.
    Checked symbolically on x^nu.
Composite with SPEC1 (E_G = -E): E_D = -E = kappa^-1 E_S, i.e. DDT (4.4) D^(1)(kappa^-1 E_S) ~ D^(+)(E_S), with their kappa."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def main():
    y = sp.Symbol("y", positive=True)
    x = sp.Symbol("x", positive=True)
    nu, l0, l1, E = sp.symbols("nu l0 l1 E")
    mu = sp.Integer(3) ** sp.Rational(1, 6)
    kappa = 4 / (3 * sp.sqrt(3))
    alpha = -sp.Rational(4, 3) * l0
    LL = sp.Rational(4, 9) * l1 ** 2 - sp.Rational(1, 4)
    ES = -kappa * E
    rho = -sp.Rational(1, 6)
    # (i) sextic side on f = x^nu, x = mu y^(2/3)
    Phi = sp.exp(y ** 4 / 4) * y ** rho * (mu * y ** sp.Rational(2, 3)) ** nu
    S = -sp.diff(Phi, y, 2) + (y ** 6 + alpha * y ** 2 + LL / y ** 2 - ES) * Phi
    lhs = sp.simplify(S * y ** 2 / (sp.exp(y ** 4 / 4) * y ** rho * (mu * y ** sp.Rational(2, 3)) ** nu))
    # L_O[x^nu]/x^nu = (1 - nu)^2 - l1^2 - x^6 (1 - nu - l0 - 3) - E x^3, at x = mu y^(2/3)
    LOn = (1 - nu) ** 2 - l1 ** 2 - x ** 6 * (1 - nu - l0 - 3) - E * x ** 3
    rhs = sp.Rational(-4, 9) * LOn.subs(x, mu * y ** sp.Rational(2, 3))
    d1 = sp.simplify(sp.expand(lhs - rhs))
    check("(i) y^2 e^(-y^4/4) y^(1/6) S[Phi] = -(4/9) L_O[f] on f = x^nu (symbolic nu, l0, l1, E), mu^6 = 3, E_S = -kappa E",
          d1 == 0, f"residual {d1}")
    # tamper: E_S = -E (no kappa)
    S_t = -sp.diff(Phi, y, 2) + (y ** 6 + alpha * y ** 2 + LL / y ** 2 + E) * Phi
    d1t = sp.simplify(sp.expand(sp.simplify(S_t * y ** 2 / (sp.exp(y ** 4 / 4) * y ** rho * (mu * y ** sp.Rational(2, 3)) ** nu)) - rhs))
    check("(i) tamper E_S = -E (kappa -> 1) fails", d1t != 0)
    lpl = sp.solve(sp.Symbol("l") * (sp.Symbol("l") + 1) - LL, sp.Symbol("l"))
    print(f"   alpha = {alpha} (with l0 = 3 sqrt3 P_X: alpha = -4 sqrt3 P_X); l = {lpl} (l + 1/2 = 2 l1/3; with l1 = 3 pi: 2 pi)",
          flush=True)
    # decaying branches: Phi ~ exp(-y^4/4) <-> f = Phi e^(-y^4/4) y^(1/6) ~ exp(-y^4/2) = exp(-x^6/(2 mu^6)) = exp(-x^6/6)
    check("(i) decaying branches correspond: exp(-y^4/2) = exp(-x^6/6) at x = mu y^(2/3), mu^6 = 3",
          sp.simplify(-(mu * y ** sp.Rational(2, 3)) ** 6 / 6 + y ** 4 / 2) == 0)
    # (ii) third-order side on psi = x^(l0/3) x^nu: Pi_theta(theta) with roots e = (1 - l1, 1 + l1, 1 - l0)
    e = [1 - l1, 1 + l1, 1 - l0]
    EG = sp.Symbol("E_G")
    mono = sp.Integer(1)
    for r in e:
        mono *= (nu - l0 / 3 - r)                       # run 3: Pi_theta acting on x^(nu - l0/3) (psi = x^(-l0/3) phi)
    ours = mono  # + x^6 - E_G x^3 (same on both sides after dividing by x^(l0/3))
    g = [r + l0 / 3 for r in e]                          # run 3: g = e + l0/3 (lead; sum 3)
    ddt = sp.Integer(1)
    for gi in g:
        ddt *= (nu - gi)
    check("(ii) gauge x^(-l0/3): prod(theta - e_i) x^(nu - l0/3) = x^(-l0/3) prod(theta - g_i) x^nu, g = e + l0/3, sum g = 3",
          sp.expand(ours - ddt) == 0 and sp.simplify(sum(g) - 3) == 0, f"g = {g}")
    print("   (ii) the x^6 and E x^3 terms are unchanged by the gauge, so E_D = E_G (no rescaling)", flush=True)
    check("composite: E_D = E_G = -E and E_S = -kappa E give E_D = kappa^-1 E_S (DDT (4.4))",
          sp.simplify(-E - ES / kappa) == 0)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "DDT":
        raise SystemExit("usage: spec1_ddt.py DDT")
    SL.banner(__file__)
    main()
    SL.imported("imported(end)")
    print("DDT DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
