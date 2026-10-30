"""B5 (i),(ii): the fusion eigenvalue as a product of two binomial polynomials; the records' W_U, W^(5) eigenvalues versus binom(2P+(p-1), 2p-1). Exact sympy."""
import sys, hashlib, os
import sympy as sp
print('script sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
PX, PI, z, PU, PV = sp.symbols('P_X pi z P_U P_V')
ok = True
def chk(name, expr):
    global ok
    z_ = sp.expand(expr) == 0; print(f'{name}: {"PASS" if z_ else "FAIL"}'); ok &= z_
def binom_poly(xx, n):  # binomial(xx, n) as the polynomial prod_{k=0}^{n-1}(xx-k)/n!
    return sp.prod([xx - k for k in range(n)]) / sp.factorial(n)
print('== (i) E_i = F_m(P_X+pi) F_m(P_X-pi), F_m(z) = prod_{j=-(m-1)}^{m-1}(z - j p) = p^(2m-1)(2m-1)! binom(z/p + m-1, 2m-1), m = i+1 = q/2, p^2 = 1/m')
for m in range(2, 9):
    i = m - 1; q = 2*m; p = 1/sp.sqrt(m)
    E = sp.prod([(PX + PI - j*p)*(PX - PI - j*p) for j in range(-i, i+1)])
    F = lambda zz: sp.prod([zz - j*p for j in range(-(m-1), m)])
    B = lambda zz: p**(2*m-1) * sp.factorial(2*m-1) * binom_poly(zz/p + m - 1, 2*m - 1)
    chk(f'q={q} (m={m}): E_i - F_m(P_X+pi)F_m(P_X-pi)', E - F(PX+PI)*F(PX-PI))
    chk(f'q={q} (m={m}): F_m(z) - p^(2m-1)(2m-1)! binom(z/p+m-1, 2m-1)', F(z) - B(z))
    chk(f'q={q} (m={m}): E_i - B(P_X+pi) B(P_X-pi)', E - B(PX+PI)*B(PX-PI))
print('== (ii) records: WUW1 (q=4, p_Ad=2): eigen(W_U) = (4/3) P_U (P_U^2 - 1/4), P_U = (P_X+pi)/sqrt2; lattice alpha = 2U, <alpha,alpha> = 4; background: T_U = (1/2)(dU)^2 + (1/2) d^2U, kappa = 1')
eWU = sp.Rational(4,3)*PU*(PU**2 - sp.Rational(1,4))
chk('eigen(W_U) == binom(2 P_U + 1, 3) exactly', eWU - binom_poly(2*PU + 1, 3))
print('   shift p-1 = 1 = 2 a kappa = 2*(1/2)*1 ->', 2*sp.Rational(1,2)*1 == 1)
print('== (ii) records: WUW2 (q=6, p_Ad=3): eigen(W^(5)) = (1/15) P_V (P_V^2-1)(4P_V^2-1), P_V = (sqrt3/2)(P_X+pi); alpha = 2V, <V,V> = kappa = 3/2, <alpha,alpha> = 6; T_U = (1/3)(dV)^2 + (2/3) d^2V')
eW5 = sp.Rational(1,15)*PV*(PV**2 - 1)*(4*PV**2 - 1)
chk('eigen(W^(5)) == binom(2 P_V + 2, 5) exactly', eW5 - binom_poly(2*PV + 2, 5))
print('   shift p-1 = 2 = 2 a kappa = 2*(2/3)*(3/2) ->', 2*sp.Rational(2,3)*sp.Rational(3,2) == 2)
print('== (ii) the records\' product constants: eigen(G_6) = -(2/9) E_1 (WUW1 P3), eigen(G_10) = -(81/4800) E_2 (WUW2 P3), versus the product of the two one-boson eigenvalues')
p1 = 1/sp.sqrt(2); E1 = sp.prod([(PX + PI - j*p1)*(PX - PI - j*p1) for j in (-1, 0, 1)])
PUx = (PX + PI)/sp.sqrt(2); PWx = (PI - PX)/sp.sqrt(2)   # WUW1 eigen(): P_U -> (P_X+pi)/sqrt2, P_W -> (pi-P_X)/sqrt2
prod1 = eWU.subs(PU, PUx) * eWU.subs(PU, PWx)
chk('q=4: eigen(W_U) eigen(W_W) == -(2/9) E_1', prod1 + sp.Rational(2,9)*E1)
p2 = 1/sp.sqrt(3); E2 = sp.prod([(PX + PI - j*p2)*(PX - PI - j*p2) for j in range(-2, 3)])
PVx = sp.sqrt(3)/2*(PX + PI); PZx = sp.sqrt(3)/2*(PI - PX)
prod2 = eW5.subs(PV, PVx) * eW5.subs(PV, PZx)
chk('q=6: eigen(W^(5)_V) eigen(W^(5)_Z) == -(81/4800) E_2', prod2 + sp.Rational(81,4800)*E2)
print('   (so the records\' G eigenvalues are exactly the products of the two binomials; consistent with the zero mode of H_U H_W acting on the grade-0 vector as (H_U)_0 (H_W)_0)')
print('== general m: binom(2P + m - 1, 2m-1) = 2^(2m-1)/(2m-1)! prod_{j=-(m-1)}^{m-1} (P - j/2): the root string P = j/2')
P = sp.Symbol('P')
for m in range(2, 9):
    chk(f'm={m}: binom(2P+m-1,2m-1) - 2^(2m-1)/(2m-1)! prod_j (P - j/2)', binom_poly(2*P + m - 1, 2*m - 1) - sp.Integer(2)**(2*m-1)/sp.factorial(2*m-1)*sp.prod([P - sp.Rational(j, 2) for j in range(-(m-1), m)]))
print('ALL PASS' if ok else 'SOME FAIL'); sys.exit(0 if ok else 1)
