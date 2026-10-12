"""EXC4 B regression: the one-pair (sigma) system from h_four_sol3.py in 'pair' mode at t = 2, bare (P,Q) = (1/2,1/3) must give the level-one
I_3 data quadratic (d1_sol3_t2.json) and I_1 = vacuum + 2."""
import json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
e, P, Q = sp.symbols('e P Q')
out = open('h_four_t2_P1-2_Q1-3_pair.out').read()
pol = sp.Poly(sp.sympify(re.search(r'Je\[1\]=(.*)', out).group(1).replace('^', '**').replace('ee', 'e'), locals={'e': e}), e)
pol = sp.Poly(pol.as_expr() / pol.LC(), e)
dd = json.load(open('d1_sol3_t2.json'))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.expand((Mx / vac.coeff_monomial(P**4) - e * sp.eye(2)).det().subs({P: sp.Rational(1, 2), Q: sp.Rational(1, 3)}))
cp = sp.Poly(cp, e); cp = sp.Poly(cp.as_expr() / cp.LC(), e)
print('pair-mode I_3 eliminant:', pol.as_expr()); print('data level-one quadratic :', cp.as_expr()); print('equal:', sp.expand(pol.as_expr() - cp.as_expr()) == 0)
