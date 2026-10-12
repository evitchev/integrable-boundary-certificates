"""EXC4 B2 (post hoc, t = 2, bare (P,Q) = (1/2,1/3)): the monic I_3 quintic over the five sigma-invariant four-point solutions (gb_elim output)
against the quintic factor of the level-two singlet block's characteristic polynomial (d2_sol3_t2_P1_2_Q1_3.json).  Exit 0 iff equal."""
import json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
TAG = 't2_P1-2_Q1-3'
out = open('gb_elim_%s.out' % TAG).read()
x, e = sp.symbols('x e')
pe_ = sp.Poly(sp.sympify(re.search(r'eliminant ee:\n(.*)', out).group(1).replace('^', '**').replace('ee', 'x'), locals={'x': x}), x)
pe_ = sp.Poly(pe_.as_expr() / pe_.LC(), x)
pi_ = sp.sympify(re.search(r'eliminant ii:\n(.*)', out).group(1).replace('ii', 'x'), locals={'x': x})
d2 = json.load(open('d2_sol3_t2_P1_2_Q1_3.json'))
facs = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in d2['factors']]
quint = [f for f in facs if f.degree() == 5][0]; quint = sp.Poly(sp.expand(quint.as_expr() / quint.LC()), e)
vdim = int(re.search(r'vdim:\n(\d+)', out).group(1))
print('t = 2, (P,Q) = (1/2,1/3): vdim %d ordered = %d unordered configurations; I_1 = %s (vacuum 31/36 + 4 = 175/36)' % (vdim, vdim // 8, sp.solve(pi_, x)))
same = sp.expand(pe_.as_expr().subs(x, e) - quint.as_expr()) == 0
print('   oper quintic : %s' % pe_.as_expr()); print('   data quintic : %s' % quint.as_expr()); print('   EQUAL: %s' % same)
sys.exit(0 if same else 2)
