"""EXC4 B2: the I_3 eliminant of the four-point sigma-invariant system (h_four_<TAG>_four.out) against the level-two singlet block at the
same point (d2_sol3_t<t>_P<P>_Q<Q>.json: characteristic polynomial = (e - level-shift) x quintic).  Also I_1 = vacuum + 4 (level two) and
the solution count.  Usage: chk_four.py <t> <P> <Q>.  Exit 0 iff the monic quintic factor of the data equals the I_3 eliminant and
I_1 = vacuum + 4 on all solutions."""
import json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
tstr, Ps, Qs = sys.argv[1:4]
TAG = 't%s_P%s_Q%s' % (tstr.replace('/', '_'), Ps.replace('/', '-'), Qs.replace('/', '-'))
out = open('h_four_%s_four.out' % TAG).read()
log = open('h_four_four_%s.log' % TAG).read() if TAG != 't2_P1-2_Q1-3' else open('h_four_four_t2.log').read()
x, e = sp.symbols('x e')
vd = int(re.search(r'vdim:\n(-?\d+)', out).group(1))
def mon(key, nm):
    p_ = sp.Poly(sp.sympify(re.search(r'%s\[1\]=(.*)' % key, out).group(1).replace('^', '**').replace(nm, 'x'), locals={'x': x}), x)
    return sp.Poly(p_.as_expr() / p_.LC(), x)
pi_, pe_ = mon('Ji', 'ii'), mon('Je', 'ee')
vac = sp.Rational(re.search(r'vacuum monic I_1 = (\S+);', log).group(1))
d2 = json.load(open('d2_sol3_t%s_P%s_Q%s.json' % (tstr.replace('/', '_'), Ps.replace('/', '_'), Qs.replace('/', '_'))))
facs = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in d2['factors']]
quint = [f for f in facs if f.degree() == 5][0]; quint = sp.Poly(sp.expand(quint.as_expr() / quint.LC()), e)
shifts = [r_ - vac for r_ in sp.solve(pi_.as_expr(), x)]
print('t = %s (P,Q) = (%s,%s): vdim %d (ordered); I_1 - vacuum over the solutions: %s (level two = 4); I_3 eliminant degree %d' % (tstr, Ps, Qs, vd, shifts, pe_.degree()))
same = sp.expand(pe_.as_expr().subs(x, e) - quint.as_expr()) == 0
div = sp.rem(sp.Poly(pe_.as_expr().subs(x, e), e), quint).is_zero if pe_.degree() >= 5 else False
print('   I_3 eliminant == data quintic: %s;  data quintic divides it: %s' % (same, div))
print('   eliminant:', pe_.as_expr()); print('   quintic  :', quint.as_expr())
sys.exit(0 if same and all(s_ == 4 for s_ in shifts) else 2)
