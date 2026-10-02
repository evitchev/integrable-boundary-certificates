"""EXC1 level 2, data side: the level-2 O(N-1)-singlet block (6 states in the N-component Y Fock space) of I_3 for Sol 1
with the record's tool; characteristic polynomial at a rational bare momentum point (P, Q).
Usage: d2_level2.py <sol> <t> <P> <Q>"""
import hashlib, json, os, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
os.environ.setdefault('IB_CODE', '<repo>/code')
os.environ.setdefault('IB_LAB', '<repo>/lab')
sys.path.insert(0, os.environ['IB_LAB']); sys.path.insert(0, os.environ['IB_CODE'])
import excited_states as XS
import vir_lib as V
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
sol = int(sys.argv[1]); t = F(sys.argv[2]); Pv = sp.Rational(sys.argv[3]); Qv = sp.Rational(sys.argv[4])
N, s_ = V.curve_point(t)
t0 = time.time()
dens = V.P4(sol, N, s_)
st = XS.singlet_states(2)
print('level-2 singlet states:', st, flush=True)
H, Gm = XS.charge_matrix(dens, N, st)
print('matrices done (%.0fs)' % (time.time() - t0), flush=True)
sub = {XS.Pm: Pv, XS.Qm: Qv}
Hn, Gn = H.subs(sub), Gm.subs(sub)
Mx = Gn.inv() * Hn
H0, G0 = XS.charge_matrix(dens, N, XS.singlet_states(0))
vac = sp.Poly(sp.expand(H0[0, 0] / G0[0, 0]), XS.Pm, XS.Qm)
lead = vac.coeff_monomial(XS.Pm**4)
e = sp.Symbol('e')
cp = sp.Poly(sp.expand((Mx / lead - e * sp.eye(len(st))).det()), e)
fac = sp.factor_list(cp.as_expr())
shift = sp.expand(vac.as_expr().subs(XS.Qm, sp.sqrt(XS.Qm**2 + 4))).subs(sub) / lead      # level-shift law at L = 2
print('characteristic polynomial factors (degree): %s' % [(sp.Poly(f, e).degree(), m) for f, m in fac[1]])
lin = [sp.solve(f, e)[0] for f, m in fac[1] if sp.Poly(f, e).degree() == 1]
print('linear roots: %s ; level-shift value (Q^2 -> Q^2 + 4): %s' % (lin, sp.simplify(shift)))
out = {'sol': sol, 't': str(t), 'P': str(Pv), 'Q': str(Qv), 'charpoly': str(cp.as_expr()), 'factors': [str(f) for f, m in fac[1]], 'gram_det': str(Gn.det())}
json.dump(out, open('d2_sol%d_t%s_P%s_Q%s.json' % (sol, str(t).replace('/', '_'), str(Pv).replace('/', '_'), str(Qv).replace('/', '_')), 'w'), indent=1)
print('done (%.0fs)' % (time.time() - t0))
