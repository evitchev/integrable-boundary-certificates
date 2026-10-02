"""EXC2 U1/U2: read the I_3 eliminant of the two-point system from the Singular output and compare with the level-1 data.
Usage: h_check.py <t> <a0> <a1>     Exit 0 = the data polynomial divides the eliminant; 2 = not."""
import hashlib, json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr, a0s, a1s = sys.argv[1], sys.argv[2], sys.argv[3]
t = sp.Rational(tstr); a0 = sp.Rational(a0s); a1 = sp.Rational(a1s)
k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1)
TAG = 't%s_%s_%s' % (tstr.replace('/', '_'), a0s.replace('/', '-'), a1s.replace('/', '-'))
out = open('h_red2_%s.out' % TAG).read()
vdim = int(re.search(r'vdim:\n(\d+)', out).group(1))
e = sp.Symbol('e')
pol = sp.Poly(sp.sympify(re.search(r'Je\[1\]=(.*)', out).group(1).replace('^', '**').replace('ee', 'e'), locals={'e': e}), e)
P, Q = sp.symbols('P Q')
dd = json.load(open('d1_sol3_t%s.json' % tstr.replace('/', '_')))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.Poly(sp.expand((Mx / vac.coeff_monomial(P**4) - e * sp.eye(2)).det()), P, Q, e)
PX2, PI2 = a0**2 / b, a1**2 / b
assert all(m[0] % 2 == 0 and m[1] % 2 == 0 for m, _ in cp.terms())
cpv = sp.Poly(sum(c * PX2 ** (m[0] // 2) * (PI2 - rho2) ** (m[1] // 2) * e ** m[2] for m, c in cp.terms()), e)
cpv = sp.Poly(sp.expand(cpv.as_expr() / cpv.LC()), e)
print('t = %s, xi-exponents (%s, %s): PX^2 = %s, pi^2 = %s' % (tstr, a0, a1, PX2, PI2))
print('two-point system: zero-dimensional, %d ordered solutions (%d unordered); I_1 = Delta + 1 on all (r21 = s21 = -b)' % (vdim, vdim // 2))
print('I_3 eliminant: degree %d; factor degrees %s' % (pol.degree(), [sp.Poly(f, e).degree() for f, m in sp.factor_list(pol.as_expr())[1]]))
print('data level-1 characteristic polynomial: %s' % cpv.as_expr())
div = sp.rem(pol, cpv).is_zero
print('data polynomial divides the eliminant: %s' % div)
sys.exit(0 if div else 2)
