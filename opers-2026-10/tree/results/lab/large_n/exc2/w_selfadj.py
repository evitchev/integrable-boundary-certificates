"""EXC2 (post-seal, labelled): which two-point solutions are the level-1 states?  Hypothesis SA: exactly the self-adjoint ones,
z2 = -z1 (the adjoint composed with xi -> -xi maps a solution at (z1, z2) to one at (-z1, -z2)); equivalently the odd-order
(even-spin) charges vanish: J_3 = -A3 = 0, J_5 = A2 A3/3 + 15 A3/2 - B2/7 = 0 [B2 ~ z1 + z2].
Test (no I_5 data): the two-point ideal + (z + z2): number of solutions and the I_3 eliminant, against the data quadratic.
Also: on the data-selected factor, are z + z2 and r11 + s11 in the ideal?
Usage: w_selfadj.py <t> <a0> <a1>.  Exit 0 iff the self-adjoint I_3 eliminant == data quadratic; 2 otherwise."""
import hashlib, json, re, subprocess, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr, a0s, a1s = sys.argv[1], sys.argv[2], sys.argv[3]
t = sp.Rational(tstr); a0 = sp.Rational(a0s); a1 = sp.Rational(a1s)
k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1)
TAG = 't%s_%s_%s' % (tstr.replace('/', '_'), a0s.replace('/', '-'), a1s.replace('/', '-'))
s2 = open('h_red2_%s.ms' % TAG).read().split('\n')
P2 = [p.strip() for p in '\n'.join(s2[2:]).split(',\n') if p.strip()]
fn = 'w_selfadj_%s.sing' % TAG
open(fn, 'w').write('ring R = 0,(ww,rho,sig,z,z2,qq,ee),dp;\noption(redSB);\nideal I =\n' + ',\n'.join(P2 + ['z + z2']) +
                    ';\nideal G = std(I);\nprint("vdim:"); vdim(G);\nideal Je = eliminate(G, ww*rho*sig*z*z2*qq);\nprint("eliminant in ee:"); Je;\n'
                    'ideal Jz = eliminate(G, ww*rho*sig*z2*qq*ee);\nprint("eliminant in z:"); Jz;\nprint("rho - sig in ideal (0 = yes):"); reduce(rho - sig, G);\nquit;\n')
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=1500).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
e = sp.Symbol('e')
vdim = int(re.search(r'vdim:\n(-?\d+)', res).group(1))
pol = sp.Poly(sp.sympify(re.search(r'Je\[1\]=(.*)', res).group(1).replace('^', '**').replace('ee', 'e'), locals={'e': e}), e)
pol = sp.Poly(pol.as_expr() / pol.LC(), e)
P, Q = sp.symbols('P Q')
dd = json.load(open('d1_sol3_t%s.json' % tstr.replace('/', '_')))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.Poly(sp.expand((Mx / vac.coeff_monomial(P**4) - e * sp.eye(2)).det()), P, Q, e)
PX2, PI2 = a0**2 / b, a1**2 / b
cpv = sp.Poly(sum(c * PX2 ** (m[0] // 2) * (PI2 - rho2) ** (m[1] // 2) * e ** m[2] for m, c in cp.terms()), e)
cpv = sp.Poly(sp.expand(cpv.as_expr() / cpv.LC()), e)
print('t = %s, xi-exponents (%s, %s): self-adjoint (z2 = -z) two-point solutions: vdim %d (ordered)' % (tstr, a0, a1, vdim))
print('  I_3 eliminant on them : %s' % pol.as_expr())
print('  data level-1 polynomial: %s' % cpv.as_expr())
print('  ' + re.search(r'eliminant in z:\n(.*)', res).group(1)[:300])
print('  rho - sig reduced: ' + res.split('rho - sig in ideal (0 = yes):\n')[1].split('\n')[0][:200])
same = sp.expand(pol.as_expr() - cpv.as_expr()) == 0
print('  self-adjoint eliminant == data polynomial: %s' % same)
sys.exit(0 if same else 2)
