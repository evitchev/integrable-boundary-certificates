"""EXC4 A1 (second pass): I_1 and I_3 of the one-point solutions of the sixth-order form (f_six_t2_<a0>_<a1>.out Groebner basis).
Effective sextic Q6 + sum_j r_j1 th^j -> symbol H_eff(T) = c^6 Q6_eff(-T/c) = T^6 + h1 T^5 + ... + h6; the Gamma-form engine (mellin_gen6,
generic h's) gives cB_2, cB_4 as polynomials in the h's; Singular eliminates ii (monic I_1) and ee (monic I_3) over the solutions.
Usage: f_six_charges.py <a0> <a1>.  Exit 0 iff Delta I_1 == 1 on every solution (the sealed half-level); reports I_3 vs the level-one data."""
import hashlib, json, re, subprocess, sys
import sympy as sp
sys.dont_write_bytecode = True
import mellin_gen6 as MG
assert hashlib.sha256(open('SEAL_EXC4.md', 'rb').read()).hexdigest() == open('SEAL_EXC4.sha256').read().split()[0]
a0s, a1s = sys.argv[1], sys.argv[2]; a0 = sp.Rational(a0s); a1 = sp.Rational(a1s); A, B = a0**2, a1**2
t = sp.Integer(2); k = sp.Integer(2); b = k / (k + 1); n = k + 4; M = 1 / k; c = n * (1 + M); f = sp.Rational(3, 2) / c; rho2 = 2 / (t * t - 1)
th, T = sp.symbols('th T')
K = 4
hs = [MG.h1, MG.h2, MG.h3, MG.h4, MG.h5, MG.h6]
hl = [MG.ONE] + [MG.P_(h) for h in hs] + [MG.ZERO] * (K + 1 - 7) if K + 1 > 7 else [MG.ONE] + [MG.P_(h) for h in hs][:K]
eng = MG.MellinWKB(n, M, [MG.ONE] + [MG.P_(h) for h in hs][:K]); eng.solve(K)
cB = {}
for i in (2, 4):
    Rv, good = eng.R(i); cB[i] = sp.expand(Rv)
# vacuum normalisation: H_vac(T) = (T^2 - l0^2)(T^2 - l1^2)(T^2 - 9/4), l0^2 = n^2 P^2/b, l1^2 = n^2 (Q^2 + rho2)/b
P, Q = sp.symbols('P Q')
Hv = sp.Poly(sp.expand((T**2 - n**2 * P**2 / b) * (T**2 - n**2 * (Q**2 + rho2) / b) * (T**2 - sp.Rational(9, 4))), T)
hv = {hs[j - 1]: Hv.coeff_monomial(T**(6 - j)) for j in range(1, 7)}
lead = {i: sp.Poly(sp.expand(cB[i].subs(hv)), P, Q).coeff_monomial(P**i) for i in (2, 4)}
vac = {i: sp.expand(cB[i].subs(hv).subs({P: sp.sqrt(A / b), Q: sp.sqrt(B / b - rho2)}) / lead[i]) for i in (2, 4)}
print('vacuum monic I_1 = %s, I_3 = %s at (a0, a1) = (%s, %s)' % (vac[2], vac[4], a0s, a1s))
# effective sextic with the residues
r = {j: sp.Symbol('r%d1' % j) for j in range(0, 5)}
Qeff = (th**2 - A) * (th**2 - B) * (th**2 - f**2) + sum(r[j] * th**j for j in range(0, 5))
He = sp.Poly(sp.expand(c**6 * Qeff.subs(th, -T / c)), T)
assert He.coeff_monomial(T**6) == 1
he = {hs[j - 1]: He.coeff_monomial(T**(6 - j)) for j in range(1, 7)}
I1 = sp.expand(cB[2].subs(he) / lead[2]); I3 = sp.expand(cB[4].subs(he) / lead[4])
print('monic I_1 of the deformed operator = vacuum + (%s)' % sp.factor(I1 - vac[2]))
out = open('f_six_t2_%s_%s.out' % (a0s.replace('/', '-'), a1s.replace('/', '-'))).read()
basis = [m.group(1) for m in re.finditer(r'G\[\d+\]=(.*)', out)]
ring_line = open('f_six_t2_%s_%s.sing' % (a0s.replace('/', '-'), a1s.replace('/', '-'))).readline()
vars_ = [v for v in re.search(r'\((.*)\),dp', ring_line).group(1).split(',') if v != 'ww']
def txt(p_):
    # integer coefficients (Singular parses 'x^2/4' as x^(2/4)): clear denominators
    pp = sp.Poly(sp.expand(p_), *sorted(sp.expand(p_).free_symbols, key=str)) if sp.expand(p_).free_symbols else sp.Poly(sp.expand(p_), sp.Symbol('dummy'))
    _, pc = pp.clear_denoms()
    return str(pc.as_expr()).replace('**', '^')
ii, ee = sp.symbols('ii ee')
lines = ['ring R = 0,(%s,ww,ii,ee),dp;' % ','.join(vars_), 'option(redSB);', 'ideal I =', ',\n'.join(basis + [txt(ii - I1), txt(ee - I3)]) + ';',
         'ideal G = std(I);', 'print("vdim:"); vdim(G);', 'ideal Ji = eliminate(G, %s*ww*ee);' % '*'.join(vars_), 'print("eliminant ii:"); Ji;',
         'ideal Je = eliminate(G, %s*ww*ii);' % '*'.join(vars_), 'print("eliminant ee:"); Je;', 'quit;']
fn = 'f_six_charges_%s_%s.sing' % (a0s.replace('/', '-'), a1s.replace('/', '-'))
open(fn, 'w').write('\n'.join(lines) + '\n')
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=3000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res); print(res[:800])
x = sp.Symbol('x')
pi_ = sp.Poly(sp.sympify(re.search(r'Ji\[1\]=(.*)', res).group(1).replace('^', '**').replace('ii', 'x'), locals={'x': x}), x)
pe_ = sp.Poly(sp.sympify(re.search(r'Je\[1\]=(.*)', res).group(1).replace('^', '**').replace('ee', 'x'), locals={'x': x}), x)
pe_ = sp.Poly(pe_.as_expr() / pe_.LC(), x)
shifts = [rt - vac[2] for rt in sp.solve(pi_.as_expr(), x)]
print('I_1 over the solutions minus vacuum: %s  (sealed: 1 = half a level)' % shifts)
print('I_3 (monic) over the solutions: %s ; numeric %s' % (pe_.as_expr(), pe_.nroots(n=10)))
dd = json.load(open('d1_sol3_t2.json'))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vacd = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
e = sp.Symbol('e')
cp = sp.Poly(sp.expand((Mx / vacd.coeff_monomial(P**4) - e * sp.eye(2)).det()), P, Q, e)
cpv = sp.Poly(sum(cf * (A / b)**(m[0] // 2) * (B / b - rho2)**(m[1] // 2) * e**m[2] for m, cf in cp.terms()), e); cpv = sp.Poly(cpv.as_expr() / cpv.LC(), e)
print('level-one data quadratic at this point: %s ; numeric %s' % (cpv.as_expr(), cpv.nroots(n=10)))
ok = all(sh == 1 for sh in shifts)
print('A1b half-level:', 'HOLDS' if ok else 'FAILS'); sys.exit(0 if ok else 2)
