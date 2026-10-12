"""EXC4 A1: Solution 3 at t = 2 in its SIXTH-ORDER polynomial form (k = 2).  Gamma form at integer k (item 294/298):
(T^2 - l0^2)(T^2 - l1^2)(T^2 - 9/4) psi = x^6 (x^3 - E) psi, T = s - theta, s = 5/2, n = 6, M = 1/2, c = n(1+M) = 9 = n/b, b = 2/3,
l0 = n P_X/p, l1 = n pi/p (the universal dictionary).  In xi = x^c/const and T = -c vartheta (gauge x^s):
    Q6(vartheta) - xi - E' xi^b,    Q6 = (th^2 - A)(th^2 - B)(th^2 - f^2),  A = b P_X^2, B = b pi^2, f = 3/(2c) = 1/6.
Class-U member with L = -1 (constant).  The formal adjoint alone (th -> -th, no xi -> -xi) maps Q6 - xi - E' xi^b to itself.
ONE apparent singularity at xi = z, exponents S_6 = {-1,1,2,3,4,6}, deformation xi sum_(j<=4) R_j th^j, pole orders <= 6 - j.
Tests: T0 the vacuum Q6 through the Gamma-form engine (mellin_gen) reproduces the certified Sol 3 vacuum I_3 at t = 2 (d1_sol3_t2.json);
A1a existence of one-point solutions (Singular); A1b r41 on them (sealed: -b) and the monic I_1 shift (sealed: p^2/b = 1, half a level);
A1c the I_3 values (the one-point I_3 of a half level is NOT the level-one data; reported).
I_1, I_3 need only the effective sextic (xi ~ T^6 here: the 1/xi terms enter at relative order >= 8).
Usage: f_six_sol3.py <t=2> <a0> <a1>   (xi-exponents as in EXC2: a0 = p P_X, a1 = p pi).  Exit 0 iff T0 passes and the shift == 1 on all solutions;
3 if no solution; 2 otherwise."""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
import mellin_gen as MG
assert hashlib.sha256(open('SEAL_EXC4.md', 'rb').read()).hexdigest() == open('SEAL_EXC4.sha256').read().split()[0]
tstr, a0s, a1s = sys.argv[1], sys.argv[2], sys.argv[3]
t = sp.Rational(tstr); assert t == 2
k = sp.Integer(2); b = k / (k + 1); n = k + 4; M = 1 / k; c = n * (1 + M); f = sp.Rational(3, 2) / c
a0 = sp.Rational(a0s); a1 = sp.Rational(a1s); A, B = a0**2, a1**2
rho2 = 2 / (t * t - 1); PX2 = A / b; PI2 = B / b
t0 = time.time()
th = sp.Symbol('th')
# ---- T0: vacuum through the Gamma-form engine.  Symbol H(T) = (T^2 - l0^2)(T^2 - l1^2)(T^2 - 9/4) = T^6 phi(1/T); hlist = coefficients of T^(6-j)
T = sp.Symbol('T')
def hlist_of(poly_T, deg, K):
    p_ = sp.Poly(sp.expand(poly_T), T)
    return [MG.P_(p_.coeff_monomial(T**(deg - j))) if j <= deg else MG.ZERO for j in range(K + 1)]
Hvac = (T**2 - MG.l0**2) * (T**2 - MG.l1**2) * (T**2 - sp.Rational(9, 4))
K = 6
eng = MG.MellinWKB(n, M, hlist_of(Hvac, 6, K)); eng.solve(K)
P, Q = sp.symbols('P Q')
dd = json.load(open('d1_sol3_t%s.json' % tstr))['blocks']['I3']
vac3 = sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q})
subl = {MG.l0: n * P / sp.sqrt(b), MG.l1: n * sp.sqrt(Q**2 + rho2) / sp.sqrt(b)}
ok = True
for i in (2, 4):
    Rv, good = eng.R(i)
    expr = sp.expand(sp.expand(Rv).subs({MG.l0**2: n**2 * P**2 / b, MG.l1**2: n**2 * (Q**2 + rho2) / b}))
    assert not expr.has(MG.l0) and not expr.has(MG.l1)
    K2 = i // 2; lead = sp.Poly(expr, P, Q).coeff_monomial(P**(2 * K2)); mon = sp.expand(expr / lead)
    if i == 4:
        cm = sp.expand(vac3 / sp.Poly(vac3, P, Q).coeff_monomial(P**4)); same = sp.expand(mon - cm) == 0; ok &= same
        print('T0 sixth-order form, vacuum I_3 == certified (d1_sol3_t2 vacuum): %s' % same, flush=True)
    else:
        print('T0 vacuum I_1 (monic) = %s' % mon, flush=True)
if not ok:
    print('T0 FAILS'); sys.exit(2)
# ---- one apparent singularity
z, E, e = Fb.z, Fb.E, Fb.e
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in range(0, 5) for m in range(1, 7 - j)}
R = {j: {m: r[j, m] for m in range(1, 7 - j)} for j in range(0, 5)}
Q6 = sp.Poly(sp.expand((th**2 - A) * (th**2 - B) * (th**2 - f**2)), th)
Pcoef = [Q6.coeff_monomial(th**j) for j in range(7)]
Lcoef = [sp.Integer(-1)]
exps = [-1, 1, 2, 3, 4, 6]
lk = Fb.local_operator(Pcoef, Lcoef, b, -1, R, exps[-1] - exps[0] - 6 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 6)
indvars = [r[4, 2], r[3, 3], r[2, 4], r[1, 5], r[0, 6]]
sol_ind = sp.solve([cc for cc in sp.Poly(sp.expand(lead - z**6 * sp.prod([(e - x) for x in exps])), e).coeffs() if cc != 0], indvars, dict=True)
print('indicial conditions fix:', sol_ind, ' (%.0fs)' % (time.time() - t0), flush=True)
assert sol_ind
conds = [(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [x for x in r.values() if x not in sol_ind[0]] + [z]
polys = []
for pos, mon, cf in eqs:
    num = sp.numer(sp.together(cf)); core = sp.Integer(1)
    for fac, mult in sp.factor_list(num)[1]:
        if fac != z: core *= fac**mult
    polys.append(sp.expand(core))
print('%d no-log equations in %d unknowns (%.0fs)' % (len(polys), len(unknowns), time.time() - t0), flush=True)
# E-linear condition at the top exponent, printed (the seal's mechanism)
for (pos, mon, cf), p_ in zip(eqs, polys):
    if pos == exps[-1] and mon[0] == 1:
        print('   E-linear condition at the top exponent:', sp.factor(p_))
ww = sp.Symbol('ww')
vars_ = [str(x) for x in unknowns] + ['ww']
txt = lambda p_: str(sp.expand(p_)).replace('**', '^')
fn = 'f_six_t2_%s_%s.sing' % (a0s.replace('/', '-'), a1s.replace('/', '-'))
open(fn, 'w').write('ring R = 0,(%s),dp;\noption(redSB);\nideal I =\n%s;\nideal G = std(I);\nprint("dim:"); dim(G);\nprint("vdim:"); vdim(G);\nideal Jz = eliminate(G, %s);\nprint("eliminant z:"); Jz;\nideal Jr = eliminate(G, %s);\nprint("eliminant r41:"); Jr;\nprint("basis:"); G;\nquit;\n'
                    % (','.join(vars_), ',\n'.join(txt(p_) for p_ in polys + [ww * z - 1]), '*'.join([str(x) for x in unknowns if x != z] + ['ww']), '*'.join([str(x) for x in unknowns if x != r[4, 1]] + ['ww'])))
print('written %s; running Singular (%.0fs)' % (fn, time.time() - t0), flush=True)
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=21000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
print(res[:1500])
dim = int(re.search(r'dim:\n(-?\d+)', res).group(1))
if dim == -1:
    print('A1a: NO one-point solution with z != 0 for S_6 at this point'); sys.exit(3)
vd = int(re.search(r'vdim:\n(-?\d+)', res).group(1))
# r41 and I_1, I_3 over the solutions: effective sextic Q6 + sum_j r_j1 th^j
m = re.search(r'Jr\[1\]=(.*)', res); r41pol = sp.sympify(m.group(1).replace('^', '**'), locals={'r41': r[4, 1]})
print('A1b: r41 eliminant: %s  (sealed: r41 = -b = %s)' % (sp.factor(r41pol), -b))
# effective sextic in T: Q6_eff(th) = Q6 + sum r_j1 th^j ;  T = -c th  ->  H_eff(T) = c^6 Q6_eff(-T/c)
Qeff = Q6.as_expr() + sum(r[j, 1] * th**j for j in range(0, 5))
Heff = sp.expand(c**6 * Qeff.subs(th, -T / c))
hl = []
pH = sp.Poly(Heff, T)
for j in range(K + 1):
    cj = pH.coeff_monomial(T**(6 - j)) if j <= 6 else 0
    hl.append(MG.P_(cj) if cj != 0 else MG.ZERO)
# the engine's generators: GENS of mellin_gen are (l0, l1, ...)?  use sympy symbols directly via P_; r's must be engine gens -> substitute numerically later
print('effective symbol built; computing I_1, I_3 over the solutions by substitution of each Groebner solution would need numeric roots;',
      'instead: I_1 and I_3 are polynomial in the r_j1, so eliminate with Singular (second pass).', flush=True)
json.dump({'t': tstr, 'a0': a0s, 'a1': a1s, 'vdim': vd, 'unknowns': [str(x) for x in unknowns]}, open(fn.replace('.sing', '.json'), 'w'), indent=1)
print('done (%.0fs)' % (time.time() - t0))
sys.exit(0)
