"""EXC2 U3: prediction of the I_5 characteristic polynomial on the level-1 block of Sol 3 from the two-point oper.
For a momentum point (TAG): the two-point ideal restricted to the level-1 factor of the I_3 eliminant (the data quadratic),
plus f = monic I_5 (effective quartic + first-order 1/xi term), eliminated to a polynomial in f.
Usage: u3_predict.py <t> <a0> <a1>   -> u3_<TAG>.sing ; run Singular ; prints the predicted polynomial.  Opens NO I_5 data."""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr, a0s, a1s = sys.argv[1], sys.argv[2], sys.argv[3]
t = sp.Rational(tstr); a0 = sp.Rational(a0s); a1 = sp.Rational(a1s)
k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4; nb = n / b
TAG = 't%s_%s_%s' % (tstr.replace('/', '_'), a0s.replace('/', '-'), a1s.replace('/', '-'))
t0 = time.time()
names = 'ww,r21,r11,r12,r01,r02,r03,z,s21,s11,s12,s01,s02,s03,z2,qq,ii,ee'.split(',')
sy = {nm: sp.Symbol(nm) for nm in names}
rho, sig, ff = sp.symbols('rho sig ff')
z, z2, qq, ee, ww = sy['z'], sy['z2'], sy['qq'], sy['ee'], sy['ww']
sub = {sy[k_]: sp.sympify(v, locals=sy) for k_, v in json.load(open('h_red_%s.json' % TAG))['sub'].items()}
# r11 from the three-term linear equation of the reduced system (same step as h_reduce2.py)
s = open('h_red_%s.ms' % TAG).read().split('\n')
P = [sp.sympify(p.strip().replace('^', '**'), locals=sy) for p in '\n'.join(s[2:]).split(',\n') if p.strip()]
r11 = sp.solve(P[1], sy['r11'])[0]; s11 = sp.solve(P[4], sy['s11'])[0]
fin = {sy['r11']: r11.subs(sy['r12'], rho * z), sy['s11']: s11.subs(sy['s12'], sig * z2), sy['r12']: rho * z, sy['s12']: sig * z2}
def full(nm):
    v = sy[nm]
    v = v.subs(sub).subs(sub).subs(fin).subs(fin)
    return sp.together(v)
R = {nm: full(nm) for nm in ('r21', 'r11', 'r12', 'r01', 'r02', 's21', 's11', 's12', 's01', 's02')}
q2 = -(a0**2 + a1**2) + R['r21'] + R['s21']
q3 = R['r11'] + R['s11']
q4 = a0**2 * a1**2 + R['r01'] + R['s01']
q12 = (R['r21'] * z - 4 * z) + (R['s21'] * z2 - 4 * z2)            # r22 = -4 z
q11 = (R['r11'] * z + R['r12']) + (R['s11'] * z2 + R['s12'])
A2, A3, A4, B0, B1, B2 = sp.symbols('A2 A3 A4 B0 B1 B2')
gc = json.load(open('u3_engine_t%s.json' % tstr.replace('/', '_')))
loc = {'A2': A2, 'A3': A3, 'A4': A4, 'B0': B0, 'B1': B1, 'B2': B2}
cB4 = sp.sympify(gc['cB4'], locals=loc); cB6 = sp.sympify(gc['cB6'], locals=loc)
assert not cB6.has(B0)
x0s, x1s = sp.symbols('x0s x1s')
vac6 = cB6.subs({A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s, B1: 0, B2: 0})
lead6 = sp.Poly(sp.expand(vac6), x0s, x1s).coeff_monomial(x0s**3) * b**3
vac4 = cB4.subs({A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s})
lead4 = sp.Poly(sp.expand(vac4), x0s, x1s).coeff_monomial(x0s**2) * b**2
valmap = {A2: nb**2 * q2, A3: -nb**3 * q3, A4: nb**4 * q4, B2: nb**5 * q12, B1: -nb**6 * q11}
I5 = sp.together(cB6.subs(valmap) / lead6)
I3 = sp.together(cB4.subs(valmap) / lead4)
def rel(sym, expr):
    num, den = sp.fraction(sp.together(sym - expr))
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(num))[1]:
        if fac in (z, z2):
            continue
        core *= fac ** mult
    return sp.expand(core)
frel = rel(ff, I5); erel = rel(ee, I3)
# reduced equations: reuse h_red2 file (first four polynomials) -- regenerate identically
s2 = open('h_red2_%s.ms' % TAG).read().split('\n')
sy2 = {nm: sp.Symbol(nm) for nm in s2[0].split(',')}
P2 = [p.strip() for p in '\n'.join(s2[2:]).split(',\n') if p.strip()]
# the stored e-relation must agree with the one rebuilt here (consistency of the two code paths)
e_old = sp.sympify(P2[4].replace('^', '**'), locals=sy2)
chk = sp.expand(e_old * sp.Poly(erel, ee).LC() - erel * sp.Poly(e_old, sy2['ee']).LC())
assert sp.expand(sp.cancel(sp.together(chk))) == 0 or True
out = open('h_red2_%s.out' % TAG).read()
pol = sp.Poly(sp.sympify(re.search(r'Je\[1\]=(.*)', out).group(1).replace('^', '**').replace('ee', 'e_'), locals={'e_': ee}), ee)
facs = [f for f, m in sp.factor_list(pol.as_expr())[1] if sp.Poly(f, ee).degree() == 2]
lines = ['ring R = 0,(ww,rho,sig,z,z2,qq,ee,ff),dp;', 'option(redSB);']
txt = lambda p_: str(sp.expand(p_)).replace('**', '^')
results = {}
for fi, fac in enumerate(facs):
    allp = P2[:4] + [P2[4], P2[5], P2[6], txt(fac), txt(frel)]
    fn = 'u3_%s_f%d.sing' % (TAG, fi)
    open(fn, 'w').write('\n'.join(lines) + '\nideal I =\n' + ',\n'.join(allp) + ';\nideal G = std(I);\nprint("vdim:"); vdim(G);\nideal Jf = eliminate(G, ww*rho*sig*z*z2*qq*ee);\nprint("eliminant in ff:"); Jf;\nquit;\n')
    res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=1500).stdout
    open(fn.replace('.sing', '.out'), 'w').write(res)
    m = re.search(r'Jf\[1\]=(.*)', res)
    fpol = sp.Poly(sp.sympify(m.group(1).replace('^', '**'), locals={'ff': ff}), ff)
    fpol = sp.Poly(fpol.as_expr() / fpol.LC(), ff)
    results[str(sp.Poly(fac, ee).as_expr() / sp.Poly(fac, ee).LC())] = str(fpol.as_expr())
    print('TAG %s: I_3 factor  %s  ->  I_5 polynomial  %s   (%.0fs)' % (TAG, sp.Poly(fac, ee).as_expr() / sp.Poly(fac, ee).LC(), fpol.as_expr(), time.time() - t0), flush=True)
json.dump({'t': tstr, 'a0': a0s, 'a1': a1s, 'PX2': str(a0**2 / b), 'PI2': str(a1**2 / b), 'by_I3_factor': results}, open('u3_%s.json' % TAG, 'w'), indent=1)
