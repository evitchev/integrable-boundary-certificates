"""EXC2 U3: I_3 and I_5 (monic) on the SELF-ADJOINT two-point solutions (z2 = -z1) at one momentum point.  Opens NO data.
Needs h_red_<TAG>.json and h_red2_<TAG>.ms (h_run.sh steps 1-3).  Writes u3p_<TAG>.json:
  vdim, I3 = [trace, det], I5 = [trace, det] of the two-state block (monic eliminants), A3 (= even-spin charge J_3) zero or not.
Usage: u3_point.py <t> <a0> <a1> [flipB1]     (flipB1: scoring tamper used by the comparison script only)"""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr, a0s, a1s = sys.argv[1], sys.argv[2], sys.argv[3]
flip = len(sys.argv) > 4
t = sp.Rational(tstr); a0 = sp.Rational(a0s); a1 = sp.Rational(a1s)
k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4; nb = n / b
TAG = 't%s_%s_%s' % (tstr.replace('/', '_'), a0s.replace('/', '-'), a1s.replace('/', '-'))
t0 = time.time()
names = 'ww,r21,r11,r12,r01,r02,r03,z,s21,s11,s12,s01,s02,s03,z2,qq,ii,ee'.split(',')
sy = {nm: sp.Symbol(nm) for nm in names}
rho, sig, ff, gg = sp.symbols('rho sig ff gg')
z, z2, qq, ee, ww = sy['z'], sy['z2'], sy['qq'], sy['ee'], sy['ww']
sub = {sy[k_]: sp.sympify(v, locals=sy) for k_, v in json.load(open('h_red_%s.json' % TAG))['sub'].items()}
s = open('h_red_%s.ms' % TAG).read().split('\n')
P = [sp.sympify(p.strip().replace('^', '**'), locals=sy) for p in '\n'.join(s[2:]).split(',\n') if p.strip()]
r11 = sp.solve(P[1], sy['r11'])[0]; s11 = sp.solve(P[4], sy['s11'])[0]
fin = {sy['r11']: r11.subs(sy['r12'], rho * z), sy['s11']: s11.subs(sy['s12'], sig * z2), sy['r12']: rho * z, sy['s12']: sig * z2}
full = lambda nm: sp.together(sy[nm].subs(sub).subs(sub).subs(fin).subs(fin))
R = {nm: full(nm) for nm in ('r21', 'r11', 'r12', 'r01', 's21', 's11', 's12', 's01')}
assert R['r21'] == -b and R['s21'] == -b
q2 = -(a0**2 + a1**2) + R['r21'] + R['s21']
q3 = R['r11'] + R['s11']
q4 = a0**2 * a1**2 + R['r01'] + R['s01']
q12 = (R['r21'] * z - 4 * z) + (R['s21'] * z2 - 4 * z2)            # r_21 z + r_22, r_22 = -4 z
q11 = (R['r11'] * z + R['r12']) + (R['s11'] * z2 + R['s12'])
A2, A3, A4, B0, B1, B2 = sp.symbols('A2 A3 A4 B0 B1 B2')
gc = json.load(open('u3_engine_t%s.json' % tstr.replace('/', '_')))
loc = {'A2': A2, 'A3': A3, 'A4': A4, 'B0': B0, 'B1': B1, 'B2': B2}
cB4 = sp.sympify(gc['cB4'], locals=loc); cB6 = sp.sympify(gc['cB6'], locals=loc)
assert not cB6.has(B0) and not cB4.has(B0) and not cB4.has(B1) and not cB4.has(B2)
x0s, x1s = sp.symbols('x0s x1s')
vacmap = {A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s, B1: 0, B2: 0}
lead6 = sp.Poly(sp.expand(cB6.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**3) * b**3
lead4 = sp.Poly(sp.expand(cB4.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**2) * b**2
valmap = {A2: nb**2 * q2, A3: -nb**3 * q3, A4: nb**4 * q4, B2: nb**5 * q12, B1: (1 if flip else -1) * nb**6 * q11}
def rel(sym, expr):
    num, den = sp.fraction(sp.together(sym - expr))
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(num))[1]:
        if fac in (z, z2):
            continue
        core *= fac ** mult
    return sp.expand(core)
erel = rel(ee, cB4.subs(valmap) / lead4); frel = rel(ff, cB6.subs(valmap) / lead6); grel = rel(gg, q3)
s2 = open('h_red2_%s.ms' % TAG).read().split('\n')
P2 = [p.strip() for p in '\n'.join(s2[2:]).split(',\n') if p.strip()]
# the e-relation written by h_reduce2.py and the one rebuilt here must be proportional (two code paths, one map)
e_old = sp.sympify(P2[4].replace('^', '**'), locals={**sy, 'rho': rho, 'sig': sig})
assert sp.expand(e_old * sp.Poly(erel, ee).LC() - erel * sp.Poly(e_old, ee).LC()) == 0, 'e-relation mismatch between h_reduce2 and u3_point'
txt = lambda p_: str(sp.expand(p_)).replace('**', '^')
fn = 'u3p_%s%s.sing' % (TAG, '_flip' if flip else '')
open(fn, 'w').write('ring R = 0,(ww,rho,sig,z,z2,qq,ee,ff,gg),dp;\noption(redSB);\nideal I =\n' + ',\n'.join(P2[:4] + P2[5:7] + ['z + z2', txt(erel), txt(frel), txt(grel)]) +
                    ';\nideal G = std(I);\nprint("vdim:"); vdim(G);\nideal Je = eliminate(G, ww*rho*sig*z*z2*qq*ff*gg);\nprint("eliminant in ee:"); Je;\n'
                    'ideal Jf = eliminate(G, ww*rho*sig*z*z2*qq*ee*gg);\nprint("eliminant in ff:"); Jf;\n'
                    'ideal Jg = eliminate(G, ww*rho*sig*z*z2*qq*ee*ff);\nprint("eliminant in gg:"); Jg;\n'
                    'ideal Jz = eliminate(G, ww*rho*sig*z2*qq*ee*ff*gg);\nprint("eliminant in z:"); Jz;\nquit;\n')
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=3000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
vdim = int(re.search(r'vdim:\n(-?\d+)', res).group(1))
def monic(tag, var):
    m = re.search(r'%s\[1\]=(.*)' % tag, res)
    pol = sp.Poly(sp.sympify(m.group(1).replace('^', '**'), locals={var: sp.Symbol(var)}), sp.Symbol(var))
    return sp.Poly(pol.as_expr() / pol.LC(), sp.Symbol(var))
pe, pf, pg = monic('Je', 'ee'), monic('Jf', 'ff'), monic('Jg', 'gg')
out = {'t': tstr, 'a0': a0s, 'a1': a1s, 'PX2': str(a0**2 / b), 'PI2': str(a1**2 / b), 'vdim_ordered': vdim,
       'I3_poly': str(pe.as_expr()), 'I5_poly': str(pf.as_expr()), 'q3_poly': str(pg.as_expr()),
       'z_eliminant': re.search(r'Jz\[1\]=(.*)', res).group(1)}
json.dump(out, open('u3p_%s%s.json' % (TAG, '_flip' if flip else ''), 'w'), indent=1)
print('%s%s: vdim %d; I_3: %s ; I_5: %s ; q3 = r11 + s11: %s   (%.0fs)' % (TAG, ' FLIP' if flip else '', vdim, pe.as_expr(), pf.as_expr(), pg.as_expr(), time.time() - t0), flush=True)
