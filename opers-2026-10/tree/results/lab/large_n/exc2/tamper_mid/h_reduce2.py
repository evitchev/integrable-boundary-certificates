"""EXC2 addendum: second reduction of the two-point system.  r11 = (3 r12 - 14 z)/(9 z) [from the 3-term equation; at t = 2],
r12 = rho z, s12 = sig z2.  Unknowns (rho, sig, z, z2) and q = 1/(z - z2); 4 equations.  Writes h_red2_t2.sing / .ms."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
t0 = time.time()
TAG = sys.argv[1]
s = open('h_red_%s.ms' % TAG).read().split('\n')
names = s[0].split(',')
sy = {n: sp.Symbol(n) for n in names}
P = [sp.sympify(p.strip().replace('^', '**'), locals=sy) for p in '\n'.join(s[2:]).split(',\n') if p.strip()]
rho, sig = sp.symbols('rho sig')
z, z2, qq, ee, ww = sy['z'], sy['z2'], sy['qq'], sy['ee'], sy['ww']
r11 = sp.solve(P[1], sy['r11'])[0]; s11 = sp.solve(P[4], sy['s11'])[0]
sub = {sy['r11']: r11.subs(sy['r12'], rho * z), sy['s11']: s11.subs(sy['s12'], sig * z2), sy['r12']: rho * z, sy['s12']: sig * z2}
print('r11 =', sp.simplify(sub[sy['r11']]), '; s11 =', sp.simplify(sub[sy['s11']]))
out = []
for i in (0, 2, 3, 5, 8):
    ex = sp.together(P[i].subs(sub))
    num = sp.numer(ex)
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(num))[1]:
        if fac in (z, z2):
            continue
        core *= fac ** mult
    core = sp.expand(core)
    out.append(core)
    print('eq %d -> degree %d, %d terms' % (i, sp.Poly(core, rho, sig, z, z2, qq, ee).total_degree(), len(sp.Poly(core, rho, sig, z, z2, qq, ee).terms())), flush=True)
allp = out + [sp.expand(qq * (z - z2) - 1), sp.expand(ww * z * z2 - 1)]
vars_ = ['ww', 'rho', 'sig', 'z', 'z2', 'qq', 'ee']
txt = lambda p_: str(p_).replace('**', '^')
open('h_red2_%s.sing' % TAG, 'w').write('ring R = 0,(%s),dp;\noption(redSB);\nideal I =\n%s;\nideal G = std(I);\nprint("dim:"); dim(G);\nprint("vdim:"); vdim(G);\nideal Je = eliminate(G, ww*rho*sig*z*z2*qq);\nprint("eliminant in ee:"); Je;\nprint("basis:"); G;\nquit;\n'
                              % (','.join(vars_), ',\n'.join(txt(p_) for p_ in allp)))
open('h_red2_%s.ms' % TAG, 'w').write(','.join(vars_) + '\n0\n' + ',\n'.join(txt(p_) for p_ in allp) + '\n')
for i, p_ in zip((0, 2, 3, 5, 8), out):
    print('--- eq', i, ':', str(sp.collect(p_, qq))[:700])
print('written (%.0fs)' % (time.time() - t0))
