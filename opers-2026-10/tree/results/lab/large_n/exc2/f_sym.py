"""EXC2: the one-point system with exponents {-1,1,2,4} for SYMBOLIC xi-exponents (A = a0^2, B = a1^2) at t = 2 and symbolic b:
equations printed and a Singular script with parameters written (to see the structure of the solution)."""
import sys
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
A, B, b = sp.symbols('A B b')
z, E, e = Fb.z, Fb.E, Fb.e
th = sp.Symbol('th')
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (2, 1, 0) for m in range(1, 5 - j)}
R = {j: {m: r[j, m] for m in range(1, 5 - j)} for j in (2, 1, 0)}
Pcoef = [A * B, 0, -(A + B), 0, 1]
exps = [-1, 1, 2, 4]
lk = Fb.local_operator(Pcoef, [sp.Rational(1, 2), 1], b, -1, R, 2)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
ind = {r[2, 2]: -4 * z, r[1, 3]: 8 * z**2, r[0, 4]: -8 * z**3}
eqs = Fb.split_conditions([(p_, sp.expand(c_.subs(ind))) for p_, c_ in conds], kaps)
unknowns = [r[2, 1], r[1, 1], r[1, 2], r[0, 1], r[0, 2], r[0, 3], z]
polys = []
for pos, mon, cf in eqs:
    kn = ' '.join('%s^%d' % (s_, p_) for s_, p_ in zip(['E'] + [str(x) for x in kaps], mon) if p_)
    print('--- resonance at %d, [%s]:' % (pos, kn or '1'))
    print('   ', cf)
    polys.append(sp.expand(sp.numer(sp.together(cf))))
w = sp.Symbol('ww')
lines = ['ring R = (0,A,B,b),(ww,r11,r12,r01,r02,r03,r21,z),lp;', 'option(redSB);', 'ideal I =',
         ',\n'.join(str(p_).replace('**', '^') for p_ in polys + [w * z - 1]) + ';', 'ideal G = std(I);', 'G;', 'quit;']
open('f_sym.sing', 'w').write('\n'.join(lines) + '\n')
