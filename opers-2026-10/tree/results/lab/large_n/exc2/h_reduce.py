"""EXC2 addendum: reduction of the two-point system (h_two_v2_t2.sing) by its linear equations.
Per point: r03 = -z (r12 + 4 z), r21 = -b (these two are the same as for one point: no background enters), then r02 and
r01 from the equations linear in them.  Remaining: 3 equations per point in (r11, r12, z | s11, s12, z2) and q = 1/(z - z2).
Writes h_red_t2.sing / .ms and the expression of e (monic I_3) on the reduced variables."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
t0 = time.time()
TAG = sys.argv[1]; b = sp.Rational(sys.argv[2])
s = open('h_two_v2_%s.sing' % TAG).read()
body = s[s.index('ideal I =') + len('ideal I ='):s.index(';\nideal G')]
names = 'ww,r21,r11,r12,r01,r02,r03,z,s21,s11,s12,s01,s02,s03,z2,qq,ii,ee'.split(',')
sy = {n: sp.Symbol(n) for n in names}
P = [sp.sympify(p.strip().replace('^', '**'), locals=sy) for p in body.split(',\n')]
sub = {}
# point 1: eqs 0..6 ; point 2: eqs 7..13
for (off, R03, R12, R21, R02, R01, Z) in ((0, 'r03', 'r12', 'r21', 'r02', 'r01', 'z'), (7, 's03', 's12', 's21', 's02', 's01', 'z2')):
    e1 = P[off + 1]; e3 = P[off + 3]
    s03 = sp.solve(e1, sy[R03])[0]
    s21 = sp.solve(e3.subs(sy[R03], s03), sy[R21])
    s21 = [v for v in s21][0]
    assert sp.simplify(s21 + b) == 0, s21
    sub[sy[R03]] = sp.expand(s03); sub[sy[R21]] = -b
for (off, R02, R01) in ((0, 'r02', 'r01'), (7, 's02', 's01')):
    e0 = sp.expand(P[off + 0].subs(sub))
    sol02 = sp.solve(e0, sy[R02])
    assert len(sol02) == 1
    sub[sy[R02]] = sp.together(sol02[0])
for (off, R02, R01) in ((0, 'r02', 'r01'), (7, 's02', 's01')):
    e2 = sp.together(P[off + 2].subs(sub))
    sol01 = sp.solve(sp.numer(e2), sy[R01])
    assert len(sol01) == 1
    sub[sy[R01]] = sp.together(sol01[0])
print('linear eliminations done (%.0fs)' % (time.time() - t0), flush=True)
red = []
keep = [4, 5, 6, 11, 12, 13]
unk = [sy[n] for n in ('r11', 'r12', 'z', 's11', 's12', 'z2', 'qq')]
for i in keep:
    ex = sp.together(P[i].subs(sub).subs(sub))
    num = sp.numer(ex)
    # strip monomial factors z^a z2^c
    f = sp.factor_list(sp.expand(num))
    core = sp.Integer(1)
    for fac, mult in f[1]:
        if fac in (sy['z'], sy['z2']):
            continue
        core *= fac ** mult
    core = sp.expand(core)
    red.append(core)
    print('eq %2d -> degree %d, %d terms  (%.0fs)' % (i, sp.Poly(core, *unk).total_degree(), len(sp.Poly(core, *unk).terms()), time.time() - t0), flush=True)
# e and its relation
eexpr = sp.solve(P[16], sy['ee'])[0]
eexpr = sp.together(eexpr.subs(sub).subs(sub))
enum, eden = sp.fraction(eexpr)
erel = sp.expand(sy['ee'] * eden - enum)
f = sp.factor_list(erel)
print('e-relation: degree %d' % sp.Poly(erel, *(unk + [sy['ee']])).total_degree())
allp = red + [sp.expand(sy['qq'] * (sy['z'] - sy['z2']) - 1), sp.expand(sy['ww'] * sy['z'] * sy['z2'] - 1), erel]
vars_ = ['ww', 'r11', 'r12', 's11', 's12', 'z', 'z2', 'qq', 'ee']
txt = lambda p_: str(p_).replace('**', '^')
open('h_red_%s.sing' % TAG, 'w').write('ring R = 0,(%s),dp;\noption(redSB);\nideal I =\n%s;\nideal G = std(I);\nprint("dim:"); dim(G);\nprint("vdim:"); vdim(G);\nideal Je = eliminate(G, ww*r11*r12*s11*s12*z*z2*qq);\nprint("eliminant in ee:"); Je;\nquit;\n'
                             % (','.join(vars_), ',\n'.join(txt(p_) for p_ in allp)))
open('h_red_%s.ms' % TAG, 'w').write(','.join(vars_) + '\n0\n' + ',\n'.join(txt(p_) for p_ in allp) + '\n')
json.dump({'sub': {str(k_): str(v) for k_, v in sub.items()}, 'e': str(eexpr)}, open('h_red_%s.json' % TAG, 'w'), indent=1)
print('written h_red_%s.sing / .ms  (%.0fs)' % (TAG, time.time() - t0))
