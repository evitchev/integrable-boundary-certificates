"""Kit verifier (v1.2.5, new; pure Python/sympy, no Singular) of EXC4 A.2: one apparent singularity in the SIXTH-order polynomial form of Solution 3
at t = 2, momenta (a0, a1) = (1/2, 1/3), exponents S_6 = {-1,1,2,3,4,6}.  Run with cwd = results/lab/large_n/exc4.
Inputs (archived): the system f_six_t2_1-2_1-3.sing (written by f_six_sol3.py; checked byte-identical by kit check exc4_regen_six), Singular's reduced
Groebner basis in f_six_t2_1-2_1-3.out, and the I_1 / I_3 relations in f_six_charges_1-2_1-3.sing (written by f_six_charges.py from the Gamma-form
engine; checked byte-identical by exc4_regen_six_charges).
Checks: (1) every input polynomial reduces to 0 modulo Singular's basis G (V(G) is contained in V(I)), and G satisfies Buchberger's criterion
(all S-polynomials reduce to 0), so G is a Groebner basis; (2) G has 2 standard monomials (vdim 2) and its 2 solutions are found exactly; (3) on both, r41 = -2/3 = -b
(the sealed residue law) and the monic I_1 is vacuum + 1 (vacuum = I_1 at zero residues): HALF a level; (4) I_3 = 41/64 on both, which is not a
root of the level-one data quadratic (d1_sol3_t2.json) at this point.
Not checked: that G generates all of I (no solution lost) -- Singular's std; (1) shows only containment V(G) in V(I)."""
import json, re, sys, time
import sympy as sp
T0 = time.time(); ok = True
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
sing = open('f_six_t2_1-2_1-3.sing').read()
VARS = re.search(r'ring R = 0,\(([^)]*)\),dp;', sing).group(1).split(',')
S = sp.symbols(VARS); L = dict(zip(VARS, S))
def parse(s): return sp.sympify(s.replace('^', '**'), locals=L)
i = sing.index('ideal I =') + 9; I = [parse(p) for p in re.split(r',\s*\n', sing[i:sing.index(';', i)]) if p.strip()]
if mode == 'tamper': I[-2] = I[-2] + 1          # negative control: perturb one no-log equation
out = open('f_six_t2_1-2_1-3.out').read()
Gs = [parse(m.group(1)) for m in re.finditer(r'G\[\d+\]=(.*)', out)]
print('input: %d polynomials in %d unknowns; Singular basis: %d elements' % (len(I), len(VARS), len(Gs)))
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from sympy.polys.orderings import grevlex
Rg, *Xg = ring(','.join(VARS), QQ, grevlex)
Ir = [Rg(p) for p in I]; Gr = [Rg(g) for g in Gs]
nz = sum(1 for f in Ir if f.rem(Gr) != 0)
print('input polynomials reducing to 0 modulo the Singular basis: %d of %d (every solution of the basis solves the system)' % (len(Ir) - nz, len(Ir))); ok &= nz == 0
LM = [g.LM for g in Gr]; nv = len(VARS)
std = []; fr = [tuple([0] * nv)]; seen = set(fr)
while fr:
    m = fr.pop()
    if any(all(a >= c for a, c in zip(m, l)) for l in LM): continue
    std.append(m)
    for j in range(nv):
        m2 = tuple(m[q] + (j == q) for q in range(nv))
        if m2 not in seen and sum(m2) < 40: seen.add(m2); fr.append(m2)
print('standard monomials of the basis (vdim): %d' % len(std)); ok &= len(std) == 2
# Buchberger criterion: every S-polynomial of the basis reduces to 0 (so it is a Groebner basis and vdim is exact)
def spoly(f, g):
    lcm = tuple(max(a, c) for a, c in zip(f.LM, g.LM))
    mf = Rg({tuple(a - c for a, c in zip(lcm, f.LM)): QQ(1)}); mg = Rg({tuple(a - c for a, c in zip(lcm, g.LM)): QQ(1)})
    return mf * f * (1 / f.LC) - mg * g * (1 / g.LC)
bad = sum(1 for a in range(len(Gr)) for b in range(a + 1, len(Gr)) if spoly(Gr[a], Gr[b]).rem(Gr) != 0)
print('Buchberger criterion on the Singular basis: %d of %d S-polynomials fail to reduce to 0' % (bad, len(Gr) * (len(Gr) - 1) // 2)); ok &= bad == 0
sols = sorted(sp.solve(Gs, S, dict=True), key=lambda s_: str(sp.nsimplify(s_[L['r32']])))
print('solutions: %d (vdim 2 expected)' % len(sols)); ok &= len(sols) == 2
ch = open('f_six_charges_1-2_1-3.sing').read(); i = ch.index('ideal I =') + 9
rel = [p for p in re.split(r',\s*\n', ch[i:ch.index(';', i)]) if p.strip()][-2:]
ii, ee = sp.symbols('ii ee'); L.update({'ii': ii, 'ee': ee})
rI1, rI3 = parse(rel[0]), parse(rel[1])
I1 = lambda sub: sp.solve(rI1.subs(sub), ii)[0]; I3 = lambda sub: sp.solve(rI3.subs(sub), ee)[0]
vac = I1({v: 0 for v in S}); print('relations from the Gamma-form engine: %s = 0, %s = 0; vacuum monic I_1 = %s' % (rel[0], rel[1], vac))
r41 = L['r41']
for s_ in sols:
    res = [sp.simplify(sp.expand(p.subs(s_))) for p in I]
    d1 = sp.nsimplify(sp.simplify(I1(s_) - vac)); e3 = sp.simplify(I3(s_))
    good = all(x == 0 for x in res) and s_[r41] == sp.Rational(-2, 3) and d1 == 1
    ok &= good
    print('  solution z = %s, r32 = %s: all %d input equations vanish exactly: %s; r41 = %s (sealed -b = -2/3); Delta I_1 = %s (half level = 1); I_3 = %s'
          % (s_[L['z']], sp.nsimplify(s_[L['r32']]), len(I), all(x == 0 for x in res), s_[r41], d1, e3))
    ok &= e3 == sp.Rational(41, 64)
dd = json.load(open('d1_sol3_t2.json'))['blocks']['I3']
P, Q, e = sp.symbols('P Q e'); b = sp.Rational(2, 3); A, B = sp.Rational(1, 4), sp.Rational(1, 9); rho2 = sp.Rational(2, 3)
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][a][c], locals={'P': P, 'Q': Q}) for a in range(2) for c in range(2)])
vacd = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.Poly(sp.expand((Mx / vacd.coeff_monomial(P**4) - e * sp.eye(2)).det()), P, Q, e)
cpv = sp.Poly(sum(cf * (A / b)**(m[0] // 2) * (B / b - rho2)**(m[1] // 2) * e**m[2] for m, cf in cp.terms()), e).monic()
notroot = cpv.eval(sp.Rational(41, 64)) != 0
print('level-one data quadratic at this point: %s; 41/64 is a root: %s' % (cpv.as_expr(), not notroot)); ok &= notroot
print('SIX VERIFY [%s]: %s' % (mode, 'PASS' if ok else 'FAIL')); sys.exit(0 if ok else 1)
