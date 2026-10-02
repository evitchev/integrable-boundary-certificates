"""EXC2 T1: existence of level-1 apparent singularities for the Sol-3 three-term operator in the sealed class
    L = (th^2 - PX^2)(th^2 - pi^2) + xi (th + 1/2) - E xi^b + xi [R_2 th^2 + R_1 th + R_0],  R_j = sum_{m=1}^{4-j} r_jm/(xi - z)^m.
Usage: f_sol3.py <t> <exps e.g. -1,1,2,4> [PX pi | sym] [tamper]
Prints the split no-log conditions (exponent pair, power of E) and the Groebner analysis of the system.
Exit 0 = solutions exist (z != 0); 3 = none."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC2.md', 'rb').read()).hexdigest() == open('SEAL_EXC2.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); exps = sorted(int(x) for x in sys.argv[2].split(','))
args = sys.argv[3:]
tamper = 'tamper' in args
args = [a for a in args if a != 'tamper']
k = 2 * (t - 1) / (3 - t); b = k / (k + 1)
PXs, PIs = sp.symbols('PX PI')
if args and args[0] != 'sym':
    PXv, PIv = sp.Rational(args[0]), sp.Rational(args[1])
else:
    PXv, PIv = PXs, PIs
z, E, e = Fb.z, Fb.E, Fb.e
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (2, 1, 0) for m in range(1, 5 - j)}
R = {j: {m: r[j, m] for m in range(1, 5 - j)} for j in (2, 1, 0)}
t0 = time.time()
P = sp.Poly(sp.expand((sp.Symbol('th')**2 - PXv**2) * (sp.Symbol('th')**2 - PIv**2)), sp.Symbol('th'))
Pcoef = [P.coeff_monomial(sp.Symbol('th')**j) for j in range(5)]
Lcoef = [sp.Rational(1, 2) if not tamper else sp.Integer(1), 1]
span = exps[-1] - exps[0]
lk = Fb.local_operator(Pcoef, Lcoef, b, -1, R, span - 4 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
target = z**4 * sp.prod([(e - x) for x in exps])
ind_eqs = [c for c in sp.Poly(sp.expand(lead - target), e).coeffs() if c != 0]
sol_ind = sp.solve(ind_eqs, [r[2, 2], r[1, 3], r[0, 4]], dict=True)
print('t = %s (b = %s), exponents %s, momenta (%s, %s)%s' % (t, b, exps, PXv, PIv, ' TAMPER' if tamper else ''))
print('indicial conditions fix:', sol_ind, flush=True)
if not sol_ind:
    print('indicial conditions inconsistent'); sys.exit(3)
conds = [(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [r[2, 1], r[1, 1], r[1, 2], r[0, 1], r[0, 2], r[0, 3], z]
print('no-log conditions: %d polynomial equations in %d unknowns  (%.0fs)' % (len(eqs), len(unknowns), time.time() - t0))
for pos, mon, cf in eqs:
    kn = ' '.join('%s^%d' % (s, p) for s, p in zip(['E'] + [str(x) for x in kaps], mon) if p)
    print('   at exponent %d, monomial [%s]: degree %s, %d terms' % (pos, kn or '1', sp.Poly(cf, *unknowns).total_degree(), len(sp.Poly(cf, *unknowns).terms())))
polys = [sp.numer(sp.together(cf)) for _, _, cf in eqs]
# saturate z != 0
w = sp.Symbol('w')
G = sp.groebner(polys + [1 - w * z], w, *unknowns, order='grevlex')
print('Groebner basis computed: %d elements  (%.0fs)' % (len(G.exprs), time.time() - t0), flush=True)
if G.exprs == [1]:
    print('RESULT: NO solution with z != 0 for exponents %s' % exps)
    # which subsets are already inconsistent?  E-graded diagnosis: keep only E^0 equations
    for label, sel in (('E^0 equations only', [p_ for (pos, mon, cf), p_ in zip(eqs, polys) if mon[0] == 0]),
                       ('equations at the first resonance only', [p_ for (pos, mon, cf), p_ in zip(eqs, polys) if pos == exps[1]]),
                       ('first two resonances', [p_ for (pos, mon, cf), p_ in zip(eqs, polys) if pos in exps[1:3]])):
        G2 = sp.groebner(sel + [1 - w * z], w, *unknowns, order='grevlex')
        print('   subsystem [%s] (%d equations): %s' % (label, len(sel), 'inconsistent' if G2.exprs == [1] else 'consistent'))
    sys.exit(3)
try:
    dim0 = G.is_zero_dimensional
except Exception:
    dim0 = None
print('RESULT: solutions exist; zero-dimensional: %s' % dim0)
json.dump({'t': str(t), 'exps': exps, 'PX': str(PXv), 'PI': str(PIv), 'basis': [str(g) for g in G.exprs], 'indicial': {str(k_): str(v) for k_, v in sol_ind[0].items()}},
          open('f_sol3_t%s_e%s%s.json' % (str(t).replace('/', '_'), '_'.join(str(x) for x in exps), '_tamper' if tamper else ''), 'w'), indent=1)
for g in G.exprs[:12]:
    print('   ', str(g)[:240])
sys.exit(0)
