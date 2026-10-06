"""EXC3 T1: ONE apparent singularity for the Sol 2 normal form (seal sec. 3):
   L = vartheta(vartheta^2 - a1^2) + xi((vartheta+1/2)^2 - a0^2) - E xi^b + xi [R_1 vartheta + R_0],  R_1 = r11/(xi-z) + r12/(xi-z)^2,
   R_0 = r01/(xi-z) + r02/(xi-z)^2 + r03/(xi-z)^3.
Usage: f_sol2.py <t> <exps e.g. -1,1,3> <a0> <a1> [tamper]   (tamper: middle term xi((vartheta+1)^2 - a0^2)).
Prints the indicial solution, the split no-log conditions, the Groebner analysis; writes f_sol2_t<t>_e<exps>[_tamper].json.
Exit 0 = solutions with z != 0 exist; 3 = none."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); exps = sorted(int(x) for x in sys.argv[2].split(',')); a0 = sp.Rational(sys.argv[3]); a1 = sp.Rational(sys.argv[4])
tamper = len(sys.argv) > 5
k = 2 * (t - 1) / (3 - t); b = k / (k + 1)
z, E, e = Fb.z, Fb.E, Fb.e
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (1, 0) for m in range(1, 4 - j)}
R = {j: {m: r[j, m] for m in range(1, 4 - j)} for j in (1, 0)}
t0 = time.time()
Pcoef = [0, -a1**2, 0, 1]                                      # vartheta^3 - a1^2 vartheta
shift = 1 if tamper else sp.Rational(1, 2)
Lcoef = [shift**2 - a0**2, 2 * shift, 1]                        # (vartheta + shift)^2 - a0^2
span = exps[-1] - exps[0]
lk = Fb.local_operator(Pcoef, Lcoef, b, -1, R, span - 3 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 3)
target = z**3 * sp.prod([(e - x) for x in exps])
ind_eqs = [c for c in sp.Poly(sp.expand(lead - target), e).coeffs() if c != 0]
sol_ind = sp.solve(ind_eqs, [r[1, 2], r[0, 3]], dict=True)
print('t = %s (b = %s), exponents %s, (a0, a1) = (%s, %s)%s' % (t, b, exps, a0, a1, ' TAMPER' if tamper else ''))
print('indicial polynomial / z^3:', sp.factor(lead / z**3), ' -> indicial conditions fix:', sol_ind, flush=True)
if not sol_ind:
    print('indicial conditions inconsistent (exponent sum must be 3)'); sys.exit(3)
conds = [(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [r[1, 1], r[0, 1], r[0, 2], z]
print('no-log conditions: %d polynomial equations in %d unknowns  (%.0fs)' % (len(eqs), len(unknowns), time.time() - t0))
polys = []
for pos, mon, cf in eqs:
    kn = ' '.join('%s^%d' % (s, p) for s, p in zip(['E'] + [str(x) for x in kaps], mon) if p)
    num = sp.numer(sp.together(cf))
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(num)[1]:
        if fac != z:
            core *= fac**mult
    polys.append(sp.expand(core))
    print('   at exponent %d, [%s]: degree %d, %d terms%s' % (pos, kn or '1', sp.Poly(core, *unknowns).total_degree(), len(sp.Poly(core, *unknowns).terms()), (': ' + str(sp.factor(core))) if len(str(core)) < 200 else ''))
w = sp.Symbol('w')
G = sp.groebner(polys + [1 - w * z], w, *unknowns, order='lex')
print('Groebner basis (lex): %d elements  (%.0fs)' % (len(G.exprs), time.time() - t0), flush=True)
tag = 'f_sol2_t%s_e%s%s' % (str(t).replace('/', '_'), '_'.join(str(x) for x in exps), '_tamper' if tamper else '')
if G.exprs == [1]:
    print('RESULT: NO solution with z != 0 for exponents %s' % exps)
    for label, sel in (('E^0 equations only', [p_ for (pos, mon, cf), p_ in zip(eqs, polys) if mon[0] == 0]),
                       ('first resonance only', [p_ for (pos, mon, cf), p_ in zip(eqs, polys) if pos == exps[1]])):
        G2 = sp.groebner(sel + [1 - w * z], w, *unknowns, order='lex')
        print('   subsystem [%s] (%d equations): %s' % (label, len(sel), 'inconsistent' if G2.exprs == [1] else 'consistent'))
    json.dump({'t': str(t), 'exps': exps, 'a0': str(a0), 'a1': str(a1), 'tamper': tamper, 'solutions': 0}, open(tag + '.json', 'w'), indent=1)
    sys.exit(3)
elz = [g for g in G.exprs if g.free_symbols <= {z}]
print('RESULT: solutions exist; zero-dimensional: %s; eliminant in z: %s' % (G.is_zero_dimensional, sp.factor(elz[0]) if elz else None))
for g in G.exprs:
    print('   ', str(g)[:200])
json.dump({'t': str(t), 'exps': exps, 'a0': str(a0), 'a1': str(a1), 'tamper': tamper, 'basis': [str(g) for g in G.exprs],
           'indicial': {str(k_): str(v) for k_, v in sol_ind[0].items()}, 'Lcoef': [str(x) for x in Lcoef]}, open(tag + '.json', 'w'), indent=1)
sys.exit(0)
