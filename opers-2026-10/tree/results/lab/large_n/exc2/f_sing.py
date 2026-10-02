"""EXC2 T1 via Singular (the sympy Groebner run for {-2,1,3,4} timed out): one apparent singularity, exponent set from argv.
Usage: f_sing.py <t> <exps> <a0> <a1>.  Exit 0 = solutions with z != 0 exist; 3 = none (basis {1}); 4 = Singular did not finish."""
import hashlib, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC2.md', 'rb').read()).hexdigest() == open('SEAL_EXC2.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); exps = sorted(int(x) for x in sys.argv[2].split(',')); a0 = sp.Rational(sys.argv[3]); a1 = sp.Rational(sys.argv[4])
k = 2 * (t - 1) / (3 - t); b = k / (k + 1)
z, E, e = Fb.z, Fb.E, Fb.e
th = sp.Symbol('th')
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (2, 1, 0) for m in range(1, 5 - j)}
R = {j: {m: r[j, m] for m in range(1, 5 - j)} for j in (2, 1, 0)}
t0 = time.time()
P = sp.Poly(sp.expand((th**2 - a0**2) * (th**2 - a1**2)), th)
Pcoef = [P.coeff_monomial(th**j) for j in range(5)]
span = exps[-1] - exps[0]
lk = Fb.local_operator(Pcoef, [sp.Rational(1, 2), 1], b, -1, R, span - 4 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
target = z**4 * sp.prod([(e - x) for x in exps])
sol_ind = sp.solve([c for c in sp.Poly(sp.expand(lead - target), e).coeffs() if c != 0], [r[2, 2], r[1, 3], r[0, 4]], dict=True)
print('t = %s, exponents %s, xi-exponents (%s, %s); indicial: %s' % (t, exps, a0, a1, sol_ind), flush=True)
conds = [(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [r[2, 1], r[1, 1], r[1, 2], r[0, 1], r[0, 2], r[0, 3], z]
polys = []
for pos, mon, cf in eqs:
    num = sp.numer(sp.together(cf))
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(num)[1]:
        if fac != z:
            core *= fac**mult
    polys.append(sp.expand(core))
    kn = ' '.join('%s^%d' % (s, p) for s, p in zip(['E'] + [str(x) for x in kaps], mon) if p)
    print('   at exponent %d, [%s]: degree %d, %d terms%s' % (pos, kn or '1', sp.Poly(core, *unknowns).total_degree(), len(sp.Poly(core, *unknowns).terms()),
                                                             (': ' + str(sp.factor(core))) if len(str(core)) < 160 else ''), flush=True)
print('%d equations, %d unknowns  (%.0fs)' % (len(polys), len(unknowns), time.time() - t0), flush=True)
fn = 'f_sing_t%s_e%s.sing' % (str(t).replace('/', '_'), '_'.join(str(x) for x in exps))
open(fn, 'w').write('ring R = 0,(ww,r21,r11,r12,r01,r02,r03,z),dp;\noption(redSB);\nideal I =\n' + ',\n'.join(str(p_).replace('**', '^') for p_ in polys + [sp.Symbol('ww') * z - 1]) +
                    ';\nideal G = std(I);\nprint("dim:"); dim(G);\nprint("vdim:"); vdim(G);\nprint("basis size:"); size(G);\nif (size(G) < 30) { G; }\nquit;\n')
try:
    res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=1500).stdout
except subprocess.TimeoutExpired:
    print('Singular did not finish in 1500 s'); sys.exit(4)
open(fn.replace('.sing', '.out'), 'w').write(res)
print(res[:1500])
dim = int(re.search(r'dim:\n(-?\d+)', res).group(1))
print('RESULT for %s: %s  (%.0fs)' % (exps, 'NO solution with z != 0' if dim == -1 else 'solutions exist, dimension %d' % dim, time.time() - t0))
sys.exit(3 if dim == -1 else 0)
