"""EXC3 T1 via Singular (the sympy Groebner run for {-2,1,4} did not finish): one apparent singularity, exponent set from argv.
Usage: f_sing_sol2.py <t> <exps> <a0> <a1>.  Exit 0 = solutions with z != 0 exist; 3 = none; 4 = Singular did not finish."""
import hashlib, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); exps = sorted(int(x) for x in sys.argv[2].split(',')); a0 = sp.Rational(sys.argv[3]); a1 = sp.Rational(sys.argv[4])
k = 2 * (t - 1) / (3 - t); b = k / (k + 1)
z, E, e = Fb.z, Fb.E, Fb.e
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (1, 0) for m in range(1, 4 - j)}
R = {j: {m: r[j, m] for m in range(1, 4 - j)} for j in (1, 0)}
t0 = time.time()
Pcoef = [0, -a1**2, 0, 1]; Lcoef = [sp.Rational(1, 4) - a0**2, 1, 1]
span = exps[-1] - exps[0]
lk = Fb.local_operator(Pcoef, Lcoef, b, -1, R, span - 3 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 3)
sol_ind = sp.solve([c for c in sp.Poly(sp.expand(lead - z**3 * sp.prod([(e - x) for x in exps])), e).coeffs() if c != 0], [r[1, 2], r[0, 3]], dict=True)
print('t = %s, exponents %s, (a0, a1) = (%s, %s); indicial: %s' % (t, exps, a0, a1, sol_ind), flush=True)
conds = [(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [r[1, 1], r[0, 1], r[0, 2], z]
polys = []
for pos, mon, cf in eqs:
    num = sp.numer(sp.together(cf)); core = sp.Integer(1)
    for fac, mult in sp.factor_list(num)[1]:
        if fac != z: core *= fac**mult
    polys.append(sp.expand(core))
    kn = ' '.join('%s^%d' % (s, p) for s, p in zip(['E'] + [str(x) for x in kaps], mon) if p)
    print('   at exponent %d, [%s]: degree %d, %d terms' % (pos, kn or '1', sp.Poly(core, *unknowns).total_degree(), len(sp.Poly(core, *unknowns).terms())), flush=True)
fn = 'f_sing_sol2_t%s_e%s.sing' % (str(t).replace('/', '_'), '_'.join(str(x) for x in exps))
open(fn, 'w').write('ring R = 0,(ww,r11,r01,r02,z),dp;\noption(redSB);\nideal I =\n' + ',\n'.join(str(p_).replace('**', '^') for p_ in polys + [sp.Symbol('ww') * z - 1]) +
                    ';\nideal G = std(I);\nprint("dim:"); dim(G);\nprint("vdim:"); vdim(G);\nif (size(G) < 30) { G; }\nquit;\n')
try:
    res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=1500).stdout
except subprocess.TimeoutExpired:
    print('Singular did not finish in 1500 s'); sys.exit(4)
open(fn.replace('.sing', '.out'), 'w').write(res); print(res[:1200])
dim = int(re.search(r'dim:\n(-?\d+)', res).group(1))
print('RESULT for %s: %s  (%.0fs)' % (exps, 'NO solution with z != 0' if dim == -1 else 'solutions exist, dimension %d' % dim, time.time() - t0))
sys.exit(3 if dim == -1 else 0)
