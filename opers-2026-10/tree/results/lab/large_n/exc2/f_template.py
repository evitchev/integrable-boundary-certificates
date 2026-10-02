"""EXC2 T0 (control, known answer): the local analysis code on the Sol-1 template class
    L = (vartheta^2 - lam^2) + xi (2 vartheta + 1 - 2 x0) + E xi^b + xi [ r01/(xi - z) + r02/(xi - z)^2 ],  exponents {-1, 2}.
Expected (EXC1): r02 = -2 z, r01 = -b, and (2 - b) z^2 + 2 x0 (1 - b) z - b lam^2 + b (1 - b)^2/4 = 0.   Exit 0 iff so."""
import hashlib, sys
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC2.md', 'rb').read()).hexdigest() == open('SEAL_EXC2.sha256').read().split()[0]
lam, x0, b, r01, r02 = sp.symbols('lam x0 b r01 r02')
z, E, e = Fb.z, Fb.E, Fb.e
lk = Fb.local_operator([-lam**2, 0, 1], [1 - 2 * x0, 2], b, +1, {0: {1: r01, 2: r02}}, 3)
lead, conds, kaps = Fb.nolog_conditions(lk, [-1, 2], 2)
ind = sp.Poly(sp.expand(lead), e)
print('indicial polynomial:', sp.factor(lead))
sol_ind = sp.solve(sp.Poly(sp.expand(lead - z**2 * (e + 1) * (e - 2)), e).coeffs(), [r02], dict=True)
print('exponents {-1, 2}  =>', sol_ind)
eqs = Fb.split_conditions([(p_, sp.expand(c_.subs(sol_ind[0]))) for p_, c_ in conds], kaps)
for pos, mon, cf in eqs:
    print('   no-log condition at exponent %d, monomial E^%d: %s = 0' % (pos, mon[0], cf))
polys = [cf for _, _, cf in eqs]
sol = sp.solve([cf for _, mon, cf in eqs if mon[0] == 1], [r01], dict=True)   # v1: r01 from the E^1 condition (v0 asked solve() for the overdetermined pair and got [])
print('solution for r01:', sol)
ok = sol_ind[0][r02] == -2 * z and any(sp.simplify(s_[r01] + b) == 0 for s_ in sol)
rest = [sp.factor(sp.expand(p_.subs(r01, -b))) for p_ in polys]
target = (2 - b) * z**2 + 2 * x0 * (1 - b) * z - b * lam**2 + b * (1 - b)**2 / 4
quad = [r for r in rest if r != 0]
prop = all(sp.simplify(sp.cancel(r / target).diff(z)) == 0 and sp.cancel(r / target).free_symbols <= {b, z} for r in quad) and len(quad) >= 1
print('remaining condition(s) with r01 = -b:', quad)
print('proportional to (2-b) z^2 + 2 x0 (1-b) z - b lam^2 + b(1-b)^2/4 (EXC1 quadratic):', prop)
print('T0', 'PASS' if (ok and prop) else 'FAIL')
sys.exit(0 if (ok and prop) else 2)
