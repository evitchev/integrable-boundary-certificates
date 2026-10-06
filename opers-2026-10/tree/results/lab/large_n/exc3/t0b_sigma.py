"""EXC3 T0(b): the symmetry sigma: L -> -(L^+)(-xi) on the Sol 2 normal form, by machine (symbolic a0, a1, b, E).
Operators are represented as polynomials in vartheta with coefficients that are finite sums of xi^r (r rational), acting on
the left: xi^r vartheta^j.  Adjoint: (xi^r vartheta^j)^+ = (-vartheta)^j xi^r = xi^r (-(vartheta + r))^j.  Then xi -> -xi: xi^r -> (-1)^r xi^r
(for r = b we keep the factor (-1)^b as a symbol s_b and allow E -> E' as the seal says).
Also T0(c): frob.py regression on the EXC2 one-point Sol 3 system (t = 2, exponents {-1,1,2,4}, (a0, a1) = (1/2, 1/3)):
4 solutions, r21 = -b.  Tamper for (b): middle term xi((vartheta+1)^2 - a0^2) must NOT be sigma-invariant.
Exit 0 iff sigma(L_0) == L_0 (with E' = -s_b E), the tamper is not invariant, and the regression holds."""
import hashlib, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
th, xi, a0, a1, b, E, sb = sp.symbols('vartheta xi a0 a1 b E s_b')
# operator = dict {r: polynomial in th}  meaning sum_r xi^r * poly_r(th)
def adjoint_flip(op):
    out = {}
    for r, poly in op.items():
        padj = sp.expand(poly.subs(th, -(th + r)))
        sign = sb if r == b else (-1) ** int(r)
        out[r] = sp.expand(-sign * padj)          # the overall minus of sigma
    return out
def same(op1, op2):
    keys = set(op1) | set(op2)
    return all(sp.expand(op1.get(k, 0) - op2.get(k, 0)) == 0 for k in keys)
L0 = {0: th * (th**2 - a1**2), 1: (th + sp.Rational(1, 2))**2 - a0**2, b: -E}
S = adjoint_flip(L0)
# E' = -s_b E  <=>  the xi^b coefficient -E maps to -s_b*(-E) = s_b E ; L_0 with E -> -s_b E has coefficient s_b E
ok_sigma = same(S, {0: L0[0], 1: L0[1], b: sb * E})
print('sigma(L_0) = ', {k: sp.factor(v) for k, v in S.items()})
print('T0(b) sigma(L_0) == L_0 with E -> -(-1)^b E:', ok_sigma)
L0t = {0: th * (th**2 - a1**2), 1: (th + 1)**2 - a0**2, b: -E}
St = adjoint_flip(L0t)
ok_tamper = not same(St, {0: L0t[0], 1: L0t[1], b: sb * E})
print('   ordering tamper xi((vartheta+1)^2 - a0^2): sigma-invariant?', not ok_tamper, '(must be False)')
# a deformation term and its image, for the record (one point at z, R_1 = r/(xi - z)): expand in powers is not finite; check
# instead on the finite Laurent pieces xi^(-m) that the rule V_1 -> -V~_1 + 2 vartheta V~_2 holds with V_2 = 0:
D = {-2: sp.Symbol('r') * th}                     # xi^(-2) r vartheta  (a piece of xi R_1 vartheta)
SD = adjoint_flip(D)
print('   sigma on xi^-2 r vartheta:', SD, ' (expected: -(xi^-2) r (-(vartheta - 2)) * (-1)^(-2) * (-1) = xi^-2 r (vartheta - 2))')
# T0(c): EXC2 regression with frob.py
import frob as Fb
z, e, Ez = Fb.z, Fb.e, Fb.E
A, B = sp.Rational(1, 4), sp.Rational(1, 9)   # a0^2, a1^2 at (1/2, 1/3)
bv = sp.Rational(2, 3)
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (2, 1, 0) for m in range(1, 5 - j)}
R = {j: {m: r[j, m] for m in range(1, 5 - j)} for j in (2, 1, 0)}
Pcoef = [A * B, 0, -(A + B), 0, 1]
lk = Fb.local_operator(Pcoef, [sp.Rational(1, 2), 1], bv, -1, R, 2)
lead, conds, kaps = Fb.nolog_conditions(lk, [-1, 1, 2, 4], 4)
ind = sp.solve([c for c in sp.Poly(sp.expand(lead - z**4 * (e + 1) * (e - 1) * (e - 2) * (e - 4)), e).coeffs() if c != 0], [r[2, 2], r[1, 3], r[0, 4]], dict=True)[0]
eqs = Fb.split_conditions([(p_, sp.expand(c_.subs(ind))) for p_, c_ in conds], kaps)
polys = [sp.numer(sp.together(cf)) for _, _, cf in eqs]
unknowns = [r[2, 1], r[1, 1], r[1, 2], r[0, 1], r[0, 2], r[0, 3], z]
w = sp.Symbol('w')
G = sp.groebner(polys + [1 - w * z], w, *unknowns, order='lex')
elz = [g for g in G.exprs if g.free_symbols <= {z}]
r21v = [g for g in G.exprs if g.free_symbols <= {r[2, 1]}]
nsol = sp.Poly(elz[0], z).degree() if elz else None
print('T0(c) EXC2 one-point regression: eliminant in z degree %s (expect 4): %s ; r21 relation: %s (expect 3 r21 + 2)' % (nsol, sp.factor(elz[0]) if elz else None, r21v))
ok_reg = nsol == 4 and any(sp.expand(g - (3 * r[2, 1] + 2)) == 0 or sp.expand(g + (3 * r[2, 1] + 2)) == 0 for g in r21v)
print('T0(b)+(c):', 'PASS' if ok_sigma and ok_tamper and ok_reg else 'FAIL')
sys.exit(0 if ok_sigma and ok_tamper and ok_reg else 2)
