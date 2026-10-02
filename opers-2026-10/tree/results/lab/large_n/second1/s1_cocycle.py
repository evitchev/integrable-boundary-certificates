"""SECOND1 test: is Sol 3's Gamma-form recurrence (item 298) a multiplier image of the tensor product of two second-order equations with common
right side (x^a - E x^b)?  Cocycle condition F(nu) == 1 (see SEAL_SECOND1.md), solved exactly for the exponents (r1, q1, q2; r2 = 0 fixes the
r <-> q translation redundancy).  Usage: s1_cocycle.py k [assign] [control]; assign 1: (a, b) = (c, n); assign 2: (a, b) = (n, c)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
k = sp.Rational(sys.argv[1]); assign = int(sys.argv[2]) if len(sys.argv) > 2 else 1; ctrl = sys.argv[3] if len(sys.argv) > 3 else ''
nu, l0, l1 = sp.symbols('nu l0 l1')
r1, q1, q2 = sp.symbols('r1 q1 q2'); r2 = sp.Integer(0)
n = k + 4; sig = n / k; c = n + sig; s = (n - 1) / 2
kG = k + sp.Rational(1, 10) if ctrl == 'wrongshift' else k
U = r1 + r2 + q1 + q2
def p0(D):
    if ctrl == 'wrongrot':     # impose the UNROTATED exponent set {s -+ l0, s -+ l1} as products r_i + q_j cannot realise: force e_ij = (s+l0, s-l0, s+l1, s-l1) with r1-r2 = 2 l0
        return sp.Mul(*[(D - e) for e in (r1 + q1, r1 + q2, r2 + q1, r2 + q2)]).subs({r1: 2 * l0, q1: s - l0, q2: s - l1 + 0 * l1})
    return sp.Mul(*[(D - e) for e in (r1 + q1, r1 + q2, r2 + q1, r2 + q2)])
a, b = (c, n) if assign == 1 else (n, c)
pa = lambda D: -(2 * D + a - U) * (2 * D + 2 * a - U)
pb = lambda D: (2 * D + b - U) * (2 * D + 2 * b - U)          # E stripped (cancels in the cocycle)
P4 = lambda T: (T**2 - l0**2) * (T**2 - l1**2)
def Gratio(T):  # G_k(T + sigma)/G_k(T), sigma = n/k
    sg = n / kG
    return (T / sg + (kG + 1) / 2) / (T / sg + (1 - kG) / 2)
# f_a(nu) = m(nu)/m(nu - a) = -pa(nu - a) H(s - nu)/p0(nu);  f_b(nu) = pb(nu - b) H(s - nu)/p0(nu);  H = P4 G.
# cocycle F = f_a(nu) f_b(nu - a) / (f_b(nu) f_a(nu - b)) == 1; the H factors reduce to H(s - nu + a)/H(s - nu + b).
Hrat = sp.Integer(1)
# H(s - nu + a)/H(s - nu + b): a - b = +-sigma
if a - b == sig:  Hrat = P4(s - nu + a) / P4(s - nu + b) * Gratio(s - nu + b)
else:             Hrat = P4(s - nu + a) / P4(s - nu + b) / Gratio(s - nu + a)
F = (-pa(nu - a)) * pb(nu - a - b) / (pb(nu - b) * (-pa(nu - b - a))) * p0(nu - b) / p0(nu - a) * Hrat
num = sp.numer(sp.together(sp.expand(F) - 1)) if False else sp.numer(sp.together(F - 1))
P = sp.Poly(sp.expand(num), nu)
eqs = [sp.expand(cf) for cf in P.coeffs()]
print(f'k = {k}, assign {assign}, control "{ctrl}": n = {n}, sigma = {sig}, c = {c}, s = {s}; cocycle numerator degree {P.degree()} in nu, {len(eqs)} coefficient equations')
unk = [r1, q1, q2] if ctrl != 'wrongrot' else []
if unk:
    sols = sp.solve(eqs, unk, dict=True)
else:
    sols = [{}] if all(sp.simplify(e) == 0 for e in eqs) else []
print('solutions (exponents may depend on l0, l1):', sols)
good = []
for so in sols:
    ok = sp.simplify(F.subs(so) - 1) == 0
    print('  check F == 1 identically:', ok, ' exponents e_ij =', [sp.simplify(e.subs(so)) for e in (r1 + q1, r1 + q2, r2 + q1, r2 + q2)])
    if ok: good.append(so)
print('RESULT:', 'CLOSES' if good else 'NONE')
