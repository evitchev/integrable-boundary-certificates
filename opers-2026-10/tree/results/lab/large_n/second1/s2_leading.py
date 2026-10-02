"""SECOND1: leading nu-coefficient of the cocycle numerator, k symbolic, exponents symbolic -- a NECESSARY condition valid for ALL k and ALL exponents.
Plus the next coefficient at the surviving k values."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
kk, nu, l0, l1, r1, q1, q2 = sp.symbols('k nu l0 l1 r1 q1 q2')
r2 = sp.Integer(0)
n = kk + 4; sig = n / kk; c = n + sig; s = (n - 1) / 2; U = r1 + r2 + q1 + q2
p0 = lambda D: sp.Mul(*[(D - e) for e in (r1 + q1, r1 + q2, r2 + q1, r2 + q2)])
P4 = lambda T: (T**2 - l0**2) * (T**2 - l1**2)
Gratio = lambda T: (T / sig + (kk + 1) / 2) / (T / sig + (1 - kk) / 2)
for assign in (1, 2):
    a, b = (c, n) if assign == 1 else (n, c)
    pa = lambda D: -(2 * D + a - U) * (2 * D + 2 * a - U)
    pb = lambda D: (2 * D + b - U) * (2 * D + 2 * b - U)
    Hrat = P4(s - nu + a) / P4(s - nu + b) * (Gratio(s - nu + b) if assign == 1 else 1 / Gratio(s - nu + a))
    F = (-pa(nu - a)) * pb(nu - a - b) / (pb(nu - b) * (-pa(nu - b - a))) * p0(nu - b) / p0(nu - a) * Hrat
    num = sp.Poly(sp.expand(sp.numer(sp.together(F - 1))), nu)
    lc = sp.factor(num.LC())
    print(f'assign {assign}: cocycle numerator degree {num.degree()}; LEADING coefficient = {lc}; depends on exponents: {lc.has(r1, q1, q2)}', flush=True)
    print(f'   roots in k: {sp.solve(lc, kk)}', flush=True)
