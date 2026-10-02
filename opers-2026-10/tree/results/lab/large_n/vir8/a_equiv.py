"""VIR8 (a): the three-term ODE forms of Sols 3 and 2 against their Gamma forms, exact.
Usage: a_equiv.py <sol> <k> <order> [tamper]     (tamper: Lc = L, i.e. the middle term x^c L(T) instead of x^(c/2) L x^(c/2))
Exit 0 = every order agrees (real) / some order disagrees (tamper); 2 otherwise."""
import hashlib, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3 as W3
import mellin_gen as MG
assert hashlib.sha256(open('SEAL_VIR8.md', 'rb').read()).hexdigest() == open('SEAL_VIR8.sha256').read().split()[0]
sol = int(sys.argv[1]); k = sp.Rational(sys.argv[2]); K = int(sys.argv[3]); tamper = len(sys.argv) > 4
M = 1 / k
l0, l1 = W3.l0, W3.l1
t0 = time.time()
if sol == 3:
    n = k + 4
    c = n * (1 + M)
    sh = 0 if tamper else c / 2
    eng = W3.WKB3(4, 1, M, [1, 0, -(l0**2 + l1**2), 0, l0**2 * l1**2], [1, -sh])
    blocks = [[1, 0, -MG.l0**2] + [0] * (K - 1), [1, 0, -MG.l1**2] + [0] * (K - 1), MG.string_series(k, 0, n / k, K + 1)]
else:
    n = 2 * k + 3
    c = n * (1 + M)
    sh = 0 if tamper else c / 2
    eng = W3.WKB3(3, 2, M, [1, 0, -l1**2, 0], [1, -2 * sh, sh**2 - l0**2])
    blocks = [[1, 0, -MG.l1**2] + [0] * (K - 1), MG.string_series(k, MG.l0, n / k, K + 1), MG.string_series(k, -MG.l0, n / k, K + 1)]
assert eng.n == n
eng.solve(K)
g = MG.MellinWKB(n, M, MG.symbol_from_blocks(blocks, K + 1)); g.solve(K)


def monic(Pp):
    Pp = sp.Poly(Pp, l0, l1)
    if Pp.is_zero:
        return None
    lead = Pp.coeff_monomial(l0 ** Pp.total_degree())
    return None if lead == 0 else sp.Poly(sp.expand(Pp.as_expr() / lead), l0, l1)


ok = True
for i in range(2, K + 1):
    cB, cG = eng.J(i)
    Rg, _ = g.R(i)
    Rg = sp.Poly(sp.expand(sp.sympify(Rg).subs({MG.l0: l0, MG.l1: l1})), l0, l1)
    mB, mG, mR = monic(cB), monic(cG), monic(Rg)
    sameB = (cB.is_zero and Rg.is_zero) or (mB is not None and mR is not None and mB == mR)
    sameG = cG.is_zero or (mG is not None and mR is not None and mG == mR)
    ok &= sameB and sameG
    print('Sol %d k = %s%s order %d (spin %d): Gamma-form R zero: %s; cB zero: %s; cG zero: %s; monic(cB) == monic(R): %s  (%.0fs)'
          % (sol, k, ' TAMPER' if tamper else '', i, i - 1, Rg.is_zero, cB.is_zero, cG.is_zero, sameB and sameG, time.time() - t0), flush=True)
print('(a) Sol %d k = %s%s:' % (sol, k, ' TAMPER' if tamper else ''), 'AGREE' if ok else 'DISAGREE')
sys.exit((0 if ok else 2) if not tamper else (0 if not ok else 2))
