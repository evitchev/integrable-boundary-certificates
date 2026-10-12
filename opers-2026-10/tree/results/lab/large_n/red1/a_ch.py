"""RED1 (A), charge level.
A3: oscillator charges from the engine (T^2 - l^2) psi = x^2 (x^2 + E) psi (n = 2, M = 1), monic, == 2^(2K) B_2K(1/2 + l/2).
A4: LIMIT k -> infinity of the monic Sol-3 operator charges (exact rationals at k = 10^4 and 10^8):
    F_(2K-1)(PX, pi; k) -> B_2K(1/2 + PX) + B_2K(1/2 + pi), K = 1..5.  Control: without the 1/2 it must NOT converge.
Exit 0 = as sealed and control fires; 2 = not as sealed; 1 = control misbehaves."""
import hashlib, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_gen as MG
import mellin_wkb as MW
assert hashlib.sha256(open('SEAL_RED1.md', 'rb').read()).hexdigest() == open('SEAL_RED1.sha256').read().split()[0]
PX, PI, L = sp.symbols('PX PI L')
t0 = time.time()
ok = True
# A3
K = 10
h = [MG.ONE, MG.ZERO, MG.P_(-MG.l0**2)] + [MG.ZERO] * (K - 1)
eng = MG.MellinWKB(2, 1, h); eng.solve(K)
for i in range(2, K + 1, 2):
    R, _ = eng.R(i)
    Pm = sp.Poly(R, MG.l0)
    mon = sp.Poly(sp.expand(Pm.as_expr() / Pm.LC()), MG.l0)
    bern = sp.Poly(sp.expand(2**i * sp.bernoulli(i, sp.Rational(1, 2) + MG.l0 / 2)), MG.l0)
    same = mon == bern
    ok &= same
    print('A3 spin %d: monic oscillator charge == 2^%d B_%d(1/2 + l/2): %s' % (i - 1, i, i, same), flush=True)
for i in range(3, K + 1, 2):
    R, _ = eng.R(i)
    ok &= sp.expand(R) == 0
# A4
errs = {}
ctrl = {}
for kexp in (4, 8):
    k = sp.Integer(10) ** kexp
    polys, odd = MW.family_monic(k, 10, PX, PI)
    for s in (1, 3, 5, 7, 9):
        Pm, lead = polys[s]
        i = s + 1
        target = sp.Poly(sp.expand(sp.bernoulli(i, sp.Rational(1, 2) + PX) + sp.bernoulli(i, sp.Rational(1, 2) + PI)), PX, PI)
        wrong = sp.Poly(sp.expand(sp.bernoulli(i, PX) + sp.bernoulli(i, PI)), PX, PI)
        mons = set(Pm.monoms()) | set(target.monoms())
        errs[kexp, s] = max(abs(Pm.coeff_monomial(m) - target.coeff_monomial(m)) for m in mons)
        mons = set(Pm.monoms()) | set(wrong.monoms())
        ctrl[kexp, s] = max(abs(Pm.coeff_monomial(m) - wrong.coeff_monomial(m)) for m in mons)
    print('   k = 1e%d done (%.0fs)' % (kexp, time.time() - t0), flush=True)
cfire = True
for s in (1, 3, 5, 7, 9):
    r = errs[4, s] / errs[8, s] if errs[8, s] != 0 else sp.oo
    good = (errs[8, s] == 0 and errs[4, s] == 0) or (errs[4, s] < sp.Rational(1, 100) and sp.Rational(10**4, 2) < r < 2 * 10**4)
    ok &= bool(good)
    fire = ctrl[8, s] > sp.Rational(1, 1000) if s > 1 else True
    cfire &= bool(fire) if s > 1 else True
    print('A4 spin %d: max|monic F - (B(1/2+PX) + B(1/2+pi))| = %.4e (k=1e4), %.4e (k=1e8), ratio %.4e -> %s;   control without 1/2: %.3e (k=1e8)'
          % (s, float(errs[4, s]), float(errs[8, s]), float(r) if r != sp.oo else float('inf'), 'O(1/k)' if good else 'NOT converging', float(ctrl[8, s])))
print('(A) charge level:', 'AS SEALED' if ok else 'NOT as sealed', '| control fires:', cfire)
sys.exit(1 if not cfire else (0 if ok else 2))
