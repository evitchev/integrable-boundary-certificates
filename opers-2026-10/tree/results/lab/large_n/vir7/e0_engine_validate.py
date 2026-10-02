"""VIR7 engine validation (ODE side only, before any Sol-1 comparison):
 (i)  mellin_gen (general, odd powers allowed) == mellin_wkb (even) on the Sol-3 operator at k = 2 and k = 10/3, orders 2..10;
 (ii) BLZ point of the class-U equivalence:  (T^2 - l^2) G_a(T), n = 2 + a, M = 1/a   versus the Schrodinger
      equation T^2 - lS^2 = x^2 (x^(2 M_S) - E), M_S = 1 + 2/a, with l = (2 + a) lS / 2:  monic charge polynomials equal.
Exit 0 iff both hold."""
import sys
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
import mellin_gen as MG
ok = True
for k in (sp.Integer(2), sp.Rational(10, 3)):
    K = 10
    n = k + 4
    h_even = MW.family_symbol(k, K // 2 + 1)
    e1 = MW.MellinWKB(n, 1 / k, h_even); e1.solve(K)
    h_gen = []
    for j in range(K + 1):
        h_gen.append(MG.P_(h_even[j // 2].as_expr()) if j % 2 == 0 else MG.ZERO)
    e2 = MG.MellinWKB(n, 1 / k, h_gen); e2.solve(K)
    for i in range(2, K + 1):
        r1, _ = e1.R(i); r2, _ = e2.R(i)
        same = sp.expand(r1 - r2) == 0
        ok &= same
    print('(i) Sol-3 operator k = %s: general engine == even engine at orders 2..%d: %s' % (k, K, ok), flush=True)
for a in (sp.Rational(1, 2), sp.Rational(5, 8), sp.Integer(2)):
    K = 8
    n = 2 + a
    MS = 1 + 2 / a
    # dual: (1 - l0^2/T^2) x frozen string of length a at 0, step sigma = n/a
    blocks = [[1, 0, -MG.l0**2] + [0] * (K - 1), MG.string_series(a, 0, n / a, K + 1)]
    hd = MG.symbol_from_blocks(blocks, K + 1)
    ed = MG.MellinWKB(n, 1 / a, hd); ed.solve(K)
    # Schrodinger: symbol T^2 - l0^2, n = 2, M = M_S
    hs = [MG.ONE, MG.ZERO, MG.P_(-MG.l0**2)] + [MG.ZERO] * (K - 1)
    es = MG.MellinWKB(2, MS, hs); es.solve(K)
    for i in range(2, K + 1, 2):
        rd, _ = ed.R(i); rs, _ = es.R(i)
        pd = sp.Poly(rd, MG.l0); ps = sp.Poly(rs, MG.l0)
        if pd.is_zero or ps.is_zero:
            print('     a = %s order %d: dual zero: %s, Schrodinger zero: %s' % (a, i, pd.is_zero, ps.is_zero)); continue
        # l_dual = (2 + a)/2 * l_S
        ps2 = sp.Poly(ps.as_expr().subs(MG.l0, 2 * MG.l0 / (2 + a)), MG.l0)
        same = sp.expand(pd.as_expr() / pd.LC() - ps2.as_expr() / ps2.LC()) == 0
        ok &= same
        print('(ii) a = %s (n = %s, M = %s) vs Schrodinger M_S = %s, spin %d: monic polynomials equal: %s' % (a, n, 1 / a, MS, i - 1, same), flush=True)
    for i in range(3, K + 1, 2):
        rd, _ = ed.R(i)
        ok &= sp.expand(rd) == 0
print('engine validation', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
