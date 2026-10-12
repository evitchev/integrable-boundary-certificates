"""RED1 (B), charge level.  Limits with exact rationals at two step sizes h = 1e-8, 1e-16.
B1: Sol-3 operator charges, k -> -4 (n -> 0): monic polynomials -> certified Sol-3 table at t = 5; the limit depends on
    PX^2 + pi^2 only; and equals the M -> 0 limit of the BLZ charges at l^2 = 4 Delta (1+M)^2, Delta = (PX^2+pi^2)/2 - 1/24.
    Spin 9 also against the polynomial printed in item 271 (Q = 0, Z = P_art^2 = -PX^2).
B2: the same for the Sol-2 operator, k -> -3/2, against the Sol-2 table at t = -5.
Exit 0 = as sealed; 2 otherwise."""
import hashlib, json, pickle, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
import mellin_gen as MG
assert hashlib.sha256(open('SEAL_RED1.md', 'rb').read()).hexdigest() == open('SEAL_RED1.sha256').read().split()[0]
import os as _os_kit; ROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..') + '/'  # KIT PATCH: results/lab/ from __file__ (was an absolute path)
PX, PI, ts, pp = sp.symbols('PX PI t p')
t0 = time.time()
Fsym = {k_: sp.Poly(sp.sympify(v), PX, PI) for k_, v in pickle.load(open('cft_F.pkl', 'rb')).items()}


def table(sol, s, tv):
    rho2 = sp.Rational(2) / (tv**2 - 1)
    if sol == 3 and s == 3:
        p2 = 2 * (tv - 1) / (tv + 1)
        out = 0
        for (i, j), cf in Fsym[4].terms():
            num, den = sp.fraction(sp.together(cf))
            pn, pd = sp.Poly(num, pp), sp.Poly(den, pp)
            out += sum(c_ * p2 ** (e // 2) for (e,), c_ in pn.terms()) / sum(c_ * p2 ** (e // 2) for (e,), c_ in pd.terms()) * PX**i * PI**j
        return sp.Poly(out, PX, PI)
    if s == 3:
        return None
    tab = json.load(open(ROOT + 'vev/vev_sol%d_w%d.json' % (sol, s + 1)))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        f = sp.sympify(val, locals={'t': ts}, rational=True)
        num, den = sp.fraction(sp.together(f))
        assert den.subs(ts, tv) != 0, 'table pole at t = %s' % tv
        expr += (num.subs(ts, tv) / den.subs(ts, tv)) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    P = sp.Poly(sp.expand(expr), PX, PI)
    return sp.Poly(P.as_expr() / P.coeff_monomial(PX ** (s + 1)), PX, PI)


def sol3_monic(k, K):
    polys, _ = MW.family_monic(k, K, PX, PI)
    return {s: polys[s][0] for s in polys}


def sol2_monic(k, K):
    n = 2 * k + 3
    blocks = [[1, 0, -MG.l1**2] + [0] * (K - 1), MG.string_series(k, MG.l0, n / k, K + 1), MG.string_series(k, -MG.l0, n / k, K + 1)]
    g = MG.MellinWKB(n, 1 / k, MG.symbol_from_blocks(blocks, K + 1)); g.solve(K)
    tv = (3 * k + 2) / (k + 2); p2 = 2 * (tv - 1) / (tv + 1)
    S1 = n * n / p2; S0 = S1 / k
    out = {}
    for i in range(2, K + 1, 2):
        R, _ = g.R(i)
        Pl = sp.Poly(R, MG.l0, MG.l1)
        Rm = sp.Poly(sp.expand(sum(cf * S0 ** (a // 2) * S1 ** (b // 2) * PX**a * PI**b for (a, b), cf in Pl.terms())), PX, PI)
        out[i - 1] = sp.Poly(sp.expand(Rm.as_expr() / Rm.coeff_monomial(PX**i)), PX, PI)
    return out


def blz_monic(M, K):
    h = [MG.ONE, MG.ZERO, MG.P_(-MG.l0**2)] + [MG.ZERO] * (K - 1)
    e = MG.MellinWKB(2, M, h); e.solve(K)
    out = {}
    for i in range(2, K + 1, 2):
        R, _ = e.R(i)
        Pm = sp.Poly(R, MG.l0)
        # substitute l^2 = (2(PX^2+PI^2) - 1/6)(1+M)^2, then make monic in PX
        ex = sum(c * ((2 * (PX**2 + PI**2) - sp.Rational(1, 6)) * (1 + M)**2) ** (d // 2) for (d,), c in Pm.terms())
        Q = sp.Poly(sp.expand(ex), PX, PI)
        out[i - 1] = sp.Poly(sp.expand(Q.as_expr() / Q.coeff_monomial(PX**i)), PX, PI)
    return out


def maxdiff(A, B):
    mons = set(A.monoms()) | set(B.monoms())
    return max(abs(A.coeff_monomial(m) - B.coeff_monomial(m)) for m in mons)


K = 10
status = 0
blz = {hexp: blz_monic(sp.Rational(1, 10**hexp), K) for hexp in (8, 16)}
print('BLZ engine at M = 1e-8, 1e-16 done (%.0fs)' % (time.time() - t0), flush=True)
Zs = sp.Symbol('Z')
printed9 = -Zs**5 - sp.Rational(15, 4) * Zs**4 - sp.Rational(51, 8) * Zs**3 - sp.Rational(1243, 224) * Zs**2 - sp.Rational(18351, 8960) * Zs - sp.Rational(60507, 394240)
for name, solfun, k0, tv0, sol in (('B1 Sol 3, k -> -4, t = 5', sol3_monic, sp.Integer(-4), sp.Integer(5), 3), ('B2 Sol 2, k -> -3/2, t = -5', sol2_monic, sp.Rational(-3, 2), sp.Integer(-5), 2)):
    ops = {hexp: solfun(k0 + sp.Rational(1, 10**hexp), K) for hexp in (8, 16)}
    print('%s: operator runs done (%.0fs)' % (name, time.time() - t0), flush=True)
    for s in (3, 5, 7, 9):
        cert = table(sol, s, tv0)
        line = '   spin %d:' % s
        if cert is not None:
            e8, e16 = maxdiff(ops[8][s], cert), maxdiff(ops[16][s], cert)
            r = e8 / e16 if e16 != 0 else sp.oo
            good = e8 < sp.Rational(1, 1000) and sp.Rational(10**8, 2) < r < 2 * 10**8
            line += ' |operator - certified table| = %.3e, %.3e (%s);' % (float(e8), float(e16), 'O(h)' if good else 'NOT converging')
            if not good:
                status = 2
            # rotation invariance of the certified limit: function of PX^2 + PI^2 only
            ex = cert.as_expr()
            rot = sp.expand(ex.subs({PX: sp.Symbol('r'), PI: 0}) .subs(sp.Symbol('r'), sp.sqrt(PX**2 + PI**2)) - ex) == 0
            line += ' table depends on PX^2 + pi^2 only: %s;' % rot
            if not rot:
                status = 2
        # operator limit vs BLZ limit (both at the same h): difference O(h)
        d8, d16 = maxdiff(ops[8][s], blz[8][s]), maxdiff(ops[16][s], blz[16][s])
        r2 = d8 / d16 if d16 != 0 else sp.oo
        good2 = d8 < sp.Rational(1, 1000) and sp.Rational(10**8, 4) < r2 < 4 * 10**8
        line += ' |operator - BLZ(M=h) at l^2 = 4 Delta (1+M)^2| = %.3e, %.3e (%s)' % (float(d8), float(d16), 'O(h)' if good2 else 'NOT converging')
        if not good2:
            status = 2
        if s == 9:
            # item 271's printed spin-9 polynomial at Q_art = 0 (pi^2 = rho^2 = 1/12), Z = -PX^2, normalised with leading -Z^5
            ex9 = sp.expand(printed9.subs(Zs, -PX**2))
            op9 = sp.Poly(ops[16][s].as_expr().subs(PI, sp.sqrt(sp.Rational(1, 12))), PX)
            d9 = max(abs(op9.coeff_monomial(m) - sp.Poly(ex9, PX).coeff_monomial(m)) for m in set(op9.monoms()) | set(sp.Poly(ex9, PX).monoms()))
            line += '; vs item 271 printed polynomial (h = 1e-16): %.3e' % float(d9)
            if d9 > sp.Rational(1, 10**10):
                status = 2
        print(line, flush=True)
print('(B) charge level:', 'AS SEALED' if status == 0 else 'NOT as sealed')
sys.exit(status)
