"""VIR8 (b): the three-term ODE forms exactly AT nZ points.  For each case prints whether cB, cG vanish, and, if not,
compares monic(cB) (and monic(cG)) with the certified charge of the named solution.  Exact.
Exit 0 = b1-b4 as sealed (cB != 0, == certified); 2 otherwise."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3 as W3
assert hashlib.sha256(open('SEAL_VIR8.md', 'rb').read()).hexdigest() == open('SEAL_VIR8.sha256').read().split()[0]
ROOT = '<repo>/results/lab/'
l0, l1 = W3.l0, W3.l1
ts, PX, PI = sp.symbols('t PX PI')


def table(sol, s, tv):
    rho2 = sp.Rational(2) / (tv**2 - 1)
    tab = json.load(open(ROOT + 'vev/vev_sol%d_w%d.json' % (sol, s + 1)))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        f = sp.sympify(val, locals={'t': ts}, rational=True)
        assert sp.fraction(sp.together(f))[1].subs(ts, tv) != 0
        expr += f.subs(ts, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    P = sp.Poly(sp.expand(expr), PX, PI)
    return sp.Poly(P.as_expr() / P.coeff_monomial(PX ** (s + 1)), PX, PI)


def to_PXPI(Pp, S0, S1):
    Pp = sp.Poly(Pp, l0, l1)
    out = sp.Poly(sp.expand(sum(c * S0 ** sp.Rational(i, 2) * S1 ** sp.Rational(j, 2) * PX**i * PI**j for (i, j), c in Pp.terms())), PX, PI)
    lead = out.coeff_monomial(PX ** out.total_degree())
    return None if lead == 0 else sp.Poly(sp.expand(out.as_expr() / lead), PX, PI)


cases = [('b1', 3, sp.Integer(3), 7, 3), ('b2', 3, sp.Integer(1), 5, 3), ('b3', 2, sp.Integer(1), 5, 2), ('b4', 3, sp.Rational(1, 2), 9, 3)]
status = 0
t0 = time.time()
for name, sol, k, spin, cert_sol in cases:
    M = 1 / k
    tv = (3 * k + 2) / (k + 2)
    p2 = 2 * (tv - 1) / (tv + 1)
    if sol == 3:
        n = k + 4; c = n * (1 + M)
        eng = W3.WKB3(4, 1, M, [1, 0, -(l0**2 + l1**2), 0, l0**2 * l1**2], [1, -c / 2])
        S1 = n * n / p2; S0 = S1
    else:
        n = 2 * k + 3; c = n * (1 + M)
        eng = W3.WKB3(3, 2, M, [1, 0, -l1**2, 0], [1, -c, c * c / 4 - l0**2])
        S1 = n * n / p2; S0 = S1 / k
    i = spin + 1
    eng.solve(i)
    try:
        cB, cG = eng.J(i)
    except ZeroDivisionError as e:
        print('%s: master reduction singular: %s' % (name, e)); status = 2; continue
    red, A0, B0 = eng.masters(i)
    cert = table(cert_sol, spin, tv)
    line = '%s  Sol-%d ODE at t = %s (k = %s, n = %s), spin %d (nu = %s; A0 = %s, B0 + 1 = %s): cB zero: %s; cG zero: %s' % (
        name, sol, tv, k, n, spin, sp.Rational(spin) / n, A0, B0 + 1, cB.is_zero, cG.is_zero)
    okc = False
    if not cB.is_zero:
        mB = to_PXPI(cB, S0, S1)
        okc = mB is not None and mB == cert
        line += '; monic(cB) == certified Sol-%d spin-%d charge: %s' % (cert_sol, spin, okc)
    if not cG.is_zero:
        mG = to_PXPI(cG, S0, S1)
        line += '; monic(cG) == certified: %s' % (mG is not None and mG == cert)
    print(line + '  (%.0fs)' % (time.time() - t0), flush=True)
    if not okc:
        status = 2
    # lower orders at the same point, for the record
    for j in range(2, i, 2):
        b_, g_ = eng.J(j)
        print('      order %d (spin %d): cB zero: %s; cG zero: %s' % (j, j - 1, b_.is_zero, g_.is_zero))
print('(b)', 'AS SEALED' if status == 0 else 'NOT as sealed')
sys.exit(status)
