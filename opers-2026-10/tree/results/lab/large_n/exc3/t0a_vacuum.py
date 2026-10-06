"""EXC3 T0(a): the normal form of the seal (sec. 1) against the certified Sol 2 vacuum: the engine (dP = 3, dL = 2) with
P_T = T(T^2 - l1^2), Lc = (T - c/2)^2 - l0^2, l1 = c a1, l0 = c a0, a1^2 = b pi^2, a0^2 = (1-b) P_X^2, must reproduce the
level-0 eigenvalues of I_1, I_3, I_5 (d0_sol2_t<t>.json; bare P = P_X, Q^2 = pi^2 - rho^2) monic in P at the fibre.
Also reported (H3): the odd-order charges J_3, J_5 of the vacuum (even spins).  Usage: t0a_vacuum.py <t> [tamper]
(tamper: a0^2 -> (1-b) P_X^2 * 11/10, must fail).  Exit 0 iff spins 1, 3, 5 agree (tamper: iff they do not)."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
import wkb3x2 as W
tstr = sys.argv[1]; tamper = len(sys.argv) > 2
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = 2 * k + 3; M = 1 / k; c = n * (1 + M)
assert c == n / b
rho2 = 2 / (t * t - 1); Nv = (t * t - 25) / (t * t - 1)
t0 = time.time()
eng = W.WKB3(3, 2, M, [1, 0, -W.l1**2, 0], [1, -c, c**2 / 4 - W.l0**2])
eng.solve(6)
P, Q = sp.symbols('P Q')
PX2 = P**2; PI2 = Q**2 + rho2
a0sq = (1 - b) * PX2 * (sp.Rational(11, 10) if tamper else 1); a1sq = b * PI2
d0 = json.load(open('d0_sol2_t%s.json' % tstr.replace('/', '_')))['vacuum']
cert = {1: -(P**2 + Q**2) + (Nv + 1) / 12,                # article P^2 + Q^2 + (N+1)/12 with article P^2 = -P_X^2, Q^2 = -(bare Q^2)
        3: sp.sympify(d0['I3'], locals={'P': P, 'Q': Q}), 5: sp.sympify(d0['I5'], locals={'P': P, 'Q': Q})}
ok = True
for i in (2, 3, 4, 5, 6):
    cB, cG = eng.J(i)
    assert cG.is_zero
    expr = sp.expand(cB.as_expr().subs({W.l1**2: c**2 * a1sq, W.l0**2: c**2 * a0sq}))
    expr = sp.expand(expr.subs({W.l1: sp.sqrt(c**2 * a1sq), W.l0: sp.sqrt(c**2 * a0sq)}))
    if i % 2 == 1:
        print('odd order %d (spin %d): vacuum charge = %s' % (i, i - 1, sp.factor(expr)), flush=True)
        continue
    K = i // 2
    lead = sp.Poly(expr, P, Q).coeff_monomial(P**(2 * K))
    mon = sp.expand(expr / lead)
    cm = sp.expand(cert[i - 1] / sp.Poly(cert[i - 1], P, Q).coeff_monomial(P**(2 * K)))
    same = sp.expand(mon - cm) == 0
    print('spin %d: engine monic = %s' % (i - 1, mon))
    print('        certified    = %s   -> equal: %s  (%.0fs)' % (cm, same, time.time() - t0), flush=True)
    ok &= same
print('T0(a) normal form vs certified Sol 2 vacuum at t = %s%s: %s' % (tstr, ' TAMPER' if tamper else '', 'PASS' if ok else 'FAIL'))
sys.exit((0 if not ok else 1) if tamper else (0 if ok else 2))
