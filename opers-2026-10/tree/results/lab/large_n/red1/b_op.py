"""RED1 (B), operator level: the double-scaling limit of the Sol-3 / Sol-2 three-term operators at N = 0, exact operator algebra.
Operators are sums zeta^j f_j(D), D = zeta d/dzeta, with (zeta^j f(D))(zeta^l g(D)) = zeta^(j+l) f(D + l) g(D).
Rescaled equations (exact for every k):
  Sol 3 (R3): P4(tau) - yt^(1+M) (tau - (1+M)/2) - Et yt,          P4 = (tau^2-a^2)(tau^2-b^2), a^2 = PX^2/p^2, b^2 = pi^2/p^2, zeta = yt^gamma, gamma = (1+M)/3
  Sol 2 (R2): tau(tau^2-b^2) - yt^(1+M) ((tau-(1+M)/2)^2 - a^2) - Et yt,   a^2 = PX^2/(p^2 k), b^2 = pi^2/p^2,        zeta = yt^gamma, gamma = 1+M
Conjugation phi = exp(-th0 zeta/gamma) zeta^(..) chi:  tau -> tau_c = th0 zeta + q - gamma D.
Modes: real | tamper (th0 = a SIMPLE root of the characteristic polynomial: the next-to-leading coefficient must not vanish).
Exit 0 = as sealed / tamper fires; 2 = not as sealed."""
import hashlib, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_RED1.md', 'rb').read()).hexdigest() == open('SEAL_RED1.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
D, q, gam, a2, b2, th, PX, PI, k = sp.symbols('D q gamma a2 b2 th PX PI k')


def omul(A, B):
    out = {}
    for j, f in A.items():
        for l, g in B.items():
            out[j + l] = sp.expand(out.get(j + l, 0) + f.subs(D, D + l) * g)
    return {j: f for j, f in out.items() if f != 0}


def oadd(*ops):
    out = {}
    for A, c in ops:
        for j, f in A.items():
            out[j] = sp.expand(out.get(j, 0) + c * f)
    return {j: f for j, f in out.items() if f != 0}


def opow(A, m):
    out = {0: sp.Integer(1)}
    for _ in range(m):
        out = omul(out, A)
    return out


ONE = {0: sp.Integer(1)}
tau_c = {1: th, 0: q - gam * D}
status = 0
for sol in (3, 2):
    if sol == 3:
        # P4(tau_c) - zeta^3 (tau_c - 3 gamma/2)        [yt^(1+M) = zeta^3, (1+M)/2 = 3 gamma/2]
        t2 = opow(tau_c, 2); t4 = omul(t2, t2)
        op = oadd((t4, 1), (t2, -(a2 + b2)), (ONE, a2 * b2), (omul({3: sp.Integer(1)}, oadd((tau_c, 1), (ONE, -sp.Rational(3, 2) * gam))), -1))
        top, char = 4, th**4 - th
        dbl = sp.Poly(4 * th**3 - 1, th)                      # double-root condition
        simple_root = sp.Integer(1)                            # th = 1 solves th^4 - th = 0 (E_0 = 0), simple
        kN0 = sp.Integer(-4)
        gamma_k = (1 + 1 / k) / 3; a2_k = PX**2 * (k + 1) / k; b2_k = PI**2 * (k + 1) / k
    else:
        # tau_c (tau_c^2 - b^2) - zeta ((tau_c - gamma/2)^2 - a^2)   [yt^(1+M) = zeta, (1+M)/2 = gamma/2]
        t2 = opow(tau_c, 2); t3 = omul(t2, tau_c)
        sh = oadd((tau_c, 1), (ONE, -gam / 2))
        op = oadd((t3, 1), (tau_c, -b2), (omul({1: sp.Integer(1)}, oadd((omul(sh, sh), 1), (ONE, -a2))), -1))
        top, char = 3, th**3 - th**2
        dbl = sp.Poly(3 * th - 2, th)
        simple_root = sp.Integer(1)                            # th = 1 solves th^3 - th^2 = 0 (E_0 = 0), simple
        kN0 = sp.Rational(-3, 2)
        gamma_k = 1 + 1 / k; a2_k = PX**2 * (k + 1) / k**2; b2_k = PI**2 * (k + 1) / k
    red = (lambda e_: sp.expand(sp.rem(sp.Poly(sp.expand(e_), th), dbl).as_expr())) if mode == 'real' else (lambda e_: sp.expand(e_.subs(th, simple_root)))
    c_top = red(op.get(top, 0) - char)
    c_next = red(op.get(top - 1, 0))
    print('Sol %d: zeta^%d coefficient minus the characteristic polynomial (%s): %s' % (sol, top, char, c_top))
    print('       zeta^%d coefficient at the %s root: %s' % (top - 1, 'double' if mode == 'real' else 'SIMPLE (tamper)', sp.factor(c_next)))
    if mode == 'tamper':
        fires = c_next != 0
        print('       tamper: next-to-leading coefficient nonzero: %s' % fires)
        if not fires:
            status = 1
        continue
    if c_top != 0 or c_next != 0:
        status = 2
    c2 = red(op.get(top - 2, 0))
    P2 = sp.Poly(c2, D)
    lead = P2.coeff_monomial(D**2)
    # c2 = lead * [ (D - d0)^2 - nu2 ]
    d0 = sp.cancel(-P2.coeff_monomial(D) / (2 * lead))
    nu2 = sp.cancel(d0**2 - P2.coeff_monomial(1) / lead)
    print('       zeta^%d coefficient = %s * [ (D - d0)^2 - nu^2 ],  d0 = %s' % (top - 2, sp.factor(lead), sp.factor(d0)))
    print('       nu^2 (general k, symbols a2, b2, gamma) = %s' % sp.factor(nu2))
    # value at N = 0
    th0val = sp.Rational(1, 4) ** sp.Rational(1, 3) if sol == 3 else sp.Rational(2, 3)
    nu2_k = nu2.subs({gam: gamma_k, a2: a2_k, b2: b2_k})
    if sol == 2:
        nu2_k = nu2_k.subs(th, sp.Rational(2, 3))
    nu2_N0 = sp.simplify(sp.limit(sp.simplify(nu2_k), k, kN0))
    nu2_N0 = sp.expand(nu2_N0)
    target = 2 * (PX**2 + PI**2) - sp.Rational(1, 6)
    okn = sp.simplify(nu2_N0 - target) == 0
    print('       at N = 0 (k -> %s): nu^2 = %s ;  == 2(PX^2 + pi^2) - 1/6 = 4 Delta: %s' % (kN0, nu2_N0, okn))
    if not okn:
        status = 2
    # the E-term: zeta^top (E_0 - Et zeta^(-eps)) with Et = E_0 (1 + eps e - (eps/2) log eps), zeta = Z/sqrt(eps)
    eps, e_, Z = sp.symbols('eps e Z', positive=True)
    expr = (1 - (1 + eps * e_ - eps / 2 * sp.log(eps)) * (Z / sp.sqrt(eps)) ** (-eps)) / eps
    lim = sp.limit(expr, eps, 0, '+')
    print('       double scaling: lim (1/eps)[1 - (Et/E_0) zeta^(-eps)] = %s   (expected log Z - e)' % sp.simplify(lim))
    if sp.simplify(lim - (sp.log(Z) - e_)) != 0:
        status = 2
    E0 = red(char)
    gcoef = sp.simplify(E0 / lead) if sol == 2 else sp.factor(sp.cancel(red(char) / lead))
    print('       limit operator: gamma^2-normalised  [ (D - d0)^2 - nu^2 ] chi + (E_0/lead) Z^2 (log Z - e) chi = 0,  E_0 = %s, E_0/lead = %s' % (E0, gcoef))
    # lower coefficients are suppressed by zeta^-1, zeta^-2 ~ sqrt(eps), eps relative to the zeta^(top-2) term
    lower = [j for j in op if j < top - 2]
    print('       remaining coefficients at zeta^j, j in %s: suppressed by eps^((%d-j)/2) -> 0' % (sorted(lower), top - 2))
print('(B) operator level (%s):' % mode, {0: 'AS SEALED' if mode == 'real' else 'tamper fires', 2: 'NOT as sealed', 1: 'tamper does not fire'}[status])
sys.exit(status)
