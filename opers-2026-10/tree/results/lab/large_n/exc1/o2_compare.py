"""EXC1 oper side, step 2: ONE apparent singularity in the xi-plane against the level-1 singlet block of I_3 (Sol 1).
(0) the trivial-monodromy condition derived directly from the potential at M = 5 (series at a root of x^6 = z):
        2M z^2 + alpha (M-1) z - 2 lam^2 + (M-1)^2/2 = 0;
(1) I_1: the one-point oper has I_1 = Delta + 1 (level 1);
(2) I_3: Res_z( condition, e - Lambda_3(z) ) against the characteristic polynomial of the data block (d1_sol1_t*.json),
    as polynomials in e with coefficients in (PX^2, pi^2) -- exact, symbolic in the momenta.
Usage: o2_compare.py [tamper]     Exit 0 = (1), (2) hold at t = 2 and 9/4 (real) / tamper fires; 2 otherwise."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
tamper = len(sys.argv) > 1
l0, l1, z1, z2, e, PX, PI = sp.symbols('l0 l1 z1 z2 e PX PI')
P, Q, lam = sp.symbols('P Q lam')
status = 0
if not tamper:
    # (0) direct derivation at M = 5
    x, x0, al, la = sp.symbols('x x0 alpha lam_', positive=True)
    M = 5
    d = sp.Symbol('d')
    z = x0**(M + 1)
    V0 = lambda xx: xx**(2 * M) + al * xx**(M - 1) + (la**2 - sp.Rational(1, 4)) / xx**2
    f = ((x0 + d)**(M + 1) - z)
    dV = -2 * sp.diff(sp.log(f), d, 2)
    reg = sp.series(sp.together(dV) - 2 / d**2, d, 0, 2).removeO()
    w1 = sp.diff(V0(x0 + d), d).subs(d, 0) + reg.coeff(d, 1)
    cond = sp.expand(sp.simplify(w1 * x0**3))
    target = sp.expand((2 * M * z**2 + al * (M - 1) * z - 2 * la**2 + sp.Rational((M - 1)**2, 2)))
    ok0 = sp.simplify(cond - target) == 0
    print('(0) M = 5: x0^3 x d/dx[regular part] at a root of x^6 = z  ==  2M z^2 + alpha(M-1) z - 2 lam^2 + (M-1)^2/2: %s' % ok0)
    if not ok0:
        status = 2
for tstr in (('2',) if tamper else ('2', '9/4')):
    t = sp.Rational(tstr)
    M = (t + 3) / (t - 1)
    rho2 = 2 / (t * t - 1)
    od = json.load(open('o1_t%s_n1%s.json' % (tstr.replace('/', '_'), '_tamper' if tamper else '')))
    loc = {'Symbol': sp.Symbol}
    R2 = sp.sympify(od['R']['2'], locals=loc); R4 = sp.sympify(od['R']['4'], locals=loc)
    c2 = sp.Poly(R2, l0, l1, z1).coeff_monomial(l1**2); c4 = sp.Poly(R4, l0, l1, z1).coeff_monomial(l1**4)
    # (1) I_1
    I1 = sp.expand((R2 / c2).subs({l0**2: 4 * M * (M + 1) * PX**2, l1**2: (M + 1) * PI**2}) / (M + 1))
    lvl = sp.expand(I1 - (PX**2 + PI**2 - sp.Rational(1, 6))) / 2
    print('t = %s (M = %s): one apparent singularity: I_1 = PX^2 + pi^2 - 1/6 + 2 x (%s)  -> level %s' % (tstr, M, lvl, lvl))
    if lvl != 1 and not tamper:
        status = 2
    # (2) I_3
    cond = 2 * M * z1**2 + l0 * (M - 1) * z1 - 2 * l1**2 + (M - 1)**2 / 2
    vac_lead = sp.Poly(sp.expand((R4 / c4).subs(z1, 0)), l0, l1).coeff_monomial(l0**4) * (4 * M * (M + 1))**2     # PX^4 coefficient of the normalised vacuum
    L3 = sp.expand(R4 / c4 / vac_lead)
    res = sp.resultant(sp.expand(cond), sp.expand(e - L3), z1)
    res = sp.expand(res)
    assert sp.Poly(res, l0).as_dict() and all(k[0] % 2 == 0 for k in sp.Poly(res, l0).as_dict()), 'resultant not even in alpha'
    resP = sp.Poly(res, l0)
    res_m = sum(c * (4 * M * (M + 1) * PX**2) ** (k[0] // 2) for k, c in resP.as_dict().items())
    res_m = sp.Poly(sp.expand(res_m.subs(l1**2, (M + 1) * PI**2)), e)
    res_m = sp.Poly(sp.expand(res_m.as_expr() / res_m.LC()), e)
    # data
    dd = json.load(open('d1_sol1_t%s.json' % tstr.replace('/', '_')))['blocks']['I3']
    Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
    vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
    lead = vac.coeff_monomial(P**4)
    cp = sp.expand((Mx / lead - e * sp.eye(2)).det())
    cp = sp.Poly(sp.expand(cp.subs({P: PX}).subs(Q**2, PI**2 - rho2)), e)
    assert not cp.as_expr().has(Q), 'odd power of Q in the characteristic polynomial'
    same = sp.expand(cp.as_expr() - res_m.as_expr()) == 0
    print('   I_3 level-1 block: characteristic polynomial (data) == Res_z(trivial monodromy, e - Lambda_3(z)) (oper): %s' % same)
    if not tamper:
        print('      data trace  = %s' % sp.factor(-cp.coeff_monomial(e)))
        print('      oper trace  = %s' % sp.factor(-res_m.coeff_monomial(e)))
        print('      discriminant (data) = %s' % sp.factor(sp.discriminant(cp.as_expr(), e)))
    if tamper:
        status = 0 if not same else 1
    elif not same:
        status = 2
print('EXC1 one-point level-1 test%s:' % (' (tamper)' if tamper else ''), {0: 'PASS' if not tamper else 'tamper fires', 2: 'FAIL', 1: 'tamper does not fire'}[status])
sys.exit(status)
