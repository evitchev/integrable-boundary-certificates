"""EXC4 A2: the DUAL Solution-2 form (item 303: M' = -1/(k+1), identical charges).  Three-term form with the same P = T(T^2 - l1^2),
L = T^2 - l0^2 but c' = n(1 + M') = n b, so beta' = n/c' = 1/b.  Normal form (xi' = x^(c')/c'):
    th(th^2 - B') + xi'((th + 1/2)^2 - A') - E xi'^(beta'),   a1' = l1/c' = pi/(p b),  a0' = l0/c' = P_X/(p b sqrt(k)),
i.e. A' = P_X^2/(k b^2) = (1-b) P_X^2 / b^2,  B' = pi^2/b^3   [b = p^2].
T0: the dual vacuum through wkb3x2 (M' = -1/(k+1)) reproduces the certified Sol 2 vacuum I_1, I_3, I_5 at the fibre (d0_sol2_t<t>.json).
Then ONE apparent singularity with S_3 = {-1,1,3}: existence, r11, and the monic I_1 shift (sealed prediction: r11 = -beta' = -1/b,
Delta I_1 = p^2/beta' = b^2 = 4/9 at t = 2).  Usage: f_dual_sol2.py <t> <P> <Q>  (bare momenta).  Exit 0 iff T0 passes and, if solutions
exist, Delta I_1 == b^2 on all of them; 3 = T0 passes but no solution; 2 = T0 fails or Delta I_1 != b^2."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
import wkb3x2 as W
assert hashlib.sha256(open('SEAL_EXC4.md', 'rb').read()).hexdigest() == open('SEAL_EXC4.sha256').read().split()[0]
tstr, Ps, Qs = sys.argv[1], sys.argv[2], sys.argv[3]
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = 2 * k + 3; rho2 = 2 / (t * t - 1); Nv = (t * t - 25) / (t * t - 1)
Mp = -1 / (k + 1); cp = n * (1 + Mp); betap = n / cp
assert betap == 1 / b and cp == n * b
Pv, Qv = sp.Rational(Ps), sp.Rational(Qs); PX2 = Pv**2; PI2 = Qv**2 + rho2
Ap = (1 - b) * PX2 / b**2; Bp = PI2 / b**3
print('t = %s: k = %s, b = %s, n = %s, dual M\' = %s, c\' = %s, beta\' = %s; A\' = %s, B\' = %s' % (tstr, k, b, n, Mp, cp, betap, Ap, Bp), flush=True)
t0 = time.time()
# ---- T0: dual vacuum through the engine (dP = 3, dL = 2, M = M')
eng = W.WKB3(3, 2, Mp, [1, 0, -W.l1**2, 0], [1, -cp, cp**2 / 4 - W.l0**2]); eng.solve(6)
P, Q = sp.symbols('P Q')
d0 = json.load(open('d0_sol2_t%s.json' % tstr.replace('/', '_')))['vacuum']
cert = {1: -(P**2 + Q**2) + (Nv + 1) / 12, 3: sp.sympify(d0['I3'], locals={'P': P, 'Q': Q}), 5: sp.sympify(d0['I5'], locals={'P': P, 'Q': Q})}
l1sq = cp**2 * (Q**2 + rho2) / b**3; l0sq = cp**2 * (1 - b) * P**2 / b**2        # l1 = n pi/p, l0 = n P_X/(p sqrt k) -- the same l's as the primal form
ok = True
for i in (2, 4, 6):
    cB, cG = eng.J(i); assert cG.is_zero
    expr = sp.expand(cB.as_expr().subs({W.l1**2: l1sq, W.l0**2: l0sq}))
    K = i // 2; lead = sp.Poly(expr, P, Q).coeff_monomial(P**(2 * K))
    mon = sp.expand(expr / lead); cm = sp.expand(cert[i - 1] / sp.Poly(cert[i - 1], P, Q).coeff_monomial(P**(2 * K)))
    same = sp.expand(mon - cm) == 0; ok &= same
    print('T0 dual vacuum spin %d == certified: %s' % (i - 1, same), flush=True)
if not ok:
    print('T0 FAILS: the dual form through this engine does not reproduce the vacuum; A2 cannot be tested this way'); sys.exit(2)
# ---- one apparent singularity, S_3 = {-1, 1, 3}
z, E, e = Fb.z, Fb.E, Fb.e
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (1, 0) for m in range(1, 4 - j)}
R = {j: {m: r[j, m] for m in range(1, 4 - j)} for j in (1, 0)}
Pcoef = [0, -Bp, 0, 1]; Lcoef = [sp.Rational(1, 4) - Ap, 1, 1]
exps = [-1, 1, 3]
lk = Fb.local_operator(Pcoef, Lcoef, betap, -1, R, exps[-1] - exps[0] - 3 + 1)
lead, conds, kaps = Fb.nolog_conditions(lk, exps, 3)
sol_ind = sp.solve([c for c in sp.Poly(sp.expand(lead - z**3 * sp.prod([(e - x) for x in exps])), e).coeffs() if c != 0], [r[1, 2], r[0, 3]], dict=True)[0]
print('indicial:', sol_ind)
conds = [(p_, sp.expand(c_.subs(sol_ind))) for p_, c_ in conds]
eqs = Fb.split_conditions(conds, kaps)
unknowns = [r[1, 1], r[0, 1], r[0, 2], z]
polys = []
for pos, mon, cf in eqs:
    num = sp.numer(sp.together(cf)); core = sp.Integer(1)
    for fac, mult in sp.factor_list(num)[1]:
        if fac != z: core *= fac**mult
    polys.append(sp.expand(core))
    kn = ' '.join('%s^%d' % (s_, p_) for s_, p_ in zip(['E'] + [str(x) for x in kaps], mon) if p_)
    print('   at exponent %d, [%s]: %s' % (pos, kn or '1', sp.factor(core) if len(str(core)) < 160 else '%d terms' % len(sp.Poly(core, *unknowns).terms())))
w = sp.Symbol('w')
G = sp.groebner(polys + [1 - w * z], w, *unknowns, order='lex')
if G.exprs == [1]:
    print('RESULT A2: no one-point solution with z != 0 in the dual form at this point (the rule is not tested here)'); sys.exit(3)
elz = [g for g in G.exprs if g.free_symbols <= {z}][0]
print('solutions: eliminant in z: %s' % sp.factor(elz))
# I_1 of the deformed operator: effective cubic th^3 + (-B' + r11) th + r01 -> engine A3 = c'^2 q1, A4 = -c'^3 q0 (I_1 needs no 1/xi term)
eng2 = W.WKB3(3, 2, Mp, [1, 0, W.A3, W.A4], [1, -cp, cp**2 / 4 - W.l0**2]); eng2.solve(2)
cB2, _ = eng2.J(2)
vac = {W.A3: -cp**2 * Bp, W.A4: 0, W.l0: cp * sp.sqrt(Ap)}
lead2 = sp.Poly(sp.expand(cB2.as_expr().subs({W.A3: -cp**2 * sp.Symbol('Bs'), W.A4: 0, W.l0**2: cp**2 * sp.Symbol('As')})), sp.Symbol('As'), sp.Symbol('Bs')).coeff_monomial(sp.Symbol('Bs')) / (1 / b**3)   # monic in pi^2 <-> P^2 coefficient: use the P_X^2 coefficient
# simpler: normalise by the P_X^2 coefficient of the vacuum I_1 (via A' = (1-b)P_X^2/b^2)
Asym, Bsym = sp.symbols('As Bs')
vacI1 = sp.expand(cB2.as_expr().subs({W.A3: -cp**2 * Bsym, W.A4: 0, W.l0**2: cp**2 * Asym}))
lead_PX2 = sp.Poly(vacI1, Asym, Bsym).coeff_monomial(Asym) * (1 - b) / b**2
I1 = sp.expand(cB2.as_expr().subs({W.A3: cp**2 * (-Bp + r[1, 1]), W.A4: -cp**3 * r[0, 1], W.l0**2: cp**2 * Ap}) / lead_PX2)
vacv = sp.expand(vacI1.subs({Asym: Ap, Bsym: Bp}) / lead_PX2)
ev = sp.Symbol('ev')
Ge = sp.groebner(G.exprs + [sp.expand(ev - I1)], w, *unknowns, ev, order='lex')
pe = [g for g in Ge.exprs if g.free_symbols <= {ev}][0]
shifts = [sp.nsimplify(0) if False else (root - vacv) for root in sp.solve(pe, ev)]
r11v = sp.solve([g for g in G.exprs if g.free_symbols <= {r[1, 1]}][0], r[1, 1]) if any(g.free_symbols <= {r[1, 1]} for g in G.exprs) else 'not univariate'
print('r11 on the solutions: %s  (sealed: -beta\' = %s)' % (r11v, -betap))
print('monic I_1 on the solutions minus the vacuum value: %s   (sealed prediction p^2/beta\' = b^2 = %s)' % (shifts, b**2))
good = all(sh == b**2 for sh in shifts)
print('A2:', 'RULE HOLDS' if good else 'RULE REFUTED', '(%.0fs)' % (time.time() - t0)); sys.exit(0 if good else 2)
