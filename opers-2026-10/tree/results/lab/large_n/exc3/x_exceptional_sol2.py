"""EXC3 (post-hoc, after the hold-out): the exceptional locus of the level-one oper of Solution 2.  From the closed form,
c7 U^2 + c5 U + c3 = 0 degenerates when c7 = 0 [4A = (1+b)^2, a0 = (1+b)/2: one root U -> infinity] or c3 = 0 [4B = (2-b)^2, a1 = 1 - b/2,
or 4B = (2b-1)^2, a1 = b - 1/2: one root U = 0].  Claim (EXC2 analogue): there the lost level-one state is the VACUUM-type oper with the
exponent shifted by a level-one amount: a0 -> (3-b)/2 [A -> A + 2(1-b)], a1 -> 1 + b/2 or b + 1/2 [B -> B + 2b]; each shift raises the
monic I_1 = A/(1-b) + B/b + const by exactly 2 (= one level).  Test on the record's blocks (I_3 and I_5, both fibres), symbolic in the
other momentum: the vacuum eigenvalue at the shifted exponent is a root of the block's characteristic polynomial on the locus.
Exit 0 iff all twelve checks hold and the control (shift + 1/10) fails."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
P, Q, lam, X, Y = sp.symbols('P Q lam X Y')
ok = True; ctrl_hits = 0
for tstr in ('2', '9/4'):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1)
    dd = json.load(open('d1_sol2_t%s_with5.json' % tstr.replace('/', '_')))
    for name in ('I3', 'I5'):
        blk = dd['blocks'][name]
        Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
        cp = sp.Poly(sp.expand((Mx - lam * sp.eye(2)).det()), P, Q, lam); vp = sp.Poly(sp.sympify(blk['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
        # in X = P^2 = A/(1-b), Y = Q^2 = B/b - rho2
        cpx = sum(c * X**(m[0] // 2) * Y**(m[1] // 2) * lam**m[2] for m, c in cp.terms()); vx = sum(c * X**(m[0] // 2) * Y**(m[1] // 2) for m, c in vp.terms())
        cases = {'a0 = (1+b)/2 -> (3-b)/2': (X, ((1 + b) / 2)**2 / (1 - b), ((3 - b) / 2)**2 / (1 - b)),
                 'a1 = 1 - b/2 -> 1 + b/2': (Y, (1 - b / 2)**2 / b - rho2, (1 + b / 2)**2 / b - rho2),
                 'a1 = b - 1/2 -> b + 1/2': (Y, (b - sp.Rational(1, 2))**2 / b - rho2, (b + sp.Rational(1, 2))**2 / b - rho2)}
        for label, (var, lo, hi) in cases.items():
            res = sp.expand(cpx.subs(var, lo).subs(lam, vx.subs(var, hi)))
            ctrl = sp.expand(cpx.subs(var, lo).subs(lam, vx.subs(var, hi + sp.Rational(1, 10))))
            print('t = %-4s %s  %s: vacuum at the shifted exponent is an eigenvalue of the level-one block on the locus: %s   (control +1/10: %s)' % (tstr, name, label, res == 0, ctrl == 0))
            ok &= res == 0; ctrl_hits += (ctrl == 0)
print('x_exceptional_sol2:', 'PASS' if ok and ctrl_hits == 0 else 'FAIL'); sys.exit(0 if ok and ctrl_hits == 0 else 2)
