"""EXC2 (post-hoc): the exceptional locus.  Closed form: F(U) has the root U = 0 (z = 0, no apparent singularity) iff
a0 = +-(1 - b/2) or a1 = +-(1 - b/2).  Claim: there the lost level-1 state is the VACUUM-type oper with the exponent
1 - b/2 replaced by 1 + b/2 (a level-1 singular vector).  Test on the record's blocks, both charges, both fibres, symbolic
in the other momentum: char. polynomial of the block at a0 = 1 - b/2 has the root vacuum(a0 -> 1 + b/2); same for a1.
Exit 0 iff all eight hold and the control (shift to 1 + b/2 + 1/10) fails."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
P, Q, lam = sp.symbols('P Q lam')
ok = True; ctrl = False
for tstr in ('2', '9/4'):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1)
    dd = json.load(open('d1_sol3_t%s_with5.json' % tstr.replace('/', '_')))
    for name in ('I3', 'I5'):
        blk = dd['blocks'][name]
        Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
        vac = sp.sympify(blk['vacuum'], locals={'P': P, 'Q': Q})
        cp = sp.Poly(sp.expand((Mx - lam * sp.eye(2)).det()), P, Q, lam)
        vp = sp.Poly(vac, P, Q)
        X, Y = sp.symbols('X Y')          # X = P^2, Y = Q^2
        cpx = sum(c * X**(m[0] // 2) * Y**(m[1] // 2) * lam**m[2] for m, c in cp.terms())
        vx = sum(c * X**(m[0] // 2) * Y**(m[1] // 2) for m, c in vp.terms())
        lo, hi = (1 - b / 2)**2 / b, (1 + b / 2)**2 / b
        # X-side: PX^2 = lo ; Y-side: pi^2 = lo, i.e. Y = lo - rho2
        resX = sp.expand(cpx.subs(X, lo).subs(lam, vx.subs(X, hi)))
        resY = sp.expand(cpx.subs(Y, lo - rho2).subs(lam, vx.subs(Y, hi - rho2)))
        cX = sp.expand(cpx.subs(X, lo).subs(lam, vx.subs(X, (1 + b / 2 + sp.Rational(1, 10))**2 / b)))
        print('t = %-4s %s: vacuum at exponent 1 + b/2 is an eigenvalue of the level-1 block at 1 - b/2:  X side %s, Y side %s   (control, shifted exponent + 1/10: %s)'
              % (tstr, name, resX == 0, resY == 0, cX == 0))
        ok &= resX == 0 and resY == 0; ctrl |= cX == 0
print('x_exceptional:', 'PASS' if ok and not ctrl else 'FAIL')
sys.exit(0 if ok and not ctrl else 2)
