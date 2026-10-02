"""EXC2 U2 (full) and U3 (prediction): interpolate trace and determinant of the level-1 I_3 and I_5 blocks, as polynomials in
the record's bare variables (P^2 = PX^2, Q^2 = pi^2 - rho^2), from the self-adjoint two-point opers at the batch points.
Degree bounds: I_3 trace 2, det 4; I_5 trace 3, det 6 (in P^2, Q^2).  Fit on the first NFIT points, CHECK on the rest.
I_3: compared with the stored data block (U2 as a polynomial identity).  I_5: written to the prediction; NO I_5 data opened.
Usage: u3_interp.py <t> [flip]    Exit 0 iff the I_3 polynomials == data and all held-back points agree; 2 otherwise."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr = sys.argv[1]; flip = len(sys.argv) > 2
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1)
tt = tstr.replace('/', '_')
P, Q, e = sp.symbols('P Q e')
pts = []
for line in open('u3_points_t%s.txt' % tt):
    _, a0s, a1s, _ = line.split()
    TAG = 't%s_%s_%s' % (tt, a0s.replace('/', '-'), a1s.replace('/', '-'))
    d = json.load(open('u3p_%s%s.json' % (TAG, '_flip' if flip else '')))
    assert d['vdim_ordered'] == 4
    X = sp.Rational(d['PX2']); Y = sp.Rational(d['PI2']) - rho2
    p3 = sp.Poly(sp.sympify(d['I3_poly'], locals={'ee': sp.Symbol('ee')}), sp.Symbol('ee')); p5 = sp.Poly(sp.sympify(d['I5_poly'], locals={'ff': sp.Symbol('ff')}), sp.Symbol('ff'))
    assert p3.degree() == 2 and p5.degree() == 2
    pts.append((X, Y, -p3.all_coeffs()[1], p3.all_coeffs()[2], -p5.all_coeffs()[1], p5.all_coeffs()[2]))
def fit(idx, deg, nfit):
    mons = [(i, j) for i in range(deg + 1) for j in range(deg + 1 - i)]
    assert nfit >= len(mons)
    A = sp.Matrix([[X**i * Y**j for (i, j) in mons] for (X, Y, *_) in pts[:nfit]])
    rhs = sp.Matrix([p[idx] for p in pts[:nfit]])
    sol, params = A.gauss_jordan_solve(rhs) if nfit == len(mons) else ((A.T * A).LUsolve(A.T * rhs), None)
    resid = A * sol - rhs
    exact_fit = all(r == 0 for r in resid)
    poly = sum(c * P**(2 * i) * Q**(2 * j) for c, (i, j) in zip(sol, mons))
    bad = [(X, Y) for (X, Y, *rest) in pts[nfit:] if sum(c * X**i * Y**j for c, (i, j) in zip(sol, mons)) != (X, Y, *rest)[idx]]
    return sp.expand(poly), exact_fit, len(pts) - nfit, len(bad)
out = {'t': tstr, 'normalisation': 'eigenvalues divided by the P^(2K) coefficient of the vacuum eigenvalue; P^2 = PX^2, Q^2 = pi^2 - rho^2 (bare variables of lab/excited_states.py)'}
ok = True
for name, idx, deg, nfit in (('I3_trace', 2, 2, 12), ('I3_det', 3, 4, 24), ('I5_trace', 4, 3, 18), ('I5_det', 5, 6, 30)):
    poly, exact_fit, nchk, nbad = fit(idx, deg, nfit)
    print('%s: degree <= %d in (P^2, Q^2); fitted on %d points (exact fit %s); held-back points %d, mismatches %d' % (name, deg, nfit, exact_fit, nchk, nbad))
    print('   ', poly)
    ok &= exact_fit and nbad == 0
    out[name] = str(poly)
# I_3 against the stored data block
dd = json.load(open('d1_sol3_t%s.json' % tt))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
Mn = Mx / vac.coeff_monomial(P**4)
tr_ok = sp.expand(Mn.trace() - sp.sympify(out['I3_trace'], locals={'P': P, 'Q': Q})) == 0
det_ok = sp.expand(Mn.det() - sp.sympify(out['I3_det'], locals={'P': P, 'Q': Q})) == 0
print('I_3 level-1 block (data) : trace identical %s, determinant identical %s' % (tr_ok, det_ok))
ok &= tr_ok and det_ok
json.dump(out, open('u3_interp_t%s%s.json' % (tt, '_flip' if flip else ''), 'w'), indent=1)
print('u3_interp t = %s%s: %s' % (tstr, ' FLIP' if flip else '', 'PASS' if ok else 'FAIL'))
sys.exit(0 if ok else 2)
