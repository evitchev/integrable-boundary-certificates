"""EXC3 U2/U3: parse the batch (u_one_sym_<TAG>.log), check I_1 = Delta + 1, J_3 = 0, sigma-invariance (kk = 0) at every point,
interpolate trace/det of the level-one I_3 and I_5 blocks as polynomials in the bare (P^2, Q^2) with held-back points, compare I_3
with the stored data block (polynomial identity), write u_interp_t<t>.json (I_5 polynomials = the prediction; NO I_5 data opened).
Usage: u_interp.py <t>.  Exit 0 iff all point checks pass, all fits exact with held-back points, and I_3 trace+det == data."""
import hashlib, json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
tstr = sys.argv[1]; t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1); Nv = (t * t - 25) / (t * t - 1)
tt = tstr.replace('/', '_')
P, Q = sp.symbols('P Q'); x = sp.Symbol('x')
pts, skipped, bad = [], [], []
for line in open('u_points_t%s.txt' % tt):
    _, _, a0s, a1s = line.split()
    TAG = 't%s_%s_%s' % (tt, a0s.replace('/', '-'), a1s.replace('/', '-'))
    out = open('u_one_sym_%s.log' % TAG).read()
    vd = int(re.search(r'vdim:\n(-?\d+)', out).group(1))
    if vd != 4:
        skipped.append((a0s, a1s, vd)); continue
    def mon(nm):
        m = re.search(r'J%s\[1\]=(.*)' % nm, out)
        p_ = sp.Poly(sp.sympify(m.group(1).replace('^', '**').replace(nm, 'x'), locals={'x': x}), x)
        return sp.Poly(p_.as_expr() / p_.LC(), x)
    pe, pf, pi, pj, pk = mon('ee'), mon('ff'), mon('ii'), mon('jj'), mon('kk')
    a0 = sp.Rational(a0s); a1 = sp.Rational(a1s)
    X = a0**2 / (1 - b); Y = a1**2 / b - rho2
    vacI1 = X + Y - (Nv + 1) / 12            # monic vacuum I_1 = P^2 + Q^2 - (N+1)/12 (bare)
    lvl = pi.degree() == 1 and sp.expand(-pi.all_coeffs()[1] - vacI1 - 2) == 0
    jz = pj.as_expr() == x; kz = pk.as_expr() == x
    if not (lvl and jz and kz and pe.degree() == 2 and pf.degree() == 2):
        bad.append((a0s, a1s, str(pi.as_expr()), str(pj.as_expr()), str(pk.as_expr()), pe.degree(), pf.degree()))
    pts.append((X, Y, -pe.all_coeffs()[1], pe.all_coeffs()[2], -pf.all_coeffs()[1], pf.all_coeffs()[2]))
print('t = %s: %d points usable, %d skipped (vdim != 4): %s; point checks failed: %s' % (tstr, len(pts), len(skipped), skipped, bad))
def fit(idx, deg, nfit):
    mons = [(i, j) for i in range(deg + 1) for j in range(deg + 1 - i)]
    assert nfit >= len(mons) and len(pts) > nfit
    A = sp.Matrix([[X**i * Y**j for (i, j) in mons] for (X, Y, *_) in pts[:nfit]]); rhs = sp.Matrix([p[idx] for p in pts[:nfit]])
    sol = (A.T * A).LUsolve(A.T * rhs)
    exact = all(r_ == 0 for r_ in (A * sol - rhs))
    poly = sp.expand(sum(c * P**(2 * i) * Q**(2 * j) for c, (i, j) in zip(sol, mons)))
    nbad = sum(1 for (X, Y, *rest) in pts[nfit:] if sum(c * X**i * Y**j for c, (i, j) in zip(sol, mons)) != (X, Y, *rest)[idx])
    return poly, exact, len(pts) - nfit, nbad
out = {'t': tstr, 'normalisation': 'eigenvalues / (P^(2K) coefficient of the vacuum eigenvalue); bare P^2 = P_X^2, Q^2 = pi^2 - rho^2', 'points_used': len(pts), 'skipped': skipped}
ok = not bad
for name, idx, deg, nfit in (('I3_trace', 2, 2, 12), ('I3_det', 3, 4, 24), ('I5_trace', 4, 3, 18), ('I5_det', 5, 6, 30)):
    poly, exact, nchk, nbad = fit(idx, deg, nfit)
    print('%s: degree <= %d in (P^2,Q^2); fit on %d points exact %s; held-back %d, mismatches %d' % (name, deg, nfit, exact, nchk, nbad)); print('    ', poly)
    ok &= exact and nbad == 0; out[name] = str(poly)
dd = json.load(open('d1_sol2_t%s.json' % tt))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
Mn = Mx / sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q).coeff_monomial(P**4)
tr_ok = sp.expand(Mn.trace() - sp.sympify(out['I3_trace'], locals={'P': P, 'Q': Q})) == 0
det_ok = sp.expand(Mn.det() - sp.sympify(out['I3_det'], locals={'P': P, 'Q': Q})) == 0
print('I_3 level-one block (data): trace identical %s, determinant identical %s' % (tr_ok, det_ok)); ok &= tr_ok and det_ok
json.dump(out, open('u_interp_t%s.json' % tt, 'w'), indent=1)
print('u_interp t = %s: %s' % (tstr, 'PASS' if ok else 'FAIL')); sys.exit(0 if ok else 2)
