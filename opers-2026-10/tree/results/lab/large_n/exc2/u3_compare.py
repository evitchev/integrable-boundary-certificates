"""EXC2 U3 (hold-out): the registered I_5 level-1 prediction against the block computed with the record's tool.
Modes: real | tamper (one data matrix entry + 1/1000, after the hash guard; must FAIL) | flip (oper-side: B1 sign flipped,
t = 2 only; must FAIL) | cross (each fibre's data against the OTHER fibre's prediction; must FAIL).
Exit: real 0 iff trace and det identical at both fibres and all 72 pointwise characteristic polynomials agree; controls 0 iff they fire."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('PREDICTION_EXC2.json', 'rb').read()).hexdigest() == open('PREDICTION_EXC2.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
pred = json.load(open('PREDICTION_EXC2.json'))
P, Q, lam = sp.symbols('P Q lam')
loc = {'P': P, 'Q': Q, 'lam': lam}
def data(tstr):
    dd = json.load(open('d1_sol3_t%s_with5.json' % tstr.replace('/', '_')))
    assert dd['commute'] is True
    blk = dd['blocks']['I5']
    Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][i][j], locals=loc) for i in range(2) for j in range(2)])
    if mode == 'tamper':
        Mx[0, 0] += sp.Rational(1, 1000)
    c6 = sp.Poly(sp.sympify(blk['vacuum'], locals=loc), P, Q).coeff_monomial(P**6)
    Mn = Mx / c6
    return sp.expand(Mn.trace()), sp.expand(Mn.det()), c6
fib = ['2', '9/4']
ok = True
for tstr in fib:
    if mode == 'flip' and tstr != '2':
        continue
    src = tstr if mode != 'cross' else [x for x in fib if x != tstr][0]
    if mode == 'flip':
        pr = json.load(open('u3_interp_t2_flip.json'))
    else:
        pr = pred['fibres'][src]
    tr, det, c6 = data(tstr)
    dtr = sp.expand(tr - sp.sympify(pr['I5_trace'], locals=loc)); ddet = sp.expand(det - sp.sympify(pr['I5_det'], locals=loc))
    ntr = len(sp.Poly(tr, P, Q).terms()); ndet = len(sp.Poly(det, P, Q).terms())
    mtr = len(sp.Poly(dtr, P, Q).terms()) if dtr != 0 else 0; mdet = len(sp.Poly(ddet, P, Q).terms()) if ddet != 0 else 0
    # pointwise characteristic polynomials (registered) against the data block
    npt = bad = 0
    if mode in ('real', 'tamper'):
        t = sp.Rational(tstr); rho2 = 2 / (t * t - 1); ff = sp.Symbol('ff')
        for p in pred['fibres'][tstr]['points']:
            X = sp.Rational(p['PX2']); Y = sp.Rational(p['PI2']) - rho2
            ev = lambda ex: sp.Poly(ex, P, Q).eval({P: sp.sqrt(X), Q: sp.sqrt(Y)}) if False else sum(c * X**(m[0] // 2) * Y**(m[1] // 2) for m, c in sp.Poly(ex, P, Q).terms())
            cp = sp.expand(ff**2 - ev(tr) * ff + ev(det))
            npt += 1; bad += sp.expand(cp - sp.sympify(p['I5_charpoly'], locals={'ff': ff})) != 0
    print('t = %-4s [%s] data vacuum P^6 coefficient %s; trace: %d coefficients, %d differ; det: %d coefficients, %d differ; pointwise char. polynomials: %d compared, %d differ'
          % (tstr, mode + ('' if mode != 'cross' else ' vs prediction for t = ' + src), c6, ntr, mtr, ndet, mdet, npt, bad))
    ok &= (dtr == 0 and ddet == 0 and bad == 0)
print('U3 hold-out [%s]: %s' % (mode, 'PASS' if ok else 'FAIL'))
if mode == 'real':
    sys.exit(0 if ok else 2)
sys.exit(0 if not ok else 1)
