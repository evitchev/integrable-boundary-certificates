"""EXC1 O3: the registered I_5 level-1 prediction (PREDICTION_EXC1.json) against the data block computed with the record's
tool (d1_sol1_t*_with5.json).  To be run after the lead's registration.  Modes: real | tamper (data trace + 1/1000).
Exit real: 0 = equal at both fibres, 2 = not; tamper: 0 = fires."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
PRED = open('PREDICTION_EXC1.sha256').read().split()[0]
assert hashlib.sha256(open('PREDICTION_EXC1.json', 'rb').read()).hexdigest() == PRED and PRED.startswith('de8366b1')
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
e, PX, PI, P, Q = sp.symbols('e PX PI P Q')
pred = json.load(open('PREDICTION_EXC1.json'))['cells']
bad = 0
for tstr in ('2', '9/4'):
    t = sp.Rational(tstr); rho2 = 2 / (t * t - 1)
    dd = json.load(open('d1_sol1_t%s_with5.json' % tstr.replace('/', '_')))
    assert dd['commute'] is True
    b = dd['blocks']['I5']
    Mx = sp.Matrix(2, 2, [sp.sympify(b['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
    vac = sp.Poly(sp.sympify(b['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
    lead = vac.coeff_monomial(P**6)
    cp = sp.expand((Mx / lead - e * sp.eye(2)).det()).subs(P, PX).subs(Q**2, PI**2 - rho2)
    cp = sp.Poly(sp.expand(cp), e)
    if mode == 'tamper':
        cp = sp.Poly(cp.as_expr() + e / 1000, e)
    pr = sp.Poly(sp.sympify(pred['t=' + tstr]['charpoly'], locals={'e': e, 'PX': PX, 'PI': PI}), e)
    same = sp.expand(cp.as_expr() - pr.as_expr()) == 0
    bad += not same
    print('t = %s: I_5 level-1 characteristic polynomial, data == registered oper prediction: %s  (transverse = level shift: %s; [I_3, I_5] = 0: %s)'
          % (tstr, same, b['transverse_is_level_shift'], dd['commute']))
print('O3 (%s):' % mode, 'PASS' if bad == 0 else 'FAIL')
sys.exit((0 if bad == 0 else 2) if mode == 'real' else (0 if bad else 1))
