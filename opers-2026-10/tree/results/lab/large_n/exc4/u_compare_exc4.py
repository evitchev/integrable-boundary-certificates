"""EXC4 B3 (hold-out): the registered level-two I_3 quintics against the blocks computed with the record's tool (d2_sol3_t<t>_P<P>_Q<Q>.json,
produced after the lead's registration and comparison).  Modes: real | tamper (quintic constant term + 1/1000 on the data side; must FAIL).
Exit: real 0 iff all three monic quintics and the level-shift roots agree; tamper 0 iff it fires."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('PREDICTION_EXC4.json', 'rb').read()).hexdigest() == open('PREDICTION_EXC4.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
pred = json.load(open('PREDICTION_EXC4.json'))
e = sp.Symbol('e'); ok = True
for p in pred['points']:
    d2 = json.load(open('d2_sol3_t%s_P%s_Q%s.json' % (p['t'].replace('/', '_'), p['P'].replace('/', '_'), p['Q'].replace('/', '_'))))
    facs = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in d2['factors']]
    quint = [f for f in facs if f.degree() == 5][0]; quint = sp.Poly(sp.expand(quint.as_expr() / quint.LC()), e)
    if mode == 'tamper':
        quint = sp.Poly(quint.as_expr() + sp.Rational(1, 1000), e)
    lin = [f for f in facs if f.degree() == 1]
    pq = sp.Poly(sp.sympify(p['I3_quintic_monic'], locals={'e': e}), e)
    diff = sp.Poly(sp.expand(pq.as_expr() - quint.as_expr()), e)
    nd = len([c for c in diff.all_coeffs() if c != 0])
    print('t = %-4s (P,Q) = (%s,%s) [%s]: data factors degrees %s; quintic coefficients differing: %d of 6; level-shift root present: %s' % (p['t'], p['P'], p['Q'], mode, [f.degree() for f in facs], nd, len(lin) == 1))
    ok &= nd == 0 and len(lin) == 1
print('U3 level-two hold-out [%s]: %s' % (mode, 'PASS' if ok else 'FAIL'))
sys.exit((0 if ok else 2) if mode == 'real' else (0 if not ok else 1))
