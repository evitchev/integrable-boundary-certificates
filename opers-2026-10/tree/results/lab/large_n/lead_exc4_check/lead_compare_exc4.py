"""Lead's blind comparison for EXC4 (B): Fable's registered PREDICTION_EXC4.json (sha256 cff1aead..., snapshot 11:20:35Z) vs the lead's own
level-two I_3 blocks (EXC1's d2_level2.py, record's excited_states tool), Solution 3.  Usage: lead_compare_exc4.py [tamper]"""
import json, sys, hashlib, sympy as sp
e = sp.Symbol('e'); tamper = len(sys.argv) > 1
pp = 'PREDICTION_EXC4.json'  # KIT PATCH: the lead's snapshot next to this script (was <home>/ib-lab/prereg_snapshots/exc4/PREDICTION_EXC4.json, outside the archive); same sha256 cff1aeadc77aab1d...
assert hashlib.sha256(open(pp, 'rb').read()).hexdigest().startswith('cff1aeadc77aab1d')
pred = json.load(open(pp)); bad = 0
for pt in pred['points']:
    tag = 't%s_P%s_Q%s' % tuple(pt[x].replace('/', '_') for x in ('t', 'P', 'Q'))
    d = json.load(open('d2_sol3_%s.json' % tag))
    quint = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in d['factors'] if sp.Poly(sp.sympify(f, locals={'e': e}), e).degree() == 5]
    assert len(quint) == 1
    q = sp.Poly(quint[0].as_expr() / quint[0].LC(), e)
    pq = sp.Poly(sp.sympify(pt['I3_quintic_monic'], locals={'e': e}), e)
    if tamper: pq = sp.Poly(pq.as_expr() + sp.Rational(1, 1000), e)
    same = sp.expand(q.as_expr() - pq.as_expr()) == 0
    print('t = %s (P,Q) = (%s,%s): monic quintic identical: %s' % (pt['t'], pt['P'], pt['Q'], same), flush=True)
    bad += not same
print('RESULT', 'MATCH' if bad == 0 else 'MISMATCH'); sys.exit(0 if bad == 0 else 2)
