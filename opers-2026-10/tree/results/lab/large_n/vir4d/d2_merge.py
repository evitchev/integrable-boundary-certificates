"""VIR4d: merge the per-fibre prediction files (same code, run in parallel) into
PREDICTION_VIR4d.json.  Opens no data file."""
import json, hashlib
out = None
for k in (4, 6, 8, 3):
    d = json.load(open('prediction_partial_%d.json' % k))
    if out is None:
        out = {'seal': d['seal'], 'normalisation': d['normalisation'], 'cells': {}}
    assert d['seal'] == out['seal']
    out['cells'].update(d['cells'])
assert sorted(out['cells']) == sorted('k=%d,spin=%d' % (k, s) for k in (3, 4, 6, 8) for s in (11, 13))
for name, c in out['cells'].items():
    assert c['coefficients'] is not None and c['odd_orders_R11_R13_vanish'] and c['symmetric_PX_pi'], name
    print(name, 't =', c['t'], 'order', c['order'], 'monomials', len(c['coefficients']), 'leading', c['leading_coefficient'])
json.dump(out, open('PREDICTION_VIR4d.json', 'w'), indent=1, sort_keys=True)
print('PREDICTION_VIR4d.json sha256', hashlib.sha256(open('PREDICTION_VIR4d.json', 'rb').read()).hexdigest())
