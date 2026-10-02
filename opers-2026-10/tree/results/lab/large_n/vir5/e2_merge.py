"""VIR5a: merge the per-fibre prediction files into PREDICTION_VIR5a.json.  Opens no data file."""
import json, hashlib
tags = ['10_3', '5_3', '2_3', '1_2', 'm6', 'm18_5', 'm6_5', 'm2_5', '4']
out = {'cells': {}}
for tg in tags:
    d = json.load(open('prediction_k%s.json' % tg))
    assert d['mode'] == 'real'
    out['seal'] = d['seal']
    out['cells']['k=%s' % d['k']] = d
    sp_ = d['spins']
    print('k = %-6s t = %-6s n = %-5s spins with a charge: %s ; without: %s ; odd orders vanish: %s'
          % (d['k'], d['t'], d['n'], [s for s in sorted(sp_, key=int) if sp_[s]['coefficients'] is not None],
             [s for s in sorted(sp_, key=int) if sp_[s]['coefficients'] is None], all(d['odd_orders_vanish'].values())))
json.dump(out, open('PREDICTION_VIR5a.json', 'w'), indent=1, sort_keys=True)
print('PREDICTION_VIR5a.json sha256', hashlib.sha256(open('PREDICTION_VIR5a.json', 'rb').read()).hexdigest())
