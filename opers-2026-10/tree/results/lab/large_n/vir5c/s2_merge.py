"""VIR5c S2: merge the per-fibre engine predictions -> PREDICTION_VIR5c.json, and the engine control:
spins <= 13 must reproduce the registered PREDICTION_VIR5a.json at the shared fibres.  Opens no data table."""
import hashlib, json, sys
sys.dont_write_bytecode = True
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
old_raw = open('<home>/fable-work/vir5/PREDICTION_VIR5a.json', 'rb').read()
assert hashlib.sha256(old_raw).hexdigest().startswith('0c12a7e9')
old = json.loads(old_raw)['cells']
files = {'10/3': 'pred_k10_3_o20.json', '2/3': 'pred_k2_3_o20.json', '-6': 'pred_km6_o20.json', '2': 'pred_k2_o20.json',
         '5/3': 'pred_k5_3_o16.json', '-18/5': 'pred_km18_5_o16.json', '-6/5': 'pred_km6_5_o16.json',
         '-2/5': 'pred_km2_5_o16.json', '4': 'pred_k4_o16.json'}
out = {'seal': SEAL, 'cells': {}}
ok = True
for k, fn in files.items():
    cell = json.load(open(fn))
    assert cell['k'] == k and cell['mode'] == 'real' and cell['seal'] == SEAL
    top = max(int(s) for s in cell['spins'])
    assert top == cell['order'] - 1, (k, top, cell['order'])
    out['cells']['k=' + k] = cell
    if 'k=' + k in old:
        same = all(old['k=' + k]['spins'][str(s)]['coefficients'] == cell['spins'][str(s)]['coefficients'] for s in range(1, 14, 2))
        ok &= same
        print('k = %-5s t = %-6s spins up to %d; spins 1..13 identical to the registered VIR5a prediction: %s' % (k, cell['t'], top, same))
    else:
        print('k = %-5s t = %-6s spins up to %d; (fibre not in VIR5a)' % (k, cell['t'], top))
json.dump(out, open('PREDICTION_VIR5c.json', 'w'), indent=1, sort_keys=True)
print('PREDICTION_VIR5c.json sha256', hashlib.sha256(open('PREDICTION_VIR5c.json', 'rb').read()).hexdigest())
print('engine control', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
