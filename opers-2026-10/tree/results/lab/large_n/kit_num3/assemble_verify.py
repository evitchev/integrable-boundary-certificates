"""Kit verifier (v1.2.5, new): every fitted NUM3 data file data_<tag>.json equals the assembly of its archived per-energy point files
(points/pt_<tag>_<i>.json, i = 0..35), field by field (assemble3.py's mapping).  Run with cwd = results/lab/large_n/num3."""
import json, sys, hashlib
ok = True
for tag in ('2_A_60', '2_B_60', '10_3_A_60', '10_3_B_60', '2_A_120', '2_A_60_wrong'):
    w = tag.endswith('_wrong'); base = tag[:-6] if w else tag
    P = [json.load(open('points/pt_%s_%d%s.json' % (base, i, '_wrong' if w else ''))) for i in range(36)]
    D = json.load(open('data_%s.json' % tag))
    exp = {'k': P[0]['k'], 'point': P[0]['point'], 'dps': P[0]['dps'], 'wrong_order': P[0]['wrong_order'], 'grid': [p['E'] for p in P],
           'logabsQ': [p['logabsQ'] for p in P], 'sign': [p['sign'] for p in P], 'thetas': P[0]['thetas'], 'seconds': [p['seconds'] for p in P]}
    same = D == exp and all(p['wrong_order'] == w and p['index'] == i for i, p in enumerate(P))
    ok &= same
    print('data_%s.json (sha256 %s...): 36 points, equal to the assembly of points/: %s' % (tag, hashlib.sha256(open('data_%s.json' % tag, 'rb').read()).hexdigest()[:16], same))
print('ASSEMBLY: ' + ('PASS' if ok else 'FAIL')); sys.exit(0 if ok else 1)
