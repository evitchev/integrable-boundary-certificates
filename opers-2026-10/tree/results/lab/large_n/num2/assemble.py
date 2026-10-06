"""NUM2: assemble points/pt_<K>_<POINT>_<DPS>_<i>.json (i = 0..35) into data_<K>_<POINT>_<DPS>.json (the format num2_fit.py reads)."""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1]
fs = {int(f.rsplit('_', 1)[1][:-5]): f for f in glob.glob(os.path.join(HERE, 'points', f'pt_{tag}_*.json'))}
missing = [i for i in range(36) if i not in fs]
if missing: print(f'{tag}: missing {missing}'); sys.exit(1)
P = [json.load(open(fs[i])) for i in range(36)]
out = {'k': P[0]['k'], 'point': P[0]['point'], 'dps': P[0]['dps'], 'wrong_order': False, 'grid': [p['E'] for p in P],
       'logabsQ': [p['logabsQ'] for p in P], 'sign': [p['sign'] for p in P], 'thetas': P[0]['thetas'], 'seconds': [p['seconds'] for p in P]}
json.dump(out, open(os.path.join(HERE, f'data_{tag}.json'), 'w'), indent=1)
print(f'{tag}: assembled 36 energies; total CPU {sum(out["seconds"]):.0f}s')
