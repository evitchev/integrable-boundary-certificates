"""NUM3: points/pt_<tag>_<i>.json (i = 0..35) -> data_<tag>.json"""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); tag = sys.argv[1]
fs = {}
for f in glob.glob(os.path.join(HERE, 'points', f'pt_{tag.replace("_wrong", "")}_*.json')):
    b = os.path.basename(f)[3:-5]; w = b.endswith('_wrong'); b = b[:-6] if w else b
    if w != tag.endswith('_wrong'): continue
    fs[int(b.rsplit('_', 1)[1])] = f
miss = [i for i in range(36) if i not in fs]
if miss: print(f'{tag}: missing {miss}'); sys.exit(1)
P = [json.load(open(fs[i])) for i in range(36)]
json.dump({'k': P[0]['k'], 'point': P[0]['point'], 'dps': P[0]['dps'], 'wrong_order': P[0]['wrong_order'], 'grid': [p['E'] for p in P],
           'logabsQ': [p['logabsQ'] for p in P], 'sign': [p['sign'] for p in P], 'thetas': P[0]['thetas'], 'seconds': [p['seconds'] for p in P]},
          open(os.path.join(HERE, f'data_{tag}.json'), 'w'), indent=1)
print(f'{tag}: assembled; CPU {sum(p["seconds"] for p in P):.0f}s')
