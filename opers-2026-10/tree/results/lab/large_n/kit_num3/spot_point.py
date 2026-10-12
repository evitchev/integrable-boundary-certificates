"""Kit verifier (v1.2.5, new): recompute ONE archived NUM3 data point with the archived script num3_point.py and compare with the archived file.
Run with cwd = results/lab/large_n/num3.  Usage: spot_point.py K POINT DPS INDEX.  The archived points/pt_<tag>.json is read first, num3_point.py is
then run unchanged (it rewrites that file in this temporary copy), and the recomputed log|Q_j| must equal the archived 60-digit values to 1e-50."""
import sys, os, json, runpy, time
import mpmath as mp
K, PT, DPS, I = sys.argv[1:5]
f = os.path.join('points', 'pt_%s_%s_%s_%s.json' % (K.replace('/', '_'), PT, DPS, I))
old = json.load(open(f)); t0 = time.time()
sys.argv = ['num3_point.py', K, PT, DPS, I]; sys.path.insert(0, os.getcwd())
runpy.run_path('num3_point.py', run_name='__main__')
new = json.load(open(f)); mp.mp.dps = int(DPS)
d = max(abs(mp.mpf(a) - mp.mpf(b)) for a, b in zip(old['logabsQ'], new['logabsQ']))
same = [old[k_] == new[k_] for k_ in ('k', 'point', 'dps', 'E', 'sign', 'thetas')]
ok = d < mp.mpf(10) ** -50 and all(same) and len(old['logabsQ']) == len(new['logabsQ']) == 3
print(f'archived point {f}: E = {old["E"]}; recomputed max |delta log|Q_j|| < 1e-50: {d < mp.mpf(10) ** -50}; E, signs, thetas identical: {all(same)} '
      f'(archived run {old["seconds"]:.0f} s, now {time.time() - t0:.0f} s)')
print('SPOT POINT: ' + ('PASS' if ok else 'FAIL')); sys.exit(0 if ok else 1)
