"""Kit verifier (v1.2.5, new; not an archived script) for NUM3's sealed T1 test (t1_zeros3.py -> t1_3.json), Solution 2 at t = 2.
Run with cwd = results/lab/large_n/num3.  (1) Recomputes the exact E = 0 constants C_j = Q3_j(0)/Q7_j(0) (Meijer G, ops3.exact_E0, as in t1_zeros3.py)
and, from the ARCHIVED Q3/Q7 values in t1_3.json, dev_j(E) = |Q3_j/Q7_map(j)/C_j - 1| at the 2 points x 5 energies x 3 j; each must equal the archived
dev to 1e-57 absolute (the archived Q values are 60-digit strings), and the worst must be < 1e-25 (sealed).  Prints the minimum and median digits (paper: 'at least 36 digits (typically 55)').
(2) Spot check: recomputes Q3 and Q7 at point B, E = 2 (the worst cell) from both ODEs with the archived code path and compares with t1_3.json."""
import sys, os, json, time
sys.path.insert(0, os.getcwd())
import mpmath as mp
from fractions import Fraction as Fr
from ops3 import sol2_three_term, seventh, exact_E0
mp.mp.dps = 60
T = json.load(open('t1_3.json'))
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
MAP = {0: 0, 1: 5, 2: 6}
ok = True; devs = []
for nm, (l0, l1) in PTS.items():
    op0, th0, F, rR = sol2_three_term(2, l0, l1, 0); e3, _ = exact_E0(th0, F['c'], mp.mpf(1), rR)
    op70, ex7 = seventh(l0, l1, 0); e7, _ = exact_E0(ex7, Fr(21, 2), mp.mpf(1), [])
    C = [e3[j] / e7[MAP[j]] for j in range(3)]
    for Es in ('-20.0', '-5.0', '2.0', '50.0', '500.0'):
        r = T[nm][Es]; Q3 = [mp.mpf(q) for q in r['Q3']]; Q7 = [mp.mpf(q) for q in r['Q7']]
        assert len(Q3) == 3 and len(Q7) == 7
        dev = [abs(Q3[j] / Q7[MAP[j]] / C[j] - 1) for j in range(3)]
        agree = all(abs(dev[j] - mp.mpf(r['dev'][j])) < mp.mpf(10) ** -57 for j in range(3))   # the archived Q's are 60-significant-digit strings
        ok &= agree; devs += dev
        print(f'point {nm} E = {Es:>6}: digits of Q3_j/Q7_j/C_j - 1 = {[mp.nstr(-mp.log10(d), 4) for d in dev]}; recomputed dev == archived dev: {agree}')
devs.sort(key=lambda d: -d); worst = devs[0]; med = devs[len(devs) // 2]
print(f'{len(devs)} cells: worst {mp.nstr(worst, 3)} = {mp.nstr(-mp.log10(worst), 4)} digits; median {mp.nstr(-mp.log10(med), 3)} digits; sealed >= 25 digits: {worst < mp.mpf(10) ** -25}')
ok &= worst < mp.mpf(10) ** -25
# (2) spot recomputation of one energy (point B, E = 2), exactly as t1_zeros3.py
l0, l1 = PTS['B']; E = mp.mpf(2); x0 = 2 / (abs(E) + 1) ** (mp.mpf(1) / 7); t0 = time.time()
op3, th3, F, rR = sol2_three_term(2, l0, l1, E); lq3, Q3 = op3.Q(th3, x0)
op7, ex7 = seventh(l0, l1, -E); lq7, Q7 = op7.Q(ex7, x0)
r = T['B']['2.0']
d3 = max(abs(Q3[j] - mp.mpf(r['Q3'][j])) / abs(mp.mpf(r['Q3'][j])) for j in range(3)); d7 = max(abs(Q7[i] - mp.mpf(r['Q7'][i])) / abs(mp.mpf(r['Q7'][i])) for i in range(7))
sp_ok = d3 < mp.mpf(10) ** -50 and d7 < mp.mpf(10) ** -50; ok &= sp_ok
print(f'spot recomputation point B, E = 2 (both ODEs): Q3 and Q7 equal the archived values to 1e-50: {sp_ok} ({time.time() - t0:.0f}s)')
print('T1 VERIFY: ' + ('PASS' if ok else 'FAIL')); sys.exit(0 if ok else 1)
