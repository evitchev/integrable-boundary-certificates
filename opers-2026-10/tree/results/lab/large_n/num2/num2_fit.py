"""NUM2 blind fits (no WKB input in the fit) and comparisons.  For each data file: log|Q_j|(E) fitted per j by least squares (3x working precision).
t = 2:   S (sealed basis): E^(1/2), log E, 1, E^(-m/2), m = 1..K2-1;   R (registered refinement): S + E^(-3/2) log E, E^(-3) log E, E^(-9/2) log E.
t = 9/4: F1: E^(p_0), log E, 1, E^(p_K) K = 2..K2;   F2: F1 + E^(-13/10).
Comparisons: grades K = 0, 2..8 vs -int s_K (wkb_pred.json); odd K vs 0; the log E coefficients vs PREDICTION_NUM2_logE.json; Weyl-pair sums (P2); tamper l0 + 1/10."""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
from fractions import Fraction as Fr
tag = sys.argv[1]
D = json.load(open(os.path.join(HERE, f'data_{tag}.json')))
dps = D['dps']; mp.mp.dps = dps
LOG = open(os.path.join(HERE, f'run_fit_{tag}.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
k = Fr(D['k']); n = k + 4; grid = [mp.mpf(e) for e in D['grid']]; NP = len(grid)
Y = [[mp.mpf(v[j]) for v in D['logabsQ']] for j in range(4)]
P = json.load(open(os.path.join(HERE, 'wkb_pred.json')))[str(k)][D['point']]
LP = json.load(open(os.path.join(HERE, 'PREDICTION_NUM2_logE.json'))).get(D['point'], {}) if k == 2 else {}
pK = lambda K: (k + 1) * (1 - K) / n
def fpow(q): return mp.mpf(q.numerator) / q.denominator
def fit(vals, K2, extra):
    with mp.workdps(3 * dps):
        cols = [('K0', lambda e: e ** fpow(pK(0))), ('log', lambda e: mp.log(e)), ('const', lambda e: mp.mpf(1))]
        for K in range(2, K2 + 1): cols.append((f'K{K}', (lambda q: (lambda e: e ** fpow(q)))(pK(K))))
        cols += extra
        A = mp.matrix([[f(e) for _, f in cols] for e in grid]); x, res = mp.qr_solve(A, mp.matrix(vals))
        return {nm: x[i] for i, (nm, _) in enumerate(cols)}, res
def digits(a, b):
    d = abs(a - b) / (abs(b) if b != 0 else 1); return mp.inf if d == 0 else -mp.log10(d)
if k == 2:
    logx = [(f'logE_{q}', (lambda q: (lambda e: e ** fpow(q) * mp.log(e)))(q)) for q in (Fr(-3, 2), Fr(-3), Fr(-9, 2))]
    VARIANTS = [('S (sealed basis, no log E terms)', []), ('R (registered refinement, + E^p log E at p = -3/2, -3, -9/2)', logx)]
else:
    VARIANTS = [('F1 (WKB lattice only)', []), ('F2 (+ E^(-13/10))', [('X1310', lambda e: e ** (-mp.mpf(13) / 10))])]
names = [f'theta = {t}' for t in D['thetas']]
RES = {}
for vname, extra in VARIANTS:
    say(f'\n=== {tag}: {vname}')
    RES[vname] = {}
    for j in range(4):
        f1, r1 = fit(Y[j], 26, extra); f2, r2 = fit(Y[j], 22, extra)
        rows = []
        for K in [0] + list(range(2, 9)):
            pv = P['pred'][str(K)]
            ext = f1[f'K{K}']; stab = digits(f1[f'K{K}'], f2[f'K{K}'])
            if pv[0] == 'pole': rows.append(f'K{K} (pole in WKB): fit {mp.nstr(ext, 12)} [stab {mp.nstr(stab, 3)}]'); continue
            w = mp.mpf(pv[1])
            if w == 0: rows.append(f'K{K}: fit {mp.nstr(ext, 3)} vs 0 [stab-abs {mp.nstr(abs(f1[f"K{K}"] - f2[f"K{K}"]), 2)}]')
            else: rows.append(f'K{K}: {mp.nstr(digits(ext, w), 4)} digits [stab {mp.nstr(stab, 3)}]')
        say(f'  Q_j {names[j]}: lsq residual {mp.nstr(r1, 3)}; ' + '; '.join(rows))
        for nm_, _ in extra:
            v = f1[nm_]; st = digits(f1[nm_], f2[nm_])
            msg = f'     extra {nm_}: {mp.nstr(v, 22)} [stab {mp.nstr(st, 3)}]'
            if nm_.startswith('logE_') and LP:
                Kmap = {'-3/2': '4', '-3': '7', '-9/2': '10'}[nm_[5:]]
                pr = mp.mpf(LP[Kmap]['logE_coefficient_pred'])
                msg += f'  vs PREDICTED {mp.nstr(pr, 22)}: ' + (f'{mp.nstr(digits(v, pr), 4)} digits' if pr != 0 else f'pred 0, |fit| {mp.nstr(abs(v), 3)}')
            say(msg)
        RES[vname][j] = {kk: str(vv) for kk, vv in f1.items()}
    # Weyl-pair sums: log|Q_(s-l0)| + log|Q_(s+l0)| and l1 pair
    for pair, (a, b) in (('l0 pair', (0, 1)), ('l1 pair', (2, 3))):
        fp, _ = fit([Y[a][i] + Y[b][i] for i in range(NP)], 26, extra); fq, _ = fit([Y[a][i] + Y[b][i] for i in range(NP)], 22, extra)
        msg = []
        for K in [0] + list(range(2, 9)):
            pv = P['pred'][str(K)]
            if pv[0] == 'pole': msg.append(f'K{K}: fit/2 {mp.nstr(fp[f"K{K}"] / 2, 14)}'); continue
            w = mp.mpf(pv[1])
            msg.append(f'K{K}: ' + (f'{mp.nstr(digits(fp[f"K{K}"] / 2, w), 4)} d' if w != 0 else f'|fit| {mp.nstr(abs(fp[f"K{K}"]), 2)}'))
        say(f'  Weyl {pair} (sum/2 vs WKB): ' + '; '.join(msg))
        for nm_, _ in extra: say(f'     {pair} sum: extra {nm_} = {mp.nstr(fp[nm_], 18)} [stab {mp.nstr(digits(fp[nm_], fq[nm_]), 3)}]')
# tamper (a): WKB with l0 + 1/10 vs the extraction (first variant, j = dominant)
f1, _ = fit(Y[0], 26, VARIANTS[-1][1])
msg = []
for K in (2, 6, 8):
    tv = P['tamper_l0'][str(K)]
    if tv is None: continue
    msg.append(f'K{K}: {mp.nstr(digits(f1[f"K{K}"], mp.mpf(tv)), 3)} digits')
say(f'\nTAMPER (a) WKB l0 + 1/10 vs extraction (Q_0): ' + '; '.join(msg) + ' (must be << sealed)')
json.dump(RES, open(os.path.join(HERE, f'fit_{tag}.json'), 'w'), indent=1)
