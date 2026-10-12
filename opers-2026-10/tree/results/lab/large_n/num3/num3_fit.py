"""NUM3 blind fits (no WKB input) and comparisons.  Sealed bases:
t = 2:   SEALED = lattice E^(3(1-K)/7) (K = 0, 2..K2), log E, 1, E^(-3/2), E^(-9/2), E^(-15/2), E^(-3) log E, E^(-6) log E;   BARE = lattice only.
t = 9/4: SEALED = lattice E^(13(1-K)/29), log E, 1, E^(-13m/10) (m = 1, 2, 3);   BARE = lattice only.
Comparisons: grades K = 0, 2..8 vs -int s_K (wkb_pred3.json); odd K vs 0; log E coefficients vs PREDICTION_NUM3_logE.json; tamper l0 + 1/10."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
from fractions import Fraction as Fr
tag = sys.argv[1]
D = json.load(open(os.path.join(HERE, f'data_{tag}.json'))); dps = D['dps']; mp.mp.dps = dps
LOG = open(os.path.join(HERE, f'run_fit_{tag}.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
k = Fr(D['k']); n = 2 * k + 3; grid = [mp.mpf(e) for e in D['grid']]; NP = len(grid); NJ = len(D['thetas'])
Y = [[mp.mpf(v[j]) for v in D['logabsQ']] for j in range(NJ)]
P = json.load(open(os.path.join(HERE, 'wkb_pred3.json')))[str(k)][D['point']]
LP = json.load(open(os.path.join(HERE, 'PREDICTION_NUM3_logE.json'))).get(D['point'], {}) if k == 2 else {}
pK = lambda K: (k + 1) * (1 - K) / n; fp = lambda q: mp.mpf(q.numerator) / q.denominator
def fit(vals, K2, extra):
    with mp.workdps(3 * dps):
        cols = [('K0', (fp(pK(0)), False)), ('log', None), ('const', (mp.mpf(0), False))] + [(f'K{K}', (fp(pK(K)), False)) for K in range(2, K2 + 1)] + extra
        def val(e, spec):
            if spec is None: return mp.log(e)
            p, lg = spec; return e ** p * (mp.log(e) if lg else 1)
        A = mp.matrix([[val(e, sp_) for _, sp_ in cols] for e in grid]); x, res = mp.qr_solve(A, mp.matrix(vals))
    return {nm: x[i] for i, (nm, _) in enumerate(cols)}, res
def digits(a, b):
    d = abs(a - b) / (abs(b) if b != 0 else 1); return mp.inf if d == 0 else -mp.log10(d)
if k == 2:
    SEALED = [('N-3/2', (-mp.mpf(3) / 2, False)), ('N-9/2', (-mp.mpf(9) / 2, False)), ('N-15/2', (-mp.mpf(15) / 2, False)),
              ('logE_-3', (mp.mpf(-3), True)), ('logE_-6', (mp.mpf(-6), True))]
else:
    SEALED = [(f'N{m}', (-mp.mpf(13) * m / 10, False)) for m in (1, 2, 3)]
VARIANTS = [('BARE (lattice only)', []), ('SEALED basis', SEALED)]
RES = {}
for vname, extra in VARIANTS:
    say(f'\n=== {tag}: {vname}'); RES[vname] = {}
    for j in range(NJ):
        f1, r1 = fit(Y[j], 28, extra); f2, r2 = fit(Y[j], 24, extra)
        row = []
        for K in [0] + list(range(2, 9)):
            pv = P['pred'][str(K)]
            if pv[0] == 'pole': row.append(f'K{K} (pole): fit {mp.nstr(f1[f"K{K}"], 12)} [stab {mp.nstr(digits(f1[f"K{K}"], f2[f"K{K}"]), 3)}]'); continue
            if pv[0] == 'zero': w = mp.mpf(0)
            else: w = mp.mpf(pv[1])
            if w == 0: row.append(f'K{K}: |fit| {mp.nstr(abs(f1[f"K{K}"]), 2)} vs 0 [stab-abs {mp.nstr(abs(f1[f"K{K}"] - f2[f"K{K}"]), 2)}]')
            else: row.append(f'K{K}: {mp.nstr(digits(f1[f"K{K}"], w), 4)} digits [stab {mp.nstr(digits(f1[f"K{K}"], f2[f"K{K}"]), 3)}]')
        say(f'  Q (theta = {mp.nstr(mp.mpf(D["thetas"][j]), 8)}): lsq residual {mp.nstr(r1, 3)}; ' + '; '.join(row))
        for nm_, _ in extra:
            v = f1[nm_]; st = digits(f1[nm_], f2[nm_]); msg = f'     extra {nm_}: {mp.nstr(v, 22)} [stab {mp.nstr(st, 3)}]'
            if nm_.startswith('logE_') and LP:
                Kmap = {'logE_-3': '8', 'logE_-6': '15'}[nm_]
                if Kmap in LP:
                    pr = mp.mpf(LP[Kmap]['logE_coefficient_pred']); msg += f'  vs PREDICTED {mp.nstr(pr, 22)}: ' + (f'{mp.nstr(digits(v, pr), 4)} digits' if pr != 0 else f'pred 0, |fit| {mp.nstr(abs(v), 3)}')
            say(msg)
        RES[vname][j] = {kk: str(vv) for kk, vv in f1.items()}
f1, _ = fit(Y[0], 28, SEALED); msg = []
for K in (2, 4, 6):
    tv = P['tamper_l0'][str(K)]
    if tv: msg.append(f'K{K}: {mp.nstr(digits(f1[f"K{K}"], mp.mpf(tv)), 3)} digits')
say('\nTAMPER (a) WKB l0 + 1/10 vs extraction (Q_0): ' + '; '.join(msg) + ' (must be << sealed)')
json.dump(RES, open(os.path.join(HERE, f'fit_{tag}.json'), 'w'), indent=1)
