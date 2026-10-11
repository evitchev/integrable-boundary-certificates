"""NUMTAB (sealed SEAL_NUMTAB.md 79f073bf...): the spectral-evidence table with stability tests, from ARCHIVED data only (no ODE integration).
Datasets D1..D5 with their base models of record; rows per fitted power with status, prediction, fit, agreement digits d; stability: windows W-/W+,
truncation K2 vs K2-4, precision 60 vs 120 digits (where archived); error estimate eps = max spread, digits e.  Writes table.json."""
import os, sys, json, hashlib, glob, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'inputs', 'num1'))
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUMTAB.md'), 'rb').read()).hexdigest().startswith('79f073bf')
for l in open(os.path.join(HERE, 'inputs', 'ORIGIN_SHA256')):
    h, f = l.split(); assert hashlib.sha256(open(os.path.join(HERE, 'inputs', f), 'rb').read()).hexdigest() == h, f
import mpmath as mp
from fractions import Fraction as Fr
LOG = open(os.path.join(HERE, 'run_numtab.log'), 'w'); T00 = time.time()
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
IN = lambda *p: os.path.join(HERE, 'inputs', *p)
fp = lambda q: mp.mpf(q.numerator) / q.denominator
# ---------------------------------------------------------------- fitting
def fit(grid, vals, cols, wdps):
    """cols: list of (name, spec); spec None = log E; (p, lg) = E^p (log E if lg)"""
    with mp.workdps(wdps):
        def v(e, s):
            if s is None: return mp.log(e)
            p, lg = s; return e ** p * (mp.log(e) if lg else 1)
        A = mp.matrix([[v(e, s) for _, s in cols] for e in grid]); x, res = mp.qr_solve(A, mp.matrix(vals))
        return {nm: +x[i] for i, (nm, _) in enumerate(cols)}
def model(pK, K2, extras):
    cols = [('K0', (pK(0), False)), ('log', None), ('const', (mp.mpf(0), False))] + [(f'K{K}', (pK(K), False)) for K in range(2, K2 + 1)]
    return cols + extras
def k2_window(NP, base, n_extra):
    return min(base, NP - 6 - n_extra)
def stability(grid, vals, pK, K2, extras, dps, prec=None):
    """returns base fit, and the spread dicts"""
    W = 3 * dps; NP = len(grid); q = NP // 4
    base = fit(grid, vals, model(pK, K2, extras), W)
    trunc = fit(grid, vals, model(pK, K2 - 4, extras), W)
    kw = k2_window(NP - q, K2, len(extras))
    wm = fit(grid[q:], vals[q:], model(pK, kw, extras), W)          # W-: drop lowest quarter
    wp = fit(grid[:NP - q], vals[:NP - q], model(pK, kw, extras), W)   # W+: drop highest quarter
    pr = None
    if prec is not None:
        g2, v2, g1, v1, kp = prec                                    # dps120 data and the dps60 data on the same window
        a = fit(g2, v2, model(pK, kp, extras), 3 * 120); b = fit(g1, v1, model(pK, kp, extras), W)
        pr = {nm: abs(a[nm] - b[nm]) for nm in a}
    sp = {}
    for nm in base:
        s = {'W-': abs(base[nm] - wm[nm]) if nm in wm else None, 'W+': abs(base[nm] - wp[nm]) if nm in wp else None,
             'trunc': abs(base[nm] - trunc[nm]) if nm in trunc else None, 'prec': pr.get(nm) if pr else None}
        sp[nm] = s
    return base, sp, kw
def digits(x, ref):
    if ref == 0: return None
    d = abs(x - ref) / abs(ref); return mp.inf if d == 0 else -mp.log10(d)
ROWS = []
def emit(ds, pt, nm, power, status, pred, val, sp, note=''):
    eps = max([v for v in sp.values() if v is not None])
    e = (-mp.log10(eps / abs(val)) if (eps > 0 and val != 0) else mp.inf)
    if pred is None: d = None; flag = 'no prediction'
    elif pred == 0: d = None; flag = 'consistent' if abs(val) <= 10 * eps else 'TENSION'
    else: d = digits(val, pred); flag = 'consistent' if abs(val - pred) <= 10 * eps else 'TENSION'
    r = dict(dataset=ds, point=pt, row=nm, power=power, status=status, pred=None if pred is None else mp.nstr(pred, 30), fit=mp.nstr(val, 30),
             d=None if d is None else float(min(d, 999)), abs_fit_if_pred0=mp.nstr(abs(val), 3) if pred == 0 else None,
             eps=mp.nstr(eps, 3), e=float(min(e, 999)), spreads={k: (None if v is None else mp.nstr(v, 3)) for k, v in sp.items()}, flag=flag, note=note)
    ROWS.append(r)
    say(f'  {ds} {pt} {nm:12s} {status:34s} pred {"-" if pred is None else mp.nstr(pred, 14):>22s} fit {mp.nstr(val, 14):>22s} '
        f'd {"-" if d is None else mp.nstr(d, 4):>6s} e {mp.nstr(e, 4):>6s} [{flag}]' + (f' |fit| {mp.nstr(abs(val), 2)}' if pred == 0 else ''))
# ---------------------------------------------------------------- D1: NUM1
mp.mp.dps = 60
import rwkb
say('== D1 NUM1 (Sol 1, M = 5)')
J60 = json.load(open(IN('num1', 'num1_dps60.json'))); J120 = json.load(open(IN('num1', 'num1_dps120.json'))); AV = json.load(open(IN('num1', 'avgpm.json')))
mu = Fr(3, 5); pK1 = lambda K: fp(mu * (1 - K))
S5 = rwkb.wkb(5, 8); eps_ = Fr(1, 10 ** 30); Sp = rwkb.wkb(Fr(5) + eps_, 6); Sm = rwkb.wkb(Fr(5) - eps_, 6)
PTS1 = {'A': (mp.mpf(1), mp.mpf(1) / 3), 'B': (mp.mpf(-3) / 2, mp.mpf(7) / 10), 'C': (mp.mpf(0), mp.mpf(1) / 3)}
def mellin(s, nu): return mp.gamma(s / 2) * mp.gamma(s / 2 + nu) * mp.gamma(mp.mpf(1) / 2 - s / 2) / (4 * mp.sqrt(mp.pi) * mp.gamma(1 + nu - s / 2))
grid1 = [mp.mpf(e) for e in J60['grid']]
for pt, (al, l) in PTS1.items():
    lam = l * (l + 1); nu = l + mp.mpf(1) / 2
    pred = {K: -rwkb.integral_coeff(S5[K], 5, al, lam)[0] for K in [0] + list(range(2, 9))}
    lim6 = -(rwkb.integral_coeff(Sp[6], Fr(5) + eps_, al, lam)[0] + rwkb.integral_coeff(Sm[6], Fr(5) - eps_, al, lam)[0]) / 2
    vals = [mp.mpf(v) for v in J60['data'][pt]]
    prec = None
    if pt == 'A':
        g2 = [mp.mpf(e) for e in J120['grid']]; v2 = [mp.mpf(v) for v in J120['data']['A']]; prec = (g2, v2, grid1, vals, 30)
    base, sp, kw = stability(grid1, vals, pK1, 30, [], 60, prec)
    for K in [0] + list(range(2, 9)):
        if K == 6 and al != 0:
            emit('D1', pt, 'e^-3 (K6)', '-3', 'non-WKB post hoc (lim WKB + a Mellin)', lim6 + al * mellin(mp.mpf(6), nu), base['K6'], sp['K6'],
                 'resonant grade; formula derived after the sealed failure; confirmed at M = 51/10 by a registered prediction')
        else:
            emit('D1', pt, f'K{K}', str(mu * (1 - K)), 'WKB predicted', pred[K], base[f'K{K}'], sp[f'K{K}'])
    emit('D1', pt, 'log e', 'log', 'fitted free', None, base['log'], sp['log'])
    emit('D1', pt, 'const', '0', 'fitted free', None, base['const'], sp['const'])
    if pt in AV:
        vav = [(mp.mpf(a) + mp.mpf(b)) / 2 for a, b in zip(J60['data'][pt], AV[pt]['minus_alpha_data'])]
        baseA, spA, _ = stability(grid1, vav, pK1, 30, [], 60)
        for K in [0] + list(range(2, 9)):
            pe = (pred[K] - rwkb.integral_coeff(S5[K], 5, -al, lam)[0]) / 2 if True else None
            pe = (-rwkb.integral_coeff(S5[K], 5, al, lam)[0] - rwkb.integral_coeff(S5[K], 5, -al, lam)[0]) / 2
            emit('D1', pt + '(+-a avg)', f'K{K}', str(mu * (1 - K)), 'WKB predicted (alpha-even)', pe, baseA[f'K{K}'], spA[f'K{K}'])
# ---------------------------------------------------------------- generic NUM2/3 datasets
def theta_dom(D):
    th = [mp.mpf(t) for t in D['thetas']]; return min(range(len(th)), key=lambda j: th[j])
def run_ds(ds, D, pK, K2, extras, P, statusmap, prec=None, pts_label=None):
    j = theta_dom(D); grid = [mp.mpf(e) for e in D['grid']]; vals = [mp.mpf(v[j]) for v in D['logabsQ']]
    base, sp, kw = stability(grid, vals, pK, K2, extras, D['dps'], prec(j) if prec else None)
    say(f'  ({ds} {D["point"]}: dominant theta = {mp.nstr(mp.mpf(D["thetas"][j]), 8)}; window K2 = {kw})')
    for nm, power, status, pred, note in statusmap(P, base):
        emit(ds, D['point'], nm, power, status, pred, base[nm], sp[nm], note)
def wkbrows(P, Ks, pK_str, pole_status):
    out = []
    for K in Ks:
        pv = P['pred'][str(K)]
        if pv[0] == 'pole': out.append((f'K{K}', pK_str(K), pole_status, None, 'Gamma-pole grade; finite part not identified'))
        elif pv[0] == 'zero': out.append((f'K{K}', pK_str(K), 'WKB predicted', mp.mpf(0), ''))
        else: out.append((f'K{K}', pK_str(K), 'WKB predicted', mp.mpf(pv[1]), ''))
    return out
# ---------------------------------------------------------------- D2: NUM2 t = 2
mp.mp.dps = 60
say('== D2 NUM2 (Sol 3, t = 2)')
WP2 = json.load(open(IN('num2', 'wkb_pred.json'))); LP2 = json.load(open(IN('num2', 'PREDICTION_NUM2_logE.json')))
pK2 = lambda K: mp.mpf(1 - K) / 2
ex2 = [('L-3/2', (-mp.mpf(3) / 2, True)), ('L-3', (mp.mpf(-3), True)), ('L-9/2', (-mp.mpf(9) / 2, True))]
D120 = json.load(open(IN('num2', 'data_2_A_120.json')))
for pt in ('A', 'B', 'C'):
    D = json.load(open(IN('num2', f'data_2_{pt}_60.json'))); P = WP2['2'][pt]; L = LP2[pt]
    def smap(P, base, L=L):
        rows = wkbrows(P, [0, 2, 3, 4, 5, 6, 7, 8], lambda K: str(Fr(1 - K, 2)), 'resonant finite part, unresolved')
        rows += [('L-3/2', '-3/2 log', 'non-WKB registered before data', mp.mpf(L['4']['logE_coefficient_pred']), '-R_4/2'),
                 ('L-3', '-3 log', 'non-WKB registered before data', mp.mpf(L['7']['logE_coefficient_pred']), 'R_7 = 0'),
                 ('L-9/2', '-9/2 log', 'non-WKB registered before data', mp.mpf(L['10']['logE_coefficient_pred']), '-(3/2) R_10'),
                 ('log', 'log', 'fitted free', None, ''), ('const', '0', 'fitted free', None, '')]
        return rows
    prec = None
    if pt == 'A':
        prec = lambda j, D=D: ([mp.mpf(e) for e in D120['grid']], [mp.mpf(v[j]) for v in D120['logabsQ']], [mp.mpf(e) for e in D['grid']], [mp.mpf(v[j]) for v in D['logabsQ']], 26)
    run_ds('D2', D, pK2, 26, ex2, P, smap, prec)
# ---------------------------------------------------------------- D3: NUM2 t = 9/4
say('== D3 NUM2 (Sol 3, t = 9/4)')
pK3 = lambda K: mp.mpf(13) * (1 - K) / 22
ex3 = [(f'N{m}', (-mp.mpf(13) * m / 10, False)) for m in (1, 2, 3)]
for pt in ('A', 'B'):
    D = json.load(open(IN('num2', f'data_10_3_{pt}_60.json'))); P = WP2['10/3'][pt]
    def smap(P, base):
        rows = wkbrows(P, [0, 2, 3, 4, 5, 6, 7, 8], lambda K: str(Fr(13 * (1 - K), 22)), 'unresolved')
        rows += [('N1', '-13/10', 'non-WKB position sealed, coefficient free', None, 'sealed F2'),
                 ('N2', '-13/5', 'non-WKB post hoc (position)', None, 'F3 post hoc'), ('N3', '-39/10', 'non-WKB post hoc (position)', None, 'F3 post hoc'),
                 ('log', 'log', 'fitted free', None, ''), ('const', '0', 'fitted free', None, '')]
        return rows
    run_ds('D3', D, pK3, 24, ex3, P, smap)
# ---------------------------------------------------------------- D4, D5: NUM3 (extended grid, E >= 1e4)
WP3 = json.load(open(IN('num3', 'wkb_pred3.json'))); LP3 = json.load(open(IN('num3', 'PREDICTION_NUM3_logE.json')))
def ext_data(k, pt):
    D = json.load(open(IN('num3', f'data_{k}_{pt}_60.json')))
    ext = sorted([json.load(open(f)) for f in glob.glob(IN('num3', 'points_ext', f'pt_{k}_{pt}_60_*.json')) if not f.endswith('_wrong.json')], key=lambda p: mp.mpf(p['E']))
    assert len(ext) == 12
    keep = [i for i, e in enumerate(D['grid']) if mp.mpf(e) >= mp.mpf(10) ** 4 - 1]
    D2 = dict(D); D2['grid'] = [D['grid'][i] for i in keep] + [p['E'] for p in ext]; D2['logabsQ'] = [D['logabsQ'][i] for i in keep] + [p['logabsQ'] for p in ext]
    return D, D2
say('== D4 NUM3 (Sol 2, t = 2; extended grid E >= 1e4)')
pK4 = lambda K: mp.mpf(3) * (1 - K) / 7
ex4 = [('N-3/2', (-mp.mpf(3) / 2, False)), ('N-9/2', (-mp.mpf(9) / 2, False)), ('N-15/2', (-mp.mpf(15) / 2, False)), ('L-3', (mp.mpf(-3), True)), ('L-6', (mp.mpf(-6), True))]
D120_3 = json.load(open(IN('num3', 'data_2_A_120.json')))
for pt in ('A', 'B'):
    D, Dx = ext_data('2', pt); P = WP3['2'][pt]
    def smap(P, base):
        rows = wkbrows(P, [0, 2, 3, 4, 5, 6, 7, 8], lambda K: str(Fr(3 * (1 - K), 7)), 'resonant finite part, unresolved')
        rows += [('N-3/2', '-3/2', 'non-WKB position sealed, coefficient free', None, 'odd m = 1'),
                 ('N-9/2', '-9/2', 'non-WKB position sealed, coefficient free', None, 'odd m = 3'),
                 ('N-15/2', '-15/2', 'non-WKB position sealed, coefficient free', None, 'odd m = 5'),
                 ('L-3', '-3 log', 'non-WKB registered before data', mp.mpf(0), 'R_8 = 0 (refined before data)'),
                 ('L-6', '-6 log', 'non-WKB registered before data', mp.mpf(0), 'R_15 = 0'),
                 ('log', 'log', 'fitted free', None, ''), ('const', '0', 'fitted free', None, '')]
        return rows
    prec = None
    if pt == 'A':
        def prec(j, D=D):
            idx120 = [i for i, e in enumerate(D120_3['grid']) if mp.mpf(e) >= mp.mpf(10) ** 4 - 1]
            idx60 = [i for i, e in enumerate(D['grid']) if mp.mpf(e) >= mp.mpf(10) ** 4 - 1]
            g2 = [mp.mpf(D120_3['grid'][i]) for i in idx120]; v2 = [mp.mpf(D120_3['logabsQ'][i][j]) for i in idx120]
            g1 = [mp.mpf(D['grid'][i]) for i in idx60]; v1 = [mp.mpf(D['logabsQ'][i][j]) for i in idx60]
            return (g2, v2, g1, v1, k2_window(len(g1), 24, len(ex4)))
    run_ds('D4', Dx, pK4, 24, ex4, P, smap, prec)
say('== D5 NUM3 (Sol 2, t = 9/4; extended grid E >= 1e4; extra)')
pK5 = lambda K: mp.mpf(13) * (1 - K) / 29
for pt in ('A', 'B'):
    D, Dx = ext_data('10_3', pt); P = WP3['10/3'][pt]
    def smap(P, base):
        rows = wkbrows(P, [0, 2, 3, 4, 5, 6, 7, 8], lambda K: str(Fr(13 * (1 - K), 29)), 'unresolved')
        rows += [('N1', '-13/10', 'non-WKB position sealed, coefficient free', None, ''), ('N2', '-13/5', 'non-WKB position sealed, coefficient free', None, ''),
                 ('N3', '-39/10', 'non-WKB position sealed, coefficient free', None, ''), ('log', 'log', 'fitted free', None, ''), ('const', '0', 'fitted free', None, '')]
        return rows
    run_ds('D5', Dx, pK5, 24, ex3, P, smap)
json.dump(ROWS, open(os.path.join(HERE, 'table.json'), 'w'), indent=1)
nt = sum(1 for r in ROWS if r['flag'] == 'TENSION')
say(f'\n{len(ROWS)} rows; TENSION rows: {nt}; done ({time.time()-T00:.0f}s)')
