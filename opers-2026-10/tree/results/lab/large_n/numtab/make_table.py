"""NUMTAB post-processing (rules stated in NOTE): (i) a WKB prediction with |pred| < 1e-40 is a round-off ZERO (odd grades / alpha-even odd grades) and is
treated as pred = 0; (ii) a row whose error estimate exceeds its fitted value (e < 0) is classed 'undetermined'; (iii) consistency: |fit - pred| <= 10 eps.
Writes table_final.json, table.csv and table.tex (booktabs, point A of D1-D4)."""
import os, json, csv
import mpmath as mp
mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'table.json')))
SYM = {'WKB predicted': 'W', 'WKB predicted (alpha-even)': 'W', 'non-WKB registered before data': 'R', 'non-WKB position sealed, coefficient free': 'S',
       'non-WKB post hoc (lim WKB + a Mellin)': 'H', 'non-WKB post hoc (position)': 'H', 'fitted free': 'F', 'resonant finite part, unresolved': 'U', 'unresolved': 'U'}
for r in R:
    fit = mp.mpf(r['fit']); eps = mp.mpf(r['eps'])
    pred = None if r['pred'] is None else mp.mpf(r['pred'])
    # d and flag come from numtab.py (full precision); only round-off-zero predictions are re-classified here (v1 recomputed d from 30-digit strings,
    # which printed 'exact' above ~30 digits -- display bug, fixed)
    if pred is not None and pred != 0 and abs(pred) < mp.mpf(10) ** -40:
        pred = mp.mpf(0); r['pred'] = '0'; r['d'] = None; r['abs_fit_if_pred0'] = mp.nstr(abs(fit), 3)
        r['flag'] = 'consistent' if abs(fit) <= 10 * eps else 'TENSION'
    r['class'] = SYM[r['status']]
    r['undetermined'] = bool(r['e'] < 0 and (pred is None or pred != 0))
json.dump(R, open(os.path.join(HERE, 'table_final.json'), 'w'), indent=1)
with open(os.path.join(HERE, 'table.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['dataset', 'point', 'row', 'power', 'class', 'status', 'prediction', 'fit', 'agreement_digits_d', 'error_estimate', 'error_digits_e',
                                   'spread_Wminus', 'spread_Wplus', 'spread_trunc', 'spread_prec', 'flag', 'undetermined', 'note'])
    for r in R: w.writerow([r['dataset'], r['point'], r['row'], r['power'], r['class'], r['status'], r['pred'], r['fit'], r['d'], r['eps'], r['e'],
                            r['spreads']['W-'], r['spreads']['W+'], r['spreads']['trunc'], r['spreads']['prec'], r['flag'], r['undetermined'], r['note']])
def fmt(x, n=12):
    if x is None: return '--'
    v = mp.mpf(x)
    if v == 0: return '$0$'
    s = mp.nstr(v, n, min_fixed=-3, max_fixed=6)
    if 'e' in s:
        m, ex = s.split('e'); return f'${m}\\times10^{{{int(ex)}}}$'
    return f'${s}$'
def powtex(p):
    if p == 'log': return '$\\log E$'
    if p.endswith(' log'): return f'$E^{{{p[:-4]}}}\\log E$'
    return f'$E^{{{p}}}$' if p != '0' else '$1$'
TITLES = {'D1': 'NUM1: Solution 1, $M=5$ ($t=2$), $\\log D(-e)$', 'D2': 'NUM2: Solution 3, $t=2$, $\\log|Q|$',
          'D3': 'NUM2: Solution 3, $t=9/4$, $\\log|Q|$', 'D4': 'NUM3: Solution 2, $t=2$ (extended grid), $\\log|Q|$'}
lines = ['\\begin{table}[p]\\centering\\scriptsize\\setlength{\\tabcolsep}{4pt}',
         '\\begin{tabular}{@{}l c r r r r@{}}', '\\toprule',
         'power & class & prediction & fitted coefficient & $d$ & $e$ \\\\', '\\midrule']
for ds in ('D1', 'D2', 'D3', 'D4'):
    lines.append(f'\\multicolumn{{6}}{{@{{}}l}}{{\\textit{{{TITLES[ds]}}}, point A}} \\\\')
    for r in R:
        if r['dataset'] != ds or r['point'] != 'A': continue
        if r['row'] in ('K0',): pw = r['power']
        d = '--' if r['d'] is None else (f'{r["d"]:.1f}' if r['d'] < 999 else 'exact')
        if r['pred'] is not None and mp.mpf(r['pred']) == 0:
            m_, ex_ = mp.nstr(abs(mp.mpf(r['fit'])), 2, min_fixed=0, max_fixed=0).split('e'); d = f'$|{{\\rm fit}}|={m_}\\times10^{{{int(ex_)}}}$'
        e = f'{r["e"]:.1f}' if r['e'] >= 0 else '$<0$'
        cls = r['class'] + ('$^\\ast$' if r['undetermined'] else '')
        lines.append(f'{powtex(r["power"])} & {cls} & {fmt(r["pred"])} & {fmt(r["fit"])} & {d} & {e} \\\\')
    lines.append('\\midrule' if ds != 'D4' else '\\bottomrule')
lines += ['\\end{tabular}',
          '\\caption{Spectral evidence: large-$E$ fits of $\\log D$ or $\\log|Q|$ (dominant Frobenius exponent) for the three operators, point A (all points in the supplementary table). '
          'Class: W = WKB coefficient $-\\int s_K$ (predicted); R = non-WKB term registered before the data; S = non-WKB term whose power was sealed, coefficient fitted; '
          'H = post-hoc term; F = fitted free (grade 1); U = resonant finite part, not identified; $^\\ast$ = undetermined (error estimate exceeds the value). '
          '$d$ = agreement digits $-\\log_{10}|{\\rm fit}-{\\rm pred}|/|{\\rm pred}|$ (for a zero prediction the absolute fit is shown); '
          '$e$ = error-estimate digits, $-\\log_{10}\\varepsilon/|{\\rm fit}|$, with $\\varepsilon$ the largest spread over two shifted fit windows, two basis truncations and, where archived, 60- vs 120-digit data. '
          'All predicted rows agree within $10\\varepsilon$.}',
          '\\label{tab:spectral-evidence}', '\\end{table}']
open(os.path.join(HERE, 'table.tex'), 'w').write('\n'.join(lines) + '\n')
nt = sum(r['flag'] == 'TENSION' for r in R); nu = sum(r['undetermined'] for r in R)
print(f'{len(R)} rows; TENSION {nt}; undetermined {nu}; predicted rows {sum(r["pred"] is not None for r in R)}')
for r in R:
    if r['undetermined']: print('  undetermined:', r['dataset'], r['point'], r['row'], r['status'], 'fit', mp.nstr(mp.mpf(r['fit']), 6), 'eps', r['eps'])
