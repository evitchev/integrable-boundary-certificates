"""VIR5a step 3: compare the REGISTERED predictions of the continued operator with the certified
Solution-3 vacuum eigenvalues.  Run only after the lead's "registered".
Exit 0 = cell A passes and controls fire; 2 = cell A fails; 1 = control failure."""
import hashlib, json, pickle, sys
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
import sympy as sp
sys.dont_write_bytecode = True
SEAL = '2efea4749de60fcf25885ce40afbfa0c60c2576f2025f42292fb277de94312c5'
PRED = '0c12a7e9e2108abcde2f3f651326a9eb3148d9aa82862f1ffb802b3c50df8751'
assert hashlib.sha256(open('SEAL_VIR5a.md', 'rb').read()).hexdigest() == SEAL
assert hashlib.sha256(open('PREDICTION_VIR5a.json', 'rb').read()).hexdigest() == PRED
print('seal and prediction hashes ok')
t, p, PX, PI, Pa, Qa = sp.symbols('t p PX PI Pa Qa')
ROOT = _KITROOT + 'results/lab/'
Fsym = {k: sp.Poly(sp.sympify(v), PX, PI) for k, v in pickle.load(open('cft_F.pkl', 'rb')).items()}
T10 = json.load(open(ROOT + 'vev/vev_sol3_w10.json'))
TAB = {11: json.load(open(ROOT + 'anchor11/vev_profile_sol3_w12.json')), 13: json.load(open(ROOT + 'anchor13/vev_profile_sol3_w14.json'))}


def certified(s, tv):
    p2 = 2 * (tv - 1) / (tv + 1)
    rho2 = sp.Rational(2) / (tv**2 - 1)
    if s == 1:
        return sp.Poly(PX**2 + PI**2 - sp.Rational(1, 6), PX, PI), 'ok'
    if s in (3, 5, 7):
        out = 0
        for (i, j), cf in Fsym[s + 1].terms():
            num, den = sp.fraction(sp.together(cf))
            pn, pd = sp.Poly(num, p), sp.Poly(den, p)
            assert all(e % 2 == 0 for (e,), _ in pn.terms()) and all(e % 2 == 0 for (e,), _ in pd.terms())
            dv = sum(c_ * p2 ** (e // 2) for (e,), c_ in pd.terms())
            if dv == 0:
                return None, 'pole'
            out += sum(c_ * p2 ** (e // 2) for (e,), c_ in pn.terms()) / dv * PX**i * PI**j
        return sp.Poly(out, PX, PI), 'ok'
    expr = 0
    if s == 9:
        for key, val in T10['coefficients'].items():
            a_, b_ = [int(x.split('^')[1]) for x in key.split()]
            f = sp.sympify(val, locals={'t': t})
            if sp.fraction(sp.together(f))[1].subs(t, tv) == 0:
                return None, 'pole'
            expr += f.subs(t, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    else:
        for key, val in TAB[s]['vev'].items():
            a_, b_ = [int(x) for x in key.split(',')]
            f = sp.sympify(val, locals={'t': t})
            if sp.fraction(sp.together(f))[1].subs(t, tv) == 0:
                return None, 'pole'
            expr += f.subs(t, tv) * (-PX**2) ** a_ * (rho2 - PI**2) ** b_
    P = sp.Poly(sp.expand(expr), PX, PI)
    lead = P.coeff_monomial(PX**(s + 1))
    if lead == 0:
        return None, 'top zero'
    return sp.Poly(sp.expand(P.as_expr() / lead), PX, PI), 'ok'


def compare_cell(d, label):
    tv = sp.Rational(d['t']); n = sp.Rational(d['n'])
    tot_bad, tot_n, rows = 0, 0, []
    for s in (1, 3, 5, 7, 9, 11, 13):
        cell = d['spins'][str(s)]
        cert, info = certified(s, tv)
        mult = (sp.Rational(s) / n).is_integer
        if cert is None:
            rows.append('spin %d: UNSCREENED (%s)' % (s, info)); continue
        if cell['coefficients'] is None:
            rows.append('spin %d: operator gives NO charge (spin/n = %s%s); certified polynomial exists' % (s, sp.Rational(s) / n, ', an integer' if mult else ''))
            if not mult:
                tot_bad += 1; tot_n += 1
            continue
        pred = {tuple(int(x) for x in key.split(',')): sp.Rational(v) for key, v in cell['coefficients'].items()}
        mons = sorted(set(pred) | set(cert.monoms()))
        bad = [m for m in mons if pred.get(m, 0) != cert.coeff_monomial(m)]
        tot_bad += len(bad); tot_n += len(mons)
        row = 'spin %d: %d / %d mismatches' % (s, len(bad), len(mons))
        if bad:
            m = bad[0]
            row += ' (first: PX^%d pi^%d predicted %s, certified %s)' % (m[0], m[1], pred.get(m, 0), cert.coeff_monomial(m))
        rows.append(row)
    print('%s k = %s (t = %s, n = %s): %d mismatches in %d coefficients' % (label, d['k'], d['t'], d['n'], tot_bad, tot_n))
    for r in rows:
        print('      ' + r)
    return tot_bad == 0, tot_bad, tot_n


pred = json.load(open('PREDICTION_VIR5a.json'))['cells']
GROUP = {'10/3': 'A', '5/3': 'A', '2/3': 'B', '1/2': 'B', '-6': 'C (regression, cc item 297)', '-18/5': 'C', '-6/5': 'C', '-2/5': 'D', '4': 'CONTROL'}
res = {}
for kk in ('10/3', '5/3', '2/3', '1/2', '-6', '-18/5', '-6/5', '-2/5', '4'):
    res[kk] = compare_cell(pred['k=' + kk], GROUP[kk])
status = 0
if not (res['10/3'][0] and res['5/3'][0]):
    status = 2
# control: k = 4 against the registered VIR4d prediction
try:
    old = json.load(open('../vir4d/PREDICTION_VIR4d.json'))['cells']
    same = all(old['k=4,spin=%d' % s]['coefficients'] == pred['k=4']['spins'][str(s)]['coefficients'] for s in (11, 13))
    print(('PASS  ' if same else 'CONTROL-FAIL  ') + 'k = 4 from the continued formula equals the registered VIR4d prediction at spins 11, 13 (string-identical coefficients)')
    if not same or not res['4'][0]:
        status = 1
except FileNotFoundError:
    print('VIR4d prediction not found'); status = 1
for fn, label in (('prediction_k10_3_tamper_gamma.json', 'TAMPER Gamma ratio -> u^k'), ('prediction_k10_3_tamper_M.json', 'TAMPER M = 1/k + 1/10')):
    try:
        d = json.load(open(fn))
    except FileNotFoundError:
        print('%s: not present' % label); status = max(status, 1); continue
    ok, nb, nt = compare_cell(d, label)
    print(('PASS  ' if not ok else 'CONTROL-FAIL  ') + '%s fails (%d of %d)' % (label, nb, nt))
    if ok:
        status = 1
print('\nSUMMARY:', {k_: '%d/%d mismatches' % (v[1], v[2]) for k_, v in res.items()})
print('RESULT:', {0: 'cell A PASSES; controls fired', 2: 'cell A FAILS', 1: 'control failure'}[status])
sys.exit(status)
