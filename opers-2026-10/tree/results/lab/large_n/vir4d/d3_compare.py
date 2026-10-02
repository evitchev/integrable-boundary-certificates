"""VIR4d step 3: compare the REGISTERED predictions with the certified spin-11 and spin-13 tables.
Guards: seal and prediction hashes.  Table format (inspected only after registration): 'vev' is
a dict "a,b" -> rational function of t = coefficient of (P_art^2)^a (Q_art^2)^b.
Exit 0 = all blind cells pass (and, if present, the tamper fails and the k = 2 cell passes);
2 = some blind cell fails; 1 = screen / control failure."""
import hashlib, json, sys
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
import sympy as sp
sys.dont_write_bytecode = True
SEAL = 'ab1bd1d98fed6f7df04134ed26acccd8bbe07f5bcb631b9d02f83325f7025dc4'
PRED = 'b8446a9f28bc875c97ea80d1c6f6d5e51282b0b00b86ba5120cb2b5b9918d67c'
assert hashlib.sha256(open('SEAL_VIR4d.md', 'rb').read()).hexdigest() == SEAL
assert hashlib.sha256(open('PREDICTION_VIR4d.json', 'rb').read()).hexdigest() == PRED
print('seal and prediction hashes ok')
t, PX, PI = sp.symbols('t PX PI')
TAB = {11: json.load(open(_KITROOT + 'results/lab/anchor11/vev_profile_sol3_w12.json')),
       13: json.load(open(_KITROOT + 'results/lab/anchor13/vev_profile_sol3_w14.json'))}
for s, d in TAB.items():
    assert d['solution'] == 3 and d['weight'] == s + 1


def certified(s, tv):
    """monic (in PX) polynomial in PX, pi at the fibre; None if the screen fails."""
    p2 = 2 * (tv - 1) / (tv + 1)
    rho2 = sp.Rational(2) / (tv**2 - 1)                 # rho = p/2 - 1/p, rho^2 = p^2/4 - 1 + 1/p^2
    assert sp.simplify(p2 / 4 - 1 + 1 / p2 - rho2) == 0
    expr = 0
    for key, val in TAB[s]['vev'].items():
        a_, b_ = [int(x) for x in key.split(',')]
        f = sp.sympify(val, locals={'t': t})
        num, den = sp.fraction(sp.together(f))
        if den.subs(t, tv) == 0:
            return None, 'pole in coefficient %s' % key
        expr += f.subs(t, tv) * (-PX**2) ** a_ * (rho2 - PI**2) ** b_
    P = sp.Poly(sp.expand(expr), PX, PI)
    lead = P.coeff_monomial(PX**(s + 1))
    if lead == 0:
        return None, 'PX^%d coefficient vanishes' % (s + 1)
    return sp.Poly(sp.expand(P.as_expr() / lead), PX, PI), 'ok'


def compare(cells, label):
    results = {}
    for name in sorted(cells, key=lambda x: (int(x.split(',')[0][2:]), int(x.split('=')[2]))):
        c = cells[name]
        k = int(name.split(',')[0][2:]); s = int(name.split('=')[2])
        tv = sp.Rational(3 * k + 2, k + 2)
        assert str(tv) == c['t']
        cert, info = certified(s, tv)
        if cert is None:
            print('%s %s: UNSCREENED (%s)' % (label, name, info)); results[name] = None; continue
        pred = {tuple(int(x) for x in key.split(',')): sp.Rational(v) for key, v in c['coefficients'].items()}
        sym = sp.expand(cert.as_expr() - cert.as_expr().subs({PX: PI, PI: PX}, simultaneous=True)) == 0
        mons = sorted(set(pred) | set(cert.monoms()))
        indep = [m for m in mons if m[0] >= m[1]]
        bad = [m for m in indep if pred.get(m, 0) != cert.coeff_monomial(m)]
        bad_all = [m for m in mons if pred.get(m, 0) != cert.coeff_monomial(m)]
        print('%s %s (t = %s): certified table symmetric under PX <-> pi: %s ; monomials %d (independent %d) ; mismatches: %d of %d independent, %d of %d in all'
              % (label, name, tv, sym, len(mons), len(indep), len(bad), len(indep), len(bad_all), len(mons)))
        if bad:
            m = bad[0]
            print('      first mismatch PX^%d pi^%d: predicted %s, certified %s' % (m[0], m[1], pred.get(m, 0), cert.coeff_monomial(m)))
        else:
            m = (0, 0)
            print('      e.g. constant term: predicted %s = certified %s' % (pred.get(m, 0), cert.coeff_monomial(m)))
        results[name] = (len(bad) == 0 and sym)
    return results


pred = json.load(open('PREDICTION_VIR4d.json'))
res = compare(pred['cells'], 'BLIND')
status = 0
if any(v is False for v in res.values()):
    status = 2
print('\nBLIND CELLS: %d pass, %d fail, %d unscreened' % (sum(1 for v in res.values() if v), sum(1 for v in res.values() if v is False), sum(1 for v in res.values() if v is None)))
for extra, label, must_pass in (('prediction_tamper_k4.json', 'TAMPER (M = 1/4 + 1/10)', False), ('prediction_partial_2.json', 'NON-BLIND k = 2', True)):
    try:
        d = json.load(open(extra))
    except FileNotFoundError:
        print('%s: file %s not present yet' % (label, extra)); continue
    r = compare(d['cells'], label)
    ok = all(v is must_pass for v in r.values())
    print(('PASS  ' if ok else 'CONTROL-FAIL  ') + '%s: %s' % (label, 'all cells pass' if must_pass else 'all cells fail'))
    if not ok:
        status = 1
print('RESULT:', {0: 'all blind cells PASS', 2: 'a blind cell FAILS', 1: 'control failure'}[status])
sys.exit(status)
