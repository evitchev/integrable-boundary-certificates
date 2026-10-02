"""VIR5c S5: BLIND comparison of the registered spin-15 predictions with the certified Solution-3 table
results/lab/anchor15/vev_profile_sol3_w16.json.  Run ONLY after the lead's "registered".
Modes: real | tamper (one table coefficient + 1/1000, after the hash guards).
Exit: real 0 = all fibres agree, 2 = mismatch; tamper 0 = fires, 1 = does not."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
PRED = open('PREDICTION_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('PREDICTION_VIR5c.json', 'rb').read()).hexdigest() == PRED
assert open('REGISTERED_BY_LEAD.txt').read().strip() != ''
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
SRC = '<repo>/results/lab/anchor15/vev_profile_sol3_w16.json'
raw = open(SRC, 'rb').read()
print('table sha256', hashlib.sha256(raw).hexdigest())
tab = json.loads(raw)['vev']
ts = sp.Symbol('t')
K = 8
TAB = {}
for key, val in tab.items():
    a, b = [int(x) for x in key.split(',')]
    num, den = sp.fraction(sp.together(sp.sympify(val, locals={'t': ts}, rational=True)))
    TAB[a, b] = (sp.Poly(num, ts), sp.Poly(den, ts))
print('table monomials:', len(TAB), 'max total degree', max(a + b for a, b in TAB))
pred = json.load(open('PREDICTION_VIR5c.json'))['cells']
allbad, alltot = 0, 0
for key in sorted(pred):
    cell = pred[key]
    t = sp.Rational(cell['t'])
    k, n, rho2 = L.fibre(t)
    co = cell['spins']['15']['coefficients']
    if co is None:
        print('t = %s: operator gives no spin-15 charge' % cell['t']); allbad += 1; alltot += 1
        continue
    P = L.ZERO
    for kk, v in co.items():
        a, b = [int(x) for x in kk.split(',')]
        P += (L.A * (-1)) ** (a // 2) * (L.ONE * L.q_(rho2) - L.B) ** (b // 2) * L.q_(sp.Rational(v))
    got = L.as_dict(P, 1 / sp.Rational(int(P.coeff(L.A ** K).numerator), int(P.coeff(L.A ** K).denominator)))
    if any(de.eval(t) == 0 for _, de in TAB.values()):
        print('t = %s: table has a pole, UNSCREENED' % cell['t']); continue
    vals = {m: sp.Rational(nu.eval(t)) / sp.Rational(de.eval(t)) for m, (nu, de) in TAB.items()}
    top = vals[K, 0]
    if top == 0:
        print('t = %s: table top vanishes, UNSCREENED' % cell['t']); continue
    vals = {m: v / top for m, v in vals.items()}
    if mode == 'tamper':
        vals[3, 2] = vals[3, 2] + sp.Rational(1, 1000)
    mons = sorted(set(vals) | set(got))
    bad = [m for m in mons if got.get(m, 0) - vals.get(m, 0) != 0]
    per_loss = {}
    for m in mons:
        per_loss.setdefault(K - m[0] - m[1], [0, 0])
        per_loss[K - m[0] - m[1]][1] += 1
        if m in bad:
            per_loss[K - m[0] - m[1]][0] += 1
    allbad += len(bad); alltot += len(mons)
    print('t = %-6s k = %-6s spin 15: %d / %d mismatches;  by loss: %s' % (cell['t'], cell['k'], len(bad), len(mons),
          ', '.join('%d: %d/%d' % (l_, v[0], v[1]) for l_, v in sorted(per_loss.items()))), flush=True)
    if bad:
        m = bad[0]
        print('      first: X^%d Y^%d predicted %s, certified %s' % (m[0], m[1], got.get(m, 0), vals.get(m, 0)))
print('S5 mode %s: %d mismatches / %d coefficients' % (mode, allbad, alltot))
if mode == 'real':
    sys.exit(0 if allbad == 0 and alltot > 0 else 2)
sys.exit(0 if allbad > 0 else 1)
