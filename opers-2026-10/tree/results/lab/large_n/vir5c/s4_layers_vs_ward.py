"""VIR5c S4: the operator's layers at losses 1..4 against the Ward-derived shift amplitudes (Sol 3:
loss 1 from the item-185 archive, losses 2-4 from the item-180 archive), at every spin of the registered
prediction file; spins 15, 17, 19 are beyond every table used so far by this seat.
Modes: real | tamper (A_(2,1) + 1/1000 after the hash guards).
Exit: real 0 = all agree, 2 = mismatch;  tamper 0 = fires, 1 = does not."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
PRED = open('PREDICTION_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('PREDICTION_VIR5c.json', 'rb').read()).hexdigest() == PRED
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
AUD = '<repo>/results/lab/audits/'
SRC = {1: AUD + 'codex_wardsol12b_27403fc_2026-09-28/AMPLITUDES_sol3_loss1.json',
       2: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol3_loss2.json',
       3: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol3_loss3.json',
       4: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol3_loss4.json'}
ks, ts = sp.symbols('k t')
def load(d):
    loc = {'k': ks, 't': ts}
    return sp.Poly(sp.sympify(d['numerator'], locals=loc, rational=True), ks, ts), sp.Poly(sp.sympify(d['denominator'], locals=loc, rational=True), ks, ts)
LAW = {}
for m, fn in SRC.items():
    raw = open(fn, 'rb').read()
    print('loss %d input sha256 %s' % (m, hashlib.sha256(raw).hexdigest()))
    law = json.loads(raw)
    assert law['solution'] == 3 and law['loss'] == m
    LAW[m] = {int(j): load(v) for j, v in law['shift_amplitudes'].items()}
pred = json.load(open('PREDICTION_VIR5c.json'))['cells']
tot, bad, new_tot, new_bad, skipped = 0, 0, 0, 0, []
for key in sorted(pred):
    cell = pred[key]
    t = sp.Rational(cell['t'])
    k, n, rho2 = L.fibre(t)
    assert k == sp.Rational(cell['k'])
    for s in sorted(int(x) for x in cell['spins']):
        K = (s + 1) // 2
        co = cell['spins'][str(s)]['coefficients']
        if K < 2:
            continue
        if co is None:
            skipped.append((cell['t'], s, 'operator gives no charge'))
            continue
        # record variables
        P = L.ZERO
        for kk, v in co.items():
            a, b = [int(x) for x in kk.split(',')]
            P += (L.A * (-1)) ** (a // 2) * (L.ONE * L.q_(rho2) - L.B) ** (b // 2) * L.q_(sp.Rational(v))
        lead = P.coeff(L.A ** K)
        P = P * (1 / lead)
        got = L.as_dict(P)
        nu = sp.Rational(2 * K - 1) / n
        p = q = -nu
        row = []
        for m in range(1, 5):
            if m > K:
                continue
            d = K - m
            vals, sing = {}, False
            for j, (nu_, de_) in LAW[m].items():
                de = de_.eval({ks: K, ts: t})
                if de == 0 or sp.rf(p + j, d) == 0:
                    sing = True
                    break
                vals[j] = sp.Rational(nu_.eval({ks: K, ts: t})) / sp.Rational(de)
            if sing:
                skipped.append((cell['t'], s, 'loss %d: amplitude or kernel pole' % m))
                continue
            if mode == 'tamper' and m == 2:
                vals[1] = vals[1] + sp.Rational(1, 1000)
            fb = 0
            for a in range(d + 1):
                expect = sum(v * L.Kd_coeff(d, p + j, q, 1, a) for j, v in vals.items())
                if got.get((a, d - a), 0) - expect != 0:
                    fb += 1
            row.append('loss %d: %d/%d' % (m, fb, d + 1))
            tot += d + 1; bad += fb
            if s >= 15:
                new_tot += d + 1; new_bad += fb
        print('t = %-6s k = %-6s spin %2d  mismatches  %s' % (cell['t'], cell['k'], s, ';  '.join(row)), flush=True)
print('skipped:', skipped)
print('S4 mode %s: %d mismatches / %d coefficients in all; at spins >= 15: %d / %d' % (mode, bad, tot, new_bad, new_tot))
if mode == 'real':
    sys.exit(0 if bad == 0 and new_tot > 0 else 2)
sys.exit(0 if bad > 0 else 1)
