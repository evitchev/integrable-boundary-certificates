"""VIR5c S3: formula F1 (record variables) against Codex's Ward-derived Sol-3 loss-1 shift amplitudes
(item 185 archive, AMPLITUDES_sol3_loss1.json), exactly, K = 2..KMAX at nine fibres.
Modes: real | neg_const (4 -> n-1) | neg_rho (rho^2 -> 0) | tamper (A_(1,1) + 1/1000, after the hash guard).
Exit: real 0 = all agree, 2 = mismatch; control modes 0 = the control FIRES (mismatches found), 1 = it does not."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 40
SRC = '<repo>/results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/AMPLITUDES_sol3_loss1.json'
raw = open(SRC, 'rb').read()
print('input sha256', hashlib.sha256(raw).hexdigest())
law = json.loads(raw)
ks, ts = sp.symbols('k t')
def load(d):
    return sp.sympify(d['numerator'], locals={'k': ks, 't': ts}, rational=True), sp.sympify(d['denominator'], locals={'k': ks, 't': ts}, rational=True)
amps = {int(j): load(v) for j, v in law['shift_amplitudes'].items()}
p_law = load(law['p']); q_law = load(law['q']); r_law = load(law['r'])
assert law['solution'] == 3 and law['loss'] == 1
# --- everything below is scoring ---
FIBRES = ['9/4', '21/11', '3/2', '4', '11/2', '-2', '1/2', '2', '7/3']
tot, bad, skipped = 0, 0, []
for tstr in FIBRES:
    t = sp.Rational(tstr)
    k, n, rho2 = L.fibre(t)
    fb, ft = 0, 0
    for K in range(2, KMAX + 1):
        sub = {ks: K, ts: t}
        vals = {}
        sing = False
        for j, (nu_, de_) in amps.items():
            de = de_.subs(sub)
            if de == 0:
                sing = True
                break
            vals[j] = sp.Rational(nu_.subs(sub)) / sp.Rational(de)
        p = sp.Rational(p_law[0].subs(sub)) / sp.Rational(p_law[1].subs(sub))
        q = sp.Rational(q_law[0].subs(sub)) / sp.Rational(q_law[1].subs(sub))
        r = sp.Rational(r_law[0].subs(sub)) / sp.Rational(r_law[1].subs(sub))
        nu = sp.Rational(2 * K - 1) / n
        assert p == -nu and q == -nu and r == 1
        d = K - 1
        if not sing and any(sp.rf(p + j, d) == 0 for j in vals):
            sing = True
        if mode == 'neg_const':
            pred = L.record_loss1_from_F1(K, t, const=n - 1)
        elif mode == 'neg_rho':
            pred = L.record_loss1_from_F1(K, t, with_rho=False)
        else:
            pred = L.record_loss1_from_F1(K, t)
        if sing or pred is None:
            skipped.append((tstr, K, 'amplitude or kernel pole' if sing else 'X^K coefficient vanishes'))
            continue
        if mode == 'tamper':
            vals[1] = vals[1] + sp.Rational(1, 1000)
        for a in range(d + 1):
            expect = sum(v * L.Kd_coeff(d, p + j, q, r, a) for j, v in vals.items())
            ft += 1
            if pred.get((a, d - a), 0) - expect != 0:
                fb += 1
    tot += ft
    bad += fb
    print('t = %-6s (k = %-6s n = %-6s): K = 2..%d, coefficients compared %d, mismatches %d' % (tstr, k, n, KMAX, ft, fb), flush=True)
print('skipped (K, t):', skipped)
print('S3 mode %s: %d mismatches / %d coefficients' % (mode, bad, tot))
if mode == 'real':
    sys.exit(0 if bad == 0 and tot > 0 else 2)
sys.exit(0 if bad > 0 else 1)
