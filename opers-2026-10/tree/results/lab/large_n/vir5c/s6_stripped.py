"""VIR5c S6: the question itself -- does the Gamma ratio produce the record's x/sinh x law?
(a) exact: record-normalised loss-1 layer of the STRIPPED operator (Gamma ratio -> u^k) minus the true one
    == dA(K) K_(K-1)(p, q; 1),  dA(K) = (n-5) nu K / (24 (nu - K + 1)),  nu = (2K-1)/n,
    a single shift amplitude growing like K, while the true A_(1,j)(K) grow like K^2.
(b) illustration (not gated): size of (stripped - true) relative to the true layer at losses 1, 2, 3.
ODE side only (engine files); the Ward amplitudes are read only to print A_(1,j)(K) next to dA(K).
Exit 0 = (a) holds at every K; 2 otherwise."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
true = json.load(open('pred_k10_3_o20.json'))
strp = json.load(open('pred_k10_3_stripped_o20.json'))
assert true['mode'] == 'real' and strp['mode'] == 'stripped' and true['k'] == strp['k'] == '10/3'
t = sp.Rational(true['t'])
k, n, rho2 = L.fibre(t)
law = json.load(open('<repo>/results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/AMPLITUDES_sol3_loss1.json'))
ks, ts = sp.symbols('k t')
A1 = {int(j): sp.sympify(v['numerator'], locals={'k': ks, 't': ts}, rational=True) / sp.sympify(v['denominator'], locals={'k': ks, 't': ts}, rational=True)
      for j, v in law['shift_amplitudes'].items()}


def record_layers(cell, s):
    K = (s + 1) // 2
    P = L.ZERO
    for kk, v in cell['spins'][str(s)]['coefficients'].items():
        a, b = [int(x) for x in kk.split(',')]
        P += (L.A * (-1)) ** (a // 2) * (L.ONE * L.q_(rho2) - L.B) ** (b // 2) * L.q_(sp.Rational(v))
    c = P.coeff(L.A ** K)
    return L.as_dict(P, 1 / sp.Rational(int(c.numerator), int(c.denominator)))


bad = 0
smax = min(max(int(x) for x in true['spins']), max(int(x) for x in strp['spins']))
print('fibre t = %s (k = %s, n = %s); spins up to %d' % (t, k, n, smax))
for s in range(3, smax + 1, 2):
    K = (s + 1) // 2
    gt, gs = record_layers(true, s), record_layers(strp, s)
    nu = sp.Rational(2 * K - 1) / n
    dA = (n - 5) * nu * K / (24 * (nu - K + 1))
    d = K - 1
    fb = sum(1 for a in range(d + 1) if gs.get((a, d - a), 0) - gt.get((a, d - a), 0) - dA * L.Kd_coeff(d, -nu, -nu, 1, a) != 0)
    bad += fb
    amps = [A1[j].subs({ks: K, ts: t}) for j in range(3)]
    rel = {}
    for m in (1, 2, 3):
        if m > K:
            continue
        dd = K - m
        num = max(abs(gs.get((a, dd - a), 0) - gt.get((a, dd - a), 0)) for a in range(dd + 1))
        den = max(abs(gt.get((a, dd - a), 0)) for a in range(dd + 1))
        rel[m] = float(num / den) if den != 0 else None
    print('K = %2d (spin %2d): (a) mismatches %d/%d;  dA = %-12s (%.4f)   A_(1,0..2) = %s;   max|stripped - true|/max|true| at losses 1,2,3: %s'
          % (K, s, fb, d + 1, dA, float(dA), ['%.3f' % float(x) for x in amps], ['%.4f' % rel[m] if rel.get(m) is not None else '-' for m in (1, 2, 3)]), flush=True)
print('S6(a)', 'PASS' if bad == 0 else 'FAIL', bad)
sys.exit(0 if bad == 0 else 2)
