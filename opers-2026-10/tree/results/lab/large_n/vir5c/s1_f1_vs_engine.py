"""VIR5c S1 (ODE side only): formula F1 against the Mellin-WKB engine's loss-1 layers, K = 1..7.
Engine values: the per-fibre VIR5a prediction files (operator output; no data table involved).
Exit 0 = F1 == engine everywhere (true operator with const = 4; stripped operator with const = n - 1)."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
V5 = '<home>/fable-work/vir5/'
cases = [('10/3', 'prediction_k10_3.json', 'true'), ('2/3', 'prediction_k2_3.json', 'true'),
         ('-6', 'prediction_km6.json', 'true'), ('10/3', 'prediction_k10_3_tamper_gamma.json', 'stripped')]
bad = 0
for kstr, fn, kind in cases:
    cell = json.load(open(V5 + fn))
    k = sp.Rational(kstr)
    n = k + 4
    assert sp.Rational(cell['k']) == k and cell['mode'] == ('real' if kind == 'true' else 'tamper_gamma')
    const = 4 if kind == 'true' else n - 1
    for K in range(1, 8):
        eng = cell['spins'][str(2 * K - 1)]['coefficients']
        nu = sp.Rational(2 * K - 1) / n
        lead = sp.binomial(nu, K) * (-1) ** K
        f1d = L.as_dict(L.loss1_operator(K, n, sgn=-1, const=const), 1 / lead)
        mism = 0
        cnt = 0
        for a in range(0, K):
            b = K - 1 - a
            e = sp.Rational(eng.get('%d,%d' % (2 * a, 2 * b), '0'))
            cnt += 1
            if sp.cancel(e - f1d.get((a, b), 0)) != 0:
                mism += 1
        # also the top layer as a sanity check of normalisation
        top = L.as_dict(L.top_operator(K, n, sgn=-1), 1 / lead)
        tmis = sum(1 for (a, b), c in top.items() if sp.Rational(eng.get('%d,%d' % (2 * a, 2 * b), '0')) - c != 0)
        bad += mism + tmis
        print('k = %-5s %-8s K = %d (spin %2d): loss-1 mismatches %d/%d, top mismatches %d/%d'
              % (kstr, kind, K, 2 * K - 1, mism, cnt, tmis, K + 1), flush=True)
print('S1', 'PASS' if bad == 0 else 'FAIL', 'total mismatches', bad)
sys.exit(0 if bad == 0 else 2)
