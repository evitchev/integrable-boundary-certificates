"""VIR5d Z0(b): F1' (polynomial form of the loss-1 layer) == F1, and == Codex's Sol-3 loss-1 amplitudes.
    loss-1 = (nu/6) [w^(K-1)] phi^(nu-1) { 1 + (K - 3/2 - nu) e1 w + (4 nu - 2K + 2) e2 w^2 }
Also: at nu = i integer the layer vanishes identically iff 2i <= K-1 (scan of lattice points).
Exit 0 iff all hold."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import f1_lib as L
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
A, B, ONE, ZERO, QQ = L.A, L.B, L.ONE, L.ZERO, L.QQ


def f1prime(K, n, sgn=-1):
    """un-normalised loss-1 layer, phi = (1 + sgn A w)(1 + sgn B w)"""
    n = sp.Rational(n)
    nu = sp.Rational(2 * K - 1) / n
    N = K - 1
    e1w = (A + B) * QQ(-sgn)           # phi = 1 - e1 w + e2 w^2  with e1 = -sgn (A + B), e2 = A B
    e2 = A * B
    br = [ONE, e1w * L.q_(K - sp.Rational(3, 2) - nu), e2 * L.q_(4 * nu - 2 * K + 2)]
    ph = L.ser_mul(L.binom_series(nu - 1, A * QQ(sgn), N), L.binom_series(nu - 1, B * QQ(sgn), N), N)
    tot = ZERO
    for j in range(min(2, N) + 1):
        tot += br[j] * ph[N - j]
    return tot * L.q_(nu / 6)


bad = 0
for kstr in ('10/3', '2/3', '-6', '5/3', '-2/5'):
    k = sp.Rational(kstr); n = k + 4
    for K in range(1, 13):
        if f1prime(K, n) != L.loss1_operator(K, n, sgn=-1, const=4):
            bad += 1
            print('F1 != F1prime at k = %s, K = %d' % (kstr, K))
print('F1prime == F1 for K = 1..12 at five fibres: mismatches', bad)

# against Codex's amplitudes, as S3 of VIR5c
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
SRC = _KITROOT + 'results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/AMPLITUDES_sol3_loss1.json'
law = json.load(open(SRC))
ks, ts = sp.symbols('k t')
amps = {int(j): (sp.sympify(v['numerator'], locals={'k': ks, 't': ts}, rational=True), sp.sympify(v['denominator'], locals={'k': ks, 't': ts}, rational=True))
        for j, v in law['shift_amplitudes'].items()}
tot = wbad = 0
for tstr in ('9/4', '3/2', '4', '-2', '2'):
    t = sp.Rational(tstr); k, n, rho2 = L.fibre(t)
    for K in range(2, 41):
        nu = sp.Rational(2 * K - 1) / n
        lead = sp.binomial(nu, K)
        l1 = f1prime(K, n, sgn=+1)
        shift = ZERO
        for a in range(K):
            shift += A**a * B**(K - 1 - a) * L.q_(-rho2 * nu * sp.binomial(nu, a) * sp.binomial(nu - 1, K - 1 - a))
        pred = L.as_dict(l1 + shift, 1 / lead)
        vals = {j: sp.Rational(nu_.subs({ks: K, ts: t})) / sp.Rational(de_.subs({ks: K, ts: t})) for j, (nu_, de_) in amps.items()}
        d = K - 1
        for a in range(d + 1):
            expect = sum(v * L.Kd_coeff(d, -nu + j, -nu, 1, a) for j, v in vals.items())
            tot += 1
            wbad += pred.get((a, d - a), 0) != expect
print('F1prime == Ward loss-1 amplitudes, K = 2..40, five fibres: %d mismatches / %d' % (wbad, tot))

# vanishing at lattice points: nu = i, K, n = (2K-1)/i
vbad = 0
rows = []
for K in range(2, 13):
    for i in range(1, K):
        n = sp.Rational(2 * K - 1, i)
        zero_l1 = f1prime(K, n) == ZERO
        zero_top = L.top_operator(K, n) == ZERO
        drop = 2 * i <= K - 1
        if zero_l1 != drop or zero_top != drop:
            vbad += 1
            rows.append((K, i, zero_top, zero_l1, drop))
print('lattice scan K = 2..12, i = 1..K-1: (top vanishes) == (loss-1 vanishes) == (2i <= K-1) violated at', vbad, 'points', rows)
ok = bad == 0 and wbad == 0 and vbad == 0
print('Z0(b)', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
