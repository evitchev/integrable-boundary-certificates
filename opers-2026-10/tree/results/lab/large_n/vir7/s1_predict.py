"""VIR7 T2: predictions of the Sol-1 operator S1(t) at one fibre, all orders 2..ORDER (odd-spin AND even-spin charges).
    (T^2 - l1^2) sigma^a R_a((T - l0)/sigma) psi = x^n (x^(n/a) - E) psi,  a = (t-1)/2, n = (t+3)/2, sigma = n/a.
Opens no data file.  Usage: s1_predict.py <t> <order> [tamper_a]"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_gen as MG
assert hashlib.sha256(open('SEAL_VIR7_addendum.md', 'rb').read()).hexdigest() == open('SEAL_VIR7_addendum.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); ORDER = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'real'
a = (t - 1) / 2
n = (t + 3) / 2
sigma = n / a
a_string = a + (sp.Rational(1, 10) if mode == 'tamper_a' else 0)
p2 = 2 * (t - 1) / (t + 1)
aY = (t - 1) / (t + 3)
t0 = time.time()
blocks = [[1, 0, -MG.l1**2] + [0] * (ORDER - 1), MG.string_series(a_string, MG.l0, sigma, ORDER + 1)]
if mode == 'tamper_a':      # keep the classical power a: multiply by T^(a - a_string) is implicit (degree n fixed); only the string length changes
    pass
h = MG.symbol_from_blocks(blocks, ORDER + 1)
eng = MG.MellinWKB(n, 1 / a, h)
PX, PI = sp.symbols('PX PI')
S1 = n * n / p2            # l1^2 = S1 pi^2
S0 = n * n / (p2 * aY)     # l0^2 = S0 PX^2
cell = {'t': str(t), 'a': str(a), 'n': str(n), 'M': str(1 / a), 'sigma': str(sigma), 'mode': mode, 'order': ORDER,
        'l1^2/pi^2': str(S1), 'l0^2/PX^2': str(S0), 'odd_spins': {}, 'even_spins': {}}
tag = str(t).replace('/', '_').replace('-', 'm') + ('' if mode == 'real' else '_' + mode)
fn = 's1_pred_t%s.json' % tag
eng.W = {}
for N in range(1, ORDER + 1):
    rest = eng.stage(N)
    eng.W[N] = MG.el_scale(rest, sp.Rational(-1) / MG.R_(eng.n))
    if N < 2:
        continue
    Ri, ok = eng.R(N)
    assert ok
    Pl = sp.Poly(Ri, MG.l0, MG.l1)
    if N % 2 == 0:
        assert Pl.is_zero or all(i % 2 == 0 and l % 2 == 0 for (i, l), _ in Pl.terms())
        Rm = sp.Poly(sp.expand(sum(cf * S0 ** (i // 2) * S1 ** (l // 2) * PX**i * PI**l for (i, l), cf in Pl.terms())), PX, PI)
        lead = Rm.coeff_monomial(PX**N)
        if lead == 0:
            cell['odd_spins'][str(N - 1)] = {'leading_coefficient': '0', 'identically_zero': Pl.is_zero, 'coefficients': None}
        else:
            Rn = sp.Poly(sp.expand(Rm.as_expr() / lead), PX, PI)
            cell['odd_spins'][str(N - 1)] = {'leading_coefficient': str(lead), 'coefficients': {'%d,%d' % m: str(c) for m, c in sorted(Rn.terms())}}
    else:
        # even spin N-1: odd in l0; store R/l0 in (PX^2, pi^2) up to the irrational factor sqrt(S0) PX, normalised by its top PX coefficient
        assert Pl.is_zero or all(i % 2 == 1 and l % 2 == 0 for (i, l), _ in Pl.terms())   # v1: a vanishing charge (spin in nZ) is allowed
        Rm = sp.Poly(sp.expand(sum(cf * S0 ** ((i - 1) // 2) * S1 ** (l // 2) * PX**(i - 1) * PI**l for (i, l), cf in Pl.terms())), PX, PI)
        lead = Rm.coeff_monomial(PX**(N - 1)) if not Rm.is_zero else 0
        cell['even_spins'][str(N - 1)] = {'identically_zero': Pl.is_zero, 'divided_by': 'PX * (top coefficient)',
                                          'coefficients': None if lead == 0 else {'%d,%d' % m: str(c / lead) for m, c in sorted(Rm.terms())}}
    print('t = %s order %d (spin %d) done  (%.0fs)' % (t, N, N - 1, time.time() - t0), flush=True)
    json.dump(cell, open(fn, 'w'), indent=1, sort_keys=True)
print('written', fn, hashlib.sha256(open(fn, 'rb').read()).hexdigest())
