"""VIR4d step 1: the C_k predictions at spins 11 and 13 -> PREDICTION_VIR4d.json.
Opens NO data file.  Seal SEAL_VIR4d.md (hash guard)."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
from wkb_lib import *
SEAL = open('SEAL_VIR4d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR4d.md', 'rb').read()).hexdigest() == SEAL
print('seal hash ok', SEAL)
PX, PI = sp.symbols('PX PI')
t0 = time.time()
out = {'seal': SEAL, 'normalisation': 'monic in PX; coefficient of PX^i pi^j under key "i,j"', 'cells': {}}
CELLS = [(int(a), [11, 13]) for a in sys.argv[1:]] if len(sys.argv) > 1 else [(4, [11, 13]), (6, [11, 13]), (8, [11, 13]), (3, [11, 13])]
for k, spins in CELLS:
    n = k + 4
    Mv = sp.Rational(1, k) + sp.Rational(1, 10)      # TAMPER
    s_ = sp.Rational(n - 1, 2)
    frozen = [s_ + sp.Rational(n, k) * sp.Rational(2 * j - (k - 1), 2) for j in range(k)]
    ex = [s_ + l0, s_ - l0, s_ + l1, s_ - l1] + frozen
    c2v = n**2 * sp.Rational(k + 1, k)
    pr = Chain([('D', ex[i] - i) for i in range(n)] + [('p', -1)], '1', -1 if n % 2 == 0 else +1, n, (F(1, n), -1))
    pr.solve(14)
    print('k = %d (t = %s, order %d) solved to order 14  (%.0fs)' % (k, sp.Rational(3 * k + 2, k + 2), n, time.time() - t0), flush=True)
    odd_ok = all(sp.cancel(pr.R(kk)[0].subs(M, Mv)) == 0 for kk in (11, 13))
    for s in spins:
        Rk, ok = pr.R(s + 1)
        Rk = sp.cancel(Rk.subs(M, Mv))
        Pl = sp.Poly(sp.expand(Rk), l0, l1)
        assert ok and all(i % 2 == 0 and j % 2 == 0 for (i, j), _ in Pl.terms())
        Rm = sp.Poly(sp.expand(sum(cf * c2v ** ((i + j) // 2) * PX**i * PI**j for (i, j), cf in Pl.terms())), PX, PI)
        lead = Rm.coeff_monomial(PX**(s + 1))
        cell = {'t': str(sp.Rational(3 * k + 2, k + 2)), 'order': n, 'M': str(Mv), 'frozen': [str(f) for f in frozen],
                'scale_c2': str(c2v), 'leading_coefficient': str(lead), 'odd_orders_R11_R13_vanish': bool(odd_ok)}
        if lead == 0:
            cell['coefficients'] = None
        else:
            Rn = sp.Poly(sp.expand(Rm.as_expr() / lead), PX, PI)
            cell['coefficients'] = {'%d,%d' % (i, j): str(c) for (i, j), c in sorted(Rn.terms())}
            sym = sp.expand(Rn.as_expr() - Rn.as_expr().subs({PX: PI, PI: PX}, simultaneous=True)) == 0
            cell['symmetric_PX_pi'] = bool(sym)
        out['cells']['k=%d,spin=%d' % (k, s)] = cell
        print('   spin %d: leading coefficient %s; %d monomials; symmetric: %s  (%.0fs)'
              % (s, lead, len(cell['coefficients'] or {}), cell.get('symmetric_PX_pi'), time.time() - t0), flush=True)
fn = 'prediction_tamper_k%s.json' % '_'.join(sys.argv[1:])
json.dump(out, open(fn, 'w'), indent=1, sort_keys=True)
print('written', fn, hashlib.sha256(open(fn, 'rb').read()).hexdigest())
