"""VIR5d Z1: the operator's monic charge at k0 + h against the certified Sol-3 charge at the drop point t0.
Exact rationals; the printed magnitudes are floats of exact differences.  Exit 0 = O(h) convergence at all four."""
import hashlib, json, pickle, sys
import sympy as sp
sys.dont_write_bytecode = True
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
ROOT = _KITROOT + 'results/lab/'
t, PX, PI = sp.symbols('t PX PI')
ok = True
for kstr, spin in (('1', 5), ('3', 7), ('5', 9), ('1/2', 9)):
    k0 = sp.Rational(kstr)
    tv = (3 * k0 + 2) / (k0 + 2)
    rho2 = sp.Rational(2) / (tv**2 - 1)
    tab = json.load(open(ROOT + 'vev/vev_sol3_w%d.json' % (spin + 1)))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        f = sp.sympify(val, locals={'t': t}, rational=True)
        assert sp.fraction(sp.together(f))[1].subs(t, tv) != 0
        expr += f.subs(t, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    P = sp.Poly(sp.expand(expr), PX, PI)
    lead = P.coeff_monomial(PX**(spin + 1))
    assert lead != 0
    cert = {m: c / lead for m, c in P.terms()}
    errs = {}
    for hexp in (8, 16):
        d = json.load(open('z1_k%s_s%d_h%d.json' % (kstr.replace('/', '_'), spin, hexp)))
        pred = {tuple(int(x) for x in key.split(',')): sp.Rational(v) for key, v in d['coefficients'].items()}
        mons = set(pred) | set(cert)
        errs[hexp] = max(abs(pred.get(m, 0) - cert.get(m, 0)) for m in mons)
    r = errs[8] / errs[16] if errs[16] != 0 else None
    good = errs[8] < sp.Rational(1, 10**4) and r is not None and sp.Rational(10**8, 2) < r < 2 * 10**8
    ok &= bool(good)
    print('k0 = %-4s t0 = %-5s spin %d: %d coefficients; max|monic(k0 + h) - certified| = %.6e (h = 1e-8), %.6e (h = 1e-16); ratio %.6e  -> %s'
          % (kstr, tv, spin, len(cert), float(errs[8]), float(errs[16]), float(r), 'O(h): PASS' if good else 'FAIL'))
print('Z1', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
