"""VIR5d gate at K = 6, 7: the certified Sol-3 table charge (anchor11 / anchor13 profiles) lies in the kernel
computed by z2_kernel.py at t = 13/7, 23/9 (K = 6) and t = 9/7 (K = 7).  Exact.  Exit 0 iff all hold."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
ROOT = _KITROOT + 'results/lab/'
t = sp.Symbol('t')
ok = True
for tstr, K, fn in (('13/7', 6, 'anchor11/vev_profile_sol3_w12.json'), ('23/9', 6, 'anchor11/vev_profile_sol3_w12.json'), ('9/7', 7, 'anchor13/vev_profile_sol3_w14.json')):
    tv = sp.Rational(tstr)
    rho2 = sp.Rational(2) / (tv**2 - 1)
    PX, PI = sp.symbols('PX PI')
    tab = json.load(open(ROOT + fn))['vev']
    expr = 0
    for key, val in tab.items():
        a_, b_ = [int(x) for x in key.split(',')]
        f = sp.sympify(val, locals={'t': t}, rational=True)
        num, den = sp.fraction(sp.together(f))
        assert den.subs(t, tv) != 0
        expr += (num.subs(t, tv) / den.subs(t, tv)) * (-PX**2) ** a_ * (rho2 - PI**2) ** b_
    P = sp.Poly(sp.expand(expr), PX, PI)
    S = {m: c for m, c in P.terms()}
    kd = json.load(open('z2_kernel_t%s_K%d.json' % (tstr.replace('/', '_'), K)))
    ker = [{(int(k.split(',')[0]), int(k.split(',')[1])): sp.Rational(v) for k, v in c['eigenvalue'].items()} for c in kd['classes']]
    mons = sorted(set(S) | set(m for e in ker for m in e))
    rk = sp.Matrix([[e.get(m, 0) for m in mons] for e in ker]).rank()
    rk2 = sp.Matrix([[e.get(m, 0) for m in mons] for e in ker] + [[S.get(m, 0) for m in mons]]).rank()
    sym = all(S.get((b, a), 0) == c for (a, b), c in S.items())
    good = rk2 == rk == len(ker)
    ok &= good
    print('t = %-5s K = %d: kernel dimension %d; Sol-3 table charge (X^K coefficient %s) in the kernel: %s; table charge symmetric: %s'
          % (tstr, K, len(ker), 'nonzero' if S.get((2 * K, 0), 0) != 0 else 'ZERO', good, sym))
print('gate', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
