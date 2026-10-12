"""VIR2 (post-seal): the commutant of {pair(s), Virasoro beta_+} in ALL two-boson densities at
WEIGHT 8, one fibre (p = 3/2): is it one-dimensional mod d and spanned by the promoted I_7?
Exit 0 iff yes for all three solutions."""
import sys, time
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
sys.path.insert(0, _KITROOT + 'code')
import hashlib, json, pickle
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, b, rho
t = sp.Symbol('t'); q = sp.Symbol('q'); p = sp.Symbol('p')
src = open('w4_main.py').read()
exec(src[src.index('I5SYM = '):src.index('# ---------------------------------------------------------------- fibres')])
s6 = open('w6_structure.py').read()
exec(s6[s6.index('def fq(x):'):s6.index('for pv in (sp.Rational(3, 2)')])
t0 = time.time()
fails = 0
fib = Fibre(sp.Rational(3, 2))
print('p = 3/2: t = %s, N = %s' % (fib.t, fib.N), flush=True)
B = W.two_boson_basis(8)
idx = {m: i for i, m in enumerate(B)}
for sol in (1, 2, 3):
    pairs = [(qv, bv) for _, bv, qv in PAIRS(sol, fib.p)]
    dens = fib.Ipromoted(sol, 7)
    dvec = [F(0)] * len(B)
    for m, c in dens.items():
        dvec[idx[m]] = fq(c)
    dim, ns, _, drows = commutant(8, pairs + [(0, fib.p)], fib.rho)
    inside = in_span_mod_d(dvec, ns, drows)
    print('Sol %d weight 8 (%d monomials): commutant of pair(s) + Virasoro: dim mod d = %d ; promoted I_7 inside: %s  (%.0fs)'
          % (sol, len(B), dim, inside, time.time() - t0), flush=True)
    fails += (dim != 1) or (not inside)
print('WEIGHT8', 'FAILED' if fails else 'OK')
sys.exit(1 if fails else 0)
