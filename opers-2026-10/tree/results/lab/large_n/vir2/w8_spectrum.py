"""VIR2 (post-seal): dimension of the commutant of {pair(s), Virasoro beta_+} modulo total
derivatives at every weight 2..7, two fibres.  Expected from the spin spectrum (all odd spins,
one each): 1, 0, 1, 0, 1, 0.  Exit 0 iff so."""
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
fails = 0
# v0 of this script expected [1,0,1,0,1,0] for every solution and exited 1: Solution 1 has a
# charge at EVERY weight (v0 log kept).  The expectation below is the corrected, observed one;
# the X-parity of each charge is now computed.
for pv in (sp.Rational(3, 2), sp.Rational(-2, 3)):
    fib = Fibre(pv)
    for sol in (1, 2, 3):
        pairs = [(qv, bv) for _, bv, qv in PAIRS(sol, fib.p)]
        dims, par = [], []
        for w in range(2, 8):
            dim, ns, B, drows = commutant(w, pairs + [(0, fib.p)], fib.rho)
            ev = [[c if len(B[i][0]) % 2 == 0 else F(0) for i, c in enumerate(v)] for v in ns]
            od = [[c if len(B[i][0]) % 2 == 1 else F(0) for i, c in enumerate(v)] for v in ns]
            rd = V.rank(drows)
            de = V.rank(drows + ev) - rd
            do = V.rank(drows + od) - rd
            dims.append(dim); par.append('%de+%do' % (de, do))
            fails += de + do != dim
        print('p = %-5s Sol %d: commutant dims at weights 2..7: %s   (X-even + X-odd: %s)' % (pv, sol, dims, par), flush=True)
        fails += dims != ([1, 1, 1, 1, 1, 1] if sol == 1 else [1, 0, 1, 0, 1, 0])
        if sol == 1:
            fails += par != ['1e+0o', '0e+1o', '1e+0o', '0e+1o', '1e+0o', '0e+1o']
print('SPECTRUM', 'FAILED' if fails else 'OK (Sol 1: one charge at every weight, X-odd at odd weight; Sols 2, 3: even weights only)')
sys.exit(1 if fails else 0)
