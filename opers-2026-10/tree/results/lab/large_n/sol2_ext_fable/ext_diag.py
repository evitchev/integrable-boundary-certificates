"""REVIEW diagnostic: the two mismatching cells of ext_fibres (Sol 2, t = 4 spin 3; t = 11/2 spin 7): print both polynomials; Beta arguments."""
import sys
from fractions import Fraction as F
import sympy as sp
from common import *
from ext_lib import run, score
for tstr, spin in (('4', 3), ('11/2', 7), ('7', 5)):
    tv = sp.Rational(tstr); kq = 2 * (tv - 1) / (3 - tv); k = F(int(kq.p), int(kq.q))
    tv_, g, Cv, ch = run(k, 1 / k)
    cert = certified(2, tv)
    n = g['D']; sig = n / k
    print(f"t = {tstr}: k = {k}, n = {n}, sigma = n/k = {sig}; spin {spin}: nu = spin/n = {F(spin)/n}; A0 = -spin/sigma = {-F(spin)/sig}; B0+1 = {F(spin)*(1+k)/n}")
    w = monic(ch[spin], (spin + 1) // 2); c = cert[spin]
    print('   WKB monic :', w)
    print('   certified :', c)
    print('   difference:', sp.factor(sp.expand(w - c)))
    print('   raw WKB charge (not normalised):', sp.factor(ch[spin]))
