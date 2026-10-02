"""SUPER3 engine controls: (a) the continued integral of a total derivative vanishes; (b) alpha = 0 reproduces BLZ vacuum I_3, I_5
(hep-th/9812247 dictionary, Delta-polynomials: I3 = D^2 - (c+2)/12 D + c(5c+22)/2880, I5 = D^3 - (c+4)/8 D^2 + (c+2)(3c+20)/576 D - c(3c+14)(7c+68)/290304)."""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge, fder
M = F(7, 3)
R = riccati(M, 6)
for n in range(2, 7, 2):
    ch = charge(fder(R[n], M), M)
    print('int d(R_%d) by half-integer class (must vanish):' % n, {str(c): all(not v for v in p.values()) for c, (b0, f0, p) in ch.items() if c[0] != 'INTEGER-f'})
for n in range(1, 7):
    ch = charge(R[n], M)
    print('R_%d classes:' % n, {str(c): {k: str(v) for k, v in p.items()} for c, (b0, f0, p) in ch.items()})
