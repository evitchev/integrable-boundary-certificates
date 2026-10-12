"""VIR2 C0 gates.  Exit 0 iff all pass."""
import sys
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, b, rho

fails = []
def req(c, msg):
    print(('PASS  ' if c else 'FAIL  ') + msg, flush=True)
    if not c:
        fails.append(msg)

# g1: Virasoro algebra of T_phi with c = 1 - 12 rho^2
bad = n = 0
states = [dict(W.VAC2), {((), (1,)): sp.Integer(1)}, {((), (1, 1)): sp.Integer(1)}, {((), (2,)): sp.Integer(1)},
          {((), (1, 2)): sp.Integer(1)}, {((), (1, 1, 1)): sp.Integer(1)}, {((), (3,)): sp.Integer(1)}]
c = 1 - 12 * rho**2
for st in states:
    for m in range(-3, 4):
        for k in range(-3, 4):
            if m == k:
                continue
            lhs = W.Lphi(m, W.Lphi(k, st))
            for mm, cc in W.Lphi(k, W.Lphi(m, st)).items():
                W.add_to(lhs, mm, -cc)
            for mm, cc in W.Lphi(m + k, st).items():
                W.add_to(lhs, mm, -(m - k) * cc)
            if m + k == 0:
                for mm, cc in st.items():
                    W.add_to(lhs, mm, -c * sp.Rational(m**3 - m, 12) * cc)
            n += 1
            if W.clean(lhs):
                bad += 1
req(bad == 0, 'g1 Virasoro algebra of T_phi with c = 1 - 12 rho^2 on %d (state, m, n) triples' % n)

# g2: h(e^{b phi}) = b^2/2 - rho b, e^{a X} has h = a^2/2 under T_X: double pole of T(z)V(w)
Rv = {}
for (xs, fs), cc in W.T_PHI.items():
    # coefficient of (z-w)^-2 with nothing left: all letters contracted
    fac, D = 1, 0
    for m in fs:
        fac = fac * b * (-1) ** (m - 1) * sp.factorial(m - 1)
        D += m
    assert D == 2
    Rv[()] = Rv.get((), 0) + cc * fac
req(sp.expand(Rv[()] - (b**2 / 2 - rho * b)) == 0, 'g2 h(e^{b phi}) = b^2/2 - rho b')

# g3: the (T_phi T_phi) gate: (11)(11) - 2(22) -> 4 (T_phi T_phi) mod d
B4 = V.basis(4)
i4 = {m: i for i, m in enumerate(B4)}
vec = [F(0)] * len(B4)
vec[i4[((), ((1, 1), (1, 1)))]] = F(1)
vec[i4[((), ((2, 2),))]] = F(-2)
img, hv = W.to_two_boson(vec, 4)
req(hv == {((), (2, 2)): 4}, 'g3 (11)(11) - 2(22) = 4 L_-2 L_-2 |0> + total derivative: HV coefficients %s' % hv)
TT = W.vir_phi((2, 2))
diff = dict(img)
for m, cc in TT.items():
    W.add_to(diff, m, -4 * cc)
req(not W.clean(diff), 'g3 its two-boson image is exactly 4 (T_phi T_phi) = %s' % {k: str(v) for k, v in img.items()})
# the same statement without relying on the generator ordering: (11)(11) - 2(22) built from
# 4 L_-2^2 - 4 L_-4 (the v0 decomposition) differs from 4 (T_phi T_phi) by a total derivative
alt = {}
for m, cc in W.vir_phi((2, 2)).items():
    W.add_to(alt, m, 4 * cc)
for m, cc in W.vir_phi((4,)).items():
    W.add_to(alt, m, -4 * cc)
for m, cc in TT.items():
    W.add_to(alt, m, -4 * cc)
req(W.is_total_derivative(alt, 4), 'g3 4 L_-2^2 - 4 L_-4 and 4 L_-2^2 have the same two-boson class mod d')
req(not W.is_total_derivative(TT, 4), 'g3 (T_phi T_phi) itself is not a total derivative (the test is not vacuous)')

# g4: I_3 image: no odd-in-X term, form A4 x1^4 + A2' x2^2-class + C x1^2 T + 4 (TT)
for sol in (1, 2, 3):
    N, s = V.curve_point(F(7, 2))
    p4 = V.to_row(V.P4(sol, N, s), i4)
    img, hv = W.to_two_boson(p4, 4)
    odd = [m for m in img if len(m[0]) % 2]
    req(not odd and hv.get(((), (2, 2))) == 4, 'g4 Sol %d I_3 image even in X, (TT) coefficient 4; HV coefficients %s'
        % (sol, {k: str(v) for k, v in hv.items()}))
print('GATES', 'FAILED' if fails else 'PASSED')
sys.exit(1 if fails else 0)
