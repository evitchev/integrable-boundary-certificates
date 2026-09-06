"""Pole-lattice certificate (2026-09-04): the certified profile poles of the three
cylindrical solutions are Beta-function lattices of a leading symbol
kappa^2 (1+u)^{m(t)} with ONE Moebius exponent function m(t) per solution.

Setting.  For the paperclip ODE of hep-th/0404195 eq. (62), the WKB
coefficient S_s is a sum of Beta functions whose poles in the exponent n
sit on the lattice  s n + 2j = 0,  j = 1..s-1  (the factors (2j + s n) of
eqs. (66)-(68)); the n-poles of the vacuum eigenvalues I_s are a subset of
that lattice.  Reading the certified cylindrical profile poles the same
way, the PINNED exponent functions (anchored at t -> infinity to the
paperclip values the paper states for N = 1: n = -1 for Solutions 1 and 3,
n = 2 for Solution 2) are

    m_1(t) = (1 - t)/(t + 3),   m_2(t) = (2t - 14)/(t + 5),   m_3(t) = (3 - t)/(t - 5),

each defined up to the lattice symmetry j -> s - j, i.e. m -> -2 - m.

Checks (all exact, Fraction arithmetic):
  P1  every certified pole t_p of the weight-W profile of Solution k (POLES
      of cyl_mo2_profiles, spin s = W - 1, weights 4, 6, 8, 10) satisfies
      s m_k(t_p) + 2 j = 0 with an INTEGER j in 1..s-1 -- 10 poles for
      Solution 1, 7 for Solution 2, 5 for Solution 3 (the t = -1 pole of the
      weight-10 profile of Solution 3, also the curve's pole, included);
  P2  the anchors: m_1, m_3 -> -1 and m_2 -> 2 as t -> infinity;
  P3  Solution 1's lattice index is j = (s+1)/2 - i for its i-th denominator
      factor (2i-1)t - (4s+3-6i), exactly: s m_1(t) + 2j = -[(2i-1)t - (4s+3-6i)]/(t+3);
  P4  Solution 3's certified jumping points t = (q+1)/(q-1), q = 4, 6, 8
      (weights 6, 10, 14) lie on the lattice with j = (q-2)/2 at spin
      s = 2q - 3, and no c_{1,q} (q = 2..12 even) is a lattice point at spin 19;
  P5  negative controls: shifting any pole by 1/1000, or replacing an
      exponent by its value plus 1/7, breaks P1 (the guard is live).
Fail-closed: any failure exits 1 without the CERTIFIED marker.  Runtime < 1 s.
"""
import sys
from fractions import Fraction as F
from cyl_mo2_profiles import POLES

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)

M = {1: lambda t: (1 - t) / (t + 3), 2: lambda t: (2 * t - 14) / (t + 5), 3: lambda t: (3 - t) / (t - 5)}
ANCHOR = {1: F(-1), 2: F(2), 3: F(-1)}
EXTRA = {(3, 10): [F(-1)]}      # the (t+1) pole of Solution 3's weight-10 middle coefficient (also the curve's pole)

def lattice_index(m_val, s):
    j = -m_val * s / 2
    return int(j) if j.denominator == 1 and 1 <= j <= s - 1 else None

def all_poles(sol):
    out = []
    for (so, W), ts in sorted(POLES.items()):
        if so == sol:
            out += [(W - 1, t) for t in ts] + [(W - 1, t) for t in EXTRA.get((so, W), [])]
    return out

def check_lattice(sol, m_fun, poles, label):
    bad = []
    for s, t in poles:
        try:
            j = lattice_index(m_fun(t), s)
        except ZeroDivisionError:
            j = None
        if j is None:
            bad.append((s, str(t)))
    return bad

if __name__ == "__main__":
    counts = {}
    for sol in (1, 2, 3):
        poles = all_poles(sol); counts[sol] = len(poles)
        bad = check_lattice(sol, M[sol], poles, f"Solution {sol}")
        require(not bad, f"P1: Solution {sol}: all {len(poles)} certified profile poles (weights 4-10) lie on the lattice s m_{sol}(t) + 2j = 0 with integer j in 1..s-1"
                + ("" if not bad else f" -- off-lattice: {bad}"))
    require(counts == {1: 10, 2: 7, 3: 5}, f"P1: pole counts per solution are 10, 7, 5 (found {counts})")
    # exact limits: each m_k is a Moebius function (a t + b)/(c t + d); its limit is a/c.  The coefficients are read off
    # EXACTLY from the pinned lambdas by two evaluations of numerator and denominator behaviour: m(t)(t + 1/2) ... instead,
    # verify the closed forms against the lambdas at four points, then take a/c from the closed forms.
    MOEB = {1: (-1, 1, 1, 3), 2: (2, -14, 1, 5), 3: (-1, 3, 1, -5)}
    ok2 = all(M[k](F(tv)) == F(a * tv + b_, c * tv + d) for k, (a, b_, c, d) in MOEB.items() for tv in (2, F(7, 3), -4, F(31, 5)))
    ok2 &= all(F(a, c) == ANCHOR[k] for k, (a, b_, c, d) in MOEB.items())
    require(ok2, "P2: anchors -- the pinned Moebius exponents equal (a t + b)/(c t + d) with (a, b, c, d) = (-1, 1, 1, 3), (2, -14, 1, 5), (-1, 3, 1, -5) at four points, and their EXACT limits a/c are -1, 2, -1 (paperclip values at N = 1)")
    # P3: Solution 1's closed-form denominators are the lattice factors
    ok3 = True
    for s in (3, 5, 7, 9):
        for i in range(1, (s - 1) // 2 + 1):
            j = (s + 1) // 2 - i
            for t in (F(2), F(7, 3), F(-11, 4), F(31, 5)):
                lhs = s * M[1](t) + 2 * j
                rhs = -((2 * i - 1) * t - (4 * s + 3 - 6 * i)) / (t + 3)
                ok3 &= (lhs == rhs)
    require(ok3, "P3: s m_1(t) + 2j = -[(2i-1)t - (4s+3-6i)]/(t+3) with j = (s+1)/2 - i, spins 3-9, four t each (Solution 1's denominators ARE lattice factors)")
    # P4: Solution 3's jumping points and the spin-19 statement
    ok4 = True
    for q in (4, 6, 8):
        t = F(q + 1, q - 1); s = 2 * q - 3
        ok4 &= (lattice_index(M[3](t), s) == (q - 2) // 2)
    require(ok4, "P4: Solution 3's jumping points t = (q+1)/(q-1), q = 4, 6, 8, are lattice points with j = (q-2)/2 at spin 2q-3 (weights 6, 10, 14)")
    none19 = all(lattice_index(M[3](F(q + 1, q - 1)), 19) is None for q in (2, 4, 6, 8, 10, 12))
    require(none19, "P4: no c_{1,q} (q = 2..12 even) is a lattice point of Solution 3 at spin 19 (weight 20)")
    # P5: negative controls
    drills_ok = True
    for sol in (1, 2, 3):
        poles = all_poles(sol)
        s0, t0 = poles[0]
        shifted = [(s0, t0 + F(1, 1000))] + poles[1:]
        drills_ok &= bool(check_lattice(sol, M[sol], shifted, "drill"))
        drills_ok &= bool(check_lattice(sol, lambda t, k=sol: M[k](t) + F(1, 7), poles, "drill"))
    require(drills_ok, "P5: negative controls -- a pole shifted by 1/1000 and an exponent shifted by 1/7 both leave the lattice (guard live)")
    if fails:
        print(f"POLE LATTICE: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("POLE LATTICE CERTIFIED (the certified profile poles of Solutions 1, 2, 3 at weights 4-10 are the Beta lattices "
          "s m(t) + 2j = 0 of the pinned Moebius exponents m_1 = (1-t)/(t+3), m_2 = (2t-14)/(t+5), m_3 = (3-t)/(t-5), "
          "anchored at the paperclip values; Solution 3's jumping points are lattice points at weight 2q-2)")
