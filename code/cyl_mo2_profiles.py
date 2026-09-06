"""M-O2 certificate (2026-09-02; weight 4 added 2026-09-04): the gauge-invariant
all-first-order profiles of the ad_{I_3} kernels of the three cylindrical
solutions at weights 4, 6, 8 and 10 are the PINNED exact rational functions
of t below.

Setting.  On the curve N = (t^2-25)/(t^2-1), s = -24t/(t^2-1) the
weight-W kernel of ad_{I_3} (Solution SOL, hep-th/0404195 eqs. 53/55/57)
on Z2 x O(N)-invariant classes is one-dimensional away from the
certified jumping points.  The coefficients a_j of the all-first-order
monomials (dX)^{W-2j} (dY.dY)^j, j = 0..W/2, are class invariants (every
total derivative carries a letter of derivative order >= 2), normalized
here to a_0 = 1.  The pins were obtained by nullspace rational
interpolation on 18-22 points and verified on 4 hold-outs (lab
~/ib-lab/mo2_profiles.py); the weight-10 profiles of Solutions 1 and 3
were PREREGISTERED from the weight-6/8 closed forms before being
computed (preregistration commit 5ebfca2) and matched exactly.

This certificate recomputes every kernel at every listed point (26 per
table; 25 where a certified jumping point is excluded) and requires exact
equality with the pinned functions; in addition EVERY denominator root of
every pinned profile (except the curve's poles t = +-1 and the certified
jumping point t = 3 of Solution 2, where the kernel is two-dimensional)
is visited and the unnormalized (dX)^W coefficient required to vanish
there exactly.  The formulas are certified on these fibers; a symbolic
Q(t) kernel identity would be needed to call them global.  Structural facts
pinned as consequences: Solution 3 is palindromic at all four weights
with coupling u_W = C(W,2)(t-3)/(t+2W-7); Solution 1 obeys the single
closed form a_j = C(W,2j)(2j-1)!!(t-1)^j / prod_{i<=j}[(2i-1)t-(4W-1-6i)];
Solution 2 has t = 3 as a pole of order j-1 (j >= 2).  Runtime ~6 min.
"""
import sys, time
from fractions import Fraction as F
from math import comb
from mixed_engine import gen_mixed_basis, genuine_kernel_mixed, cyl_P4, canon_mixed, canon_pair

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)

def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), F(-24) * t / (t * t - 1)

def first_order(W):
    return [canon_mixed([1] * (W - 2 * j), tuple(canon_pair(1, 1) for _ in range(j))) for j in range(W // 2 + 1)]

R = F   # rational functions of t are given as python lambdas over Fraction
def dfact2(n):
    out = 1
    while n > 1:
        out *= n; n -= 2
    return out

def sol1_closed_form(W):
    def a(j, t):
        den = F(1)
        for i in range(1, j + 1):
            den *= (2 * i - 1) * t - (4 * W - 1 - 6 * i)
        return comb(W, 2 * j) * dfact2(2 * j - 1) * (t - 1) ** j / den
    return [lambda t, j=j: (F(1) if j == 0 else a(j, t)) for j in range(W // 2 + 1)]

def sol3_profile(W):
    u = lambda t: comb(W, 2) * (t - 3) / (t + 2 * W - 7)
    if W == 4:
        return [lambda t: F(1), u, lambda t: F(1)]
    if W == 6:
        return [lambda t: F(1), u, u, lambda t: F(1)]
    if W == 8:
        m = lambda t: 42 * (t - 3) * (5 * t - 11) / ((t + 9) * (3 * t - 1))
        return [lambda t: F(1), u, m, u, lambda t: F(1)]
    if W == 10:
        m = lambda t: 30 * (t - 3) * (7 * t - 17) / ((t + 1) * (t + 13))
        return [lambda t: F(1), u, m, m, u, lambda t: F(1)]
    raise ValueError(W)

def sol2_profile(W):
    if W == 4:   # weight 4 (2026-09-04): exact fit on 18 points, t = 3 (the certified jump) a simple pole
        return [lambda t: F(1),
                lambda t: 12 * (t - 1) / (5 * t - 11),
                lambda t: 8 * (t - 1) ** 2 / ((t - 3) * (5 * t - 11))]
    if W == 6:
        return [lambda t: F(1),
                lambda t: 15 * (t - 1) / (2 * (2 * t - 5)),
                lambda t: 10 * (t - 1) ** 2 / ((t - 3) * (2 * t - 5)),
                lambda t: 2 * (t - 1) ** 2 * (7 * t - 5) / (3 * (t - 3) ** 2 * (2 * t - 5))]
    if W == 8:
        return [lambda t: F(1),
                lambda t: 56 * (t - 1) / (11 * t - 29),
                lambda t: 112 * (t - 1) ** 2 / ((t - 3) * (11 * t - 29)),
                lambda t: 448 * (t - 1) ** 3 * (9 * t - 11) / (3 * (t - 3) ** 2 * (11 * t - 29) * (13 * t - 19)),
                lambda t: 32 * (t - 1) ** 3 * (5 * t - 3) * (9 * t - 11) / (3 * (t - 3) ** 3 * (11 * t - 29) * (13 * t - 19))]
    if W == 10:
        return [lambda t: F(1),
                lambda t: 45 * (t - 1) / (7 * t - 19),
                lambda t: 120 * (t - 1) ** 2 / ((t - 3) * (7 * t - 19)),
                lambda t: 60 * (t - 1) ** 3 * (11 * t - 17) / ((t - 3) ** 2 * (4 * t - 7) * (7 * t - 19)),
                lambda t: 720 * (t - 1) ** 5 * (11 * t - 17) / ((t - 3) ** 3 * (4 * t - 7) * (7 * t - 19) * (17 * t - 23)),
                lambda t: 16 * (t - 1) ** 5 * (11 * t - 17) * (13 * t - 7) / ((t - 3) ** 4 * (4 * t - 7) * (7 * t - 19) * (17 * t - 23))]
    raise ValueError(W)

PINS = {1: sol1_closed_form, 2: sol2_profile, 3: sol3_profile}
# sample points: generic (t = +-1 and the certified jumping points excluded per (sol, W)); where a
# pinned profile has a pole (t = 19/3 at weight 8, t = 9 at weight 10 for Solution 1) the (dX)^W
# coefficient is required to vanish -- the pole is verified, not skipped
POINTS = [F(a, b) for a, b in [(2,1),(7,5),(3,2),(9,7),(11,9),(4,1),(6,1),(7,1),(8,1),(9,1),(11,2),(13,3),(-2,1),(-7,5),(-4,1),(-6,1),(-9,2),(17,5),(19,3),(-13,4),(25,3),(-15,2),(21,11),(10,1),(-11,3),(23,7)]]
EXCLUDE = {(3, 6): {F(5, 3)}, (3, 10): {F(7, 5)}}   # certified jumping points (kernel dim 2)
# every denominator root of the pinned profiles (t = +-1 and Solution 2's certified jump t = 3 excluded)
POLES = {(1, 4): [F(9)], (2, 4): [F(11, 5)], (3, 4): [],   # weight 4: t = 1 (Solution 1) and t = -1 (Solution 3) are the curve's poles
         (1, 6): [F(17), F(11, 3)], (1, 8): [F(25), F(19, 3), F(13, 5)], (1, 10): [F(33), F(9), F(21, 5), F(15, 7)],
         (3, 6): [F(-5)], (3, 8): [F(-9), F(1, 3)], (3, 10): [F(-13)],
         (2, 6): [F(5, 2)], (2, 8): [F(29, 11), F(19, 13)], (2, 10): [F(19, 7), F(7, 4), F(23, 17)]}

if __name__ == "__main__":
    t0 = time.time()
    for W in (4, 6, 8, 10):
        BASIS = gen_mixed_basis(W)
        MONOS = first_order(W)
        for sol in (1, 2, 3):
            pins = PINS[sol](W)
            pts = [t for t in POINTS if t not in EXCLUDE.get((sol, W), set())]
            bad, poles_seen = [], []
            for t in pts:
                N, s = point(t)
                ker = genuine_kernel_mixed(cyl_P4(sol, N, s), BASIS, N)
                if len(ker) != 1:
                    bad.append(f"t={t}: kernel dim {len(ker)}"); continue
                v = ker[0]
                top = v.get(MONOS[0], F(0))
                try:
                    exp = [p(t) for p in pins]
                except ZeroDivisionError:
                    # a pole of the pinned profile: the (dX)^W coefficient must vanish there
                    if top != 0:
                        bad.append(f"t={t}: pinned pole but (dX)^W coefficient {top} != 0")
                    else:
                        poles_seen.append(t)
                    continue
                if top == 0:
                    bad.append(f"t={t}: zero (dX)^W coefficient off any pinned pole"); continue
                prof = [v.get(m, F(0)) / top for m in MONOS]
                if prof != exp:
                    bad.append(f"t={t}: profile {prof} != pinned {exp}")
            require(not bad, f"P1: Solution {sol}, weight {W}: pinned profile reproduced exactly at {len(pts)} points"
                    + (f" (pinned poles met at sample points t = {[str(x) for x in poles_seen]})" if poles_seen else "")
                    + ("" if not bad else f" -- {bad[:2]}"))
            # P3: every denominator root of the pinned profile is a zero of the unnormalized (dX)^W coefficient
            badp = []
            for t in POLES[(sol, W)]:
                pins_t = PINS[sol](W)
                try:
                    [pp(t) for pp in pins_t]
                    badp.append(f"t={t}: not a pole of the pinned profile"); continue
                except ZeroDivisionError:
                    pass
                N, s = point(t)
                ker = genuine_kernel_mixed(cyl_P4(sol, N, s), BASIS, N)
                if len(ker) != 1:
                    badp.append(f"t={t}: kernel dim {len(ker)}"); continue
                if ker[0].get(MONOS[0], F(0)) != 0:
                    badp.append(f"t={t}: (dX)^W coefficient {ker[0].get(MONOS[0])} != 0")
            require(not badp, f"P3: Solution {sol}, weight {W}: all {len(POLES[(sol, W)])} pinned poles are exact zeros of the (dX)^W coefficient"
                    + ("" if not badp else f" -- {badp[:2]}"))
        # structural consequences, checked symbolically at the pin level
        u = lambda t, W=W: comb(W, 2) * (t - 3) / (t + 2 * W - 7)
        p3 = PINS[3](W)
        require(all(p3[j](F(k, 7)) == p3[W // 2 - j](F(k, 7)) for j in range(W // 2 + 1) for k in (13, 29, -37)),
                f"P2: Solution 3 palindromic at weight {W}")
        require(all(p3[1](F(k, 7)) == u(F(k, 7)) for k in (13, 29, -37)), f"P2: Solution 3 coupling u_{W} = C({W},2)(t-3)/(t+{2*W-7})")
    print(f"\n{time.time()-t0:.0f}s")
    if fails:
        print(f"MO2 PROFILES: {len(fails)} PIN(S) NOT MET"); sys.exit(1)
    print("MO2 PROFILES CERTIFIED (twelve exact invariant profiles at weights 4, 6, 8, 10; Solution 3 palindromic with "
          "u_W = C(W,2)(t-3)/(t+2W-7); Solution 1 one closed form with Kac-collision denominators; Solution 2 with "
          "t = 3 a pole of order j-1; weight-10 forms of Solutions 1 and 3 preregistered and matched)")
