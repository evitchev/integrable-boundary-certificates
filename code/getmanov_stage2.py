"""Getmanov stage 2: encode the cigar W-infinity data of
Kotousov-Lukyanov (1910.05947) and test the identification of the
enhanced c=1 tower.

Pipeline (each stage validated against KL's stated OPEs before the
next is trusted):
 1. convention map: their (phi, theta) with <dphi dphi> = -1/2 u^-2
    onto repo bosons X1, X2 (<dX dX> = +u^-2) via dphi = p dX1,
    dtheta = t dX2 with p = t = i/sqrt2 (solved from W2W2 pole-2
    = 2 W2); background B = -p/sqrt(k); central charge c(k) computed
    exactly from the pole-4 term.
 2. W3 transcribed (KL eq. Weq1); validated: W2W3 has no pole-4/3,
    pole-2 = 3 W3, pole-1 = dW3.
 3. W4 SOLVED from the W3W3 pole-2 using KL's OPE
    (W4 = [pole2 - 2:W2^2: + (2 mu - C4/2) d^2 W2] / 2), after
    validating pole-6 = C6, pole-5 = 0, pole-4 = 2 C4 W2,
    pole-3 = C4 dW2 with C6, C4, mu the stated c-dependent constants.
 4. W2W4 OPE check: pole-4 mixing = (c+10)(17c+2)/(15(c-2)) W2,
    pole-2 = 4 W4, pole-1 = dW4.
 5. I3^G ~ (k+6) W4 + 4(2k+1) :W2^2:; validated by its OWN screening
    data at a numeric k: the Getmanov vertex operators (the
    cos-pair with momenta (-i sqrt(k/2), -+sqrt((k+2)/2)) in repo
    variables -- the paperclip form with n = k! -- and the mu'
    vertex (i 4/sqrt(2k), 0)) must annihilate I3^G.
 6. IDENTIFICATION TESTS against the enhanced c=1 tower
    fingerprint (one boson, momenta {+-1/sqrt2, +-2 sqrt2}).
"""
import sys

import sympy as sp

from ope_poles import ope_pole
from screening_fit import screening_conditions
from super_engine import cur_add, normal_prod

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def scur(pairs):
    """current from {super-mono: sympy expr}: canonicalize keys
    (sorted letters; bosonic, no signs) and drop zeros."""
    out = {}
    for m, c in pairs.items():
        key = tuple(sorted(m))
        out[key] = sp.expand(out.get(key, 0) + c)
    return {m: c for m, c in out.items() if c != 0}


def add(A, B, coeff=1):
    out = dict(A)
    for m, c in B.items():
        out[m] = sp.expand(out.get(m, 0) + coeff * c)
    return {m: c for m, c in out.items() if c != 0}


def dcur(A):
    from super_engine import d_cur
    return {m: sp.expand(c) for m, c in d_cur(A).items()
            if sp.expand(c) != 0}


def bos(cur):
    """super-engine monomials -> ope_engine (m, j) monomials."""
    out = {}
    for mono, c in cur.items():
        out[tuple((m, j) for _, m, j in mono)] = c
    return out


def X(m, j):
    return (0, m, j)


rk = sp.Symbol('rk', positive=True)          # sqrt(k)
k = rk**2
p = sp.I / sp.sqrt(2)                        # dphi = p dX1
t = sp.I / sp.sqrt(2)                        # dtheta = t dX2

# ---- stage 1: W2 ---------------------------------------------------
W2 = scur({(X(1, 0), X(1, 0)): -p**2,
           (X(1, 1), X(1, 1)): -t**2,
           (X(2, 0),): -p / rk})

p2 = ope_pole(W2, W2, 2)
check("W2 W2 pole-2 = 2 W2",
      scur(p2) == scur({m: 2 * c for m, c in W2.items()}))
check("W2 W2 pole-3 = 0", not scur(ope_pole(W2, W2, 3)))
p4 = ope_pole(W2, W2, 4)
cval = sp.simplify(2 * p4.get((), 0))
print(f"  central charge c(k) = {sp.factor(cval)}", flush=True)
check("c(k) = 2 + 6/k", sp.simplify(cval - 2 - 6 / k) == 0)

# ---- stage 2: W3 (KL eq. Weq1) -------------------------------------
W3 = scur({(X(1, 1), X(1, 1), X(1, 1)): -sp.I * (6 * k + 4) / (3 * k) * t**3,
           (X(1, 0), X(1, 0), X(1, 1)): -sp.I * 2 * p**2 * t,
           (X(2, 0), X(1, 1)): sp.I * rk * p * t,
           (X(1, 0), X(2, 1)): -sp.I * (k + 2) / rk * p * t,
           (X(3, 1),): -sp.I * (k + 2) / (6 * k) * t})

check("W2 W3 pole-4 = 0", not scur(ope_pole(W2, W3, 4)))
check("W2 W3 pole-3 = 0", not scur(ope_pole(W2, W3, 3)))
check("W2 W3 pole-2 = 3 W3",
      scur(ope_pole(W2, W3, 2))
      == scur({m: 3 * c for m, c in W3.items()}))
check("W2 W3 pole-1 = dW3",
      scur(ope_pole(W2, W3, 1)) == dcur(W3))

# ---- stage 3: W3 W3 -> W4 ------------------------------------------
C6 = cval * (cval + 7) * (2 * cval - 1) / (9 * (cval - 2))
C4 = (cval + 7) * (2 * cval - 1) / (3 * (cval - 2))
mu = (2 * cval**2 + 22 * cval - 25) / (30 * (cval - 2))

p6 = ope_pole(W3, W3, 6)
check("W3 W3 pole-6 = C6",
      sp.simplify(p6.get((), 0) - C6) == 0 and len(scur(p6)) <= 1)
check("W3 W3 pole-5 = 0", not scur(ope_pole(W3, W3, 5)))
check("W3 W3 pole-4 = 2 C4 W2",
      scur(ope_pole(W3, W3, 4))
      == scur({m: 2 * C4 * c for m, c in W2.items()}))
check("W3 W3 pole-3 = C4 dW2",
      scur(ope_pole(W3, W3, 3))
      == scur({m: C4 * c for m, c in dcur(W2).items()}))

W2sq = {m: sp.expand(c) for m, c in normal_prod(W2, W2).items()
        if sp.expand(c) != 0}
d2W2 = dcur(dcur(W2))
pole2 = ope_pole(W3, W3, 2)
W4 = add(add(pole2, W2sq, -2), d2W2, 2 * mu - C4 / 2)
W4 = {m: sp.expand(c / 2) for m, c in W4.items()
      if sp.expand(c) != 0}

# ---- stage 4: W2 W4 checks -----------------------------------------
mix = (cval + 10) * (17 * cval + 2) / (15 * (cval - 2))
q4 = ope_pole(W2, W4, 4)
check("W2 W4 pole-4 = (c+10)(17c+2)/(15(c-2)) W2",
      scur(q4) == scur({m: mix * c for m, c in W2.items()}))
check("W2 W4 pole-2 = 4 W4",
      scur(ope_pole(W2, W4, 2))
      == scur({m: 4 * c for m, c in W4.items()}))
check("W2 W4 pole-1 = dW4",
      scur(ope_pole(W2, W4, 1)) == dcur(W4))

# ---- stage 5: I3^G and its own screenings --------------------------
I3G = add({m: (k + 6) * c for m, c in W4.items()},
           W2sq, 4 * (2 * k + 1))
I3G = {m: sp.expand(c) for m, c in I3G.items() if sp.expand(c) != 0}

kv = sp.Rational(1, 4)                       # numeric validation point
I3Gk = {m: sp.nsimplify(sp.simplify(c.subs(rk, sp.sqrt(kv))))
        for m, c in I3G.items()}
I3Gk = bos({m: c for m, c in I3Gk.items() if c != 0})
alph_cos_p = (-sp.I * sp.sqrt(kv / 2), sp.sqrt((kv + 2) / 2))
alph_cos_m = (-sp.I * sp.sqrt(kv / 2), -sp.sqrt((kv + 2) / 2))
alph_mup = (sp.I * 4 / sp.sqrt(2 * kv), sp.Integer(0))
for name, al in (("cos+", alph_cos_p), ("cos-", alph_cos_m),
                 ("mu'", alph_mup)):
    conds = screening_conditions(I3Gk, al, 2)
    res = [x for x in conds if sp.simplify(x) != 0]
    check(f"Getmanov screening {name} annihilates I3^G at k=1/4",
          not res)

print()
if failures:
    print(f"FAILURES: {failures}")
    print("GETMANOV STAGE-2 TRANSCRIPTION FAILED")
    sys.exit(1)
print("GETMANOV TRANSCRIPTION VALIDATED (W2, W3, W4, I3^G, own "
      "screenings at k=1/4)")

# ---- stage 6: identification tests (REWRITTEN after audit) ---------
# A first version tested whether the FULL two-boson charge commutes
# with axis-supported vertex operators and wrongly concluded
# exclusion.  The identification question is about the AXIS PULLBACK:
# keep only the monomials supported on one boson, compare the
# resulting one-boson charge class with I3' mod d.  Lift structure
# (audit round 2, certified below): at each matching point exactly
# ONE Liouville-dual pair of the fingerprint lifts to the full
# two-boson charge, branch-correlated with sqrt(k); the
# complementary pair is emergent, and nothing lifts on theta.

kk = sp.Symbol('kk')          # unrestricted k (rk positivity would
                              # suppress the k<0 matches: audit catch)


def pullback(cur, axis):
    out = {}
    for mono, c in cur.items():
        if all(j == axis for _, m, j in mono):
            key = tuple(sorted(m for _, m, j in mono))
            out[key] = sp.expand(out.get(key, 0) + c)
    return {m: c for m, c in out.items() if c != 0}


def ratio_classes(pb):
    """weight-4 one-boson classes mod d ((112)=0, (4)=0, (13)=-(22)):
    returns the (22)-class / (1111)-class ratio."""
    quart = pb.get((1, 1, 1, 1), 0)
    c22 = pb.get((2, 2), 0) - pb.get((1, 3), 0)
    return sp.factor(sp.cancel(c22 / quart))


I3Gsym = add({m: (k + 6) * c for m, c in W4.items()},
             W2sq, 4 * (2 * k + 1))
I3Gsym = {m: sp.expand(c) for m, c in I3Gsym.items()
          if sp.expand(c) != 0}

Rphi = ratio_classes(pullback(I3Gsym, 0)).subs(rk**2, kk)
Rtheta = sp.factor(sp.cancel(
    ratio_classes(pullback(I3Gsym, 1)).subs(rk**2, kk)))
check("phi-axis ratio = -(k^2+12k+16)/(2k)",
      sp.simplify(Rphi + (kk**2 + 12 * kk + 16) / (2 * kk)) == 0)
check("theta-axis ratio = -(k^3+13k^2+20k+12)/(2(k+2)(2k+3))",
      sp.simplify(Rtheta + (kk**3 + 13 * kk**2 + 20 * kk + 12)
                  / (2 * (kk + 2) * (2 * kk + 3))) == 0)

phi_match = sp.factor(sp.numer(sp.together(Rphi - sp.Rational(5, 2))))
check("phi-axis matching locus = (k+1)(k+16)",
      sp.simplify(phi_match / ((kk + 1) * (kk + 16))).is_constant())
theta_match = sp.factor(sp.numer(sp.together(Rtheta
                                             - sp.Rational(5, 2))))
check("theta-axis matching locus = k^3+23k^2+55k+42",
      sp.simplify(theta_match
                  / (kk**3 + 23 * kk**2 + 55 * kk + 42)).is_constant())

# exact ray matches at k = -1 and k = -16 (phi axis): the reduced
# charge equals I3' = (1111) - 5/2 (13) modulo d and normalization
for kval, rkval in ((-1, sp.I), (-16, 4 * sp.I)):
    I3Gv = {m: sp.nsimplify(sp.simplify(c.subs(rk, rkval)))
            for m, c in I3Gsym.items()}
    pb = pullback({m: c for m, c in I3Gv.items() if c != 0}, 0)
    quart = pb.get((1, 1, 1, 1), 0)
    c22 = pb.get((2, 2), 0) - pb.get((1, 3), 0)
    check(f"k={kval}: phi-pullback class = I3' ray (ratio 5/2)",
          quart != 0 and sp.simplify(c22 / quart
                                     - sp.Rational(5, 2)) == 0)

# lift structure at the matching points (audit round 2): exactly one
# dual pair screens the FULL two-boson charge, selected by the branch
# of sqrt(k); nothing lifts on theta
r2 = sp.sqrt(2)


def full_bos(rkval):
    cur = {m: sp.nsimplify(sp.simplify(c.subs(rk, rkval)))
           for m, c in I3Gsym.items()}
    return bos({m: c for m, c in cur.items() if c != 0})


def screens_full(I3b, al):
    conds = screening_conditions(I3b, al, 2)
    return not [x for x in conds if sp.simplify(x) != 0]


LIFT = {(-1, sp.I): (-1, 1), (-1, -sp.I): (1, -1),
        (-16, 4 * sp.I): (1, -1), (-16, -4 * sp.I): (-1, 1)}
for (kval, rkval), (s_half, s_big) in LIFT.items():
    I3b = full_bos(rkval)
    ok = (screens_full(I3b, (s_half / r2, sp.Integer(0)))
          and screens_full(I3b, (s_big * 2 * r2, sp.Integer(0)))
          and not screens_full(I3b, (-s_half / r2, sp.Integer(0)))
          and not screens_full(I3b, (-s_big * 2 * r2, sp.Integer(0))))
    check(f"k={kval}, rk={rkval}: exactly the dual pair "
          f"({s_half}/sqrt2, {s_big}*2sqrt2) lifts on phi", ok)
    th = any(screens_full(I3b, (sp.Integer(0), al))
             for al in (1 / r2, -1 / r2, 2 * r2, -2 * r2))
    check(f"k={kval}, rk={rkval}: no fingerprint momentum lifts on "
          f"theta", not th)

print()
if failures:
    print(f"FAILURES: {failures}")
    print("GETMANOV STAGE-2 FAILED")
    sys.exit(1)
print("GETMANOV AXIS-PULLBACK MATCHES FOUND (audit reversal): the "
      "phi-axis reduction of I3^G equals the enhanced I3' ray at "
      "k = -1 and k = -16, and the theta-axis matching locus is the "
      "cubic k^3+23k^2+55k+42 (real root ~ -20.41).  Getmanov is a "
      "FORMAL ANALYTIC-CONTINUATION CANDIDATE at the I3-ray level "
      "(k=-1,-16 lie outside the positive-level cigar regime, and an "
      "isolated ratio match in a 2-dim quotient is expected of a "
      "one-parameter family -- the I5 comparison is decisive).  Lift "
      "structure: exactly one branch-correlated dual pair of the "
      "fingerprint lifts to the full two-boson charge; the "
      "complementary pair is emergent; nothing lifts on theta.")
print("GETMANOV STAGE-2 COMPLETE")
