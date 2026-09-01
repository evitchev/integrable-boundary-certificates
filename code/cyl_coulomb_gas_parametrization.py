"""The cylindrical curve is the Coulomb-gas parametrization of a
Virasoro central charge (N = c_Y; the theory's full free-boson charge
is N + 1), and the detected/predicted jumping points are c_{1,q}.

Exact (sympy rationals), seconds.  Fail-closed.

(1) With x = (t-1)/(t+1):  N(t) = 13 - 6 (x + 1/x),  s(t) = 6 (x - 1/x),
    identically in t; hence s^2 = (N-25)(N-1) and the sheet exchange
    s -> -s is x -> 1/x (beta -> 1/beta).  So rational t <=> rational
    beta^2 = x = p/q, where N = c_{p,q} = 1 - 6 (p-q)^2 / (p q).
(2) Every detected or predicted jumping point is a (1,q) central charge
    with q EVEN:
    t_k = (2k+3)/(2k+1) = (q+1)/(q-1) with q = 2k+2, N(t_k) = c_{1,q},
    k = 1..5 (q = 4, 6, 8, 10, 12), and the all-spin decoupling point
    t = 3 is c_{1,2} = -2.  In terms of q the certified/predicted jump
    spins are 2q-3 on sheet 3 and (3q-2)/2 on sheet 2 (integer, odd, iff
    q = 2 mod 4).
(3) The odd-q points t = (q+1)/(q-1) for q = 3, 5, 7, 9, 11 (t = 2, 3/2,
    4/3, 5/4, 6/5) are rational points of the curve too -- listed here as
    the points that the weight-14 symbolic stratification and the
    weight-16/18 exhaustive residue scans (both primes) did NOT flag.
    This script records the arithmetic only; silence at those points
    is the cited scans' statement, not this script's.
Marker: COULOMB-GAS PARAMETRIZATION CERTIFIED."""
import sys

import sympy as sp

fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)


t, x = sp.symbols("t x")
N = (t**2 - 25) / (t**2 - 1)
S = -24 * t / (t**2 - 1)
X = (t - 1) / (t + 1)
check(sp.simplify(N - (13 - 6 * (X + 1 / X))) == 0, "N(t) = 13 - 6(x + 1/x), x = (t-1)/(t+1)")
check(sp.simplify(S - 6 * (X - 1 / X)) == 0, "s(t) = 6(x - 1/x)")
check(sp.simplify(S**2 - (N - 25) * (N - 1)) == 0, "s^2 = (N-25)(N-1) recovered")
t2 = sp.Symbol("t2")
sol = sp.solve(sp.Eq((t2 - 1) / (t2 + 1), 1 / X), t2)  # t2 with x(t2) = 1/x(t)
check(len(sol) == 1 and sp.simplify(sol[0] + t) == 0 and sp.simplify(N.subs(t, -t) - N) == 0
      and sp.simplify(S.subs(t, -t) + S) == 0,
      "sheet exchange x -> 1/x (beta -> 1/beta) is t -> -t: N fixed, s -> -s")


def c_pq(p, q):
    return sp.Integer(1) - sp.Integer(6) * (p - q) ** 2 / (p * q)


named_points = {  # k: (t_k, (N, s)) -- the detected/predicted points of the exact record
    1: (sp.Rational(5, 3), (sp.Rational(-25, 2), sp.Rational(-45, 2))),
    2: (sp.Rational(7, 5), (sp.Integer(-24), sp.Integer(-35))),
    3: (sp.Rational(9, 7), (sp.Rational(-143, 4), sp.Rational(-189, 4))),
    4: (sp.Rational(11, 9), (sp.Rational(-238, 5), sp.Rational(-297, 5))),
    5: (sp.Rational(13, 11), (sp.Rational(-119, 2), sp.Rational(-143, 2))),
}
for k, (tk, (Nk, sk)) in named_points.items():
    q = 2 * k + 2
    check(tk == sp.Rational(q + 1, q - 1), f"t_{k} = (q+1)/(q-1) with q = {q}")
    check(N.subs(t, tk) == Nk and S.subs(t, tk) == sk, f"(N,s)(t_{k}) = ({Nk}, {sk}) as recorded")
    check(Nk == c_pq(1, q), f"N(t_{k}) = c_(1,{q}) = {c_pq(1, q)}")
    check(X.subs(t, tk) == sp.Rational(1, q), f"beta^2(t_{k}) = 1/{q}")
    check(4 * k + 1 == 2 * q - 3, f"sheet-3 spin 4k+1 = 2q-3 = {2*q-3}")
    check((3 * k + 2) * 2 == 3 * q - 2, f"sheet-2 spin 3k+2 = (3q-2)/2 = {sp.Rational(3*q-2, 2)}")
check(N.subs(t, 3) == c_pq(1, 2) == -2 and S.subs(t, 3) == -9, "t = 3 (all-spin decoupling) is c_(1,2) = -2")
check(sp.simplify(N.subs(t, (x + 1) / (x - 1)) - (13 - 6 * (1 / x + x))) == 0,
      "inverse map t = (q+1)/(q-1) <-> x = 1/q, symbolic")
odd = {q: sp.Rational(q + 1, q - 1) for q in (3, 5, 7, 9, 11)}
for q, tq in odd.items():
    check(N.subs(t, tq) == c_pq(1, q), f"odd q = {q}: t = {tq} is rational with N = c_(1,{q}) = {c_pq(1, q)} "
          "(arithmetic only; scan silence is the cited scans' statement, not this script's)")
kk = sp.Symbol('k')
gen = c_pq(1, 2 * kk + 2)
check(sp.simplify(gen + (4 * kk + 1) * (3 * kk + 2) / (kk + 1)) == 0,
      "c_(1,2k+2) = -(4k+1)(3k+2)/(k+1): the jump spins are the factors of the numerator of -c")
print()
if fails:
    print(f"COULOMB-GAS PARAMETRIZATION: {len(fails)} FAILURE(S)")
    sys.exit(1)
print("COULOMB-GAS PARAMETRIZATION CERTIFIED (N = 13 - 6(x + 1/x), s = 6(x - 1/x), x = (t-1)/(t+1); "
      "t_1..t_5 = c_(1,4), c_(1,6), c_(1,8), c_(1,10), c_(1,12); t = 3 = c_(1,2))")
