"""NUM2: the operators (sealed definitions).  three_term(t-data, l0, l1, E, ordering): [P(T) - x^(c/2) T x^(c/2) - E x^n] phi = 0;
sixth(l0, l1, EG): (theta-1)(theta-4) prod (theta - 5/2 -+ l0)(theta - 5/2 -+ l1) psi = x^6 (x^3 - EG) psi.  Exact E = 0 Q-functions (Meijer G)."""
import mpmath as mp
from fractions import Fraction as Fr
from ode import Op
def polymul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out
def fib(k):
    k = Fr(k); n = k + 4; b = k / (k + 1); c = n / b; s = (n - 1) / 2
    return dict(k=k, n=n, c=c, s=s, b=b)
def three_term(k, l0, l1, E, wrong_order=False):
    F = fib(k); s = mp.mpf(F['s'].numerator) / F['s'].denominator; c = F['c']
    cf = mp.mpf(c.numerator) / c.denominator
    # (s - theta)^2 - l^2 = theta^2 - 2 s theta + s^2 - l^2
    q0 = [s * s - l0 * l0, -2 * s, mp.mpf(1)]; q1 = [s * s - l1 * l1, -2 * s, mp.mpf(1)]
    P = polymul(q0, q1)
    shift = 0 if wrong_order else cf / 2           # x^(c/2) T x^(c/2) = x^c (T - c/2); wrong ordering x^c T
    R1 = [-(s - shift), mp.mpf(1)]                 # -x^c (s - shift - theta)
    terms = [(c, R1)]
    terms.append((F['n'], [-mp.mpf(E)]))
    thetas = [s - l0, s + l0, s - l1, s + l1]
    return Op(P, terms), thetas, F
def sixth(l0, l1, EG):
    ex = [mp.mpf(1), mp.mpf(4), mp.mpf(5) / 2 - l0, mp.mpf(5) / 2 + l0, mp.mpf(5) / 2 - l1, mp.mpf(5) / 2 + l1]
    P = [mp.mpf(1)]
    for e in ex: P = polymul(P, [-e, mp.mpf(1)])
    return Op(P, [(Fr(9), [mp.mpf(-1)]), (Fr(6), [mp.mpf(EG)])]), ex
def exact_three_term_E0(k, l0, l1):
    F = fib(k); s = mp.mpf(F['s'].numerator) / F['s'].denominator; c = mp.mpf(F['c'].numerator) / F['c'].denominator
    th = [s - l0, s + l0, s - l1, s + l1]; a = [x / c for x in th]; b1 = (s + c / 2) / c
    vt = (-1 + sum(a) - b1) / 3
    out = []
    for i in range(4):
        v = mp.mpf(1)
        for j in range(4):
            if j != i: v *= mp.gamma(a[j] - a[i])
        out.append(v / mp.gamma(b1 - a[i]) * c ** (-3 * a[i]) * mp.sqrt(3) / (2 * mp.pi) * c ** (3 * vt))
    return out, c * vt
def exact_sixth_E0(l0, l1):
    ex = [mp.mpf(1), mp.mpf(4), mp.mpf(5) / 2 - l0, mp.mpf(5) / 2 + l0, mp.mpf(5) / 2 - l1, mp.mpf(5) / 2 + l1]
    a = [e / 9 for e in ex]; vt = (mp.mpf(-5) / 2 + sum(a)) / 6
    out = []
    for i in range(6):
        v = mp.mpf(1)
        for j in range(6):
            if j != i: v *= mp.gamma(a[j] - a[i])
        out.append(v * mp.mpf(9) ** (-6 * a[i]) * mp.sqrt(6) / (2 * mp.pi) ** (mp.mpf(5) / 2) * mp.mpf(9) ** (6 * vt))
    return out, 9 * vt
