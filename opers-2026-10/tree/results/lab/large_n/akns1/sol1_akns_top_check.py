"""Sol 1's F1 top layer (record item 178) vs Fateev-Lukyanov's printed top layer of the quantum AKNS
vacuum eigenvalues (hep-th/0510271, the display following their WKB integral; s = 2K-1, 2l + m = s + 1),
with n = -p^2/2 = -(t-1)/(t+1).  Prints the number of distinct ratio scale factors per assignment
(1 = agreement up to normalisation; nan = a Gamma pole at that (t, K))."""
import sympy as sp
def rf(x, n): return sp.prod([x + i for i in range(n)]) if n > 0 else 1
def ib_top(t, K):
    aY = (t - 1)/(t + 3); aX = 1; p = -aY*(K - sp.Rational(1, 2)); q = -aX*(K - sp.Rational(1, 2)); z = aY/aX
    return [sp.binomial(K, j)*rf(p, K - j)*rf(q, j)*z**j/rf(p, K) for j in range(K + 1)]
def fl_top(n, K):
    s = 2*K - 1
    return {l: 2**(s + 1 - 2*l)*(n + 2)**sp.Rational(s + 1 - 2*l, 2)*sp.gamma(sp.Rational(1, 2) - l + (n + 1)*s/(n + 2))
               / (n**(s - l)*sp.gamma(-s/(n + 2))*sp.factorial(s + 1 - 2*l)*sp.factorial(l)) for l in range(K + 1)}
for t in (sp.Rational(7), sp.Rational(11, 3), sp.Rational(-2)):
    n = -(t - 1)/(t + 1)
    for assign in ('FL P^2 ~ record Y, Q^2 ~ X', 'FL P^2 ~ record X, Q^2 ~ Y'):
        sc = set()
        for K in (2, 3, 4, 5):
            ib, fl = ib_top(t, K), fl_top(n, K)
            for j in range(K):
                r_fl = fl[j + 1]/fl[j] if assign.startswith('FL P^2 ~ record Y') else fl[K - j - 1]/fl[K - j]
                sc.add(sp.nsimplify(sp.simplify((ib[j + 1]/ib[j])/r_fl)))
        print(t, '|', assign, '| distinct scale factors:', sorted(map(str, sc)))
