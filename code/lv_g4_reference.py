"""Exact transcription of Litvinov-Vilkoviskiy G4 (arXiv:2105.04018,
eq. (2.8)) at the level of its structural coefficients, and the
free-CFT-limit analysis over all nine boundary-label pairs.

Conventions: epsilon_1 epsilon_2 = 1 (their footnote), epsilon_3 =
-epsilon_1 - epsilon_2; the free point (all background charges off)
is epsilon_3 = 0, i.e. epsilon_1 = t -> i with epsilon_2 = 1/t.

Audit-corrected (2026-08-07): a first reconnaissance claimed only
(D,D) has a finite free limit.  FALSE -- with the epsilon constraint,
(epsilon_alpha + epsilon_beta)/epsilon_3 = -1 identically for
(B,C)/(C,B), and every apparent 1/epsilon_3 divergence in the
site-dependent and double-sum terms is cancelled by explicit
epsilon_3 prefactors.  This script asserts, symbolically:

  - finite free limits exist for exactly (D,D), (B,C), (C,B);
  - the axis-quartic sum_k (dphi_k)^4 coefficient at the free point
    is -(2n-2)/3 for (D,D) -- nonzero for n >= 2, VANISHING at the
    rank-one point n = 1 -- and -(2n+1)/3 for (B,C)/(C,B), nonzero
    for all n >= 1.

SCOPE (audit round 2): what is machine-checked here is the finite
label set and the coefficient formulas.  The downstream no-go -- that
a nonzero axis-quartic admits no O(N)-invariant pair-contraction
expression under any linear change of boson frame -- is a
MATHEMATICAL ARGUMENT (rank-four tensor symmetry), not a machine
check, and it requires the would-be O(N) block to have dimension at
least two: for a single boson, Y^4 = (Y^2)^2 IS a pair invariant.
The BCD no-go is therefore a generic/nondegenerate-rank statement,
with the (D,D) case additionally degenerate at n = 1.
"""
import sys

import sympy as sp

t, n, j_, i_ = sp.symbols('t n j i', positive=False)
eps = {1: t, 2: 1 / t, 3: -t - 1 / t}
FREE = sp.I  # t -> i makes eps_3 = 0

failures = []


def coeffs(alpha, beta):
    """Structural coefficients of G4 as functions of t (eps_1 = t)."""
    e1, e2, e3 = eps[1], eps[2], eps[3]
    ea, eb = eps[alpha], eps[beta]
    quartic = -sp.Rational(1, 3) * (2 * n - (ea + eb) / e3)
    lin_left = 4 * e3 * (j_ - 1 + (e3 - ea) / (2 * e3))   # sqrt(e1 e2)=1
    lin_right = -4 * e3 * (n - j_ + (e3 - eb) / (2 * e3))
    dd = (2 * n + 4 * (n - 1) * (e1**2 + e2**2) / 3
          + (1 - 2 * e3**2) * (ea + eb - 2 * e3) / (3 * e3))
    dsum = -4 * e3**2 * ((i_ - 1 + (e3 - ea) / (2 * e3))
                         * (n - j_ + (e3 - eb) / (2 * e3)))
    return [quartic, lin_left, lin_right, dd, dsum]


def free_limit(expr):
    ex = sp.cancel(sp.together(expr))
    den = sp.denom(ex).subs(t, FREE)
    if sp.simplify(den) == 0:
        return None
    return sp.simplify(ex.subs(t, FREE))


LABELS = {1: "B", 2: "C", 3: "D"}
finite_pairs = {}
for a in (1, 2, 3):
    for b in (1, 2, 3):
        lims = [free_limit(c) for c in coeffs(a, b)]
        finite = all(l is not None for l in lims)
        if finite:
            finite_pairs[(a, b)] = lims
        print(f"({LABELS[a]},{LABELS[b]}): "
              + ("finite free limit; quartic -> "
                 f"{sp.factor(lims[0])}" if finite
                 else "divergent at the free point"))

expect = {(3, 3), (1, 2), (2, 1)}
if set(finite_pairs) != expect:
    failures.append(f"finite set {set(finite_pairs)} != {expect}")

qDD = sp.factor(finite_pairs[(3, 3)][0])
qBC = sp.factor(finite_pairs[(1, 2)][0])
if sp.simplify(qDD + sp.Rational(2, 3) * (n - 1)) != 0:
    failures.append(f"DD quartic {qDD}")
if sp.simplify(qBC + sp.Rational(1, 3) * (2 * n + 1)) != 0:
    failures.append(f"BC quartic {qBC}")
if sp.simplify(finite_pairs[(1, 2)][0] - finite_pairs[(2, 1)][0]) != 0:
    failures.append("BC != CB quartic")

if sp.simplify(qDD.subs(n, 1)) != 0:
    failures.append("DD quartic should vanish at the rank-one point n=1")
if sp.simplify(qBC.subs(n, 1)) == 0:
    failures.append("BC quartic unexpectedly vanishes at n=1")

if failures:
    print("FAILURES:", failures)
    sys.exit(1)
print("\nLV G4 FREE-LIMIT ANALYSIS PASSED: finite exactly for "
      "(D,D), (B,C), (C,B); axis-quartic -(2/3)(n-1) [DD, zero at "
      "n=1] and -(2n+1)/3 [BC/CB, nonzero].")
