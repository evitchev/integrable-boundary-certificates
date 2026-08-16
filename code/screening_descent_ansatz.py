"""Cheap character tests for naive two-dimensional screening sums.

This is not a construction of the missing Z2^2/Galois descent complex.
It tests two denominator-style guesses made from the exact one-ray numerator

    P(q) = sum_{m>=0} (-1)^m q^(m(m+1)/2)

against the augmented target 1 + q^2/(1-q^2), where q records density spin.
It also checks why the most literal Gram-weighted Z^2 theta series is not an
ordinary formal power series at physical paperclip parameters.
"""

from math import isqrt

import sympy as sp


ORDER = 14


def multiply(left, right):
    """Multiply truncated integer power series."""
    out = [0] * ORDER
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right[:ORDER - i]):
            out[i + j] += left_coefficient * right_coefficient
    return out


def two_boson_partitions():
    """Coefficients of prod_{k>=1} (1-q^k)^(-2)."""
    series = [0] * ORDER
    series[0] = 1
    for degree in range(1, ORDER):
        factor = [0] * ORDER
        # (1-x)^(-2) = sum_{r>=0} (r+1)x^r.
        for power in range(0, ORDER, degree):
            factor[power] = power // degree + 1
        series = multiply(series, factor)
    return series


def one_ray_numerator():
    series = [0] * ORDER
    m = 0
    while True:
        degree = m * (m + 1) // 2
        if degree >= ORDER:
            return series
        series[degree] = (-1) ** m
        m += 1


def times_one_minus_q(series):
    return [coefficient - (series[degree - 1] if degree else 0)
            for degree, coefficient in enumerate(series)]


def has_opposite_branch_partner(m, ell):
    """Whether (m,ell) has an integral s'=-s-1 exponent partner."""
    discriminant = (m - ell)**2 + 2 * (m + ell) + 1
    if discriminant < 0:
        return False
    root = isqrt(discriminant)
    return root * root == discriminant


def has_factor_pair_parameterization(m, ell):
    """Positive-quadrant factor-pair form of the partner condition."""
    smaller = min(m, ell)
    difference = abs(m - ell)
    for a in range(1, isqrt(smaller) + 1):
        if smaller % a == 0:
            b = smaller // a
            if difference == b - a - 1:
                return True
    return False


def positive_partner_count(bound):
    return sum(has_opposite_branch_partner(m, ell)
               for m in range(1, bound + 1)
               for ell in range(1, bound + 1))


def main():
    bosons = two_boson_partitions()
    local_functionals = times_one_minus_q(bosons)
    ray = one_ray_numerator()
    target = [int(degree == 0 or degree >= 2 and degree % 2 == 0)
              for degree in range(ORDER)]

    # Treating the two positive screening directions as independent rays.
    quadrant = multiply(local_functionals, multiply(ray, ray))

    # Inclusion-exclusion for a union of two rays sharing their origin.
    ray_union = [2 * coefficient - int(degree == 0)
                 for degree, coefficient in enumerate(ray)]
    union = multiply(local_functionals, ray_union)

    print(f"target:               {target}")
    print(f"independent quadrant: {quadrant}")
    print(f"two-ray union:        {union}")
    assert quadrant[2] == 2 != target[2]
    assert union[2] == target[2] and union[3] == 1 != target[3]

    n, m, ell = sp.symbols("n m ell")
    gram = sp.Matrix([[1, -(n + 1)], [-(n + 1), 1]])
    determinant = sp.factor(gram.det())
    grading = sp.factor(
        (m**2 + ell**2 - 2 * (n + 1) * m * ell + m + ell) / 2)
    diagonal = sp.factor(grading.subs({m: ell}))
    print(f"det Gram(alpha,beta): {determinant}")
    print(f"Gram grading Q_n(m,ell): {grading}")
    print(f"diagonal Q_n(ell,ell): {diagonal}")
    assert determinant == -n * (n + 2)
    assert sp.simplify(diagonal - ell * (-ell * n + 1)) == 0
    # Swap symmetry Q_n(m,ell) = Q_n(ell,m): the swap-antisymmetric
    # part of any fixed-coefficient kernel therefore sums to the
    # identically zero series, reducing the cone-kernel no-go to the
    # swap-symmetric case (together with its other hypotheses: fixed
    # coefficients, nonzero effective orbit coefficients on a
    # full-dimensional interior) -- review finding.
    swapped = grading.subs({m: ell, ell: m}, simultaneous=True)
    assert sp.simplify(grading - swapped) == 0

    # In the nonnegative quadrant, equality of two Q_n exponents for every n
    # forces the pairs to agree up to m <-> ell.  This is false on all of Z^2:
    # the alternative sum relation s'=-s-1 produces cross-chamber collisions.
    positive_value = sp.factor(grading.subs({m: 2, ell: 2}))
    negative_value = sp.factor(grading.subs({m: -4, ell: -1}))
    print(f"full-lattice collision: Q_n(2,2) = {positive_value}, "
          f"Q_n(-4,-1) = {negative_value}")
    assert sp.simplify(positive_value - negative_value) == 0

    # If m>=ell>=1 and d=m-ell, the square-discriminant equation is
    #
    #     r^2 = (d+1)^2 + 4*ell.
    #
    # Factoring the difference of squares gives ell=a*b and
    # d=b-a-1.  Hence there are at most divisor-many partner-admitting
    # points for each smaller coordinate: O(N log N), rather than N^2,
    # in an N by N positive box.
    for m_value in range(1, 101):
        for ell_value in range(1, 101):
            assert has_opposite_branch_partner(m_value, ell_value) \
                == has_factor_pair_parameterization(m_value, ell_value)
            assert not has_opposite_branch_partner(m_value, -ell_value)
            assert not has_opposite_branch_partner(-m_value, ell_value)
    assert not has_opposite_branch_partner(1, 1)
    assert not has_opposite_branch_partner(2, 1)
    assert has_opposite_branch_partner(2, 2)
    counts = [(bound, positive_partner_count(bound))
              for bound in (10, 100, 1000)]
    print(f"positive-quadrant partner counts: {counts}")
    assert counts == [(10, 12), (100, 301), (1000, 5131)]


if __name__ == "__main__":
    main()
