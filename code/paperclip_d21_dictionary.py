"""Exact paperclip / affine D(2,1;alpha) screening dictionary.

Feigin--Jimbo--Mukhin use three fermionic screening vectors a_1,a_2,a_3
with

    a_i^2 = 1,
    a_1.a_2 = k_3 + 1,  a_2.a_3 = k_1 + 1,
    a_3.a_1 = k_2 + 1,
    k_1 + k_2 + k_3 = -4,

and add the affine vector a_0 = -a_1-a_2-a_3.  This script verifies
symbolically that the four paperclip screening vectors are the rank-two,
critical-level specialization

    (k_1,k_2,k_3) = (-2,n,-n-2)

up to permutation.  No square-root branch choices enter the check.
"""

import sympy as sp


n = sp.symbols("n")

# A sign pair (epsilon, delta) denotes
# (epsilon*sqrt(-n/2), delta*sqrt((n+2)/2)).
a = {
    0: (-1, -1),
    1: (+1, +1),
    2: (+1, -1),
    3: (-1, +1),
}


def dot(left, right):
    """Paperclip bilinear form without choosing either square root."""
    eps, delta = left
    eps2, delta2 = right
    return sp.factor(-n * eps * eps2 / 2 + (n + 2) * delta * delta2 / 2)


def main():
    gram = sp.Matrix([[dot(a[i], a[j]) for j in range(4)]
                      for i in range(4)])

    assert all(sp.simplify(gram[i, i] - 1) == 0 for i in range(4))
    assert tuple(a[0][j] + a[1][j] + a[2][j] + a[3][j]
                 for j in range(2)) == (0, 0)

    k1 = sp.factor(dot(a[2], a[3]) - 1)
    k2 = sp.factor(dot(a[3], a[1]) - 1)
    k3 = sp.factor(dot(a[1], a[2]) - 1)
    levels = (k1, k2, k3)

    assert levels == (-2, n, -n - 2)
    assert sp.simplify(sum(levels) + 4) == 0
    assert sp.simplify(gram.det()) == 0
    assert sp.simplify(gram.extract([1, 2, 3], [1, 2, 3]).det()) == 0
    # Generically the datum really has rank two (it drops further only at
    # the excluded degenerations n=0,-2).
    assert sp.factor(gram.extract([1, 2], [1, 2]).det()) == -n * (n + 2)

    # For the generic three-vector D(2,1;alpha) datum, rank drop is
    # equivalent to one of the three sl_2 levels being critical.  Verify the
    # identity after imposing k_1+k_2+k_3=-4.
    u, v = sp.symbols("k_1 k_2")
    w = -4 - u - v
    generic_gram = sp.Matrix([
        [1, w + 1, v + 1],
        [w + 1, 1, u + 1],
        [v + 1, u + 1, 1],
    ])
    generic_det = sp.factor(generic_gram.det())
    critical_product = sp.factor(2 * (u + 2) * (v + 2) * (w + 2))
    assert sp.simplify(generic_det - critical_product) == 0

    # The derivative-free quartic part of the paperclip I_3 is not the
    # square of one quadratic form.  This rules out the most naive map
    # from a single critical sl_2 Segal--Sugawara generator.
    quartic_a = n / (6 * (3 * n + 2))
    quartic_b = n * (n + 2) / ((3 * n + 2) * (3 * n + 4))
    quartic_c = (n + 2) / (6 * (3 * n + 4))
    square_discriminant = sp.factor(quartic_b**2
                                    - 4 * quartic_a * quartic_c)
    assert square_discriminant == (
        -8 * n * (n + 2)
        / (9 * (3 * n + 2)**2 * (3 * n + 4)**2)
    )

    print("Gram matrix in the order (a_0,a_1,a_2,a_3):")
    sp.pprint(gram)
    print(f"levels (k_1,k_2,k_3) = {levels}")
    print("sum of levels =", sp.simplify(sum(levels)))
    print("generic rank = 2; k_1 = -2 is the critical sl_2 level")
    print("generic det Gram = 2*(k_1+2)*(k_2+2)*(k_3+2)")
    print("quartic B^2-4AC =", square_discriminant, "(not a square)")


if __name__ == "__main__":
    main()
