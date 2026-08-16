"""Low-grade joint-screening kernel pilot for the all-spin problem.

This is a diagnostic, not a character proof.  On either invariant sector the
sign-symmetry group permutes the screenings transitively and fixes candidate
densities, so the joint kernel equals the kernel of any one screening.  The
script therefore imposes one screening, introduces its twisted-derivative
witness, and finally quotients ordinary total derivatives.

It samples four paperclip parameters and three exact fibers on one
rationalized pillow-to-paperclip degeneration path:

* the projected rank-two paperclip fiber n=9/8;
* the nondegenerate three-boson pillow fiber n=9/8, nu=16/25;
* the unreduced three-boson critical fiber at nu=0.

All normalized screening components are rational or rational multiples of
``I``.  Avoiding an unnecessary algebraic-number extension keeps the exact
three-field blocks small.  The critical fiber is intentionally unreduced.  Its
extra kernel dimensions are checked to equal the independently computed free
one-boson charge dimensions through the displayed range.

The calculation is deliberately kept at low density spins.  Its purpose is to
validate conventions for a future invariant one-screening descent
calculation, not to extrapolate an all-spin invariant result.  It also checks
the Euler character of the *full*, undescended one-screening ray through
charge spin seven.  Mixing an invariant degree-zero term with unrestricted
nonzero-momentum terms would not define a complex, so every term in that ray
check is deliberately unrestricted.
"""

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from ope_engine import Sector, d_mono, gen_basis
from screening_fit import residue_screen, twisted_d


def _exact_rank(matrix):
    return DomainMatrix.from_Matrix(matrix).rank()


def _basis(spin, nfields):
    if spin == 0:
        return [()]
    return gen_basis(spin, nfields, z2=False)


def _invariant_basis(spin, nfields, symmetry):
    basis = _basis(spin, nfields)
    if symmetry == "none":
        return basis
    if symmetry == "even_each":
        return [mono for mono in basis if all(
            sum(field == j for _, field in mono) % 2 == 0
            for j in range(nfields)
        )]
    if symmetry == "tetrahedral":
        # Pillow sign symmetries flip two of the three fields.  Hence the
        # three letter-count parities must agree: all even or all odd.
        return [mono for mono in basis if len({
            sum(field == j for _, field in mono) % 2
            for j in range(nfields)
        }) == 1]
    raise ValueError(f"unknown symmetry {symmetry!r}")


def joint_kernel_dimension(density_spin, nfields, screenings, symmetry):
    """Dimension of the joint kernel on invariant local functionals."""
    candidates = _invariant_basis(density_spin, nfields, symmetry)
    target = _basis(density_spin - 1, nfields)
    witnesses = _basis(density_spin - 2, nfields)
    target_index = {mono: i for i, mono in enumerate(target)}

    n_cand = len(candidates)
    n_wit = len(witnesses)
    n_rows = len(screenings) * len(target)
    n_cols = n_cand + len(screenings) * n_wit
    matrix = sp.zeros(n_rows, n_cols)

    for screen_no, alpha in enumerate(screenings):
        row0 = screen_no * len(target)
        for col, mono in enumerate(candidates):
            for out_mono, coeff in residue_screen({mono: 1}, alpha).items():
                matrix[row0 + target_index[out_mono], col] = coeff
        witness0 = n_cand + screen_no * n_wit
        for local_col, mono in enumerate(witnesses):
            for out_mono, coeff in twisted_d(mono, alpha,
                                               nfields).items():
                matrix[row0 + target_index[out_mono],
                       witness0 + local_col] = -coeff

    # DomainMatrix keeps these examples over QQ_I or a small algebraic
    # extension.  Its exact fraction-free nullspace is dramatically faster
    # than expression-level Matrix.nullspace for the pillow block.
    null_rows = DomainMatrix.from_Matrix(matrix).nullspace().to_Matrix()
    projected = null_rows[:, :n_cand].T
    projected_rank = _exact_rank(projected)

    derivative_cols = []
    candidate_index = {mono: i for i, mono in enumerate(candidates)}
    for mono in _invariant_basis(density_spin - 1, nfields, symmetry):
        col = sp.zeros(n_cand, 1)
        for out_mono, coeff in d_mono(mono).items():
            col[candidate_index[out_mono], 0] = coeff
        if any(col):
            derivative_cols.append(col)
    derivatives = sp.Matrix.hstack(*derivative_cols) \
        if derivative_cols else sp.zeros(n_cand, 0)
    derivative_rank = _exact_rank(derivatives)

    # Every total derivative must occur among the screened densities.  Test
    # that inclusion explicitly rather than assuming it in the dimension.
    combined_rank = _exact_rank(sp.Matrix.hstack(projected, derivatives))
    assert combined_rank == projected_rank
    return projected_rank - derivative_rank


def paperclip_screenings(n):
    n = sp.sympify(n)
    a = sp.sqrt(-n / 2)
    b = sp.sqrt((n + 2) / 2)
    return [(eps * a, delta * b)
            for eps in (-1, 1) for delta in (-1, 1)]


def pillow_screenings(n, nu):
    """Pillow vectors of arXiv:1208.5259 divided by i*sqrt(2).

    This normalization makes every vector have square one and its nu=0
    rank-two projection equal to the paperclip square of screenings.
    """
    n = sp.sympify(n)
    nu = sp.sympify(nu)
    x = sp.I * sp.sqrt(n * nu / 2)
    y = sp.I * sp.sqrt(n * (1 - nu) / 2)
    z = sp.sqrt((n + 2) / 2)
    return [
        (+x, +y, +z),
        (+x, -y, -z),
        (-x, +y, -z),
        (-x, -y, +z),
    ]


def pillow_path_screenings(t):
    """A rationalized fixed-n path from pillow to paperclip.

    Here n=9/8 and

        x = (3/4) * 2t/(1+t^2),
        y = (3/4) * (1-t^2)/(1+t^2),
        z = 5/4,

    so x^2+y^2=9/16=n/2 and z^2=(n+2)/2.  The pillow parameter is
    nu=x^2/(n/2); t=0 is the rank-drop fiber.
    """
    t = sp.sympify(t)
    x = 3 * t / (2 * (1 + t**2))
    y = 3 * (1 - t**2) / (4 * (1 + t**2))
    z = sp.Rational(5, 4)
    return [
        (sp.I * x, sp.I * y, +z),
        (sp.I * x, -sp.I * y, -z),
        (-sp.I * x, sp.I * y, -z),
        (-sp.I * x, -sp.I * y, +z),
    ]


def kernel_row(nfields, screening, symmetry):
    return [joint_kernel_dimension(density_spin, nfields, [screening],
                                   symmetry)
            for density_spin in range(2, 7)]


def twisted_quotient_dimension(poly_spin, nfields, beta):
    """Dimension of polynomial-spin `poly_spin` modulo d_beta."""
    target = _basis(poly_spin, nfields)
    if poly_spin == 0:
        return 1
    lower = _basis(poly_spin - 1, nfields)
    target_index = {mono: i for i, mono in enumerate(target)}
    matrix = sp.zeros(len(target), len(lower))
    for col, mono in enumerate(lower):
        for out_mono, coeff in twisted_d(mono, beta, nfields).items():
            matrix[target_index[out_mono], col] = coeff
    return len(target) - _exact_rank(matrix)


def full_ray_dimensions(charge_spin, alpha):
    """Dimensions of the unrestricted fixed-grading screening ray.

    At cochain degree m the momentum is m*alpha and the polynomial spin is

        p_m = charge_spin + 1 - m*(m+1)/2.

    These are full two-boson spaces.  The sign-invariant charge sector is not
    a subcomplex because the first screening changes the momentum.
    """
    terms = []
    m = 0
    while True:
        poly_spin = charge_spin + 1 - m * (m + 1) // 2
        if poly_spin < 0:
            break
        momentum = tuple(m * component for component in alpha)
        terms.append(twisted_quotient_dimension(poly_spin, 2, momentum))
        m += 1
    return terms


def main():
    target = [1, 0, 1, 0, 1]
    n_values = [sp.Rational(9, 8), sp.Integer(2), sp.Rational(5, 3),
                sp.Rational(-7, 3)]
    for n_value in n_values:
        screening = paperclip_screenings(n_value)[0]
        actual = kernel_row(2, screening, "even_each")
        print(f"paperclip n={n_value}: {actual}", flush=True)
        assert actual == target

    pillow = kernel_row(3, pillow_path_screenings(sp.Rational(1, 2))[0],
                        "tetrahedral")
    print(f"pillow n=9/8, nu=16/25: {pillow}", flush=True)
    assert pillow == target

    critical = kernel_row(3, pillow_path_screenings(0)[0],
                          "tetrahedral")
    print(f"critical three-field nu=0 (unreduced): {critical}",
          flush=True)
    assert critical == [2, 0, 3, 0, 4]

    free = [Sector(charge_spin + 1, 1, z2=True).q_dim()
            for charge_spin in range(1, 6)]
    print(f"free Z2-even one-boson charges: {free}", flush=True)
    assert free == [1, 0, 2, 0, 3]
    assert critical == [pc + decoupled for pc, decoupled
                        in zip(target, free)]
    print("critical row = paperclip row + decoupled free-boson row",
          flush=True)

    alpha = paperclip_screenings(sp.Rational(9, 8))[0]
    expected_ray_terms = [
        [3, 1],
        [5, 3, 1],
        [10, 5, 1],
        [16, 10, 3],
        [29, 16, 5, 1],
        [45, 29, 10, 1],
        [75, 45, 16, 3],
    ]
    expected_full_kernels = [2, 3, 6, 9, 17, 25, 43]
    for charge_spin, expected_terms, expected_kernel in zip(
            range(1, 8), expected_ray_terms, expected_full_kernels):
        terms = full_ray_dimensions(charge_spin, alpha)
        euler = sum((-1) ** m * dimension
                    for m, dimension in enumerate(terms))
        direct = joint_kernel_dimension(charge_spin + 1, 2, [alpha],
                                        "none")
        print(f"full ray s={charge_spin}: {terms}, "
              f"Euler={euler}, direct kernel={direct}", flush=True)
        assert terms == expected_terms
        assert euler == direct == expected_kernel


if __name__ == "__main__":
    main()
