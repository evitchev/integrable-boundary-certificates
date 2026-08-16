"""Exact dictionary certificate for the all-spin A_2^(2) hierarchy.

Feigin--Frenkel, arXiv:hep-th/9310022, Theorem 4.7.3 proves that the
generic quantum affine-Toda hierarchy for A_2^(2) has mutually commuting
local integrals in every positive degree coprime to 6.  To apply that
theorem to Vitchev's equation (34), one still has to identify his two
branches with the two Coulomb-gas (Phi_1,2 / Phi_2,1) realizations.

This script supplies that missing normalization check.  It builds the
background-charge Virasoro field with the validated point-split normal
product, solves the exact twisted-total-derivative screening condition at
density spin six, and compares the result with Vitchev's equation (36).

The calculation is generic in the Coulomb-gas parameter ``a``.  Isolated
resonant values can have larger kernels, so they are deliberately not used
to weaken the generic theorem into an all-specializations statement.
"""

import math

import sympy as sp

from screening_fit import screening_conditions
from super_engine import d_cur, normal_prod


def _bosonic_current(super_current):
    """Drop the parity tag from a purely bosonic super-engine current."""
    out = {}
    for mono, coefficient in super_current.items():
        assert all(parity == 0 for parity, _, _ in mono)
        plain = tuple((order, field) for _, order, field in mono)
        out[plain] = sp.factor(coefficient)
    return out


def _screened_coefficient(density, momentum, coefficient):
    """Solve generically for the coefficient making ``density`` screened."""
    conditions = screening_conditions(density, (momentum,), 1)
    solutions = sp.solve(conditions, coefficient, dict=True)
    assert len(solutions) == 1 and coefficient in solutions[0]
    solution = sp.factor(solutions[0][coefficient])
    assert all(sp.factor(condition.subs(coefficient, solution)) == 0
               for condition in conditions)
    return solution, conditions


def main():
    a, coefficient = sp.symbols("a coefficient", nonzero=True)

    # In the convention <dX(z)dX(w)>=(z-w)^-2, the Coulomb-gas stress
    # tensor below has c=1-12 Q^2.  The two ordinary Virasoro screenings
    # have momenta a and -2/a.  The two degenerate perturbing momenta are
    # gamma_12=1/a and gamma_21=-a/2.
    background_charge = a / 2 - 1 / a
    current = ((0, 1, 0),)
    derivative_current = ((0, 2, 0),)
    stress_tensor = {
        tuple(sorted(current + current)): sp.Rational(1, 2),
        derivative_current: background_charge,
    }

    stress_cubed = normal_prod(normal_prod(stress_tensor, stress_tensor),
                               stress_tensor)
    d_stress_squared = normal_prod(d_cur(stress_tensor),
                                   d_cur(stress_tensor))
    trial = dict(stress_cubed)
    for mono, value in d_stress_squared.items():
        trial[mono] = sp.expand(trial.get(mono, 0) + coefficient * value)
    trial = _bosonic_current(trial)

    # Joint-kernel lemma (added after external review): the ordinary
    # Virasoro screenings S_a and S_{-2/a} impose NO condition on the
    # pencil -- every :T^3: + A :(dT)^2: lies in their kernels
    # automatically.  Hence the single perturbing condition below
    # places each sheet's charge in the FULL twisted affine-Toda joint
    # kernel: {S_a, S_{-a/2}} for one sheet, the dual {S_{-2/a},
    # S_{1/a}} for the other (root ratio -1/2 in both), which is the
    # screening datum Feigin--Frenkel's Theorem 4.7.3 refers to.
    for virasoro_momentum in (a, -2 / a):
        residual = screening_conditions(trial, (virasoro_momentum,), 1)
        assert all(sp.factor(condition) == 0 for condition in residual)

    phi12, conditions12 = _screened_coefficient(
        trial, 1 / a, coefficient)
    phi21, conditions21 = _screened_coefficient(
        trial, -a / 2, coefficient)

    expected12 = sp.factor((12 * a**4 - 5 * a**2 + 3) / (8 * a**2))
    expected21 = sp.factor((3 * a**4 - 20 * a**2 + 192)
                           / (32 * a**2))
    assert sp.factor(phi12 - expected12) == 0
    assert sp.factor(phi21 - expected21) == 0

    # Vitchev uses c=N and s^2=(N-25)(N-1).  Put beta^2=a^2/2.
    # Then c=1-6(beta-beta^-1)^2 and
    # s=6(beta^2-beta^-2).  His P6 is
    #
    #   8 :T^3: + K_s :(dT)^2:,
    #
    # so the coefficient in the normalized trial above is K_s/8.
    central_charge = sp.factor(1 - 12 * background_charge**2)
    sheet_coordinate = sp.factor(3 * a**2 - 12 / a**2)
    assert sp.factor(central_charge - (13 - 3*a**2 - 12/a**2)) == 0
    assert sp.factor(sheet_coordinate**2
                     - (central_charge - 25)*(central_charge - 1)) == 0

    vitchev_plus = sp.factor(
        (181 - 17*central_charge - 15*sheet_coordinate) / 64)
    vitchev_minus = sp.factor(
        (181 - 17*central_charge + 15*sheet_coordinate) / 64)
    assert sp.factor(vitchev_plus - phi21) == 0
    assert sp.factor(vitchev_minus - phi12) == 0

    # Also check equation (34) against equations (36), (41), and (42)
    # coefficient by coefficient in Vitchev's invariant P6 basis
    # [(33), (22)(11), (12)(12), (11)^3].
    k_plus = sp.factor(8 * vitchev_plus)
    composite_plus = [3 - k_plus/12, -6, -12 + k_plus, 1]
    printed_plus = [
        (107 + 17*central_charge + 15*sheet_coordinate) / 96,
        -6,
        -(-85 + 17*central_charge + 15*sheet_coordinate) / 8,
        1,
    ]
    assert all(sp.factor(left - right) == 0
               for left, right in zip(composite_plus, printed_plus))

    k_minus = sp.factor(8 * vitchev_minus)
    composite_minus = [3 - k_minus/12, -6, -12 + k_minus, 1]
    printed_minus = [
        (107 + 17*central_charge - 15*sheet_coordinate) / 96,
        -6,
        -(-85 + 17*central_charge - 15*sheet_coordinate) / 8,
        1,
    ]
    assert all(sp.factor(left - right) == 0
               for left, right in zip(composite_minus, printed_minus))

    degrees = [degree for degree in range(1, 50)
               if math.gcd(degree, 6) == 1]
    assert degrees[:8] == [1, 5, 7, 11, 13, 17, 19, 23]

    print(f"c(a) = {central_charge}")
    print(f"s(a) = {sheet_coordinate}")
    print(f"Phi_1,2 coefficient of :(dT)^2: = {phi12}")
    print(f"Phi_2,1 coefficient of :(dT)^2: = {phi21}")
    print(f"generic A_2^(2) IM degrees < 50: {degrees}")
    print("joint-kernel lemma: Virasoro screenings S_a, S_-2/a impose "
          "no condition on the pencil (asserted)")
    print(f"screening equations checked: {len(conditions12)} + "
          f"{len(conditions21)}")


if __name__ == "__main__":
    main()
