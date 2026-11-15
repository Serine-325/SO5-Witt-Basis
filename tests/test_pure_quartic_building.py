
from fractions import Fraction

import pytest
import sympy as sp

from building import compute_pure_quartic_building

varpi, epsilon = sp.symbols(
    "varpi epsilon", nonzero=True
)

Q = Fraction


@pytest.mark.parametrize(
    "alpha, r, expected_depths, expected_original, expected_alcove",
    [
        (
            epsilon, 1,
            [Q(0), Q(0)],
            (Q(0), Q(0)),
            (Q(0), Q(0)),
        ),
        (
            epsilon, varpi,
            [Q(1, 2), Q(1, 2)],
            (Q(1, 2), Q(0)),
            (Q(1, 2), Q(0)),
        ),
        (
            varpi, 1,
            [Q(1, 4), Q(3, 4)],
            (Q(1, 4), Q(1, 2)),
            (Q(1, 2), Q(1, 4)),
        ),
        (
            varpi, epsilon,
            [Q(1, 4), Q(3, 4)],
            (Q(1, 4), Q(1, 2)),
            (Q(1, 2), Q(1, 4)),
        ),
        (
            epsilon * varpi, 1,
            [Q(1, 4), Q(3, 4)],
            (Q(1, 4), Q(1, 2)),
            (Q(1, 2), Q(1, 4)),
        ),
        (
            epsilon * varpi, epsilon,
            [Q(1, 4), Q(3, 4)],
            (Q(1, 4), Q(1, 2)),
            (Q(1, 2), Q(1, 4)),
        ),
    ],
)
def test_pure_quartic_building(
    alpha, r, expected_depths,
    expected_original, expected_alcove
):
    result = compute_pure_quartic_building(
        alpha, sp.sympify(r), varpi, epsilon
    )

    assert result["depths"] == expected_depths
    assert result["original_point"] == expected_original
    assert result["alcove_point"] == expected_alcove

from building import pure_quartic_generator_matrices
from quartic import TARGET_GRAM


def test_pure_quartic_lie_algebra_condition():
    alpha, r = sp.symbols(
        "alpha r", nonzero=True
    )

    M1, M2 = pure_quartic_generator_matrices(alpha, r)

    J = TARGET_GRAM

    for M in [M1, M2]:
        assert (
            M.T * J + J * M
        ).applyfunc(sp.simplify) == sp.zeros(5)


def test_pure_quartic_generators_commute():
    alpha, r = sp.symbols(
        "alpha r", nonzero=True
    )

    M1, M2 = pure_quartic_generator_matrices(alpha, r)

    assert (
        M1 * M2 - M2 * M1
    ).applyfunc(sp.simplify) == sp.zeros(5)
