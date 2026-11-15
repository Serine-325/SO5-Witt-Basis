
import sympy as sp

from building import sum_quadratic_coordinate_matrices


def test_sum_quadratic_first_generator():
    alpha, beta = sp.symbols("alpha beta")
    a, b, c, d, z = sp.symbols("a b c d z")

    X1, X2 = sum_quadratic_coordinate_matrices(
        alpha, beta
    )

    v = sp.Matrix([a, b, c, d, z])

    expected = sp.Matrix([
        alpha * b, a, 0, 0, 0
    ])

    assert X1 * v == expected


def test_sum_quadratic_second_generator():
    alpha, beta = sp.symbols("alpha beta")
    a, b, c, d, z = sp.symbols("a b c d z")

    X1, X2 = sum_quadratic_coordinate_matrices(
        alpha, beta
    )

    v = sp.Matrix([a, b, c, d, z])

    expected = sp.Matrix([
        0, 0, beta * d, c, 0
    ])

    assert X2 * v == expected


def test_equal_quadratic_generators_in_witt_basis():
    from web_runner import run_constructor
    from building import sum_quadratic_generator_matrices
    from verification_panel import EXPECTED

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon",
        beta="epsilon",
        mu1="1",
        mu2="1"
    )

    assert result.valid
    assert result.case == "4.1"

    M1, M2 = sum_quadratic_generator_matrices(result)

    J = EXPECTED

    for M in (M1, M2):
        assert M.shape == (5, 5)

        check = (M.T * J + J * M).applyfunc(
            sp.simplify
        )

        assert check == sp.zeros(5)


def test_equal_quadratic_generator_coefficients():
    from web_runner import run_constructor
    from building import sum_quadratic_generator_matrices
    from verification_panel import EPS

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon",
        beta="epsilon",
        mu1="1",
        mu2="1"
    )

    M1, M2 = sum_quadratic_generator_matrices(result)

    expected = sp.Matrix([
        [0, EPS/2, 0, -sp.Rational(1, 4), 0],
        [sp.Rational(1, 2), 0, sp.Rational(1, 4), 0, 0],
        [0, EPS, 0, -sp.Rational(1, 2), 0],
        [-EPS, 0, -EPS/2, 0, 0],
        [0, 0, 0, 0, 0]
    ])

    assert (M1 - expected).applyfunc(
        sp.simplify
    ) == sp.zeros(5)


def test_equal_quadratic_second_generator_coefficients():
    from web_runner import run_constructor
    from building import sum_quadratic_generator_matrices
    from verification_panel import EPS

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon",
        beta="epsilon",
        mu1="1",
        mu2="1"
    )

    M1, M2 = sum_quadratic_generator_matrices(result)

    expected = sp.Matrix([
        [0, EPS/2, 0, sp.Rational(1, 4), 0],
        [sp.Rational(1, 2), 0, -sp.Rational(1, 4), 0, 0],
        [0, -EPS, 0, -sp.Rational(1, 2), 0],
        [EPS, 0, -EPS/2, 0, 0],
        [0, 0, 0, 0, 0]
    ])

    assert (M2 - expected).applyfunc(
        sp.simplify
    ) == sp.zeros(5)


def test_equal_unramified_building_point():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon",
        beta="epsilon",
        mu1="1",
        mu2="1"
    )

    assert result.valid

    building_result = compute_sum_quadratic_building(
        result
    )

    assert building_result["original_point"] == (0, 0)
    assert building_result["alcove_point"] == (0, 0)
    assert building_result["depths"] == [0, 0]


def test_equal_unramified_mu_pi_building_point():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon",
        beta="epsilon",
        mu1="pi",
        mu2="pi"
    )

    assert result.valid

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == (
        Fraction(1, 2), Fraction(1, 2)
    )


def test_equal_ramified_pi_building_point():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square=True,
        alpha="pi",
        beta="pi",
        mu1="1",
        mu2="1"
    )

    assert result.valid

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )


def test_equal_ramified_epsilon_pi_building_point():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon*pi",
        beta="epsilon*pi",
        mu1="1",
        mu2="1"
    )

    assert result.valid

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )



def test_equal_quadratic_minus_one_nonsquare():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    half = Fraction(1, 2)

    cases = [
        ("epsilon", "1", "1", (0, 0)),
        ("epsilon", "pi", "pi", (half, half)),
        ("pi", "1", "epsilon", (half, 0)),
        ("pi", "epsilon", "1", (half, 0)),
        ("epsilon*pi", "1", "epsilon", (half, 0)),
        ("epsilon*pi", "epsilon", "1", (half, 0)),
    ]

    for alpha, mu1, mu2, expected in cases:
        result = run_constructor(
            minus_one_square=False,
            alpha=alpha,
            beta=alpha,
            mu1=mu1,
            mu2=mu2
        )

        assert result.valid

        data = compute_sum_quadratic_building(result)

        assert data["alcove_point"] == expected, (
            f"Failed for alpha={alpha}, "
            f"mu1={mu1}, mu2={mu2}"
        )


def test_equal_ramified_pi_mu_epsilon():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square=True,
        alpha="pi",
        beta="pi",
        mu1="epsilon",
        mu2="epsilon"
    )

    assert result.valid

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )


def test_equal_ramified_epsilon_pi_mu_epsilon():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square=True,
        alpha="epsilon*pi",
        beta="epsilon*pi",
        mu1="epsilon",
        mu2="epsilon"
    )

    assert result.valid

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )


def test_case_42_mu_one():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        False, "pi", "epsilon*pi", "1", "1"
    )

    assert result.valid
    assert result.case == "4.2"

    data = compute_sum_quadratic_building(result)

    assert data["original_point"] == (
        0, Fraction(1, 2)
    )

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )

    assert data["diagonal_verified"]


def test_case_42_mu_epsilon():
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        False, "pi", "epsilon*pi",
        "epsilon", "epsilon"
    )

    assert result.valid
    assert result.case == "4.2"

    data = compute_sum_quadratic_building(result)

    assert data["original_point"] == (
        Fraction(1, 2), 0
    )

    assert data["alcove_point"] == (
        Fraction(1, 2), 0
    )

    assert data["diagonal_verified"]

import pytest


@pytest.mark.parametrize(
    "minus_one_square, alpha, beta, mu1, mu2, expected",
    [
        (True, "epsilon", "pi", "1", "epsilon", (0.5, 0)),
        (True, "epsilon", "pi", "pi", "1", (0.5, 0.5)),
        (True, "epsilon", "epsilon*pi", "1", "epsilon", (0.5, 0)),
        (True, "epsilon", "epsilon*pi", "pi", "1", (0.5, 0.5)),
        (True, "pi", "epsilon*pi", "1", "epsilon", (0.5, 0)),
        (True, "pi", "epsilon*pi", "epsilon", "1", (0.5, 0)),

        (False, "epsilon", "pi", "1", "epsilon", (0.5, 0)),
        (False, "epsilon", "pi", "pi", "1", (0.5, 0.5)),
        (False, "epsilon", "epsilon*pi", "1", "epsilon", (0.5, 0)),
        (False, "epsilon", "epsilon*pi", "pi", "1", (0.5, 0.5)),
    ],
)
def test_case_43_building(
    minus_one_square, alpha, beta, mu1, mu2, expected
):
    from web_runner import run_constructor
    from building import compute_sum_quadratic_building
    from fractions import Fraction

    result = run_constructor(
        minus_one_square, alpha, beta, mu1, mu2
    )

    assert result.valid
    assert result.case == "4.3"

    data = compute_sum_quadratic_building(result)

    assert data["alcove_point"] == tuple(
        Fraction(x) for x in expected
    )
    assert data["diagonal_verified"]
