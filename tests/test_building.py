
from fractions import Fraction

import pytest

from building import in_fundamental_alcove, fundamental_alcove


@pytest.mark.parametrize("x1, x2", [
    (0, 0),
    ("1/2", 0),
    ("1/4", "1/4"),
    ("1/2", "1/2"),
    (1, 0),
])
def test_points_inside_alcove(x1, x2):
    assert in_fundamental_alcove(x1, x2)


@pytest.mark.parametrize("x1, x2", [
    (0, "1/2"),
    ("1/2", "-1/4"),
    (1, "1/2"),
])
def test_points_outside_alcove(x1, x2):
    assert not in_fundamental_alcove(x1, x2)


@pytest.mark.parametrize("x1, x2, expected", [
    (0, 0, (0, 0)),
    (0, "1/2", ("1/2", 0)),
    ("1/4", "1/4", ("1/4", "1/4")),
    ("1/2", "1/2", ("1/2", "1/2")),
    (0, 1, (1, 0)),
    ("1/2", "-1/4", ("1/2", "1/4")),
    (1, "1/2", ("1/2", 0)),
    (-1, 0, (1, 0)),
])
def test_fundamental_alcove(x1, x2, expected):
    point, transformations = fundamental_alcove(x1, x2)

    expected = tuple(Fraction(t) for t in expected)

    assert point == expected
    assert in_fundamental_alcove(*point)


def test_exchange_hyperbolic_pairs():
    point, transformations = fundamental_alcove(0, "1/2")

    assert point == (Fraction(1, 2), Fraction(0))
    assert transformations == ["Exchange hyperbolic pairs"]

def test_biquadratic_root_coefficients():
    import sympy as sp
    from building import root_coefficients

    varpi = sp.Symbol("varpi", nonzero=True)

    # Matrix of X1 = B in the exceptional Witt basis
    M1 = sp.Matrix([
        [0, varpi, 0, 0, 0],
        [1, 0,     0, 0, 0],
        [0, 0,     0, -1, 0],
        [0, 0, -varpi, 0, 0],
        [0, 0,     0, 0, 0],
    ])

    coefficients = root_coefficients(M1)

    assert coefficients["epsilon1"] == {
        "positive": 0,
        "negative": 0,
    }

    assert coefficients["epsilon2"] == {
        "positive": 0,
        "negative": 0,
    }

    assert coefficients["epsilon1-epsilon2"] == {
        "positive": -varpi,
        "negative": 1,
    }

    assert coefficients["epsilon1+epsilon2"] == {
        "positive": 0,
        "negative": 0,
    }

def test_biquadratic_second_generator():
    import sympy as sp
    from building import root_coefficients

    varpi = sp.Symbol("varpi", nonzero=True)
    half = sp.Rational(1, 2)

    # Matrix of X2 = iB in the exceptional Witt basis
    M2 = sp.Matrix([
        [0,        0,    0, -half, 0],
        [0,        0, half,     0, 0],
        [0, -2*varpi,    0,     0, 0],
        [2*varpi,  0,    0,     0, 0],
        [0,        0,    0,     0, 0],
    ])

    coefficients = root_coefficients(M2)

    assert coefficients["epsilon1"] == {
        "positive": 0,
        "negative": 0,
    }

    assert coefficients["epsilon2"] == {
        "positive": 0,
        "negative": 0,
    }

    assert coefficients["epsilon1-epsilon2"] == {
        "positive": 0,
        "negative": 0,
    }

    assert coefficients["epsilon1+epsilon2"] == {
        "positive": half,
        "negative": 2*varpi,
    }

def test_monomial_valuation():
    import sympy as sp
    from building import monomial_valuation

    varpi = sp.Symbol("varpi", nonzero=True)

    assert monomial_valuation(0, varpi) == sp.oo
    assert monomial_valuation(1, varpi) == 0
    assert monomial_valuation(-varpi, varpi) == 1
    assert monomial_valuation(2*varpi, varpi) == 1
    assert monomial_valuation(sp.Rational(1, 2), varpi) == 0
    assert monomial_valuation(varpi**(-2), varpi) == -2


def test_monomial_valuation_p_adic():
    import sympy as sp
    from building import monomial_valuation

    varpi = sp.Symbol("varpi", nonzero=True)

    assert monomial_valuation(3, varpi, p=3) == 1
    assert monomial_valuation(9, varpi, p=3) == 2
    assert monomial_valuation(sp.Rational(1, 3), varpi, p=3) == -1
    assert monomial_valuation(2*varpi, varpi, p=3) == 1

def test_root_inequalities_first_generator():
    import sympy as sp
    from building import (
        root_coefficients,
        monomial_valuation,
        root_inequalities,
    )

    varpi = sp.Symbol("varpi", nonzero=True)

    M1 = sp.Matrix([
        [0, varpi, 0, 0, 0],
        [1, 0,     0, 0, 0],
        [0, 0,     0, -1, 0],
        [0, 0, -varpi, 0, 0],
        [0, 0,     0, 0, 0],
    ])

    coefficients = root_coefficients(M1)

    inequalities = root_inequalities(
        coefficients,
        Fraction(1, 2),
        lambda c: monomial_valuation(c, varpi),
    )

    assert inequalities["epsilon1-epsilon2"] == {
        "lower": Fraction(-1, 2),
        "upper": Fraction(-1, 2),
    }

    for root in [
        "epsilon1",
        "epsilon2",
        "epsilon1+epsilon2",
    ]:
        assert inequalities[root] == {
            "lower": -sp.oo,
            "upper": sp.oo,
        }

def test_combine_biquadratic_inequalities():
    import sympy as sp
    from building import (
        root_coefficients,
        monomial_valuation,
        root_inequalities,
        combine_inequalities,
    )

    varpi = sp.Symbol("varpi", nonzero=True)
    half = sp.Rational(1, 2)

    M1 = sp.Matrix([
        [0, varpi, 0, 0, 0],
        [1, 0,     0, 0, 0],
        [0, 0,     0, -1, 0],
        [0, 0, -varpi, 0, 0],
        [0, 0,     0, 0, 0],
    ])

    M2 = sp.Matrix([
        [0, 0, 0, -half, 0],
        [0, 0, half, 0, 0],
        [0, -2*varpi, 0, 0, 0],
        [2*varpi, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ])

    valuation = lambda c: monomial_valuation(c, varpi)

    inequalities1 = root_inequalities(
        root_coefficients(M1),
        Fraction(1, 2),
        valuation,
    )

    inequalities2 = root_inequalities(
        root_coefficients(M2),
        Fraction(1, 2),
        valuation,
    )

    combined = combine_inequalities(
        inequalities1,
        inequalities2,
    )

    assert combined["epsilon1-epsilon2"] == {
        "lower": Fraction(-1, 2),
        "upper": Fraction(-1, 2),
    }

    assert combined["epsilon1+epsilon2"] == {
        "lower": Fraction(1, 2),
        "upper": Fraction(1, 2),
    }

    for root in ["epsilon1", "epsilon2"]:
        assert combined[root] == {
            "lower": -sp.oo,
            "upper": sp.oo,
        }

def test_solve_biquadratic_building_point():
    import sympy as sp
    from building import solve_building_point

    inequalities = {
        "epsilon1": {
            "lower": -sp.oo,
            "upper": sp.oo,
        },
        "epsilon2": {
            "lower": -sp.oo,
            "upper": sp.oo,
        },
        "epsilon1-epsilon2": {
            "lower": Fraction(-1, 2),
            "upper": Fraction(-1, 2),
        },
        "epsilon1+epsilon2": {
            "lower": Fraction(1, 2),
            "upper": Fraction(1, 2),
        },
    }

    point = solve_building_point(inequalities)

    assert point == (Fraction(0), Fraction(1, 2))

    alcove_point, transformations = fundamental_alcove(
        *point
    )

    assert alcove_point == (
        Fraction(1, 2),
        Fraction(0),
    )


def test_solver_rejects_nonunique_point():
    from building import solve_building_point

    inequalities = {
        "epsilon1": {
            "lower": 0,
            "upper": 1,
        },
        "epsilon2": {
            "lower": 0,
            "upper": 1,
        },
    }

    with pytest.raises(ValueError):
        solve_building_point(inequalities)


def test_solver_rejects_unbounded_region():
    import sympy as sp
    from building import solve_building_point

    inequalities = {
        "epsilon1": {
            "lower": 0,
            "upper": sp.oo,
        },
        "epsilon2": {
            "lower": 0,
            "upper": sp.oo,
        },
    }

    with pytest.raises(ValueError):
        solve_building_point(inequalities)

def test_complete_exceptional_biquadratic_building():
    import sympy as sp
    from building import (
        compute_building_point,
        monomial_valuation,
    )

    varpi = sp.Symbol("varpi", nonzero=True)
    half = sp.Rational(1, 2)

    # First generator: X1 = B
    M1 = sp.Matrix([
        [0, varpi, 0, 0, 0],
        [1, 0,     0, 0, 0],
        [0, 0,     0, -1, 0],
        [0, 0, -varpi, 0, 0],
        [0, 0,     0, 0, 0],
    ])

    # Second generator: X2 = iB
    M2 = sp.Matrix([
        [0, 0, 0, -half, 0],
        [0, 0, half, 0, 0],
        [0, -2*varpi, 0, 0, 0],
        [2*varpi, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ])

    # Depths of the two generators
    depths = [
        Fraction(1, 2),
        Fraction(1, 2),
    ]

    valuation = lambda c: monomial_valuation(c, varpi)

    # Complete building computation
    result = compute_building_point(
        [M1, M2],
        depths,
        valuation,
    )

    # Verify the original building point
    assert result["original_point"] == (
        Fraction(0),
        Fraction(1, 2),
    )

    # Verify the fundamental alcove point
    assert result["alcove_point"] == (
        Fraction(1, 2),
        Fraction(0),
    )

    # Verify the Weyl transformation
    assert result["transformations"] == [
        "Exchange hyperbolic pairs"
    ]

    # Verify the diagonal condition
    assert result["diagonal_verified"] is True

    # Verify the coefficient valuations
    assert result["valuations"][0][
        "epsilon1-epsilon2"
    ] == {
        "positive": 1,
        "negative": 0,
    }

    assert result["valuations"][1][
        "epsilon1+epsilon2"
    ] == {
        "positive": 0,
        "negative": 1,
    }


def test_diagonal_condition_failure():
    import sympy as sp
    from building import (
        check_diagonal_condition,
        monomial_valuation,
    )

    varpi = sp.Symbol("varpi", nonzero=True)

    M = sp.diag(1, 0, -1, 0, 0)

    valuation = lambda c: monomial_valuation(c, varpi)

    with pytest.raises(ValueError):
        check_diagonal_condition(
            M,
            Fraction(1, 2),
            valuation,
        )

def test_matrix_in_exceptional_witt_basis():
    import sympy as sp
    from building import matrix_in_witt_basis

    varpi = sp.Symbol("varpi", nonzero=True)
    half = sp.Rational(1, 2)

    # Columns are e1, e2, f1, f2, g1
    # Original basis: (1, i, B, iB, g1)
    P = sp.Matrix([
        [1, 0,  half,             0, 0],
        [1, 0, -half,             0, 0],
        [0, 1,     0, -1/(2*varpi), 0],
        [0, 1,     0,  1/(2*varpi), 0],
        [0, 0,     0,             0, 1],
    ])

    # Multiplication by X1 = B
    L1 = sp.Matrix([
        [0, 0, varpi,     0, 0],
        [0, 0,     0, varpi, 0],
        [1, 0,     0,     0, 0],
        [0, 1,     0,     0, 0],
        [0, 0,     0,     0, 0],
    ])

    # Multiplication by X2 = iB
    L2 = sp.Matrix([
        [0,  0,     0, -varpi, 0],
        [0,  0, varpi,      0, 0],
        [0, -1,     0,      0, 0],
        [1,  0,     0,      0, 0],
        [0,  0,     0,      0, 0],
    ])

    # Compute matrices in the Witt basis
    M1 = matrix_in_witt_basis(L1, P)
    M2 = matrix_in_witt_basis(L2, P)

    # Expected matrices from our paper
    expected_M1 = sp.Matrix([
        [0, varpi, 0,  0, 0],
        [1,     0, 0,  0, 0],
        [0,     0, 0, -1, 0],
        [0,     0, -varpi, 0, 0],
        [0,     0, 0,  0, 0],
    ])

    expected_M2 = sp.Matrix([
        [0,        0,    0, -half, 0],
        [0,        0, half,     0, 0],
        [0, -2*varpi,    0,     0, 0],
        [2*varpi,  0,    0,     0, 0],
        [0,        0,    0,     0, 0],
    ])

    assert M1 == expected_M1
    assert M2 == expected_M2

def test_symbolic_valuation():
    import sympy as sp
    from building import symbolic_valuation

    epsilon = sp.Symbol("epsilon", nonzero=True)
    varpi = sp.Symbol("pi", nonzero=True)

    def valuation(expr):
        return symbolic_valuation(
            expr, varpi, epsilon
        )

    assert valuation(0) == sp.oo
    assert valuation(epsilon) == 0
    assert valuation(varpi) == 1
    assert valuation(epsilon*varpi) == 1
    assert valuation(-epsilon*varpi) == 1
    assert valuation(1/epsilon) == 0
    assert valuation(1/(epsilon*varpi)) == -1
    assert valuation(epsilon**2*varpi**3) == 3
    assert valuation(-varpi/epsilon) == 1
    assert valuation(2*epsilon*varpi) == 1

@pytest.mark.parametrize(
    "alpha_name, beta_name, expected",
    [
        ("epsilon", "pi", ("1/2", "0")),
        ("pi", "epsilon", ("1/4", "1/4")),
        ("epsilon*pi", "epsilon", ("1/4", "1/4")),
    ],
)
def test_automatic_biquadratic_building(
    alpha_name, beta_name, expected
):
    from biquadratic import construct_biquadratic
    from building import compute_biquadratic_building

    # Construct the Witt basis
    result = construct_biquadratic(
        minus_one_square=True,
        alpha_name=alpha_name,
        beta_name=beta_name,
    )

    # Compute the building point automatically
    building_result = compute_biquadratic_building(
        result
    )

    expected_point = tuple(
        Fraction(t) for t in expected
    )

    assert building_result["alcove_point"] == expected_point

    assert building_result["diagonal_verified"] is True

    assert in_fundamental_alcove(
        *building_result["alcove_point"]
    )

@pytest.mark.parametrize(
    "alpha_name, beta_name",
    [
        ("pi", "epsilon"),
        ("epsilon*pi", "epsilon"),
    ],
)
def test_biquadratic_ramified_minus_one_nonsquare(
    alpha_name, beta_name
):
    from biquadratic import construct_biquadratic
    from building import compute_biquadratic_building

    # Construct the Witt basis
    result = construct_biquadratic(
        minus_one_square=False,
        alpha_name=alpha_name,
        beta_name=beta_name,
    )

    # Compute the building point
    building_result = compute_biquadratic_building(
        result
    )

    # Verify the parameter valuations
    assert building_result["parameter_valuations"] == {
        "alpha": Fraction(1),
        "beta": Fraction(0),
    }

    # Verify the generator depths
    assert building_result["depths"] == [
        Fraction(0),
        Fraction(1, 2),
    ]

    # Verify the fundamental alcove point
    assert building_result["alcove_point"] == (
        Fraction(1, 4),
        Fraction(1, 4),
    )

    assert building_result["diagonal_verified"] is True

@pytest.mark.parametrize(
    "order_by_valuation",
    [False, True],
)
def test_automatic_exceptional_biquadratic_building(
    order_by_valuation
):
    from biquadratic import construct_biquadratic
    from building import compute_biquadratic_building

    # Construct the exceptional Witt basis
    result = construct_biquadratic(
        minus_one_square=False,
        alpha_name="epsilon",
        beta_name="pi",
        mode="special",
        order_by_valuation=order_by_valuation,
    )

    # Verify mu = 1
    assert result.mu0 == 1
    assert result.mu1 == 0

    # Verify the Witt basis
    assert result.verified

    # Compute the building point
    b = compute_biquadratic_building(result)

    # Verify parameter valuations
    assert b["parameter_valuations"] == {
        "alpha": Fraction(0),
        "beta": Fraction(1),
    }

    # Verify generator depths
    assert b["depths"] == [
        Fraction(1, 2),
        Fraction(1, 2),
    ]

    # Verify the fundamental alcove point
    assert b["alcove_point"] == (
        Fraction(1, 2),
        Fraction(0),
    )

    # Verify the original point
    if result.swapped_pairs:
        assert b["original_point"] == (
            Fraction(1, 2),
            Fraction(0),
        )
    else:
        assert b["original_point"] == (
            Fraction(0),
            Fraction(1, 2),
        )

    assert b["diagonal_verified"] is True
