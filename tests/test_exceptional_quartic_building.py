import sympy as sp
from fractions import Fraction

from building import (
    exceptional_quartic_generator_matrices,
    exceptional_quartic_valuation_bound,
    compute_exceptional_quartic_building,
)

from quartic import TARGET_GRAM as J


def test_exceptional_quartic_lie_algebra():
    s, t = sp.symbols("s t", nonzero=True)

    M1, M2 = exceptional_quartic_generator_matrices(s, t)

    for M in [M1, M2]:
        assert (
            M.T * J + J * M
        ).applyfunc(sp.simplify) == sp.zeros(5)


def test_exceptional_quartic_depths():
    result = compute_exceptional_quartic_building()

    assert result["depths"] == [
        Fraction(0),
        Fraction(0),
    ]


def test_exceptional_quartic_original_point():
    result = compute_exceptional_quartic_building()

    assert result["original_point"] == (
        Fraction(0),
        Fraction(0),
    )


def test_exceptional_quartic_alcove_point():
    result = compute_exceptional_quartic_building()

    assert result["alcove_point"] == (
        Fraction(0),
        Fraction(0),
    )


def test_exceptional_quartic_valuation_bounds():
    s, t = sp.symbols("s t", nonzero=True)

    valuation = exceptional_quartic_valuation_bound

    assert valuation(0, s, t) == sp.oo
    assert valuation(1 / (2*s), s, t) == 0
    assert valuation(s**2 + t**2, s, t) == 0
    assert valuation(s**2 - t**2, s, t) == 0


def test_exceptional_quartic_diagonal_condition():
    result = compute_exceptional_quartic_building()

    assert result["diagonal_verified"] is True
    assert result["valuation_bounds"] is True
