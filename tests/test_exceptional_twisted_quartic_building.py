from fractions import Fraction

import sympy as sp

from quartic import TARGET_GRAM
from building import (
    exceptional_twisted_quartic_generator_matrices,
    exceptional_twisted_quartic_valuation_bound,
    compute_exceptional_twisted_quartic_building,
)


s, t, pi = sp.symbols("s t pi", nonzero=True)


def test_generator_matrices():
    M1, M2 = exceptional_twisted_quartic_generator_matrices(
        s, t, pi
    )

    for M in [M1, M2]:
        assert M.shape == (5, 5)
        assert (
            M.T * TARGET_GRAM + TARGET_GRAM * M
        ).applyfunc(sp.simplify) == sp.zeros(5)



def test_generator_commutativity():
    M1, M2 = exceptional_twisted_quartic_generator_matrices(
        s, t, pi
    )

    def reduce_witness(expr):
        expr = sp.cancel(expr)
        numerator, denominator = sp.fraction(expr)

        numerator = sp.rem(
            sp.expand(numerator),
            t**4 + s**4 + 1,
            t,
        )

        return sp.cancel(numerator / denominator)

    commutator = (M1 * M2 - M2 * M1).applyfunc(
        reduce_witness
    )

    assert commutator == sp.zeros(5)

def test_generator_depths():
    result = compute_exceptional_twisted_quartic_building()

    assert result["depths"] == [
        Fraction(1, 2),
        Fraction(1, 2),
    ]


def test_original_coordinates():
    result = compute_exceptional_twisted_quartic_building()

    assert result["original_point"] == (
        Fraction(1, 2),
        Fraction(0),
    )


def test_fundamental_alcove():
    result = compute_exceptional_twisted_quartic_building()

    assert result["alcove_point"] == (
        Fraction(1, 2),
        Fraction(0),
    )


def test_diagonal_condition():
    result = compute_exceptional_twisted_quartic_building()

    assert result["diagonal_verified"] is True


def test_valuation_bounds():
    result = compute_exceptional_twisted_quartic_building()

    assert result["valuation_bounds"] is True


def test_uniformizer_valuation():
    assert exceptional_twisted_quartic_valuation_bound(
        pi / s, s, t, pi
    ) == Fraction(1)


def test_unit_valuation():
    assert exceptional_twisted_quartic_valuation_bound(
        1 / (2 * s), s, t, pi
    ) == Fraction(0)


def test_zero_valuation():
    assert exceptional_twisted_quartic_valuation_bound(
        sp.Integer(0), s, t, pi
    ) == sp.oo
