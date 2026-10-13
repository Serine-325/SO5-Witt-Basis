"""Verify the actual displayed sum-of-quadratics Witt vectors in every valid case."""
import pytest

from web_runner import CLASSES, MU_OPTIONS, run_constructor, valid_parameters
from verification_panel import EXPECTED, verify_sum_result


def admissible():
    return [
        (square, a, b, m1, m2)
        for square in (True, False)
        for a in CLASSES
        for b in CLASSES
        for m1 in MU_OPTIONS[a]
        for m2 in MU_OPTIONS[b]
        if valid_parameters(square, a, b, m1, m2)
    ]


CASES = admissible()


def test_expected_number_of_cases():
    assert len(CASES) == 36


@pytest.mark.parametrize('params', CASES)
def test_actual_returned_witt_basis_gram(params):
    result = run_constructor(*params)
    verification = verify_sum_result(result)
    assert verification.verified, (params, verification.gram)
    assert verification.gram == EXPECTED
    assert verification.norms == (0, 0, 0, 0, 1)


def test_fail_closed_on_corrupted_vector():
    from dataclasses import replace
    result = run_constructor(True, 'epsilon', 'epsilon', '1', '1')
    broken = replace(result, basis={**result.basis, 'g1': '(0, 0, 2)'})
    verification = verify_sum_result(broken)
    assert not verification.verified


def test_expected_gram_in_agreed_order():
    import sympy as sp
    assert EXPECTED == sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])


@pytest.mark.parametrize('mode,square,alpha,beta', [
    ('canonical', True, 'epsilon', 'pi'),
    ('special', False, 'epsilon', 'pi'),
])
def test_biquadratic_display_gram_order(mode, square, alpha, beta):
    from biquadratic import construct_biquadratic
    from verification_panel import biquad_gram_in_display_order
    result = construct_biquadratic(
        minus_one_square=square,
        alpha_name=alpha,
        beta_name=beta,
        mode=mode,
    )
    assert result.verified
    assert biquad_gram_in_display_order(result) == EXPECTED
