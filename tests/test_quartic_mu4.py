import pytest

from quartic_mu4 import mu4_options_for_case, is_admissible_mu4


CASES = [
    ("Pure quartic extension", "F(sqrt(epsilon))", "1", "varpi"),
    ("Pure quartic extension", "F(sqrt(varpi))", "1", "epsilon"),
    ("Pure quartic extension", "F(sqrt(epsilon*varpi))", "1", "epsilon"),
    ("Exceptional unramified: d = eta", "F(sqrt(-1))", "1", "varpi"),
    ("Exceptional twisted: d = eta*pi", "F(sqrt(-1))", "eta", "1"),
]


@pytest.mark.parametrize("quartic_case,subfield,accepted,rejected", CASES)
def test_mu4_options(quartic_case, subfield, accepted, rejected):
    assert mu4_options_for_case(quartic_case, subfield) == (
        accepted, rejected
    )


@pytest.mark.parametrize("quartic_case,subfield,accepted,rejected", CASES)
def test_admissible_mu4(quartic_case, subfield, accepted, rejected):
    assert is_admissible_mu4(quartic_case, subfield, accepted)


@pytest.mark.parametrize("quartic_case,subfield,accepted,rejected", CASES)
def test_rejected_mu4(quartic_case, subfield, accepted, rejected):
    assert not is_admissible_mu4(quartic_case, subfield, rejected)
