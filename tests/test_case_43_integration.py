"""Integration tests for all admissible Case 4.3 configurations."""

import pytest

from web_runner import (
    CLASSES,
    MU_OPTIONS,
    valid_parameters,
    run_constructor,
)


def case_43_parameters():
    cases = []

    for minus_one_square in (True, False):
        for alpha in CLASSES:
            for beta in CLASSES:
                for mu1 in MU_OPTIONS[alpha]:
                    for mu2 in MU_OPTIONS[beta]:
                        params = (
                            minus_one_square,
                            alpha,
                            beta,
                            mu1,
                            mu2,
                        )

                        if not valid_parameters(*params):
                            continue

                        result = run_constructor(*params)

                        if result.case == "4.3":
                            cases.append(params)

    return cases


CASES = case_43_parameters()


def test_case_43_configuration_count():
    assert len(CASES) == 20


@pytest.mark.parametrize("params", CASES)
def test_case_43_constructor_output(params):
    result = run_constructor(*params)

    assert result.valid
    assert result.case == "4.3"

    assert set(result.basis) == {
        "e1", "e2", "f1", "f2", "g1"
    }

    for name in ("e1", "e2", "f1", "f2", "g1"):
        assert result.basis[name].strip()

    assert "r" in result.witnesses
    assert "s" in result.witnesses

    assert "Use Proposition 4.8" in result.raw_output
    assert "Unsupported sign combination" not in result.raw_output
