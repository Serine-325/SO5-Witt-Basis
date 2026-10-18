
from admissibility import (
    pure_quartic_norm,
    pure_quartic_discriminant_check,
    hilbert_symbol_qp,
)


def test_pure_quartic_admissibility_q3():
    alpha = 3
    r = 1

    # First condition: lambda = D
    assert pure_quartic_norm(alpha, r) == -3
    assert pure_quartic_discriminant_check(alpha, r)

    # Second condition of Theorem 3.6
    # (mu_4, d) = 1 because mu_4 = 1.
    left_side = 1

    hasse_q1 = hilbert_symbol_qp(alpha, -1, 3)
    lambda_value = pure_quartic_norm(alpha, r)
    right_side = (
        hasse_q1
        * hilbert_symbol_qp(lambda_value, -1, 3)
    )

    assert left_side == right_side
