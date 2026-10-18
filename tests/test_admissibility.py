
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
def test_exceptional_quartic_admissibility_q3():
    """Verify the first exceptional quartic case over Q_3."""
    import sympy as sp
    from admissibility import (
        exceptional_quartic_norm,
        hilbert_symbol_qp,
    )

    s = sp.Symbol("s", nonzero=True)

    # Choose t = 1 and impose s**4 = -2.
    # Then s**4 + t**4 = -1.
    norm = exceptional_quartic_norm(s, 1)

    assert sp.rem(norm + 1, s**4 + 2, s) == 0

    # Theorem 3.6: lambda = -1 and D = -1.
    lambda_value = -1
    assert lambda_value == -1

    # Since mu_4 = 1, (mu_4, d) over E_1 equals 1.
    left_side = 1

    # For this case epsilon_H(q_1) = (-1, -1) over Q_3.
    hasse_q1 = hilbert_symbol_qp(-1, -1, 3)
    lambda_symbol = hilbert_symbol_qp(lambda_value, -1, 3)

    right_side = hasse_q1 * lambda_symbol

    assert left_side == right_side

def test_exceptional_twisted_quartic_admissibility_q3():
    """Verify the twisted exceptional quartic case over Q_3."""
    import sympy as sp
    from admissibility import (
        exceptional_quartic_norm,
        hilbert_symbol_qp,
    )

    s = sp.Symbol("s", nonzero=True)
    i = sp.Symbol("i")

    # Set t = 1 and impose s**4 = -2.
    # eta = s**2 + i, with N(eta) = -1.
    norm_eta = exceptional_quartic_norm(s, 1)
    norm_d = 9 * norm_eta

    # d = 3*eta, so N(d) = -9.
    assert sp.rem(norm_d + 9, s**4 + 2, s) == 0

    # In F_9, eta reduces to 1+i.
    # Its fourth power is -1, so eta is a nonsquare unit.
    fourth_power = sp.rem((1 + i)**4, i**2 + 1, i)
    assert sp.Poly(fourth_power + 1, i, modulus=3).is_zero

    # Local-field result: (eta, 3*eta) over E_1 equals -1.
    left_side = -1

    # Right-hand side of Theorem 3.6.
    hasse_q1 = hilbert_symbol_qp(-3, -1, 3)
    lambda_symbol = hilbert_symbol_qp(-1, -1, 3)

    right_side = hasse_q1 * lambda_symbol

    assert hasse_q1 == -1
    assert lambda_symbol == 1
    assert left_side == right_side

def test_hensel_lift_fourth_root_minus_two():
    from admissibility import hensel_lift_fourth_root_minus_two

    expected = [1, 7, 7, 34]

    for precision, root in enumerate(expected, start=1):
        result = hensel_lift_fourth_root_minus_two(precision)

        assert result == root
        assert (result**4 + 2) % (3**precision) == 0

    for precision in range(5, 11):
        result = hensel_lift_fourth_root_minus_two(precision)

        assert (result**4 + 2) % (3**precision) == 0