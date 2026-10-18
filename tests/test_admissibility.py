
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

def test_quadratic_element_arithmetic():
    from admissibility import QuadraticElement

    x = QuadraticElement(1, 2, 3)
    y = QuadraticElement(2, 1, 3)

    # Addition
    assert x + y == QuadraticElement(3, 3, 3)

    # Multiplication
    assert x * y == QuadraticElement(8, 5, 3)

    # Conjugation
    assert x.conjugate() == QuadraticElement(1, -2, 3)

    # Norm
    assert x.norm() == -11

    # Verify x * conjugate(x) = N(x)
    assert x * x.conjugate() == QuadraticElement(-11, 0, 3)


def test_quadratic_element_extension_mismatch():
    import pytest
    from admissibility import QuadraticElement

    x = QuadraticElement(1, 2, 3)
    y = QuadraticElement(1, 2, 5)

    with pytest.raises(ValueError):
        x + y

    with pytest.raises(ValueError):
        x * y

def test_quadratic_residue_multiplication():
    from admissibility import quadratic_residue_multiply

    mul = quadratic_residue_multiply

    # In F_9 = F_3[i], where i^2 = -1
    x = (1, 1)
    x_squared = mul(x, x, -1, 3)
    x_fourth = mul(x_squared, x_squared, -1, 3)

    assert x_squared == (0, 2)
    assert x_fourth == (2, 0)

    # Multiplicative identity
    assert mul(x, (1, 0), -1, 3) == x

def test_hilbert_symbol_q3_i():
    """Verify Hilbert symbols over the unramified field Q_3(i)."""
    from admissibility import (
        QuadraticElement,
        hilbert_symbol_q3_i,
    )

    one = QuadraticElement(1, 0, -1)
    pi = QuadraticElement(3, 0, -1)
    eta = QuadraticElement(1, 1, -1)

    assert hilbert_symbol_q3_i(eta, pi) == -1
    assert hilbert_symbol_q3_i(pi, pi) == 1
    assert hilbert_symbol_q3_i(one, eta) == 1

    # Symmetry of the Hilbert symbol
    assert hilbert_symbol_q3_i(pi, eta) == -1

    # The twisted exceptional parameter d = 3*eta
    d = eta * pi
    assert hilbert_symbol_q3_i(eta, d) == -1
def test_exceptional_quartic_extension_hilbert_symbols():
    """
    Check the extension-field Hilbert symbols for
    the two exceptional quartic constructions.
    """
    from admissibility import (
        QuadraticElement,
        hilbert_symbol_q3_i,
    )

    # s = 7 approximates a root of s^4 = -2 in Q_3
    s = 7
    assert (s**4 + 2) % 27 == 0

    one = QuadraticElement(1, 0, -1)
    eta = QuadraticElement(s**2, 1, -1)
    pi = QuadraticElement(3, 0, -1)

    # Case 1: mu = 1, d = eta
    assert hilbert_symbol_q3_i(one, eta) == 1

    # Case 2: mu = eta, d = 3*eta
    d_twisted = eta * pi
    assert hilbert_symbol_q3_i(eta, d_twisted) == -1

    # The norm of eta has the square class of -1 in Q_3
    assert (-eta.norm()) % 3 == 1
def test_exceptional_quartic_full_admissibility():
    """
    Verify both admissibility conditions in Theorem 3.6
    for the two exceptional quartic cases over Q_3.
    """
    from admissibility import (
        QuadraticElement,
        hilbert_symbol_qp,
        hilbert_symbol_q3_i,
    )

    s = 7  # 3-adic approximation: s^4 = -2 mod 27
    assert (s**4 + 2) % 27 == 0

    alpha = -1
    lam = -1

    one = QuadraticElement(1, 0, alpha)
    eta = QuadraticElement(s**2, 1, alpha)

    cases = [
        ("unramified", one, 1, 1),
        ("twisted", eta, 3, -1),
    ]

    for name, mu, scale, expected in cases:
        d = eta * QuadraticElement(scale, 0, alpha)

        a = scale * s**2
        b = scale
        D = a**2 - alpha * b**2

        # Condition 1: lambda and D have the same
        # square class in Q_3.
        assert (-D) % (scale**2) == 0
        assert ((-D) // (scale**2)) % 3 == 1

        # Condition 2: extension-field Hilbert symbol
        lhs = hilbert_symbol_q3_i(mu, d)

        # Hasse invariant of q_1
        epsilon = (
            hilbert_symbol_qp(alpha, alpha * D, 3)
            * hilbert_symbol_qp(-a, -alpha * D, 3)
        )

        rhs = epsilon * hilbert_symbol_qp(lam, -1, 3)

        assert lhs == expected, name
        assert rhs == expected, name
        assert lhs == rhs, name