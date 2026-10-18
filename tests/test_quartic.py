
import sympy as sp

from quartic import (
    TARGET_GRAM,
    pure_quartic_basis,
    pure_quartic_form,
    pure_quartic_gram,
)


def test_pure_quartic_symbolic_gram():
    alpha, r = sp.symbols("alpha r", nonzero=True)

    gram = pure_quartic_gram(alpha, r)

    assert gram == TARGET_GRAM


def test_pure_quartic_isotropic_vectors():
    alpha, r = sp.symbols("alpha r", nonzero=True)

    e1, e2, f1, f2, g1 = pure_quartic_basis(alpha, r)

    for vector in (e1, e2, f1, f2):
        assert sp.simplify(
            pure_quartic_form(vector, alpha, r)
        ) == 0

    assert sp.simplify(
        pure_quartic_form(g1, alpha, r)
    ) == 1

from quartic import exceptional_quartic_gram


def test_exceptional_quartic_gram():
    """Verify the exceptional case d = eta."""
    s, t = sp.symbols("s t", nonzero=True)

    gram = exceptional_quartic_gram(s, t)

    # Apply the defining relation s**4 + t**4 = -1.
    gram = gram.applyfunc(
        lambda entry: sp.simplify(
            entry.subs(t**4, -1 - s**4)
        )
    )

    assert gram == TARGET_GRAM

from quartic import exceptional_twisted_quartic_gram


def test_exceptional_twisted_quartic_gram():
    """
    Verify the Witt Gram matrix for d = eta*pi,
    with mu4 = eta and s**4 + t**4 = -1.
    """
    s, t, pi = sp.symbols(
        "s t pi", nonzero=True
    )

    gram = exceptional_twisted_quartic_gram(s, t, pi)

    # Apply the defining relation.
    gram = gram.applyfunc(
        lambda entry: sp.simplify(
            entry.subs(t**4, -1 - s**4)
        )
    )

    assert gram == TARGET_GRAM


def test_exceptional_quartic_isotropic_vectors():
    """Verify the four isotropic vectors in the exceptional quartic basis."""
    import sympy as sp
    from quartic import (
        exceptional_quartic_basis,
        exceptional_quartic_form,
    )

    s = sp.Symbol("s", nonzero=True)
    t = sp.Integer(1)

    basis = exceptional_quartic_basis(s, t)

    for vector in basis[:4]:
        value = exceptional_quartic_form(vector, s, t)
        numerator, denominator = sp.fraction(sp.cancel(value))

        assert sp.rem(numerator, s**4 + 2, s) == 0
        assert sp.rem(denominator, s**4 + 2, s) != 0


def test_exceptional_twisted_quartic_isotropic_vectors():
    """Verify isotropy in the exceptional twisted quartic basis."""
    from quartic import (
        exceptional_twisted_quartic_basis,
        exceptional_twisted_quartic_form,
    )

    s = sp.Symbol("s", nonzero=True)
    t = sp.Integer(1)
    pi = sp.Integer(3)

    basis = exceptional_twisted_quartic_basis(s, t, pi)

    for vector in basis[:4]:
        value = exceptional_twisted_quartic_form(vector, s, t, pi)
        numerator, denominator = sp.fraction(sp.cancel(value))

        assert sp.rem(numerator, s**4 + 2, s) == 0
        assert sp.rem(denominator, s**4 + 2, s) != 0