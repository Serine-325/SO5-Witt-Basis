"""Exact Gram matrix tests for the exceptional Case 4.2."""

import pytest
import sympy as sp


@pytest.mark.parametrize("reversed_order", [False, True])
def test_case_42_exceptional_gram(reversed_order):
    p, s, t = sp.symbols("p s t", nonzero=True)
    half = sp.Rational(1, 2)

    # Coordinates (a,b,c,d,z):
    # x = a + b*sqrt(alpha)
    # y = c + d*sqrt(beta)
    # q = -N_alpha(x) - N_beta(y) - z**2

    if not reversed_order:
        alpha, beta = p, -p

        e1 = sp.Matrix([0, 1, 0, 1, 0])
        e2 = sp.Matrix([t**2, 0, -s**2, 0, 1])
        f1 = sp.Matrix([0, 1/(2*p), 0, -1/(2*p), 0])
        f2 = sp.Matrix([t**2*half, 0, -s**2*half, 0, -half])
        g1 = sp.Matrix([s**2, 0, t**2, 0, 0])

    else:
        alpha, beta = -p, p

        e1 = sp.Matrix([0, 1, 0, 1, 0])
        e2 = sp.Matrix([-s**2, 0, t**2, 0, 1])
        f1 = sp.Matrix([0, -1/(2*p), 0, 1/(2*p), 0])
        f2 = sp.Matrix([-s**2*half, 0, t**2*half, 0, -half])
        g1 = sp.Matrix([t**2, 0, s**2, 0, 0])

    # Matrix of the symmetric bilinear form
    # B(u,v) = (q(u+v)-q(u)-q(v))/2
    D = sp.diag(-1, alpha, -1, beta, -1)

    basis = [e1, e2, f1, f2, g1]

    gram = sp.Matrix(
        5, 5,
        lambda i, j: (basis[i].T * D * basis[j])[0]
    )

    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    # Apply the mathematical relation s**4 + t**4 = -1.
    for i in range(5):
        for j in range(5):
            difference = sp.expand(gram[i, j] - expected[i, j])
            difference = difference.subs(t**4, -1 - s**4)
            assert sp.simplify(difference) == 0

    print(
        "PASSED:",
        "reversed" if reversed_order else "canonical",
        "ordering"
    )


@pytest.mark.parametrize("reversed_order", [False, True])
def test_case_42_plus_plus_gram(reversed_order):
    """Verify Proposition 4.3, Case 4.2 (+1,+1)."""

    p = sp.symbols("p", nonzero=True)
    half = sp.Rational(1, 2)

    # Normalize epsilon = -1.
    if reversed_order:
        alpha, beta = -p, p
    else:
        alpha, beta = p, -p

    # Coordinates (a,b,c,d,z):
    # x = a + b*sqrt(alpha)
    # y = c + d*sqrt(beta)
    # q = N_alpha(x) + N_beta(y) - z**2

    e1 = sp.Matrix([1, 0, 0, 0, 1])
    e2 = sp.Matrix([0, 1, 0, 1, 0])
    f1 = sp.Matrix([half, 0, 0, 0, -half])
    f2 = sp.Matrix([
        0, -1/(2*alpha), 0, 1/(2*alpha), 0
    ])
    g1 = sp.Matrix([0, 0, 1, 0, 0])

    D = sp.diag(1, -alpha, 1, -beta, -1)
    basis = [e1, e2, f1, f2, g1]

    gram = sp.Matrix(
        5, 5,
        lambda i, j: (basis[i].T * D * basis[j])[0]
    )

    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    assert all(
        sp.simplify(gram[i, j] - expected[i, j]) == 0
        for i in range(5)
        for j in range(5)
    )
