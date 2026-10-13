"""Symbolic Gram matrix verification for Case 4.1."""

import sympy as sp


def test_case_41_gram_matrix():
    alpha, mu1, mu2 = sp.symbols(
        "alpha mu1 mu2", nonzero=True
    )
    u0, u1 = sp.symbols("u0 u1")

    # Coordinates: (a, b, c, d, z)
    # x = a + b*sqrt(alpha)
    # y = c + d*sqrt(alpha)

    Q = sp.diag(
        mu1, -alpha * mu1,
        mu2, -alpha * mu2, 1
    )

    e1 = sp.Matrix([1, 0, u0, u1, 0])
    e2 = sp.Matrix([0, 1, alpha*u1, u0, 0])

    f1 = sp.Matrix([1, 0, -u0, -u1, 0]) / (2*mu1)
    f2 = sp.Matrix([0, -1, alpha*u1, u0, 0]) / (2*alpha*mu1)

    g1 = sp.Matrix([0, 0, 0, 0, 1])

    # Required Witt basis order
    P = sp.Matrix.hstack(e1, e2, f1, f2, g1)

    gram = P.T * Q * P

    # N_alpha(u) = u0**2 - alpha*u1**2 = -mu1/mu2
    gram = gram.applyfunc(
        lambda x: sp.simplify(
            sp.expand(x).subs(
                u0**2, alpha*u1**2 - mu1/mu2
            )
        )
    )

    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    print("\nCalculated Gram matrix:")
    sp.pprint(gram)

    assert gram == expected
