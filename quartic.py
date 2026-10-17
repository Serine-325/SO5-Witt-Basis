
"""
Constructive Witt bases for quartic elliptic tori in SO(5).

Based on Section 4.3 of our mathematical paper.

Cases:
    1. Pure quartic extensions
    2. Exceptional unramified case: d = eta
    3. Exceptional unramified case: d = eta * pi

Target Gram matrix:
    [0, 0, 1, 0, 0]
    [0, 0, 0, 1, 0]
    [1, 0, 0, 0, 0]
    [0, 1, 0, 0, 0]
    [0, 0, 0, 0, 1]
"""

import sympy as sp


TARGET_GRAM = sp.Matrix([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1],
])

def pure_quartic_form(v, alpha, r):
    """
    Quadratic form for the pure quartic case.

    v = (x0, x1, x2, x3, z)
    d = r * sqrt(alpha)
    mu4 = 1
    lambda = -alpha * r**2
    """
    x0, x1, x2, x3, z = v

    return sp.expand(
        x0**2
        + alpha * x1**2
        - 2 * alpha * r * x2 * x3
        - alpha * r**2 * z**2
    )

def pure_quartic_first_pair(alpha, r):
    """Construct the first hyperbolic pair (e1, f1)."""
    e1 = sp.Matrix([0, 0, 1, 0, 0])
    f1 = sp.Matrix([0, 0, 0, -sp.S.One / (alpha * r), 0])
    return e1, f1


def polar_form(u, v, alpha, r):
    """Compute B_q(u,v) from the quadratic form."""
    return sp.simplify(
        (
            pure_quartic_form(u + v, alpha, r)
            - pure_quartic_form(u, alpha, r)
            - pure_quartic_form(v, alpha, r)
        ) / 2
    )

def pure_quartic_basis(alpha, r):
    """Construct the complete pure quartic Witt basis."""

    e1, f1 = pure_quartic_first_pair(alpha, r)

    e2 = sp.Matrix([
        0, 1, 0, 0, sp.S.One / r
    ])

    f2 = sp.Matrix([
        0,
        sp.S.One / (2 * alpha),
        0,
        0,
        -sp.S.One / (2 * alpha * r)
    ])

    g1 = sp.Matrix([1, 0, 0, 0, 0])

    return [e1, e2, f1, f2, g1]


def pure_quartic_gram(alpha, r):
    """Compute the Gram matrix of the constructed basis."""

    basis = pure_quartic_basis(alpha, r)

    return sp.Matrix([
        [
            polar_form(u, v, alpha, r)
            for v in basis
        ]
        for u in basis
    ])



def exceptional_quartic_form(v, s, t):
    """
    Quadratic form for the exceptional quartic case d = eta.

    Coordinates:
        v = (a0, a1, b0, b1, z)

    where:
        a = a0 + i*a1
        b = b0 + i*b1
        eta = s**2 + i*t**2

    Assumption:
        s**4 + t**4 = -1
    """
    a0, a1, b0, b1, z = v

    return sp.expand(
        a0**2 - a1**2
        - s**2 * (b0**2 - b1**2)
        + 2 * t**2 * b0 * b1
        - z**2
    )

def exceptional_quartic_basis(s, t):
    """Construct a Witt basis for d = eta."""

    e1 = sp.Matrix([1, 1, 0, 0, 0])

    e2 = sp.Matrix([
        0, 0, 0, sp.S.One / s, 1
    ])

    f1 = sp.Matrix([
        sp.Rational(1, 2),
        -sp.Rational(1, 2),
        0, 0, 0
    ])

    f2 = sp.Matrix([
        0, 0, 0,
        sp.S.One / (2 * s),
        -sp.Rational(1, 2)
    ])

    g1 = sp.Matrix([
        0, 0, s, -t**2 / s, 0
    ])

    return [e1, e2, f1, f2, g1]


def exceptional_quartic_gram(s, t):
    """Compute the Gram matrix for d = eta."""

    basis = exceptional_quartic_basis(s, t)

    def B(u, v):
        return sp.simplify(
            (
                exceptional_quartic_form(u + v, s, t)
                - exceptional_quartic_form(u, s, t)
                - exceptional_quartic_form(v, s, t)
            ) / 2
        )

    return sp.Matrix([
        [B(u, v) for v in basis]
        for u in basis
    ])

def exceptional_twisted_quartic_form(v, s, t, pi):
    """
    Quadratic form for d = eta*pi, mu4 = eta.

    Coordinates:
        v = (a0, a1, b0, b1, z)

    Assumption:
        s**4 + t**4 = -1
    """
    a0, a1, b0, b1, z = v

    return sp.expand(
        s**2 * (a0**2 - a1**2)
        - 2 * t**2 * a0 * a1
        - pi * (
            (s**4 - t**4) * (b0**2 - b1**2)
            - 4 * s**2 * t**2 * b0 * b1
        )
        - z**2
    )

def exceptional_twisted_quartic_basis(s, t, pi):
    """Construct the Witt basis for d = eta*pi."""

    A = s**2 + t**2
    C = s**2 - t**2
    D = s**4 + t**4

    e1 = sp.Matrix([
        0, 0, A / D, C / D, 0
    ])

    e2 = sp.Matrix([
        -sp.S.One / s, 0, 0, 0, 1
    ])

    f1 = sp.Matrix([
        0, 0,
        -C / (2 * pi * D),
        A / (2 * pi * D),
        0
    ])

    f2 = sp.Matrix([
        -sp.S.One / (2 * s),
        0, 0, 0,
        -sp.Rational(1, 2)
    ])

    g1 = sp.Matrix([
        t**2 / s, s, 0, 0, 0
    ])

    return [e1, e2, f1, f2, g1]


def exceptional_twisted_quartic_gram(s, t, pi):
    """Compute the Gram matrix for d = eta*pi."""

    basis = exceptional_twisted_quartic_basis(s, t, pi)

    def B(u, v):
        return sp.factor(
            (
                exceptional_twisted_quartic_form(u + v, s, t, pi)
                - exceptional_twisted_quartic_form(u, s, t, pi)
                - exceptional_twisted_quartic_form(v, s, t, pi)
            ) / 2
        )

    return sp.Matrix([
        [B(u, v) for v in basis]
        for u in basis
    ])
