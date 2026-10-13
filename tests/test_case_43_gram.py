"""Independent symbolic Gram matrix verification for Case 4.3."""

import sympy as sp


def test_case_43_general_gram():
    A, B, M1, M2, c = sp.symbols(
        "A B M1 M2 c", nonzero=True
    )
    r0, r1, s0, s1 = sp.symbols("r0 r1 s0 s1")

    # r belongs to F(sqrt(A)); s belongs to F(sqrt(B)).
    r = (r0, r1)
    s = (s0, s1)

    def multiply_by_root(x, D):
        return (D * x[1], x[0])

    def scale(x, factor):
        return (factor * x[0], factor * x[1])

    def vector(x, y, z):
        return sp.Matrix([x[0], x[1], y[0], y[1], z])

    # The displayed construction uses the B-summand first
    # and the A-summand second.
    def q(v):
        x0, x1, y0, y1, z = v
        return (
            M2 * (x0**2 - B * x1**2)
            + M1 * (y0**2 - A * y1**2)
            + A * B * c**2 * z**2
        )

    def bilinear(v, w):
        return sp.expand((q(v + w) - q(v) - q(w)) / 2)

    e1 = vector(s, multiply_by_root(r, A), 0)
    e2 = vector(multiply_by_root(s, B), (0, 0), 1 / c)

    f1 = vector(
        scale(s, 1 / (2 * A)),
        scale(multiply_by_root(r, A), -1 / (2 * A)),
        0,
    )
    f2 = vector(
        scale(multiply_by_root(s, B), -1 / (2 * A * B)),
        (0, 0),
        1 / (2 * A * B * c),
    )
    g1 = vector((0, 0), r, 0)

    basis = [e1, e2, f1, f2, g1]

    # Norm conditions:
    # M1*N_A(r) = 1
    # M2*N_B(s) = A
    norm_relations = {
        r0**2: 1 / M1 + A * r1**2,
        s0**2: A / M2 + B * s1**2,
    }

    def reduce_norms(expression):
        return sp.simplify(
            sp.expand(expression).subs(norm_relations)
        )

    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    gram = sp.Matrix(
        5, 5,
        lambda i, j: reduce_norms(
            bilinear(basis[i], basis[j])
        ),
    )

    assert gram == expected

    # Independent check of all five quadratic values.
    assert [reduce_norms(q(v)) for v in basis] == [
        0, 0, 0, 0, 1
    ]


def test_case_43_non_swapped_general_gram():
    """Verify the original, non-swapped Case 4.3 vectors."""
    A, B, M1, M2, c = sp.symbols(
        "A B M1 M2 c", nonzero=True
    )
    r0, r1, s0, s1 = sp.symbols("r0 r1 s0 s1")

    # r = r0 + r1*sqrt(A)
    # s = s0 + s1*sqrt(B)
    r = (r0, r1)
    s = (s0, s1)

    def root(x, D):
        return (D*x[1], x[0])

    def scale(x, k):
        return (k*x[0], k*x[1])

    def vec(x, y, z):
        return sp.Matrix([x[0], x[1], y[0], y[1], z])

    # Original coordinate order:
    # F(sqrt(A)) + F(sqrt(B)) + F
    def q(v):
        a, b, d, e, z = v
        return (
            M1*(a**2 - A*b**2)
            + M2*(d**2 - B*e**2)
            + A*B*c**2*z**2
        )

    def pairing(v, w):
        return sp.expand((q(v+w) - q(v) - q(w))/2)

    # t = 1/c, so c^2*t^2 = 1.
    t = 1/c

    # Exact vectors printed by the non-swapped branch.
    e1 = vec(root(r, A), s, 0)
    e2 = vec((0, 0), root(s, B), t)

    f1 = vec(
        scale(root(r, A), -1/(2*A)),
        scale(s, 1/(2*A)),
        0,
    )
    f2 = vec(
        (0, 0),
        scale(root(s, B), -1/(2*A*B)),
        t/(2*A*B),
    )
    g1 = vec(r, (0, 0), 0)

    basis = [e1, e2, f1, f2, g1]

    # M1*N_A(r) = 1
    # M2*N_B(s) = A
    norm_relations = {
        r0**2: 1/M1 + A*r1**2,
        s0**2: A/M2 + B*s1**2,
    }

    def reduce_norms(expr):
        return sp.simplify(
            sp.expand(expr).subs(norm_relations)
        )

    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    gram = sp.Matrix(
        5, 5,
        lambda i, j: reduce_norms(
            pairing(basis[i], basis[j])
        ),
    )

    assert gram == expected

    assert [reduce_norms(q(v)) for v in basis] == [
        0, 0, 0, 0, 1
    ]
