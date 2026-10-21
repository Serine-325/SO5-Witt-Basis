import sympy as sp

from biquadratic import construct_biquadratic


def test_biquadratic_special_matches_published_example():
    r = construct_biquadratic(
        False, "epsilon", "pi", mode="special"
    )

    assert r.verified
    assert (r.mu0, r.mu1) == (1, 0)

    pi = r.beta
    half = sp.Rational(1, 2)

    expected = {
        "e1": (0, 0, 1, 1, 0),
        "e2": (1, 1, 0, 0, 0),
        "f1": (0, 0, -half / pi, half / pi, 0),
        "f2": (half, -half, 0, 0, 0),
        "g1": (0, 0, 0, 0, 1),
    }

    for name, coordinates in expected.items():
        v = r.basis[name]
        actual = (v.x.a, v.x.b, v.x.c, v.x.d, v.z)

        assert all(
            sp.simplify(a - b) == 0
            for a, b in zip(actual, coordinates)
        ), f"Incorrect vector: {name}"
