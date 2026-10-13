"""Independent norm-witness checks for Case 4.3."""

import pytest
import sympy as sp

epsilon, pi = sp.symbols("epsilon pi", nonzero=True)


@pytest.mark.parametrize(
    "name,A,B,M2,coefficient",
    [
        (
            "exact",
            pi,
            epsilon*pi,
            epsilon,
            sp.I/epsilon,
        ),
        (
            "reverse",
            epsilon*pi,
            pi,
            epsilon,
            sp.I,
        ),
        (
            "epsilon",
            epsilon*pi,
            epsilon,
            pi,
            sp.I,
        ),
    ],
)
def test_pure_quadratic_witnesses(
    name, A, B, M2, coefficient
):
    # Here -1 is a square in F.
    # s = coefficient * sqrt(B).
    # Its quadratic norm is -B*coefficient**2.
    norm_s = -B*coefficient**2

    assert sp.simplify(M2*norm_s - A) == 0, name


def test_unit_witness():
    # A=epsilon, B=epsilon*pi, M2=epsilon.
    # s=1, so N_B(s)=1.
    A = epsilon
    M2 = epsilon
    norm_s = sp.Integer(1)

    assert sp.simplify(M2*norm_s - A) == 0


def test_quartic_witness():
    # Here epsilon=-1 and -1 is not a square in F.
    # s=s0**2+i*t0**2, with i**2=-1.
    # Lemma 2.1.5 supplies s0**4+t0**4=-1.
    s0, t0 = sp.symbols("s0 t0")

    norm_s = s0**4 + t0**4
    norm_s = norm_s.subs(s0**4, -1-t0**4)

    A = -pi
    M2 = pi

    assert sp.simplify(M2*norm_s - A) == 0
