"""Independent exact checks for all biquadratic intermediate-field choices."""
import itertools
import pytest
import sympy as sp

from biquadratic import (
    CLASSES, KElement, Vector, construct_biquadratic, parse_rational,
)


def original_trace_quadratic(v, alpha, beta, mu0, mu1, rho):
    """Compute q directly from Tr_{K/F}(mu*x*sigma(x)/4) using K multiplication.

    Since Tr_{K/F}(u0+u1*A+(u2+u3*A)*B)=4*u0,
    the trace term is the constant coefficient of mu*x*sigma(x).
    """
    x = v.x
    sigma_x = KElement(x.a, x.b, -x.c, -x.d)
    prod = KElement(mu0, mu1, 0, 0).multiply(
        x.multiply(sigma_x, alpha, beta), alpha, beta
    )
    # mu*x*sigma(x) lies in E1; its A-component is killed by Tr_{E1/F}.
    assert sp.cancel(prod.c) == 0
    assert sp.cancel(prod.d) == 0
    return sp.cancel(prod.a + rho**2*v.z**2)


def add_vectors(x, y):
    return Vector(KElement(*[u+v for u,v in zip(x.x.coords(),y.x.coords())]), x.z+y.z)


@pytest.mark.parametrize("square,alpha,beta", [
    (square, alpha, beta)
    for square in (True, False)
    for alpha, beta in itertools.permutations(CLASSES, 2)
])
@pytest.mark.parametrize("mode,t", [
    ("canonical", ("1", "0", "0", "0")),
    ("witness", ("1", "2", "1", "-1")),
    ("witness", ("1/2", "0", "3", "0")),
])
@pytest.mark.parametrize("order", [True, False])
def test_all_biquadratic_cases(square, alpha, beta, mode, t, order):
    result = construct_biquadratic(square, alpha, beta, mode, t, "2", order)
    assert result.verified
    vectors = [result.basis[n] for n in ("e1", "e2", "f1", "f2", "g1")]
    for i, v in enumerate(vectors):
        qv = original_trace_quadratic(v, result.alpha, result.beta, result.mu0, result.mu1, result.rho)
        assert sp.cancel(qv - (1 if i == 4 else 0)) == 0
        for j, w in enumerate(vectors):
            qw = original_trace_quadratic(w, result.alpha, result.beta, result.mu0, result.mu1, result.rho)
            qsum = original_trace_quadratic(add_vectors(v,w), result.alpha, result.beta,
                                             result.mu0, result.mu1, result.rho)
            assert sp.cancel((qsum-qv-qw)/2 - result.gram[i,j]) == 0


@pytest.mark.parametrize("beta", ["pi", "epsilon*pi"])
@pytest.mark.parametrize("order", [True, False])
def test_unramified_minus_one_nonsquare_special(beta, order):
    r = construct_biquadratic(False, "epsilon", beta, "special", rho_input="1", order_by_valuation=order)
    assert r.verified
    assert (r.mu0, r.mu1) == (1, 0)
    assert sp.simplify(r.alpha+1) == 0
    # The special vectors are independently checked using the trace expression.
    vectors = [r.basis[n] for n in ("e1", "e2", "f1", "f2", "g1")]
    for i, v in enumerate(vectors):
        for j, w in enumerate(vectors):
            qsum = original_trace_quadratic(add_vectors(v,w), r.alpha, r.beta, r.mu0, r.mu1, r.rho)
            qv = original_trace_quadratic(v, r.alpha, r.beta, r.mu0, r.mu1, r.rho)
            qw = original_trace_quadratic(w, r.alpha, r.beta, r.mu0, r.mu1, r.rho)
            assert sp.cancel((qsum-qv-qw)/2-r.gram[i,j]) == 0


def test_ordering_uses_normalized_square_classes():
    assert construct_biquadratic(True, "epsilon", "pi").swapped_pairs
    assert not construct_biquadratic(True, "pi", "epsilon*pi").swapped_pairs
    assert not construct_biquadratic(True, "epsilon*pi", "pi").swapped_pairs


def test_validation():
    with pytest.raises(ValueError):
        construct_biquadratic(True, "epsilon", "epsilon")
    with pytest.raises(ValueError):
        construct_biquadratic(True, "epsilon", "pi", mode="special")
    with pytest.raises(ValueError):
        construct_biquadratic(False, "pi", "epsilon", mode="special")
    with pytest.raises(ValueError):
        construct_biquadratic(True, "epsilon", "pi", rho_input="0")
    with pytest.raises(ValueError):
        construct_biquadratic(True, "epsilon", "pi", mode="witness", t_coefficients=("0","0","0","0"))
    with pytest.raises(ValueError):
        parse_rational("__import__('os').system('echo hacked')")
    with pytest.raises(ValueError):
        parse_rational("1/0")
