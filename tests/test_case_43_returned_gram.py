"""Independent Gram-matrix verification of all returned Case 4.3 vectors.

Parse the CLI's printed triples as elements of the two quadratic summands,
not as arbitrary SymPy radicals. This verifies the ACTUAL constructor output
in the original coordinate order, including swapped and special branches.

B_q(x,y) = (q(x+y)-q(x)-q(y))/2 is the convention used here.
"""

import ast

import pytest
import sympy as sp

from test_case_43_integration import CASES
from web_runner import hilbert, run_constructor

EPS = sp.Symbol("epsilon", nonzero=True)
PI = sp.Symbol("pi", nonzero=True)
R0, R1, S0, S1, T0 = sp.symbols("r0 r1 s0 s1 t0")


def _parse(expr, D, epsilon, names=None):
    """Parse a limited arithmetic grammar into F(sqrt(D)) coefficient pairs.

    For D=None, return a scalar in F. No Python eval/sympify of CLI text.
    """
    names = dict(names or {})
    names.update({"epsilon": epsilon, "pi": PI, "s0": S0, "t0": T0})

    def scalar(x):
        return (x, sp.S.Zero) if D is not None else x

    def pair(x):
        return x if isinstance(x, tuple) else scalar(x)

    def add(x, y):
        if D is None:
            return x + y
        x, y = pair(x), pair(y)
        return x[0] + y[0], x[1] + y[1]

    def mul(x, y):
        if D is None:
            return x * y
        x, y = pair(x), pair(y)
        return x[0]*y[0] + D*x[1]*y[1], x[0]*y[1] + x[1]*y[0]

    def div(x, y):
        if D is None:
            return x/y
        x, y = pair(x), pair(y)
        assert y[1] == 0, f"Non-scalar denominator in {expr!r}"
        return x[0]/y[0], x[1]/y[0]

    def negate(x):
        return mul(scalar(-1), x)

    def parse(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return scalar(sp.sympify(node.value))
        if isinstance(node, ast.Name):
            if node.id == "i":
                if epsilon == -1:
                    assert D == -1, "i=sqrt(epsilon) outside its quadratic summand"
                    return (sp.S.Zero, sp.S.One)
                return scalar(sp.I)
            if node.id in names:
                return pair(names[node.id]) if D is not None else names[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return negate(parse(node.operand))
            if isinstance(node.op, ast.UAdd):
                return parse(node.operand)
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow) and isinstance(node.right, ast.Constant) and node.right.value == 2:
                x = parse(node.left)
                return mul(x, x)
            x, y = parse(node.left), parse(node.right)
            if isinstance(node.op, ast.Add):
                return add(x, y)
            if isinstance(node.op, ast.Sub):
                return add(x, negate(y))
            if isinstance(node.op, ast.Mult):
                return mul(x, y)
            if isinstance(node.op, ast.Div):
                return div(x, y)
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "sqrt" and len(node.args) == 1
                and not node.keywords):
            argument = parse(node.args[0])
            assert D is not None and argument == (D, 0), (expr, argument, D)
            return sp.S.Zero, sp.S.One
        raise AssertionError(f"Unsupported constructor expression: {expr!r}")

    return parse(ast.parse(expr, mode="eval").body)


def _triple(text):
    parsed = ast.parse(text, mode="eval").body
    assert isinstance(parsed, ast.Tuple) and len(parsed.elts) == 3, text
    return [ast.unparse(node) for node in parsed.elts]


def _witness(expression, which, D, epsilon):
    if expression == which:
        return (R0, R1) if which == "r" else (S0, S1)
    return _parse(expression, D, epsilon)


def _vector(text, alpha, beta, mu1, mu2, minus_one_square, witnesses):
    epsilon = EPS if minus_one_square else sp.Integer(-1)
    classes = {"epsilon": epsilon, "pi": PI, "epsilon*pi": epsilon*PI}
    A, B = classes[alpha], classes[beta]
    swapped = hilbert(mu1, alpha, minus_one_square) == -1

    def coordinate(expr, D, which):
        w = _witness(witnesses[which], which, D, epsilon)
        return _parse(expr, D, epsilon, {which: w})

    x, y, z = _triple(text)
    # r belongs to the +1 Hilbert-sign summand, s to the -1 summand.
    first = coordinate(x, A, "s" if swapped else "r")
    second = coordinate(y, B, "r" if swapped else "s")
    last = _parse(z, None, epsilon)
    return first + second + (last,)


def _verify(params):
    minus_one_square, alpha, beta, mu1, mu2 = params
    result = run_constructor(*params)
    assert result.valid and result.case == "4.3"

    epsilon = EPS if minus_one_square else sp.Integer(-1)
    classes = {"1": sp.S.One, "epsilon": epsilon, "pi": PI, "epsilon*pi": epsilon*PI}
    A, B = classes[alpha], classes[beta]
    M1, M2 = classes[mu1], classes[mu2]
    lam = classes[result.lam]
    swapped = hilbert(mu1, alpha, minus_one_square) == -1
    norm_A = B if swapped else A
    norm_B = A if swapped else B
    norm_M1 = M2 if swapped else M1
    norm_M2 = M1 if swapped else M2

    basis = {name: _vector(result.basis[name], alpha, beta, mu1, mu2,
                           minus_one_square, result.witnesses)
             for name in ("e1", "e2", "f1", "f2", "g1")}

    substitutions = {}
    if result.witnesses["r"] == "r":
        substitutions[R0**2] = 1/norm_M1 + norm_A*R1**2
    if result.witnesses["s"] == "s":
        substitutions[S0**2] = norm_A/norm_M2 + norm_B*S1**2
    if "s0**2" in result.witnesses["s"]:
        assert not minus_one_square
        substitutions[S0**4] = -1 - T0**4

    def reduce(expr):
        expr = sp.expand(expr)
        expr = sp.expand(expr.subs(substitutions))
        return sp.expand(expr)

    def pairing(v, w):
        return reduce(
            M1*(v[0]*w[0] - A*v[1]*w[1])
            + M2*(v[2]*w[2] - B*v[3]*w[3])
            + lam*v[4]*w[4]
        )

    def q(v):
        return reduce(
            M1*(v[0]**2 - A*v[1]**2)
            + M2*(v[2]**2 - B*v[3]**2)
            + lam*v[4]**2
        )

    order = ("e1", "e2", "f1", "f2", "g1")
    expected = sp.Matrix([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ])
    gram = sp.Matrix(5, 5, lambda i, j: pairing(basis[order[i]], basis[order[j]]))
    assert gram == expected, f"{params}: Gram mismatch: {gram - expected}"
    assert [q(basis[name]) for name in order] == [0, 0, 0, 0, 1]


@pytest.mark.parametrize("params", CASES)
def test_case_43_actual_returned_gram(params):
    _verify(params)
