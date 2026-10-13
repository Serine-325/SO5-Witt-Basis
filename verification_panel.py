"""Exact, fail-closed verification of the displayed SO(5) Witt bases.

The sum-of-quadratics verifier parses the *returned* constructor expressions
using a restricted AST grammar, represents both quadratic extensions as
coefficient pairs, and reduces the stated norm-witness relations. No eval.

The bilinear form is the half-polarization B_q(x,y)=(q(x+y)-q(x)-q(y))/2.
The expected Gram matrix has diagonal (0,0,0,0,1) in the order
(e1,e2,f1,f2,g1). These computations are formal; they do not solve norm
problems over a concrete local field.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any

import sympy as sp

from web_runner import WittResult, hilbert

EPS, PI = sp.symbols('epsilon pi', nonzero=True)
S0, T0, R0, R1, S1 = sp.symbols('s0 t0 r0 r1 s1')
S, T = sp.symbols('s t')
ORDER = ('e1', 'e2', 'f1', 'f2', 'g1')
EXPECTED = sp.Matrix([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1],
])


@dataclass(frozen=True)
class Verification:
    verified: bool
    gram: sp.Matrix
    norms: tuple[sp.Expr, ...]
    assumptions: tuple[str, ...]
    message: str


def _class(name: str, epsilon: sp.Expr) -> sp.Expr:
    return {'1': sp.S.One, 'epsilon': epsilon, 'pi': PI,
            'epsilon*pi': epsilon*PI}[name]


def _parse(expr: str, D: sp.Expr | None, epsilon: sp.Expr,
           names: dict[str, Any] | None = None) -> Any:
    """Return an F-scalar (D=None) or coefficients (a,b) for a+b*sqrt(D)."""
    names = dict(names or {})
    names.update({'epsilon': epsilon, 'pi': PI, 's0': S0, 't0': T0,
                  's': S, 't': T})

    def scalar(v):
        return v if D is None else (v, sp.S.Zero)

    def pair(v):
        return v if isinstance(v, tuple) else scalar(v)

    def add(x, y):
        if D is None:
            return x + y
        a, b = pair(x)
        c, d = pair(y)
        return a+c, b+d

    def mul(x, y):
        if D is None:
            return x*y
        a, b = pair(x)
        c, d = pair(y)
        return a*c+D*b*d, a*d+b*c

    def div(x, y):
        if D is None:
            return x/y
        a, b = pair(x)
        c, d = pair(y)
        denominator = sp.expand(c*c-D*d*d)
        if denominator == 0:
            raise ValueError(f'Zero denominator in {expr!r}')
        return ((a*c-D*b*d)/denominator, (b*c-a*d)/denominator)

    def parse(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return scalar(sp.Integer(node.value))
        if isinstance(node, ast.Name):
            if node.id == 'i':
                if epsilon == -1:
                    if D == -1:
                        return sp.S.Zero, sp.S.One
                    raise ValueError('i only permitted inside F(sqrt(-1)) or i*sqrt(pi)')
                return scalar(sp.I)
            if node.id in names:
                return pair(names[node.id]) if D is not None else names[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return mul(scalar(-1), parse(node.operand))
            if isinstance(node.op, ast.UAdd):
                return parse(node.operand)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'sqrt' and len(node.args) == 1 and not node.keywords:
            if D is None:
                raise ValueError('Square root not allowed in F-scalar coordinate')
            argument = parse(node.args[0])
            if argument != (D, sp.S.Zero):
                raise ValueError(f'Unexpected square root in {expr!r}: {argument}, expected {D}')
            return sp.S.Zero, sp.S.One
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow) and isinstance(node.right, ast.Constant) and type(node.right.value) is int and node.right.value >= 0:
                n = node.right.value
                base = parse(node.left)
                acc = scalar(1)
                for _ in range(n):
                    acc = mul(acc, base)
                return acc
            # For epsilon=-1, i*sqrt(pi) denotes sqrt(epsilon*pi) in the
            # exceptional Case 4.2 (choice of compatible square roots).
            if (epsilon == -1 and D == -PI and isinstance(node.op, ast.Mult)
                    and isinstance(node.left, ast.Name) and node.left.id == 'i'
                    and isinstance(node.right, ast.Call) and isinstance(node.right.func, ast.Name)
                    and node.right.func.id == 'sqrt' and len(node.right.args) == 1
                    and isinstance(node.right.args[0], ast.Name) and node.right.args[0].id == 'pi'):
                return sp.S.Zero, sp.S.One
            # In the exceptional (-1 nonsquare) case, i/sqrt(pi)
            # = sqrt(epsilon*pi)/pi for compatible square roots.
            if (epsilon == -1 and D == -PI and isinstance(node.op, ast.Div)
                    and ((isinstance(node.left, ast.Name) and node.left.id == 'i')
                         or (isinstance(node.left, ast.UnaryOp) and isinstance(node.left.op, ast.USub)
                             and isinstance(node.left.operand, ast.Name) and node.left.operand.id == 'i'))):
                sign = -1 if isinstance(node.left, ast.UnaryOp) else 1
                right = node.right
                if (isinstance(right, ast.BinOp) and isinstance(right.op, ast.Mult)
                        and isinstance(right.right, ast.Call)
                        and isinstance(right.right.func, ast.Name)
                        and right.right.func.id == 'sqrt' and len(right.right.args) == 1
                        and isinstance(right.right.args[0], ast.Name)
                        and right.right.args[0].id == 'pi'):
                    k = _parse(ast.unparse(right.left), None, epsilon, names)
                    return sp.S.Zero, sign/(k*PI)
            x, y = parse(node.left), parse(node.right)
            if isinstance(node.op, ast.Add):
                return add(x, y)
            if isinstance(node.op, ast.Sub):
                return add(x, mul(scalar(-1), y))
            if isinstance(node.op, ast.Mult):
                return mul(x, y)
            if isinstance(node.op, ast.Div):
                return div(x, y)
        raise ValueError(f'Unsupported constructor expression: {expr!r}')

    return parse(ast.parse(expr.replace('^', '**'), mode='eval').body)


def _triple(text: str) -> tuple[str, str, str]:
    """Extract three coordinates, allowing an outer scalar multiplier."""
    node = ast.parse(text.replace('^', '**'), mode='eval').body
    def triplet(n):
        if isinstance(n, ast.Tuple) and len(n.elts) == 3:
            return [ast.unparse(x) for x in n.elts]
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
            if isinstance(n.left, ast.Tuple):
                scale = ast.unparse(n.right)
                return [f'({scale})*({x})' for x in triplet(n.left)]
            if isinstance(n.right, ast.Tuple):
                scale = ast.unparse(n.left)
                return [f'({scale})*({x})' for x in triplet(n.right)]
        raise ValueError(f'Not a scalar times a three-coordinate vector: {text!r}')
    return tuple(triplet(node))


def verify_sum_result(result: WittResult) -> Verification:
    """Verify each of the five *returned* vectors for Cases 4.1, 4.2, 4.3.

    Raises ValueError on an unsupported printed expression; the UI must not
    claim verification if parsing or checking fails.
    """
    if not result.valid or result.case not in ('4.1', '4.2', '4.3'):
        raise ValueError('Only valid, identified sum-of-quadratics results are supported')
    epsilon = EPS if result.minus_one_square else sp.Integer(-1)
    A, B = _class(result.alpha, epsilon), _class(result.beta, epsilon)
    M1, M2 = _class(result.mu1, epsilon), _class(result.mu2, epsilon)
    lam = _class(result.lam, epsilon)
    swapped = (result.case == '4.3' and
               hilbert(result.mu1, result.alpha, result.minus_one_square) == -1)

    def coord(text: str, D: sp.Expr, summand: int):
        names = {}
        if result.case == '4.1':
            if 'u' not in result.witnesses:
                raise ValueError('Missing Case 4.1 norm witness')
            names['u'] = _parse(result.witnesses['u'], D, epsilon)
        if result.case == '4.3':
            which = ('s' if swapped else 'r') if summand == 0 else ('r' if swapped else 's')
            witness = result.witnesses.get(which)
            if witness is None:
                raise ValueError(f'Missing {which} witness')
            names[which] = ((R0, R1) if which == 'r' else (S0, S1)) if witness == which else _parse(witness, D, epsilon)
        return _parse(text, D, epsilon, names)

    vectors = {}
    for name in ORDER:
        x, y, z = _triple(result.basis[name])
        vectors[name] = (*coord(x, A, 0), *coord(y, B, 1), _parse(z, None, epsilon))

    substitutions: dict[sp.Expr, sp.Expr] = {}
    assumptions = []
    if result.case == '4.1':
        if 's^2' in result.witnesses['u']:
            substitutions[S**4] = -1-T**4
            assumptions.append('s⁴ + t⁴ = −1 (quartic norm witness)')
        assumptions.append('N_α(u) = −μ₁/μ₂ (explicit witness)')
    elif result.case == '4.2' and any('s^2' in v for v in result.basis.values()):
        substitutions[S**4] = -1-T**4
        assumptions.append('s⁴ + t⁴ = −1 (quartic witness)')
    elif result.case == '4.3':
        r = result.witnesses['r']
        s = result.witnesses['s']
        norm_A, norm_B = (B, A) if swapped else (A, B)
        norm_M1, norm_M2 = (M2, M1) if swapped else (M1, M2)
        if r == 'r':
            substitutions[R0**2] = 1/norm_M1 + norm_A*R1**2
            assumptions.append('μ₊ N(r) = 1 (symbolic norm witness)')
        if s == 's':
            substitutions[S0**2] = norm_A/norm_M2 + norm_B*S1**2
            assumptions.append('μ₋ N(s) = α₊ (symbolic norm witness)')
        if 's0**2' in s:
            substitutions[S0**4] = -1-T0**4
            assumptions.append('s₀⁴ + t₀⁴ = −1 (quartic witness)')
        if not assumptions:
            assumptions.append('Exact displayed norm witnesses')

    def reduce(expr):
        expr = sp.expand(expr)
        # Substitution may need more than one pass when square reductions
        # generate fourth powers in exceptional cases.
        for _ in range(3):
            expr = sp.expand(expr.subs(substitutions))
        return sp.cancel(expr)

    def Bq(v, w):
        return reduce(M1*(v[0]*w[0]-A*v[1]*w[1]) +
                      M2*(v[2]*w[2]-B*v[3]*w[3]) + lam*v[4]*w[4])

    gram = sp.Matrix(5, 5, lambda i, j: Bq(vectors[ORDER[i]], vectors[ORDER[j]]))
    norms = tuple(gram[i, i] for i in range(5))
    verified = all(reduce(gram[i, j]-EXPECTED[i, j]) == 0 for i in range(5) for j in range(5))
    return Verification(
        verified=verified, gram=gram, norms=norms,
        assumptions=tuple(assumptions),
        message=('All 25 half-polarized Gram entries agree with 2H ⊕ ⟨1⟩.'
                 if verified else 'The returned basis did not pass the symbolic Gram check.'),
    )


def render_sum_verification(result: WittResult) -> None:
    """Render a Streamlit panel after a sum-of-quadratics result."""
    import streamlit as st
    with st.expander('Mathematical verification of this Witt basis', expanded=True):
        st.caption('Formal verification of the actual returned vectors using the stated norm relations.')
        try:
            verification = verify_sum_result(result)
        except (ValueError, TypeError, KeyError, SyntaxError) as error:
            st.warning('Verification unavailable for this expression. No success claim is made.')
            st.code(str(error), language='text')
            return
        if verification.verified:
            st.success(verification.message)
        else:
            st.error(verification.message)
        st.latex(r'B_q(v,w)=\frac{q(v+w)-q(v)-q(w)}{2}')
        st.markdown('**Quadratic values** (order: e₁, e₂, f₁, f₂, g₁)')
        st.latex(r'\left(q(e_1),q(e_2),q(f_1),q(f_2),q(g_1)\right)='
                 + sp.latex(sp.Tuple(*verification.norms)))
        st.markdown('**Computed Gram matrix** (same order)')
        st.latex(sp.latex(verification.gram))
        st.markdown('**Expected Gram matrix**')
        st.latex(sp.latex(EXPECTED))
        st.markdown('**Norm conditions used**')
        for assumption in verification.assumptions:
            st.write('• ' + assumption)
        st.caption('A symbolic check under the stated assumptions, not a numerical '
                   'proof that a norm equation is solvable over a particular local field.')


def biquad_gram_in_display_order(result) -> sp.Matrix:
    """Permute the biquadratic constructor's internal Gram order for the UI.

    The established biquadratic constructor stores its Gram matrix in
    (e1,f1,e2,f2,g1) order. We deliberately leave that constructor and its
    existing tests unchanged, and present (e1,e2,f1,f2,g1) consistently.
    """
    return result.gram


def render_biquad_verification(result) -> None:
    """Show the existing exact biquadratic verification in the agreed order."""
    import streamlit as st
    gram = biquad_gram_in_display_order(result)
    correct = result.verified and all(sp.cancel(a - b) == 0 for a, b in zip(gram, EXPECTED))
    with st.expander('Mathematical verification of this Witt basis', expanded=True):
        if correct:
            st.success('All 25 Gram entries were checked symbolically against the original trace form.')
        else:
            st.error('Symbolic verification did not pass.')
        st.latex(r'B_q(v,w)=\frac{q(v+w)-q(v)-q(w)}{2}')
        st.markdown('**Quadratic values** (order: e₁, e₂, f₁, f₂, g₁)')
        st.latex(sp.latex(sp.Tuple(*(gram[i, i] for i in range(5)))))
        st.markdown('**Computed Gram matrix** (order: e₁, e₂, f₁, f₂, g₁)')
        st.latex(sp.latex(gram))
        st.markdown('**Expected Gram matrix** (same order)')
        st.latex(sp.latex(EXPECTED))
        if result.t is not None:
            st.latex(r'\mu N_{K/E_1}(t)=A')
            st.caption('This identity is checked coefficientwise by the biquadratic constructor.')
        st.caption('These are formal identities in the specified quadratic algebra; '
                   'no arbitrary local-field norm equation is solved here.')
