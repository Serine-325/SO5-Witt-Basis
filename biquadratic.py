"""Exact symbolic Witt-basis constructor for the biquadratic SO(5) case.

Mathematical setting (odd residual characteristic):
  K = F(A, B), A^2 = alpha, B^2 = beta, alpha != beta in square classes;
  E1 = F(A), sigma(A)=A, sigma(B)=-B;
  x = a + b*A + (c + d*A)*B;
  q_mu,lambda(x,z) = Tr_{K/F}((mu/4)*x*sigma(x)) + lambda*z^2.

This program works in the formal quadratic algebra using SymPy. It does NOT
solve arbitrary norm equations over an unspecified p-adic field. Instead:
  * canonical mode takes mu=A and t=1;
  * witness mode takes a supplied t, computes mu=A/N_{K/E1}(t);
  * special mode uses mu=1, alpha=-1 and a direct Witt basis.

The Gram matrix is computed from the ORIGINAL q_mu,lambda (not merely asserted).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re
from typing import Mapping

import sympy as sp

CLASSES = ("epsilon", "pi", "epsilon*pi")
A_NAME, B_NAME = sp.Symbol("A"), sp.Symbol("B")
EPSILON, PI = sp.symbols("epsilon pi", nonzero=True)
RHO = sp.Symbol("rho", nonzero=True)


def parse_rational(value: str | int | Fraction) -> sp.Rational:
    """Parse integers and fractions safely; never eval arbitrary user input."""
    if isinstance(value, (int, Fraction)):
        f = Fraction(value)
    else:
        s = str(value).strip()
        if not re.fullmatch(r"[+-]?\d+(?:\s*/\s*[+-]?\d+)?", s):
            raise ValueError("Enter an integer or fraction such as 1, -2, or 3/5.")
        try:
            f = Fraction(s.replace(" ", ""))
        except ZeroDivisionError as error:
            raise ValueError("Denominator must be nonzero.") from error
    return sp.Rational(f.numerator, f.denominator)


def class_expr(name: str, minus_one_square: bool) -> sp.Expr:
    if name not in CLASSES:
        raise ValueError(f"Unknown square class: {name}")
    epsilon = EPSILON if minus_one_square else sp.Integer(-1)
    return {"epsilon": epsilon, "pi": PI, "epsilon*pi": epsilon * PI}[name]


def class_product(a: str, b: str) -> str:
    """Product in F*/F*2 (only for the three nontrivial classes)."""
    if a == b:
        return "1"
    return ({*CLASSES} - {a, b}).pop()


def valuation_parity(square_class: str) -> int:
    """Valuation of canonical square-class representative modulo 2."""
    return 0 if square_class in ("1", "epsilon") else 1


@dataclass(frozen=True)
class KElement:
    """Coefficients (a,b,c,d) for a+b*A+(c+d*A)*B."""
    a: sp.Expr = sp.Integer(0)
    b: sp.Expr = sp.Integer(0)
    c: sp.Expr = sp.Integer(0)
    d: sp.Expr = sp.Integer(0)

    def coords(self) -> tuple[sp.Expr, ...]:
        return self.a, self.b, self.c, self.d

    def multiply(self, other: "KElement", alpha: sp.Expr, beta: sp.Expr) -> "KElement":
        a, b, c, d = self.coords()
        e, f, g, h = other.coords()
        return KElement(
            sp.expand(a * e + alpha * b * f + beta * (c * g + alpha * d * h)),
            sp.expand(a * f + b * e + beta * (c * h + d * g)),
            sp.expand(a * g + alpha * b * h + c * e + alpha * d * f),
            sp.expand(a * h + b * g + c * f + d * e),
        )

    def latex(self) -> str:
        a, b, c, d = self.coords()
        return sp.latex(a + b * A_NAME + (c + d * A_NAME) * B_NAME)


@dataclass(frozen=True)
class Vector:
    x: KElement
    z: sp.Expr = sp.Integer(0)

    def latex(self) -> str:
        return rf"\left({self.x.latex()},\,{sp.latex(self.z)}\right)"

    def coordinate_latex(self) -> str:
        return r"\left(" + r",\;".join(sp.latex(sp.factor(c)) for c in (*self.x.coords(), self.z)) + r"\right)"


@dataclass(frozen=True)
class BiquadraticResult:
    alpha_name: str
    beta_name: str
    alpha: sp.Expr
    beta: sp.Expr
    minus_one_square: bool
    mode: str
    mu0: sp.Expr
    mu1: sp.Expr
    t: KElement | None
    rho: sp.Expr
    basis: Mapping[str, Vector]
    gram: sp.Matrix
    norm_check: tuple[sp.Expr, sp.Expr] | None
    ordered_by_valuation: bool
    swapped_pairs: bool
    explanation: str

    @property
    def verified(self) -> bool:
        expected = sp.Matrix([
            [0, 1, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1],
        ])
        expected = expected.extract([0, 2, 1, 3, 4], [0, 2, 1, 3, 4])
        return all(sp.cancel(a - b) == 0 for a, b in zip(self.gram, expected)) and (
            self.norm_check is None or all(sp.cancel(x) == 0 for x in self.norm_check)
        )

    def mu_latex(self) -> str:
        return sp.latex(sp.factor(self.mu0) + sp.factor(self.mu1) * A_NAME)

    def q_latex(self) -> str:
        a, b, c, d, z = sp.symbols("a b c d z")
        return sp.latex(
            self.mu0 * (a*a + self.alpha*b*b - self.beta*c*c - self.alpha*self.beta*d*d)
            + 2*self.alpha*self.mu1*(a*b - self.beta*c*d)
            + self.rho**2*z*z
        )

    def report(self) -> str:
        lines = [
            "BIQUADRATIC SO(5) — EXACT WITT BASIS",
            f"alpha class: {self.alpha_name}  |  beta class: {self.beta_name}",
            f"-1 is a square in F: {self.minus_one_square}",
            f"mode: {self.mode}",
            "Coordinates: x=a+b*A+(c+d*A)*B; A^2=alpha, B^2=beta",
            f"alpha={self.alpha}; beta={self.beta}",
            f"mu = ({self.mu0}) + ({self.mu1})*A",
            f"lambda = ({self.rho})^2",
        ]
        if self.t is not None:
            lines.append(f"t coordinates (a,b,c,d): {self.t.coords()}")
            lines.append("Norm condition: mu * N_{K/E1}(t) = A")
        lines.extend([
            "",
            "WITT BASIS in the K-basis (1,A,B,AB), followed by z:",
        ])
        for name in ("e1", "f1", "e2", "f2", "g1"):
            v = self.basis[name]
            lines.append(f"{name} = ({', '.join(map(str, (*v.x.coords(), v.z)))})")
        lines.extend([
            "",
            "GRAM MATRIX in the order (e1,f1,e2,f2,g1):",
            str(self.gram),
            f"Verified: {self.verified}",
            f"Valuation ordering applied: {self.ordered_by_valuation}; pairs swapped: {self.swapped_pairs}",
            self.explanation,
            "",
            "NOTE: Formal symbolic identities; not an implementation of p-adic norm solving.",
        ])
        return "\n".join(lines)


def norm_components(t: KElement, alpha: sp.Expr, beta: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """N_{K/E1}(t) = U + V*A."""
    a, b, c, d = t.coords()
    U = sp.expand(a*a + alpha*b*b - beta*(c*c + alpha*d*d))
    V = sp.expand(2*(a*b - beta*c*d))
    return U, V


def bilinear(x: Vector, y: Vector, alpha: sp.Expr, beta: sp.Expr,
             mu0: sp.Expr, mu1: sp.Expr, rho: sp.Expr) -> sp.Expr:
    """B=(q(v+w)-q(v)-q(w))/2, evaluated from the ORIGINAL trace form."""
    a, b, c, d = x.x.coords()
    e, f, g, h = y.x.coords()
    return sp.cancel(
        mu0 * (a*e + alpha*b*f - beta*c*g - alpha*beta*d*h)
        + alpha*mu1 * (a*f + b*e - beta*(c*h + d*g))
        + rho*rho*x.z*y.z
    )


def construct_biquadratic(
    minus_one_square: bool,
    alpha_name: str,
    beta_name: str,
    mode: str = "canonical",
    t_coefficients: tuple[str, str, str, str] = ("1", "0", "0", "0"),
    rho_input: str = "1",
    order_by_valuation: bool = True,
) -> BiquadraticResult:
    """Construct and verify a Witt basis with exact symbolic arithmetic.

    Modes:
    canonical: mu=A, t=1, for all three intermediate fields.
    witness:   user supplies t in K; mu=A/N(t) is derived (never guessed).
    special:   -1 nonsquare, alpha=epsilon=-1, mu=1; direct simple basis.

    Only square lambda=rho^2 is supported, as required for split SO(5).
    """
    if alpha_name not in CLASSES or beta_name not in CLASSES or alpha_name == beta_name:
        raise ValueError("Choose distinct nontrivial square classes alpha and beta.")
    if mode not in ("canonical", "witness", "special"):
        raise ValueError("Unknown construction mode.")
    alpha = class_expr(alpha_name, minus_one_square)
    beta = class_expr(beta_name, minus_one_square)
    rho = parse_rational(rho_input)
    if rho == 0:
        raise ValueError("rho must be nonzero, since lambda=rho^2 is nonzero.")
    z1 = 1 / rho
    zero = KElement()
    if mode == "special":
        if minus_one_square or alpha_name != "epsilon":
            raise ValueError("The direct mu=1 special case requires -1 nonsquare and alpha=epsilon=-1.")
        mu0, mu1, t, norm_check = sp.Integer(1), sp.Integer(0), None, None
        basis = {
            "e1": Vector(KElement(1, 1, 0, 0)),
            "f1": Vector(KElement(sp.Rational(1,2), -sp.Rational(1,2), 0, 0)),
            "e2": Vector(KElement(0, 0, 1, 1)),
            "f2": Vector(KElement(0, 0, -1/(2*beta), 1/(2*beta))),
            "g1": Vector(zero, z1),
        }
        explanation = (
            "Special unramified case: alpha=-1, mu=1; "
            "q=a^2-b^2-beta*c^2+beta*d^2+rho^2*z^2. "
            "The simple basis is verified directly, without rescaling by an E1-scalar."
        )
    else:
        if mode == "canonical":
            t = KElement(1, 0, 0, 0)
        else:
            t = KElement(*(parse_rational(v) for v in t_coefficients))
        if all(v == 0 for v in t.coords()):
            raise ValueError("t must be nonzero.")
        U, V = norm_components(t, alpha, beta)
        D = sp.factor(U*U - alpha*V*V)
        if sp.simplify(D) == 0:
            raise ValueError("N_{K/E1}(t) is not invertible for the selected parameters.")
        # A/(U+V*A) = (-alpha*V + U*A)/(U^2-alpha*V^2).
        mu0 = sp.cancel(-alpha*V/D)
        mu1 = sp.cancel(U/D)
        # mu * N(t) = A, checked coefficientwise in E1.
        norm_check = (
            sp.cancel(mu0*U + alpha*mu1*V),
            sp.cancel(mu0*V + mu1*U - 1),
        )
        def tx(a: sp.Expr = 0, b: sp.Expr = 0, c: sp.Expr = 0, d: sp.Expr = 0) -> Vector:
            return Vector(t.multiply(KElement(a, b, c, d), alpha, beta))
        basis = {
            "e1": tx(a=1),
            "f1": tx(b=1/alpha),
            "e2": tx(c=1),
            "f2": tx(d=-1/(alpha*beta)),
            "g1": Vector(zero, z1),
        }
        explanation = (
            "General proposition: mu*N_{K/E1}(t)=A. "
            "Multiplication by t transports the Witt basis of "
            "q_{A,rho^2} to q_{mu,rho^2}. "
            "All pairings are also checked directly in the original trace form."
        )

    # Ordering compares canonical square classes of alpha and alpha*beta,
    # NOT the raw valuation of alpha*beta (which may include a square factor).
    swapped = False
    if order_by_valuation:
        second_class = class_product(alpha_name, beta_name)
        if valuation_parity(second_class) > valuation_parity(alpha_name):
            basis["e1"], basis["e2"] = basis["e2"], basis["e1"]
            basis["f1"], basis["f2"] = basis["f2"], basis["f1"]
            swapped = True

    order = ("e1", "e2", "f1", "f2", "g1")
    gram = sp.Matrix(5, 5, lambda i, j: bilinear(
        basis[order[i]], basis[order[j]], alpha, beta, mu0, mu1, rho
    ))
    result = BiquadraticResult(
        alpha_name, beta_name, alpha, beta, minus_one_square, mode,
        mu0, mu1, t, rho, basis, gram, norm_check,
        order_by_valuation, swapped, explanation,
    )
    if not result.verified:
        raise ArithmeticError("Symbolic Witt identities failed; no result will be returned.")
    return result
