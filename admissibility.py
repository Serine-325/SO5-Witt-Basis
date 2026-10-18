
import sympy as sp


def pure_quartic_norm(alpha, r):
    """
    Compute N_{E1/F}(d) for d = r*sqrt(alpha).

    The nontrivial conjugation sends sqrt(alpha)
    to -sqrt(alpha).
    """
    return sp.simplify(-alpha * r**2)


def pure_quartic_discriminant_check(alpha, r):
    """
    Verify the chosen representative lambda = -alpha*r**2
    agrees with the norm D.

    This is an algebraic check, not a computation
    of local-field Hilbert symbols.
    """
    D = pure_quartic_norm(alpha, r)
    lam = -alpha * r**2

    return sp.simplify(D - lam) == 0


def hilbert_symbol_qp(a, b, p):
    """
    Compute the Hilbert symbol (a,b)_p over Q_p
    for nonzero rational a,b and an odd prime p.
    """
    a, b = sp.Rational(a), sp.Rational(b)
    p = int(p)

    if p == 2 or not sp.isprime(p):
        raise ValueError("p must be an odd prime")
    if a == 0 or b == 0:
        raise ValueError("a and b must be nonzero")

    def valuation(x):
        num, den = x.as_numer_denom()
        return (
            sp.factorint(abs(int(num))).get(p, 0)
            - sp.factorint(int(den)).get(p, 0)
        )

    m, n = valuation(a), valuation(b)
    u = a / sp.Integer(p)**m
    v = b / sp.Integer(p)**n

    def residue(x):
        num, den = x.as_numer_denom()
        return (int(num) * pow(int(den), -1, p)) % p

    return int(
        (-1)**(m * n * ((p - 1) // 2))
        * sp.legendre_symbol(residue(u), p)**n
        * sp.legendre_symbol(residue(v), p)**m
    )
def exceptional_quartic_norm(s, t):
    """Compute the norm of eta = s**2 + i*t**2, where i**2 = -1."""
    s, t = sp.sympify(s), sp.sympify(t)
    return sp.expand(s**4 + t**4)

def hensel_lift_fourth_root_minus_two(precision):
    """
    Find s modulo 3**precision such that
    s**4 + 2 == 0 modulo 3**precision.

    Start with s = 1 modulo 3.
    """
    if precision < 1:
        raise ValueError("Precision must be positive")

    s = 1

    for n in range(1, precision):
        modulus = 3**n
        next_modulus = 3**(n + 1)

        for digit in range(3):
            candidate = s + digit * modulus

            if (candidate**4 + 2) % next_modulus == 0:
                s = candidate
                break
        else:
            raise ArithmeticError("Hensel lifting failed")

    return s
from dataclasses import dataclass


@dataclass(frozen=True)
class QuadraticElement:
    """
    Represent a + b*sqrt(alpha) in F(sqrt(alpha)).
    """

    a: object
    b: object
    alpha: object

    def __add__(self, other):
        if self.alpha != other.alpha:
            raise ValueError("Different quadratic extensions")

        return QuadraticElement(
            sp.simplify(self.a + other.a),
            sp.simplify(self.b + other.b),
            self.alpha,
        )

    def __mul__(self, other):
        if self.alpha != other.alpha:
            raise ValueError("Different quadratic extensions")

        return QuadraticElement(
            sp.simplify(
                self.a * other.a
                + self.alpha * self.b * other.b
            ),
            sp.simplify(
                self.a * other.b + self.b * other.a
            ),
            self.alpha,
        )

    def conjugate(self):
        return QuadraticElement(
            self.a,
            -self.b,
            self.alpha,
        )

    def norm(self):
        return sp.simplify(
            self.a**2 - self.alpha * self.b**2
        )
def quadratic_residue_multiply(x, y, alpha, p):
    """
    Multiply elements in F_p[sqrt(alpha)].

    x = (a, b) represents a + b*sqrt(alpha).
    """
    a, b = x
    c, d = y

    return (
        (a*c + alpha*b*d) % p,
        (a*d + b*c) % p,
    )
def quadratic_character_f9(x):
    """
    Return +1 or -1 according to whether
    x is a square in F_9 = F_3[i].

    x = (a, b) represents a + b*i.
    """
    a, b = x
    x = (a % 3, b % 3)

    if x == (0, 0):
        raise ValueError("Quadratic character requires a nonzero element")

    result = (1, 0)

    for _ in range(4):
        result = quadratic_residue_multiply(
            result, x, -1, 3
        )

    if result == (1, 0):
        return 1

    if result == (2, 0):
        return -1

    raise ArithmeticError("Unexpected quadratic character value")
def hilbert_symbol_q3_i(x, y):
    """
    Compute the Hilbert symbol (x, y) over Q_3(i),
    where i**2 = -1.

    Inputs are nonzero QuadraticElement objects
    with rational coefficients and alpha = -1.
    """

    if x.alpha != -1 or y.alpha != -1:
        raise ValueError("This function requires Q_3(i)")

    def valuation_3(a):
        a = sp.Rational(a)

        if a == 0:
            return None

        numerator, denominator = a.as_numer_denom()
        numerator = abs(int(numerator))
        denominator = int(denominator)

        valuation = 0

        while numerator % 3 == 0:
            numerator //= 3
            valuation += 1

        while denominator % 3 == 0:
            denominator //= 3
            valuation -= 1

        return valuation

    def residue_mod_3(a):
        a = sp.Rational(a)
        numerator, denominator = a.as_numer_denom()

        return (
            int(numerator) % 3
            * pow(int(denominator) % 3, -1, 3)
        ) % 3

    def valuation_and_character(z):
        a = sp.Rational(z.a)
        b = sp.Rational(z.b)

        valuations = [
            v for v in (valuation_3(a), valuation_3(b))
            if v is not None
        ]

        if not valuations:
            raise ValueError("Hilbert symbol requires nonzero elements")

        m = min(valuations)

        scale = sp.Rational(3) ** m

        residue = (
            residue_mod_3(a / scale),
            residue_mod_3(b / scale),
        )

        return m, quadratic_character_f9(residue)

    m, chi_x = valuation_and_character(x)
    n, chi_y = valuation_and_character(y)

    return (
        (chi_x if n % 2 else 1)
        * (chi_y if m % 2 else 1)
    )
