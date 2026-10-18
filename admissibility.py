
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