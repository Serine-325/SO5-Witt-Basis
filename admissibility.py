
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
