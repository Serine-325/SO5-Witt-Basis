"""
Worked Example: Biquadratic Extension

Base field: F = Q_3
Extension: K = F(A, B), where A^2 = -1 and B^2 = 3

Parameters:
    alpha = -1
    beta = 3
    mu = 1
    lambda = 1

Quadratic form:
    q(x, z) = Tr_{K/F}(x * sigma(x) / 4) + z^2

Writing x = a + b*A + (c + d*A)*B, we obtain:
    q(x, z) = a^2 - b^2 - 3*c^2 + 3*d^2 + z^2

The Witt basis is written in two-coordinate notation:
    (x, z), where x belongs to K and z belongs to F.
"""

import sympy as sp

A, B = sp.symbols("A B")

# Gram matrix of the original quadratic form
# Coordinates are (a, b, c, d, z).
Q = sp.diag(1, -1, -3, 3, 1)

# Witt basis from the biquadratic construction
# e1 = (B(A+1), 0)
# e2 = (A+1, 0)
# f1 = (B(A-1)/6, 0)
# f2 = ((1-A)/2, 0)
# g1 = (0, 1)

e1 = sp.Matrix([0, 0, 1, 1, 0])
e2 = sp.Matrix([1, 1, 0, 0, 0])

f1 = sp.Matrix([
    0, 0, -sp.Rational(1, 6), sp.Rational(1, 6), 0
])

f2 = sp.Matrix([
    sp.Rational(1, 2), -sp.Rational(1, 2), 0, 0, 0
])

g1 = sp.Matrix([0, 0, 0, 0, 1])

# Construct the change-of-basis matrix
P = sp.Matrix.hstack(e1, e2, f1, f2, g1)

# Compute the Gram matrix
gram = sp.simplify(P.T * Q * P)

expected = sp.Matrix([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1],
])

assert P.det() != 0
assert gram == expected

# Norm-coset verification:
# -2 is a square in Q_3.
# If t^2 = -2, then v = (1-A)/t belongs to K.
# Since v is fixed by sigma, N_{K/F(A)}(v) = v^2 = A.

# Verify (1-A)^2 = -2*A using A^2 = -1.
# Represent elements of F(A) by coefficient pairs (a,b).
alpha = -1

def multiply_E1(x, y):
    a, b = x
    c, d = y
    return (a*c + alpha*b*d, a*d + b*c)

square = multiply_E1((1, -1), (1, -1))
norm_v = tuple(sp.Rational(c, -2) for c in square)

assert norm_v == (0, 1)  # Represents A

print("Biquadratic extension: K = Q_3(A, B)")
print("A^2 = -1, B^2 = 3")
print("mu = 1, lambda = 1")

print("\nWitt basis in two-coordinate notation:")
print("e1 = (B(A+1), 0)")
print("e2 = (A+1, 0)")
print("f1 = (B(A-1)/6, 0)")
print("f2 = ((1-A)/2, 0)")
print("g1 = (0, 1)")

print("\nCalculated Gram matrix:")
sp.pprint(gram)

print("\nNorm-coset verification:")
print("N_{K/F(A)}((1-A)/t) = A, where t^2 = -2")
print("Therefore [A] = [1] in the norm quotient.")

print("\nSUCCESS: The biquadratic Witt basis is verified.")