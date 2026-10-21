"""
Worked Example: Sum of Two Equal Quadratic Extensions

Base field: F = Q_3
Quadratic extension: E = F(sqrt(-1))

Parameters:
    alpha = -1
    mu1 = 1
    mu2 = -1/2
    u = 1 + sqrt(-1)

Norm condition:
    N(u) = 1^2 - (-1)(1^2) = 2 = -mu1/mu2

Quadratic form:
    q = <1, 1, -1/2, -1/2, 1>

Construct the Witt basis (e1, e2, f1, f2, g1)
and verify its Gram matrix using exact SymPy arithmetic.
"""

import sympy as sp

alpha = sp.Integer(-1)
mu1 = sp.Integer(1)
mu2 = -sp.Rational(1, 2)

u0 = sp.Integer(1)
u1 = sp.Integer(1)

# Verify the required norm condition.
norm_u = u0**2 - alpha * u1**2
assert norm_u == -mu1 / mu2

# Gram matrix of the original quadratic form.
Q = sp.diag(
    mu1,
    -alpha * mu1,
    mu2,
    -alpha * mu2,
    1,
)

# Witt basis, using the formulas from the existing test.
e1 = sp.Matrix([1, 0, u0, u1, 0])
e2 = sp.Matrix([0, 1, alpha * u1, u0, 0])

f1 = sp.Matrix([1, 0, -u0, -u1, 0]) / (2 * mu1)

f2 = sp.Matrix(
    [0, -1, alpha * u1, u0, 0]
) / (2 * alpha * mu1)

g1 = sp.Matrix([0, 0, 0, 0, 1])

# Change-of-basis matrix.
P = sp.Matrix.hstack(e1, e2, f1, f2, g1)

# Calculate the Gram matrix.
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

print("Norm condition verified:", norm_u)
print("\nWitt basis:")
for name, vector in zip(
    ["e1", "e2", "f1", "f2", "g1"],
    [e1, e2, f1, f2, g1],
):
    print(f"{name} = {list(vector)}")

print("\nCalculated Gram matrix:")
sp.pprint(gram)

print("\nSUCCESS: The Witt basis is verified.")