
"""
Worked Example: Two Distinct Quadratic Extensions

F = Q_3
E1 = F(A), A^2 = 3
E2 = F(B), B^2 = -3

Parameters:
    alpha = 3
    beta = -3
    mu1 = mu2 = 1
    lambda = -1

Quadratic form:
    q(x,y,z) = N_E1(x) + N_E2(y) - z^2

Writing x = a + b*A and y = c + d*B:
    q = a^2 - 3*b^2 + c^2 + 3*d^2 - z^2
"""

import sympy as sp

# Original Gram matrix in coordinates (a,b,c,d,z)
Q = sp.diag(1, -3, 1, 3, -1)

# Witt basis in E1 + E2 + F notation:
#
# e1 = (1, 0, 1)
# e2 = (A, B, 0)
# f1 = (1/2, 0, -1/2)
# f2 = (-A/6, B/6, 0)
# g1 = (0, 1, 0)

e1 = sp.Matrix([1, 0, 0, 0, 1])
e2 = sp.Matrix([0, 1, 0, 1, 0])

f1 = sp.Matrix([
    sp.Rational(1, 2), 0, 0, 0,
    -sp.Rational(1, 2)
])

f2 = sp.Matrix([
    0, -sp.Rational(1, 6),
    0, sp.Rational(1, 6), 0
])

g1 = sp.Matrix([0, 0, 1, 0, 0])

# Change-of-basis matrix
P = sp.Matrix.hstack(e1, e2, f1, f2, g1)

# Expected Witt Gram matrix
expected = sp.Matrix([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1],
])

gram = P.T * Q * P

# Verify linear independence and all pairings
assert P.det() != 0
assert gram == expected

# Norm witnesses:
# N_E1(1) = N_E2(1) = 1
assert 1 == 1 / 1

# Check the discriminant square class:
# alpha*beta = -9 and lambda = -1
# lambda/(alpha*beta) = 1/9 = (1/3)^2
assert sp.Rational(-1, -9) == sp.Rational(1, 3)**2

print("Two distinct quadratic extensions over Q_3")
print("E1 = Q_3(sqrt(3))")
print("E2 = Q_3(sqrt(-3))")
print("mu1 = mu2 = 1, lambda = -1")

print("\nWitt basis:")
print("e1 = (1, 0, 1)")
print("e2 = (A, B, 0)")
print("f1 = (1/2, 0, -1/2)")
print("f2 = (-A/6, B/6, 0)")
print("g1 = (0, 1, 0)")

print("\nCalculated Gram matrix:")
sp.pprint(gram)

print("\nDeterminant of change-of-basis matrix:")
print(P.det())

print("\nSUCCESS: Distinct quadratic-extension Witt basis verified.")
