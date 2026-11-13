
from fractions import Fraction


def in_fundamental_alcove(x1, x2):
    """
    Check whether (x1, x2) lies in the closure of
    the fundamental alcove of split SO(5).

    Conditions:
        x1 >= x2 >= 0
        x1 + x2 <= 1
    """
    x1 = Fraction(x1)
    x2 = Fraction(x2)

    return x1 >= x2 >= 0 and x1 + x2 <= 1


def fundamental_alcove(x1, x2):
    """
    Move a point into the fundamental alcove of SO(5)
    using affine Weyl reflections.

    Returns:
        point: (x1, x2) in the fundamental alcove
        transformations: list of reflections used
    """
    x1 = Fraction(x1)
    x2 = Fraction(x2)

    transformations = []

    while not in_fundamental_alcove(x1, x2):

        # Reflection in the wall x2 = 0
        if x2 < 0:
            x2 = -x2
            transformations.append("Reflect in x2 = 0")

        # Reflection in the wall x1 = x2
        elif x1 < x2:
            x1, x2 = x2, x1
            transformations.append("Exchange hyperbolic pairs")

        # Affine reflection in the wall x1 + x2 = 1
        elif x1 + x2 > 1:
            x1, x2 = 1 - x2, 1 - x1
            transformations.append("Reflect in x1 + x2 = 1")

    return (x1, x2), transformations


# -------------------------------------------------
# Root coefficients in the normalized Witt basis
# Basis order: (e1, e2, f1, f2, g1)
# -------------------------------------------------

def root_coefficients(M):
    """
    Extract the positive and negative root coefficients
    of X from its matrix M in the normalized Witt basis.

    Convention:
        M[j, i] is the coefficient of b_j in X(b_i).

    Returns:
        Dictionary of positive and negative root
        coefficients for the four positive roots of B2.
    """

    return {
        "epsilon1": {
            "positive": M[4, 2],
            "negative": M[4, 0],
        },
        "epsilon2": {
            "positive": M[4, 3],
            "negative": M[4, 1],
        },
        "epsilon1-epsilon2": {
            "positive": M[3, 2],
            "negative": M[1, 0],
        },
        "epsilon1+epsilon2": {
            "positive": M[1, 2],
            "negative": M[3, 0],
        }
    }

# -------------------------------------------------
# Valuation of symbolic monomial coefficients
# -------------------------------------------------

def monomial_valuation(expr, varpi, p=None):
    """
    Compute the valuation of a symbolic monomial
    of the form c * varpi**n.

    Assumptions:
        nu(varpi) = 1
        nu(0) = infinity

    If p is supplied, rational coefficients are
    evaluated using the p-adic valuation.

    If p is None, nonzero rational coefficients
    are assumed to be units.
    """
    import sympy as sp

    expr = sp.cancel(expr)

    if expr == 0:
        return sp.oo

    powers = expr.as_powers_dict()
    exponent = powers.get(varpi, 0)

    coefficient = sp.cancel(expr / varpi**exponent)

    if not coefficient.is_Rational:
        raise ValueError(
            f"Unsupported monomial coefficient: {expr}"
        )

    valuation = Fraction(exponent)

    if p is not None:
        numerator, denominator = coefficient.as_numer_denom()

        numerator = abs(int(numerator))
        denominator = abs(int(denominator))

        while numerator % p == 0:
            valuation += 1
            numerator //= p

        while denominator % p == 0:
            valuation -= 1
            denominator //= p

    return valuation

# -------------------------------------------------
# Moy-Prasad root inequalities
# -------------------------------------------------

def root_inequalities(coefficients, depth, valuation):
    """
    Compute the Moy-Prasad bounds for each positive root.

    Formula:
        r - nu(c_positive) <= gamma(x)
                          <= nu(c_negative) - r

    Parameters:
        coefficients: dictionary of root coefficients
        depth: depth r of the toral element
        valuation: function computing nu of a coefficient

    Returns:
        Dictionary of lower and upper bounds.
    """
    import sympy as sp

    r = Fraction(depth)
    inequalities = {}

    for root, values in coefficients.items():

        positive = values["positive"]
        negative = values["negative"]

        v_positive = valuation(positive)
        v_negative = valuation(negative)

        if v_positive == sp.oo:
            lower = -sp.oo
        else:
            lower = r - v_positive

        if v_negative == sp.oo:
            upper = sp.oo
        else:
            upper = v_negative - r

        inequalities[root] = {
            "lower": lower,
            "upper": upper,
        }

    return inequalities

# -------------------------------------------------
# Combine inequalities from toral generators
# -------------------------------------------------

def combine_inequalities(*systems):
    """
    Combine root inequalities from multiple generators.

    For each root:
        lower = maximum of all lower bounds
        upper = minimum of all upper bounds

    Raises ValueError if the bounds are inconsistent.
    """

    if not systems:
        raise ValueError("No inequality systems supplied.")

    combined = {}

    for root in systems[0]:

        lower = max(
            system[root]["lower"]
            for system in systems
        )

        upper = min(
            system[root]["upper"]
            for system in systems
        )

        if lower > upper:
            raise ValueError(
                f"Inconsistent inequalities for {root}"
            )

        combined[root] = {
            "lower": lower,
            "upper": upper,
        }

    return combined

# -------------------------------------------------
# Solve the building root inequalities
# -------------------------------------------------

def solve_building_point(inequalities):
    """
    Solve the root inequalities in the apartment of SO(5).

    Returns:
        (x1, x2) as exact rational numbers.

    Raises ValueError if the inequalities do not
    determine a unique point.
    """
    import sympy as sp
    from itertools import combinations

    roots = {
        "epsilon1": (1, 0),
        "epsilon2": (0, 1),
        "epsilon1-epsilon2": (1, -1),
        "epsilon1+epsilon2": (1, 1),
    }

    x1, x2 = sp.symbols("x1 x2")

    boundaries = []

    # Collect all finite boundary equations
    for root, bounds in inequalities.items():

        a, b = roots[root]
        expression = a*x1 + b*x2

        for bound in [bounds["lower"], bounds["upper"]]:

            if bound not in [sp.oo, -sp.oo]:
                equation = sp.Eq(expression, bound)

                if equation not in boundaries:
                    boundaries.append(equation)

    candidates = set()

    # Intersect every pair of boundary lines
    for eq1, eq2 in combinations(boundaries, 2):

        solution = sp.solve(
            [eq1, eq2],
            [x1, x2],
            dict=True,
        )

        for sol in solution:

            if x1 not in sol or x2 not in sol:
                continue

            u = Fraction(sol[x1])
            v = Fraction(sol[x2])

            # Check all root inequalities
            valid = True

            for root, bounds in inequalities.items():

                a, b = roots[root]
                value = a*u + b*v

                if not (
                    bounds["lower"]
                    <= value
                    <= bounds["upper"]
                ):
                    valid = False
                    break

            if valid:
                candidates.add((u, v))

    if len(candidates) != 1:
        raise ValueError(
            "The inequalities do not determine "
            "a unique building point."
        )

    point = candidates.pop()

    # Check that no nonzero feasible direction exists
    directions = [
        (1, 0), (-1, 0),
        (0, 1), (0, -1),
        (1, 1), (-1, -1),
        (1, -1), (-1, 1),
    ]

    for d1, d2 in directions:

        feasible = True

        for root, bounds in inequalities.items():

            a, b = roots[root]

            value = a*point[0] + b*point[1]
            change = a*d1 + b*d2

            if value == bounds["lower"] and change < 0:
                feasible = False
                break

            if value == bounds["upper"] and change > 0:
                feasible = False
                break

        if feasible:
            raise ValueError(
                "The building point is not unique."
            )

    return point

# -------------------------------------------------
# Verify the diagonal coefficient condition
# -------------------------------------------------

def check_diagonal_condition(M, depth, valuation):
    """
    Check that all diagonal coefficients of X
    have valuation at least depth r.
    """

    r = Fraction(depth)

    for i in range(5):

        coefficient = M[i, i]
        v = valuation(coefficient)

        if v < r:
            raise ValueError(
                f"Diagonal condition failed at index {i}: "
                f"valuation {v} < depth {r}"
            )

    return True

# -------------------------------------------------
# Complete building point computation
# -------------------------------------------------

def compute_building_point(generators, depths, valuation):
    """
    Compute the building point from toral generators.

    Parameters:
        generators: list of matrices in the Witt basis
        depths: list of generator depths
        valuation: valuation function on F

    Returns:
        Dictionary containing the full computation.
    """

    if len(generators) != len(depths):
        raise ValueError(
            "Each generator must have a depth."
        )

    if not generators:
        raise ValueError("No generators supplied.")

    all_coefficients = []
    all_valuations = []
    all_inequalities = []

    for M, r in zip(generators, depths):

        # Step 1: Verify diagonal coefficients
        check_diagonal_condition(M, r, valuation)

        # Step 2: Extract root coefficients
        coefficients = root_coefficients(M)

        # Step 3: Compute their valuations
        valuations = {
            root: {
                "positive": valuation(values["positive"]),
                "negative": valuation(values["negative"]),
            }
            for root, values in coefficients.items()
        }

        # Step 4: Construct root inequalities
        inequalities = root_inequalities(
            coefficients,
            r,
            valuation,
        )

        all_coefficients.append(coefficients)
        all_valuations.append(valuations)
        all_inequalities.append(inequalities)

    # Step 5: Combine all inequalities
    combined = combine_inequalities(
        *all_inequalities
    )

    # Step 6: Solve for the building point
    point = solve_building_point(combined)

    # Step 7: Move into the fundamental alcove
    alcove_point, transformations = fundamental_alcove(
        *point
    )

    return {
        "depths": [Fraction(r) for r in depths],
        "root_coefficients": all_coefficients,
        "valuations": all_valuations,
        "generator_inequalities": all_inequalities,
        "combined_inequalities": combined,
        "original_point": point,
        "alcove_point": alcove_point,
        "transformations": transformations,
        "diagonal_verified": True,
    }

# -------------------------------------------------
# Matrix of a toral generator in the Witt basis
# -------------------------------------------------

def matrix_in_witt_basis(L_X, P):
    """
    Compute the matrix of a toral generator X
    in the normalized Witt basis.

    Parameters:
        L_X: matrix of multiplication by X
             in the original algebra basis.

        P: matrix whose columns are the Witt
           basis vectors in the order
           (e1, e2, f1, f2, g1).

    Formula:
        M_X = P^(-1) * L_X * P
    """
    import sympy as sp

    L_X = sp.Matrix(L_X)
    P = sp.Matrix(P)

    if L_X.shape != (5, 5):
        raise ValueError("L_X must be a 5 x 5 matrix.")

    if P.shape != (5, 5):
        raise ValueError("P must be a 5 x 5 matrix.")

    if sp.simplify(P.det()) == 0:
        raise ValueError("The Witt basis matrix is singular.")

    M_X = P.inv() * L_X * P

    return M_X.applyfunc(sp.simplify)

# -------------------------------------------------
# Extract the biquadratic Witt basis matrix
# -------------------------------------------------

def biquadratic_witt_matrix(result):
    """
    Construct the change-of-basis matrix P from
    an existing BiquadraticResult.

    Original algebra basis:
        (1, A, B, AB, g1)

    Witt basis order:
        (e1, e2, f1, f2, g1)
    """
    import sympy as sp

    basis_order = ["e1", "e2", "f1", "f2", "g1"]

    columns = []

    for name in basis_order:

        vector = result.basis[name]

        column = sp.Matrix([
            vector.x.a,
            vector.x.b,
            vector.x.c,
            vector.x.d,
            vector.z,
        ])

        columns.append(column)

    P = sp.Matrix.hstack(*columns)

    return P

# -------------------------------------------------
# Biquadratic toral generator matrices
# -------------------------------------------------

def biquadratic_generator_matrices(result):
    """
    Construct the matrices of X1 = B and X2 = AB
    in the Witt basis returned by the existing
    biquadratic constructor.

    Original algebra basis:
        (1, A, B, AB, g1)

    Relations:
        A**2 = alpha
        B**2 = beta
    """
    import sympy as sp

    alpha = result.alpha
    beta = result.beta

    # Extract the actual Witt basis
    P = biquadratic_witt_matrix(result)

    # Multiplication by X1 = B
    L1 = sp.Matrix([
        [0, 0, beta,    0, 0],
        [0, 0,    0, beta, 0],
        [1, 0,    0,    0, 0],
        [0, 1,    0,    0, 0],
        [0, 0,    0,    0, 0],
    ])

    # Multiplication by X2 = AB
    L2 = sp.Matrix([
        [0,     0,    0, alpha*beta, 0],
        [0,     0, beta,          0, 0],
        [0, alpha,    0,          0, 0],
        [1,     0,    0,          0, 0],
        [0,     0,    0,          0, 0],
    ])

    # Change to the actual Witt basis
    M1 = matrix_in_witt_basis(L1, P)
    M2 = matrix_in_witt_basis(L2, P)

    return M1, M2

# -------------------------------------------------
# Valuation with the nonsquare unit epsilon
# -------------------------------------------------

def symbolic_valuation(expr, varpi, epsilon, p=None):
    """
    Compute the valuation of a symbolic monomial
    c * epsilon**m * varpi**n.

    Normalization:
        nu(epsilon) = 0
        nu(varpi) = 1
        nu(0) = infinity

    Rational coefficients are handled by
    monomial_valuation.
    """
    import sympy as sp

    expr = sp.cancel(expr)

    if expr == 0:
        return sp.oo

    # Extract the power of the unit epsilon
    powers = expr.as_powers_dict()
    epsilon_power = powers.get(epsilon, 0)

    # Remove the unit factor
    reduced = sp.cancel(
        expr / epsilon**epsilon_power
    )

    # Compute the remaining valuation
    return monomial_valuation(
        reduced,
        varpi,
        p=p,
    )

# -------------------------------------------------
# Complete biquadratic building computation
# -------------------------------------------------

def compute_biquadratic_building(result):
    """
    Compute the building point directly from an
    existing BiquadraticResult.

    Generators:
        X1 = B
        X2 = AB

    Depths:
        r1 = nu(beta) / 2
        r2 = nu(alpha * beta) / 2
    """
    from biquadratic import EPSILON, PI

    alpha = result.alpha
    beta = result.beta

    # Valuation on the base field
    def valuation(expr):
        return symbolic_valuation(
            expr, PI, EPSILON
        )

    # Compute parameter valuations
    v_alpha = valuation(alpha)
    v_beta = valuation(beta)

    # Compute generator depths
    r1 = v_beta / 2
    r2 = (v_alpha + v_beta) / 2

    # Construct matrices from the actual Witt basis
    M1, M2 = biquadratic_generator_matrices(result)

    # Compute the building point
    building_result = compute_building_point(
        [M1, M2],
        [r1, r2],
        valuation,
    )

    # Store additional mathematical information
    building_result["family"] = "biquadratic"
    building_result["alpha"] = alpha
    building_result["beta"] = beta

    building_result["parameter_valuations"] = {
        "alpha": v_alpha,
        "beta": v_beta,
    }

    building_result["generator_names"] = [
        "B",
        "AB",
    ]

    building_result["generator_matrices"] = [
        M1,
        M2,
    ]

    building_result["swapped_pairs"] = (
        result.swapped_pairs
    )

    return building_result
