

import streamlit as st
import sympy as sp
from fractions import Fraction


def math_latex(value):
    if isinstance(value, Fraction):
        value = sp.Rational(
            value.numerator,
            value.denominator
        )

    return sp.latex(value)


# -------------------------------------------------
# Display the building computation
# -------------------------------------------------

def display_building(result):
    """
    Display the mathematical details of a building
    computation in Streamlit.
    """

    st.divider()
    st.subheader("Bruhat–Tits Building")

    # 1. Parameter valuations
    st.markdown("### 1. Parameter valuations")

    parameter_valuations = result.get(
        "parameter_valuations", {}
    )

    for name, value in parameter_valuations.items():
        st.latex(
            r"\nu(" + name + ")="
            + math_latex(value)
        )

    # 2. Generator depths
    st.markdown("### 2. Toral generators and depths")

    names = result.get("generator_names", [])

    for j, depth in enumerate(result["depths"]):
        name = names[j] if j < len(names) else str(j+1)

        st.latex(
            "X_{" + str(j+1) + "}=" + name
            + r",\qquad r_{" + str(j+1) + "}="
            + math_latex(depth)
        )

    # 3. Root coefficient valuations
    st.markdown("### 3. Root coefficients and valuations")

    for j, coefficients in enumerate(
        result["root_coefficients"]
    ):

        with st.expander(
            f"Generator X{j+1}: root coefficients"
        ):

            rows = []

            for root, values in coefficients.items():

                valuations = result["valuations"][j][root]

                rows.append({
                    "Root": root,
                    "Positive coefficient":
                        str(values["positive"]),
                    "Negative coefficient":
                        str(values["negative"]),
                    "Positive valuation":
                        str(valuations["positive"]),
                    "Negative valuation":
                        str(valuations["negative"]),
                })

            st.table(rows)

    # 4. Moy-Prasad inequalities
    st.markdown("### 4. Moy–Prasad inequalities")

    st.latex(
        r"r-\nu(c_\gamma(X))"
        r"\leq\gamma(x)\leq"
        r"\nu(c_{-\gamma}(X))-r"
    )

    root_labels = {
        "epsilon1": "x_1",
        "epsilon2": "x_2",
        "epsilon1-epsilon2": "x_1-x_2",
        "epsilon1+epsilon2": "x_1+x_2",
    }

    for root, bounds in result[
        "combined_inequalities"
    ].items():

        lower = bounds["lower"]
        upper = bounds["upper"]
        label = root_labels[root]

        if lower == -sp.oo and upper == sp.oo:
            continue

        if lower == upper:
            st.latex(
                label + "=" + math_latex(lower)
            )
        else:
            st.latex(
                math_latex(lower)
                + r"\leq" + label
                + r"\leq" + math_latex(upper)
            )

    if result["diagonal_verified"]:
        st.success("Diagonal coefficient condition verified.")

    # 5. Building point
    st.markdown("### 5. Building point")

    x1, x2 = result["original_point"]

    st.latex(
        r"x=\left("
        + math_latex(x1) + ","
        + math_latex(x2)
        + r"\right)"
    )

    # 6. Fundamental alcove
    st.markdown("### 6. Fundamental alcove")

    st.latex(
        r"x_1\geq x_2\geq 0,"
        r"\qquad x_1+x_2\leq 1"
    )

    transformations = result["transformations"]

    if transformations:
        for transformation in transformations:
            st.write(transformation)
    else:
        st.info("The point is already in the fundamental alcove.")

    y1, y2 = result["alcove_point"]

    st.latex(
        r"\boxed{x_{\mathrm{alcove}}=\left("
        + math_latex(y1) + ","
        + math_latex(y2)
        + r"\right)}"
    )

    # 7. Fundamental alcove diagram
    st.markdown("### 7. Building diagram")

    fig = draw_fundamental_alcove(
        result["alcove_point"]
    )

    st.pyplot(fig)

    import matplotlib.pyplot as plt
    plt.close(fig)

# -------------------------------------------------
# Fundamental alcove diagram
# -------------------------------------------------

def draw_fundamental_alcove(point):
    import matplotlib.pyplot as plt

    # Computed building point
    x1, x2 = map(float, point)

    fig, ax = plt.subplots(figsize=(6, 5))

    # Vertices of the fundamental alcove
    vertices_x = [0, 1, 0.5, 0]
    vertices_y = [0, 0, 0.5, 0]

    # Draw the fundamental alcove
    ax.fill(
        vertices_x,
        vertices_y,
        color="lightblue",
        alpha=0.4
    )

    ax.plot(
        vertices_x,
        vertices_y,
        color="navy",
        linewidth=2
    )

    # Label the three vertices
    ax.text(-0.03, -0.05, "(0,0)")
    ax.text(1.01, -0.05, "(1,0)")
    ax.text(0.52, 0.52, "(1/2,1/2)")

    # Label the three walls
    ax.text(0.25, 0.29, r"$x_1=x_2$",
            rotation=45, color="navy")

    ax.text(0.77, 0.29, r"$x_1+x_2=1$",
            rotation=-45, color="navy")

    ax.text(0.4, -0.08, r"$x_2=0$",
            color="navy")

    # Plot the computed building point
    ax.scatter(
        x1, x2,
        color="red",
        s=100,
        zorder=5
    )

    # Display exact coordinates
    label = f"({point[0]}, {point[1]})"

    ax.annotate(
        label,
        (x1, x2),
        xytext=(10, 12),
        textcoords="offset points",
        color="red",
        fontsize=12,
        fontweight="bold"
    )

    # Axes and layout
    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")

    ax.set_title(
        r"Fundamental Alcove of Split $\mathrm{SO}(5)$"
    )

    ax.set_xlim(-0.15, 1.2)
    ax.set_ylim(-0.12, 0.7)
    ax.set_aspect("equal")

    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()

    return fig
