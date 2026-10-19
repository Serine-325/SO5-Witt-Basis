"""Streamlit interface for both established SO(5) Witt basis families."""
import streamlit as st
import sympy as sp
from verification_panel import render_sum_verification
from web_runner import CLASSES, MU_OPTIONS, run_constructor, valid_parameters
from biquadratic import construct_biquadratic

st.set_page_config(
    page_title="SO(5) | Explicit Witt Basis Constructor",
    page_icon="📐",
    layout="wide",
)

st.title("Explicit Witt Bases for Elliptic Tori in Split SO(5)")
st.caption("Research software · Two quadratic extensions and biquadratic extensions")
st.markdown(
    "Construct Witt bases from elliptic torus parameters and inspect the "
    "pairing identities. The biquadratic module performs exact symbolic "
    "verification of the Gram matrix."
)

tab_sum, tab_biquad, tab_quartic = st.tabs([
    "Sums of two quadratic extensions",
    "Biquadratic extension",
    "Quartic extensions",
])

with tab_sum:
    st.subheader("Sums of two quadratic extensions")
    st.latex(r"V=F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F")
    st.latex(r"q(x,y,z)=\mu_1N_{\alpha}(x)+\mu_2N_{\beta}(y)+\lambda z^2")
    st.caption("The existing Stage 1 constructor is used unchanged.")

    minus_one_square = st.radio(
        "Is −1 a square in F?", options=(True, False),
        format_func=lambda v: "Yes" if v else "No",
        horizontal=True, key="sum_minus_one_square",
    )
    left, right = st.columns(2)
    with left:
        alpha = st.selectbox("First quadratic class α", CLASSES, key="sum_alpha")
        mu1 = st.selectbox("First coefficient μ₁", MU_OPTIONS[alpha], key="sum_mu1")
    with right:
        beta = st.selectbox("Second quadratic class β", CLASSES, key="sum_beta")
        mu2 = st.selectbox("Second coefficient μ₂", MU_OPTIONS[beta], key="sum_mu2")
    params = (minus_one_square, alpha, beta, mu1, mu2)
    try:
        admissible = valid_parameters(*params)
    except ValueError as error:
        admissible = False
        st.error(str(error))
    if admissible:
        st.success("The split SO(5) Hilbert-symbol criterion is satisfied.")
    else:
        st.warning(
            "These choices do not satisfy "
            r"$(\mu_1,\alpha)(\mu_2,\beta)=(\alpha,\beta)$."
        )
    if st.button("Construct Stage 1 Witt Basis", type="primary", disabled=not admissible):
        try:
            st.session_state["sum_result"] = run_constructor(*params)
            st.session_state["sum_params"] = params
        except Exception as error:
            st.session_state.pop("sum_result", None)
            st.error(f"Construction failed: {error}")
    if st.session_state.get("sum_params") == params and "sum_result" in st.session_state:
        result = st.session_state["sum_result"]
        st.markdown(f"**Case {result.case}** · Square class of λ: `{result.lam}`")
        for name in ("e1", "e2", "f1", "f2", "g1"):
            st.code(f"{name} = {result.basis[name]}", language="text")
        if result.witnesses:
            st.markdown("**Norm witnesses**")
            for key, value in result.witnesses.items():
                st.code(f"{key} = {value}", language="text")
        with st.expander("Complete output and norm conditions"):
            st.code(result.raw_output, language="text")
        st.download_button(
            "Download Stage 1 result (.txt)", data=result.raw_output,
            file_name="SO5_Sum_Quadratics_Witt_Basis.txt", mime="text/plain",
        )
        render_sum_verification(result)
with tab_biquad:
    st.subheader("Biquadratic extension")
    st.latex(r"K=F(A,B),\quad A^2=\alpha,\quad B^2=\beta,\quad E_1=F(A)")
    st.latex(r"x=a+bA+(c+dA)B")
    st.latex(r"q_{\mu,\lambda}(x,z)=\operatorname{Tr}_{K/F}\!\left(\frac{\mu}{4}x\sigma(x)\right)+\lambda z^2")
    st.caption(
        "Choose distinct nontrivial square classes. Here σ fixes A and sends B to −B. "
        "All displayed vectors use the original coordinates (1, A, B, AB; z)."
    )

    square = st.radio(
        "Is −1 a square in F?", (True, False),
        format_func=lambda v: "Yes" if v else "No",
        horizontal=True, key="biq_square",
    )
    col1, col2 = st.columns(2)
    with col1:
        a_class = st.selectbox("Fixed field E₁ = F(√α)", CLASSES, key="biq_alpha")
    with col2:
        b_choices = [v for v in CLASSES if v != a_class]
        b_class = st.selectbox("Choose β (different square class)", b_choices, key="biq_beta")

    exceptional = not square and a_class == "epsilon"

    if exceptional:
        mu_options = {
            "one": "mu = 1",
            "eta": "mu = eta = s0^2 + i*t0^2",
        }
    else:
        mu_options = {
            "sqrt_alpha": "mu = sqrt(alpha)",
            "one": "mu = 1",
        }

    mu_choice = st.radio(
        "Choose mu",
        list(mu_options),
        format_func=lambda key: mu_options[key],
        key="biq_mu_choice",
    )

    is_split = mu_choice == ("one" if exceptional else "sqrt_alpha")
    mode = "special" if exceptional else "canonical"
    if mode == "witness":
        st.markdown(
            "Choose the four **F-coordinates** of a nonzero witness "
            r"$t=t_0+t_1A+(t_2+t_3A)B$. "
            "The program derives μ from this witness."
        )
        c0, c1, c2, c3 = st.columns(4)
        with c0:
            t0 = st.text_input("t₀", "1", key="biq_t0")
        with c1:
            t1 = st.text_input("t₁", "0", key="biq_t1")
        with c2:
            t2 = st.text_input("t₂", "0", key="biq_t2")
        with c3:
            t3 = st.text_input("t₃", "0", key="biq_t3")
        t_coordinates = (t0, t1, t2, t3)
    else:
        t_coordinates = ("1", "0", "0", "0")

    st.markdown(r"**Split line:** choose $\rho\in F^\times$ and set $\lambda=\rho^2$.")
    rho = st.text_input("ρ (integer or rational number)", "1", key="biq_rho")
    order_pairs = st.checkbox(
        "Order hyperbolic pairs by canonical valuation class (thesis convention)",
        value=True, key="biq_order",
    )

    if st.button("Construct and verify biquadratic Witt basis", type="primary"):
        try:
            if not is_split:
                raise ValueError(
                    "The selected mu gives the non-split quadratic space. Please choose the other mu."
                )
            result = construct_biquadratic(
                minus_one_square=square,
                alpha_name=a_class,
                beta_name=b_class,
                mode=mode,
                t_coefficients=t_coordinates,
                rho_input=rho,
                order_by_valuation=order_pairs,
            )
            st.session_state["biq_result"] = result
            st.session_state["biq_settings"] = (
                square, a_class, b_class, mu_choice, mode, t_coordinates, rho, order_pairs,
            )
        except (ValueError, ArithmeticError) as error:
            st.session_state.pop("biq_result", None)
            st.error(str(error))

    current_settings = (square, a_class, b_class, mu_choice, mode, t_coordinates, rho, order_pairs)
    if st.session_state.get("biq_settings") == current_settings and "biq_result" in st.session_state:
        result = st.session_state["biq_result"]
        st.divider()
        st.success("Exact symbolic verification passed: the Gram matrix is 2H ⊕ ⟨1⟩.")
        st.markdown("**Chosen parameters**")
        st.latex(r"\alpha=" + sp.latex(result.alpha) +
                 r",\quad \beta=" + sp.latex(result.beta))
        st.latex(r"\mu=" + result.mu_latex() +
                 r",\quad \lambda=" + sp.latex(result.rho**2))
        if result.t is not None:
            st.markdown("**Norm witness**")
            st.latex(r"t=" + result.t.latex())
            st.latex(r"\mu N_{K/E_1}(t)=A")
        st.markdown("**Witt basis**")
        for name in ("e1", "e2", "f1", "f2", "g1"):
            v = result.basis[name]
            st.latex(name.replace("1", "_1").replace("2", "_2") + "=" + v.latex())
            st.caption(f"{name} coordinates (a,b,c,d,z): " + str((*v.x.coords(), v.z)))
        if result.swapped_pairs:
            st.info("The two complete hyperbolic pairs were interchanged to put the ramified square class first.")
        st.markdown("**Verified Gram matrix** (order: e₁, e₂, f₁, f₂, g₁)")
        st.latex(sp.latex(result.gram))
        with st.expander("Show quadratic form and verification notes"):
            st.latex(r"q_{\mu,\lambda}(x,z)=" + result.q_latex())
            st.write(result.explanation)
            st.caption(
                "All computations are formal identities. The application does not "
                "implement arithmetic or norm-equation solving in a specified local field."
            )
        st.download_button(
            "Download biquadratic result (.txt)", data=result.report(),
            file_name="SO5_Biquadratic_Witt_Basis.txt", mime="text/plain",
        )

st.divider()
st.markdown(
    "**Research scope:** explicit Witt bases for elliptic maximal tori "
    "in split SO(5), including sums of two quadratic extensions, "
    "biquadratic extensions, and quartic extensions with a unique "
    "quadratic subfield. The application provides symbolic Gram-matrix "
    "verification. Bruhat–Tits building calculations are not implemented here."
)
st.markdown(
    "Source repository: [SO5-Witt-Basis](https://github.com/Serine-325/SO5-Witt-Basis)."
)

with tab_quartic:
    st.subheader("Quartic Extensions")

    st.markdown(
        "Constructive Witt bases for quartic elliptic tori "
        "in the split group SO(5)."
    )

    
    import sympy as sp

    minus_one_square = st.radio(
        "Is -1 a square in F?",
        ["Yes", "No"],
        horizontal=True,
        key="quartic_minus_one_square",
    )

    epsilon = sp.Symbol("epsilon", nonzero=True)
    varpi = sp.Symbol("varpi", nonzero=True)

    if minus_one_square == "Yes":
        subfield_options = [
            "F(sqrt(epsilon))",
            "F(sqrt(varpi))",
            "F(sqrt(epsilon*varpi))",
        ]
    else:
        subfield_options = [
            "F(sqrt(-1))",
            "F(sqrt(varpi))",
            "F(sqrt(epsilon*varpi))",
        ]

    subfield = st.selectbox(
        "Choose the quadratic subfield E1",
        subfield_options,
        key="quartic_subfield",
    )

    alpha_options = {
        "F(sqrt(epsilon))": epsilon,
        "F(sqrt(varpi))": varpi,
        "F(sqrt(epsilon*varpi))": epsilon * varpi,
    }
    alpha_label = {
        "F(sqrt(epsilon))": "epsilon",
        "F(sqrt(varpi))": "varpi",
        "F(sqrt(epsilon*varpi))": "epsilon*varpi",
    }.get(subfield, "alpha")      
    if subfield == "F(sqrt(-1))":
        d_options = ["eta", "eta*varpi"]
    elif minus_one_square == "Yes":
        if subfield == "F(sqrt(epsilon))":
            d_options = ["sqrt(epsilon)", "varpi*sqrt(epsilon)"]
        else:
            d_options = ["sqrt(alpha)", "epsilon*sqrt(alpha)"]
    else:
        d_options = ["sqrt(alpha)"]

    d_choice = st.selectbox(
        "Choose the quartic parameter d",
        d_options,
        format_func=lambda d: d.replace(
        "sqrt(alpha)", f"sqrt({alpha_label})"
        ),
        key="quartic_d",
    )

    if d_choice == "eta":
        quartic_case = "Exceptional unramified: d = eta"
    elif d_choice == "eta*varpi":
        quartic_case = "Exceptional twisted: d = eta*pi"
    else:
        quartic_case = "Pure quartic extension"
        alpha_value = alpha_options[subfield]

        # d = b * sqrt(alpha)
        # The current quartic constructor uses r_value as b.

        d_coefficients = {
            "sqrt(epsilon)": sp.Integer(1),
            "sqrt(alpha)": sp.Integer(1),
            "varpi*sqrt(epsilon)": varpi,
            "epsilon*sqrt(alpha)": epsilon,
        }

        r_value = d_coefficients[d_choice]

    st.info(
        "The quartic constructions will be displayed here, "
        "together with their symbolic Gram-matrix verification."
    )

    if quartic_case == "Pure quartic extension":
        from quartic import (
            pure_quartic_basis,
            pure_quartic_gram,
            TARGET_GRAM,
        )
        import sympy as sp

        st.markdown("### Pure quartic Witt basis")

       

        alpha = alpha_options[subfield]
        r = r_value

        basis = pure_quartic_basis(alpha, r)

        if minus_one_square == "Yes":
            # i^2 = -1 in F
            # Change lambda from -alpha*r^2 to alpha.
            corrected_basis = []

            for vector in basis:
                v = vector.copy()
                v[4] = sp.I * r * v[4]
                corrected_basis.append(v)

            basis = corrected_basis
            lambda_value = alpha

        else:
            lambda_value = -alpha * r**2

        # Gram matrix for the selected quadratic form
        Q = sp.Matrix([
            [1, 0, 0, 0, 0],
            [0, alpha, 0, 0, 0],
            [0, 0, 0, -alpha*r, 0],
            [0, 0, -alpha*r, 0, 0],
            [0, 0, 0, 0, lambda_value],
        ])

        P = sp.Matrix.hstack(*basis)
        gram = sp.simplify(P.T * Q * P)
        st.latex(
            r"\alpha=" + sp.latex(alpha)
            + r",\quad d=" + sp.latex(r) + r"\sqrt{" + sp.latex(alpha) + "}"
            + r",\quad \mu_4=1"
            + r",\quad \lambda=" + sp.latex(lambda_value)
        )
        st.markdown("**Definition of the quartic extension**")

        st.latex(
            r"E_1=F(\sqrt{\alpha}),\qquad "
            r"K=E_1(\delta),\qquad \delta^2=d"
        )

        st.latex(
            r"\delta^2="
            + sp.latex(r * sp.sqrt(alpha))
        )
        if minus_one_square == "Yes":
            st.latex(r"i\in F,\qquad i^2=-1")
        st.markdown("**Constructed Witt basis in the original basis**")

        delta = sp.Symbol("delta")
        sqrt_alpha = sp.sqrt(alpha)

        original_basis = [
            sp.S.One,
            sqrt_alpha,
            delta,
            sqrt_alpha * delta,
        ]

        for name, vector in zip(
            ["e_1", "e_2", "f_1", "f_2", "g_1"],
            basis,
        ):
            x = sum(
                (vector[i] * original_basis[i] for i in range(4)),
                sp.S.Zero,
            )
            z = vector[4]

            st.latex(
                name
                + r"=\left("
                + sp.latex(x)
                + r",\,"
                + sp.latex(z)
                + r"\right)"
            )


        st.markdown("**Computed Gram matrix**")
        st.latex(sp.latex(gram))

        if gram == TARGET_GRAM:
            st.success("Symbolic verification successful.")
        else:
            st.error("Gram matrix verification failed.")
    elif quartic_case == "Exceptional unramified: d = eta":

        from quartic import (
            exceptional_quartic_basis,
            exceptional_quartic_gram,
            TARGET_GRAM,
        )

        import sympy as sp

        st.markdown("### Exceptional unramified quartic Witt basis")

        s, t = sp.symbols("s t", nonzero=True)
        i = sp.Symbol("i")
        delta = sp.Symbol("delta")

        eta = s**2 + i * t**2

        st.markdown("**Definition of the quartic extension**")

        st.latex(
            r"E_1=F(i),\qquad i^2=-1,\qquad "
            r"K=E_1(\delta)"
        )

        st.latex(
            r"\eta=s^2+it^2,\qquad "
            r"s^4+t^4=-1,\qquad "
            r"\delta^2=\eta"
        )

        st.latex(
            r"\mu_4=1,\qquad \lambda=-1"
        )

        # Construct the Witt basis
        basis = exceptional_quartic_basis(s, t)

        # Original F-basis of K
        original_basis = [
            sp.Integer(1),
            i,
            delta,
            i * delta,
        ]

        st.markdown(
            "**Constructed Witt basis in the original basis**"
        )

        for name, vector in zip(
            ["e_1", "e_2", "f_1", "f_2", "g_1"],
            basis,
        ):
            x = sp.expand(
                sum(
                    vector[j] * original_basis[j]
                    for j in range(4)
                )
            )

            z = vector[4]

            st.latex(
                name
                + r"=\left("
                + sp.latex(x)
                + r",\,"
                + sp.latex(z)
                + r"\right)"
            )

        # Compute the Gram matrix
        gram = exceptional_quartic_gram(s, t)

        # Apply the defining relation s^4 + t^4 = -1
        gram_reduced = gram.applyfunc(
            lambda entry: sp.simplify(
                sp.expand(entry).subs(
                    t**4, -1 - s**4
                )
            )
        )

        st.markdown("**Computed Gram matrix**")
        st.latex(sp.latex(gram_reduced))

        if gram_reduced == TARGET_GRAM:
            st.success("Symbolic verification successful.")
        else:
            st.error("Gram matrix verification failed.")
    elif quartic_case == "Exceptional twisted: d = eta*pi":

        from quartic import (
            exceptional_twisted_quartic_basis,
            exceptional_twisted_quartic_gram,
            TARGET_GRAM,
        )

        import sympy as sp

        st.markdown("### Exceptional twisted quartic Witt basis")

        s, t = sp.symbols("s t", nonzero=True)
        varpi = sp.Symbol("varpi", nonzero=True)

        i = sp.Symbol("i")
        delta = sp.Symbol("delta")

        st.markdown("**Definition of the quartic extension**")

        st.latex(
            r"E_1=F(i),\qquad i^2=-1,\qquad "
            r"K=E_1(\delta)"
        )

        st.latex(
            r"\eta=s^2+it^2,\qquad s^4+t^4=-1"
        )

        st.latex(
            r"\delta^2=\eta\varpi,\qquad "
            r"\mu_4=\eta,\qquad \lambda=-1"
        )

        # Construct the Witt basis
        basis = exceptional_twisted_quartic_basis(
            s, t, varpi
        )

        # Original F-basis of K
        original_basis = [
            sp.Integer(1),
            i,
            delta,
            i * delta,
        ]

        st.markdown(
            "**Constructed Witt basis in the original basis**"
        )

                # Display simplified Witt basis
        st.latex(r"e_1=\left(\frac{(1+i)\delta}{\eta},\,0\right)")

        st.latex(r"e_2=\left(-\frac{1}{s},\,1\right)")

        st.latex(
            r"f_1=\left(-\frac{(1-i)\delta}{2\varpi\eta},\,0\right)"
        )

        st.latex(r"f_2=\left(-\frac{1}{2s},\,-\frac12\right)")

        st.latex(
            r"g_1=\left(\frac{t^2}{s}+is,\,0\right)"
        )

        # Compute the Gram matrix
        gram = exceptional_twisted_quartic_gram(
            s, t, varpi
        )

        # Apply the relation s^4 + t^4 = -1
        gram_reduced = gram.applyfunc(
            lambda entry: sp.simplify(
                sp.cancel(entry).subs(
                    t**4, -1 - s**4
                )
            )
        )

        st.markdown("**Computed Gram matrix**")
        st.latex(sp.latex(gram_reduced))

        if gram_reduced == TARGET_GRAM:
            st.success("Symbolic verification successful.")
        else:
            st.error("Gram matrix verification failed.")