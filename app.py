"""Streamlit interface for both established SO(5) Witt basis families."""
import streamlit as st
import sympy as sp

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

tab_sum, tab_biquad = st.tabs([
    "Sums of two quadratic extensions",
    "Biquadratic extension",
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

    options = {
        "canonical": "Canonical split representative: μ = √α, t = 1",
        "witness": "Supply a norm witness t; compute μ = √α / N(t)",
    }
    if not square and a_class == "epsilon":
        options["special"] = "Unramified case, −1 nonsquare: μ = 1 (direct basis)"
    mode = st.radio(
        "Construction", list(options), format_func=lambda key: options[key],
        key="biq_mode",
    )

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
                square, a_class, b_class, mode, t_coordinates, rho, order_pairs,
            )
        except (ValueError, ArithmeticError) as error:
            st.session_state.pop("biq_result", None)
            st.error(str(error))

    current_settings = (square, a_class, b_class, mode, t_coordinates, rho, order_pairs)
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
    "**Research scope:** explicit Witt bases for sums of two quadratic extensions "
    "and for biquadratic extensions. Quartic extensions and Bruhat–Tits "
    "building calculations are not implemented here."
)
st.markdown(
    "Source repository: [SO5-Witt-Basis](https://github.com/Serine-325/SO5-Witt-Basis)."
)
