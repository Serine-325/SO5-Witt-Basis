"""Public-facing Streamlit interface for the Stage 1 SO(5) constructor."""
import streamlit as st

from web_runner import CLASSES, MU_OPTIONS, run_constructor, valid_parameters

st.set_page_config(
    page_title="SO(5) | Explicit Witt Basis Constructor",
    page_icon="📐",
    layout="wide",
)

st.title("Explicit Witt Bases for Elliptic Tori in Split SO(5)")
st.caption("Stage 1 · Sums of two quadratic extensions · Research software")
st.markdown(
    "This interface constructs a Witt basis using the verified Python program "
    "for elliptic maximal tori arising from sums of two quadratic extensions. "
    "No Python installation or GitHub account is required to use this page."
)
st.latex(r"V=F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F")
st.latex(r"q(x,y,z)=\mu_1N_{\alpha}(x)+\mu_2N_{\beta}(y)+\lambda z^2")

st.divider()
st.subheader("Choose the torus parameters")
minus_one_square = st.radio(
    "Is −1 a square in the base field F?",
    options=(True, False),
    format_func=lambda value: "Yes" if value else "No",
    horizontal=True,
    key="minus_one_square",
)
left, right = st.columns(2)
with left:
    alpha = st.selectbox(r"First quadratic class $\alpha$", CLASSES, key="alpha")
    mu1 = st.selectbox(r"First coefficient $\mu_1$", MU_OPTIONS[alpha], key="mu1")
with right:
    beta = st.selectbox(r"Second quadratic class $\beta$", CLASSES, key="beta")
    mu2 = st.selectbox(r"Second coefficient $\mu_2$", MU_OPTIONS[beta], key="mu2")

params = (minus_one_square, alpha, beta, mu1, mu2)
try:
    admissible = valid_parameters(*params)
except ValueError as error:
    admissible = False
    st.error(str(error))

if admissible:
    st.success("These parameters satisfy the split SO(5) Hilbert-symbol criterion.")
else:
    st.warning(
        "The selected parameters do not satisfy "
        r"$(\mu_1,\alpha)(\mu_2,\beta)=(\alpha,\beta)$. "
        "Change one or both coefficients to continue."
    )

if st.button("Construct Witt Basis", type="primary", disabled=not admissible):
    try:
        st.session_state["last_result"] = run_constructor(*params)
        st.session_state["last_params"] = params
    except Exception as error:
        st.session_state.pop("last_result", None)
        st.error(f"The construction could not be completed: {error}")

if st.session_state.get("last_params") == params and "last_result" in st.session_state:
    result = st.session_state["last_result"]
    st.divider()
    st.subheader("Constructed Witt basis")
    st.markdown(f"**Case {result.case}** · Square class of λ: `{result.lam}`")
    st.caption(
        "Coordinates are in the original order "
        "F(√α) ⊕ F(√β) ⊕ F, as printed by the research program."
    )
    for vector in ("e1", "f1", "e2", "f2", "g1"):
        st.code(f"{vector} = {result.basis[vector]}", language="text")

    if result.witnesses:
        st.markdown("**Norm witnesses**")
        for key, value in result.witnesses.items():
            st.code(f"{key} = {value}", language="text")

    st.info(
        "Some constructions use symbolic elements such as u, s, r, or s0 and t0. "
        "They are subject to the exact norm equations stated in the output; "
        "the interface does not numerically solve norm equations over an unspecified local field."
    )
    with st.expander("View full mathematical output and norm conditions"):
        st.code(result.raw_output, language="text")
    st.download_button(
        "Download this result (.txt)",
        data=result.raw_output,
        file_name="SO5_Witt_Basis_Result.txt",
        mime="text/plain",
    )

st.divider()
st.markdown(
    "**Research scope.** This interface currently covers sums of two quadratic "
    "extensions only. Biquadratic and quartic extensions, and the associated "
    "Bruhat–Tits building computations, are planned as separate stages."
)
st.markdown(
    "Source repository: [SO5-Witt-Basis](https://github.com/Serine-325/SO5-Witt-Basis) "
    "(repository access depends on its visibility settings)."
)
