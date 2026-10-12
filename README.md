# Explicit Witt Bases and Bruhat–Tits Buildings for Elliptic Tori in Split SO(5)

## Run the Interactive Witt Basis Constructor

[Launch the Witt Basis Constructor](https://so5-witt-basis.onrender.com/)

This interactive Python application constructs and verifies explicit Witt bases associated with elliptic maximal tori in split SO(5).

The application currently supports two families:

1. **Sums of two quadratic extensions:** Construct Witt bases from the parameters of two quadratic field extensions, including the relevant norm conditions.

2. **Biquadratic extensions:** Construct Witt bases using a biquadratic field extension and a chosen quadratic intermediate field. The program performs exact symbolic verification of the resulting Gram matrix.

For both families, users can select mathematical parameters and inspect the resulting Witt basis.

The application is developed using Python, Streamlit, and SymPy.

No Python installation or GitHub account is required to use the online application.

**Current scope:** Explicit Witt basis constructions for sums of two quadratic extensions and biquadratic extensions.
## Project Overview

This project develops computational methods for studying elliptic maximal tori in the split special orthogonal group

$$
G=\mathrm{SO}(2\mathbb H\oplus\langle1\rangle)\simeq\mathrm{SO}(5)
$$

over a non-archimedean local field of odd residue characteristic.

The ultimate objective is to explicitly determine the points in the Bruhat–Tits building associated with elliptic maximal tori, using the parameters defining the tori and their corresponding quadratic forms.

The first stage constructs explicit Witt bases, providing coordinates suitable for studying the associated lattices and their positions in the building.

The project builds on my doctoral research in the representation theory of p-adic groups at the University of Ottawa.

## Mathematical Framework

Let \(F\) be a non-archimedean local field of odd residue characteristic. We study elliptic maximal tori in the split special orthogonal group

$$
G=\mathrm{SO}(2\mathbb H\oplus\langle1\rangle)\simeq\mathrm{SO}(5).
$$

The project considers three families of algebra structures:

1. **Sums of two quadratic extensions:** \(F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F\).
2. **Biquadratic extensions:** \(K\oplus F\), where \(K/F\) is biquadratic.
3. **Quartic field extensions:** \(E\oplus F\), where \(E/F\) is a quartic field extension.

For each family, the objective is to construct an explicit Witt basis directly from the torus parameters and subsequently use these coordinates to study the corresponding points in the Bruhat–Tits building.

## Mathematical Implementation

### 1. Sums of Two Quadratic Extensions

**Python file:** `Witt_Basis_Sum_Quadratics.py`

Consider

$$
V=F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F
$$

equipped with the quadratic form

$$
q(x,y,z)=\mu_1N_\alpha(x)+\mu_2N_\beta(y)+\lambda z^2,
$$

where \(N_\alpha\) and \(N_\beta\) denote the corresponding quadratic field norms.

The algorithm constructs an explicit Witt basis

$$
\mathcal B=(e_1,f_1,e_2,f_2,g_1)
$$

under the relevant admissibility and norm conditions.

The implementation distinguishes two cases:

- **Equal quadratic extensions:** \(\alpha=\beta\), with explicit constructions using solutions of the corresponding norm equation.
- **Distinct quadratic extensions:** \(\alpha\ne\beta\), using norm witnesses and the appropriate Hilbert-symbol conditions.

The original implementation has been tested on 36 admissible parameter configurations.

### 2. Biquadratic Extensions

**Python file:** `biquadratic.py`

Let

$$
K=F(\sqrt{\epsilon},\sqrt{\varpi})
$$

be a biquadratic extension of \(F\), and let

$$
E_1=F(\sqrt{\alpha})
$$

be one of its three quadratic intermediate fields.

Writing $\sigma$ for the nontrivial automorphism of $K/E_1$, consider the quadratic form


```math
q_{\mu,\lambda}(x,z)=\mathrm{Tr}_{K/F}\left(\frac{\mu}{4}x\sigma(x)\right)+\lambda z^2
```

The construction uses the relative norm map

$$
N_{K/E_1}:K^\times\longrightarrow E_1^\times
$$

and the relevant norm-class criterion to determine suitable parameters for the split quadratic space.

The algorithm:

1. Allows the user to select a quadratic intermediate field and the associated square classes.
2. Distinguishes the cases where \(-1\) is or is not a square in \(F\).
3. Constructs explicit isotropic vectors and their normalized dual vectors.
4. Orders the hyperbolic pairs according to the adopted ramification convention.
5. Performs exact symbolic verification of the resulting Witt basis.

The implementation includes automated tests for the biquadratic construction.

### 3. Verification of the Witt Basis

The constructed vectors are required to satisfy

$$
q(e_i)=q(f_i)=0,\qquad(e_i,f_j)=\delta_{ij},
$$

together with vanishing cross-pairings and the normalization of \(g_1\).

Here the symmetric bilinear pairing is defined by

$$
(u,v)=q(u+v)-q(u)-q(v).
$$

In the ordered basis

$$
(e_1,f_1,e_2,f_2,g_1),
$$

the expected Gram matrix is

$$
\begin{pmatrix}
0&1&0&0&0\\
1&0&0&0&0\\
0&0&0&1&0\\
0&0&1&0&0\\
0&0&0&0&1
\end{pmatrix}.
$$

The biquadratic implementation uses SymPy to verify these identities symbolically, avoiding floating-point approximations.

Symbolic verification confirms the algebraic identities of the constructed vectors. Existence of the required norm witnesses over the local field is justified by the mathematical norm conditions.

### 4. Software Organization

| File | Description |
|---|---|
| `Witt_Basis_Sum_Quadratics.py` | Construction for sums of two quadratic extensions |
| `biquadratic.py` | Construction for biquadratic extensions |
| `web_runner.py` | Connects the original construction to the web application |
| `app.py` | Streamlit interface for both constructions |
| `tests/test_biquadratic.py` | Automated tests for the biquadratic construction |
| `requirements.txt` | Python dependencies |

The interactive application is available through the link at the beginning of this README.

## Research Development Plan

### Stage 1 — Explicit Witt Basis Constructions

- [x] Sums of two quadratic extensions — implemented and tested.
- [x] Biquadratic extensions — implemented with exact symbolic verification and automated tests.
- [ ] Quartic field extensions — planned.

### Stage 2 — Bruhat–Tits Building Computations

- [ ] Express the relevant lattice data in the constructed Witt bases.
- [ ] Determine the corresponding points in the Bruhat–Tits building.
- [ ] Implement valuation and coordinate calculations.
- [ ] Compare the computed results with the theoretical classification.

### Stage 3 — Verification and Documentation

- [x] Implement automated checks for the biquadratic Witt basis construction.
- [x] Provide an interactive interface for the first two families.
- [ ] Add detailed worked examples for both constructions.
- [ ] Document the mathematical proofs and their relationship to the algorithms.
- [ ] Extend verification to future quartic constructions.

## Author

**Serine Bairakji**

PhD in Mathematics, University of Ottawa

Research interests: representation theory of p-adic groups, algebraic groups, elliptic tori, and Bruhat–Tits buildings.
