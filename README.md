# Explicit Witt Bases and Bruhat–Tits Buildings for Elliptic Tori in Split SO(5)

## Run the Interactive Witt Basis Constructor

[Launch the Witt Basis Constructor](https://so5-witt-basis.onrender.com/)

This interactive Python application constructs and symbolically verifies explicit Witt bases associated with elliptic maximal tori in the split special orthogonal group SO(5).

The application currently supports three families:

1. **Sums of two quadratic extensions:** Construct Witt bases from the parameters of two quadratic field extensions, subject to the relevant norm and admissibility conditions.

2. **Biquadratic extensions:** Construct Witt bases using a biquadratic field extension and a selected quadratic intermediate field.

3. **Quartic field extensions with a unique quadratic subfield:** Construct Witt bases for the pure quartic and two exceptional unramified cases, using explicit formulas and their defining algebraic relations.

For each implemented family, users can explore the mathematical parameters, inspect the constructed Witt basis, and verify the corresponding Gram matrix symbolically.

The application is developed using **Python, SymPy, and Streamlit**, with automated mathematical tests implemented using **pytest and GitHub Actions**.

No Python installation or GitHub account is required to use the online application.

**Current scope:** Explicit Witt basis constructions and symbolic Gram-matrix verification for all three families. Bruhat–Tits building computations remain a future research objective.

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

1. **Sums of two quadratic extensions:** $F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F$.
2. **Biquadratic extensions:** $K\oplus F$, where $K/F$ is biquadratic.
3. **Quartic field extensions:** $E\oplus F$, where $E/F$ is a quartic field extension.

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

where $N_\alpha$ and $N_\beta$ denote the corresponding quadratic field norms.

The algorithm constructs an explicit Witt basis

$$
\mathcal B=(e_1,f_1,e_2,f_2,g_1)
$$

under the relevant admissibility and norm conditions.

The implementation distinguishes two cases:

- **Equal quadratic extensions:** $\alpha=\beta$, with explicit constructions using solutions of the corresponding norm equation.
- **Distinct quadratic extensions:** $\alpha\ne\beta$, using norm witnesses and the appropriate Hilbert-symbol conditions.

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
### 3. Quartic Witt Basis Constructions

The module `quartic.py` implements symbolic Witt basis constructions for three families of quartic elliptic tori in the split group SO(5).

#### Implemented cases

1. **Pure quartic case:** `d = r*sqrt(alpha)`, with `mu4 = 1`.
   When `-1` is a square in `F`, the application uses `lambda = alpha`.
   Otherwise, it uses `lambda = -alpha*r**2`.
2. **Exceptional unramified case:** `d = eta`, with `mu4 = 1` and `lambda = -1`.

3. **Exceptional twisted case:** `d = eta*pi`, with `mu4 = eta` and `lambda = -1`.

For the exceptional cases, `eta = s**2 + i*t**2`, where `i**2 = -1` and `s**4 + t**4 = -1`.

#### Symbolic verification

Each construction produces five vectors in the order `(e1, e2, f1, f2, g1)`.

SymPy computes the associated Gram matrix using the polar form

`B(u,v) = (q(u+v) - q(u) - q(v))/2`.

The target Gram matrix is

```text
[0 0 1 0 0]
[0 0 0 1 0]
[1 0 0 0 0]
[0 1 0 0 0]
[0 0 0 0 1]
```

For the exceptional cases, the final entry simplifies to 1 under the assumption `s**4 + t**4 = -1`.

The symbolic checks verify the constructed bases under their stated assumptions. They do not independently solve norm equations or establish local-field admissibility.

#### Running the tests

```bash
python -m pytest -q
```

The quartic tests are located in `tests/test_quartic.py`.

#### Worked Example: Pure Quartic Extension

Consider the pure quartic case with parameters

\[
\alpha=3,\qquad r=1,\qquad \mu_4=1,\qquad \lambda=-3.
\]

These parameters may be interpreted over the local field \(\mathbb Q_5\).

The associated quadratic form is

\[
q(x_0,x_1,x_2,x_3,z)
=x_0^2+3x_1^2-6x_2x_3-3z^2.
\]

The Python implementation constructs the following Witt basis:

\[
\begin{aligned}
e_1&=(0,0,1,0,0),\\
e_2&=(0,1,0,0,1),\\
f_1&=(0,0,0,-1/3,0),\\
f_2&=(0,1/6,0,0,-1/6),\\
g_1&=(1,0,0,0,0).
\end{aligned}
\]

With respect to the normalized polar form

\[
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2},
\]

the Gram matrix in the ordered basis \((e_1,e_2,f_1,f_2,g_1)\) is

\[
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}.
\]

Thus, the construction produces two hyperbolic pairs and an anisotropic vector of norm one, as required for the split five-dimensional quadratic space.

**Reproduce the calculation:**

```python
from quartic import pure_quartic_basis, pure_quartic_gram

basis = pure_quartic_basis(3, 1)

for vector in basis:
    print(list(vector))

print(pure_quartic_gram(3, 1))
```

All calculations use exact SymPy expressions rather than floating-point approximations.

### 4. Verification of the Witt Basis

The constructed vectors are required to satisfy

$$
q(e_i)=q(f_i)=0,\qquad B(e_i,f_j)=\delta_{ij},
$$

together with vanishing cross-pairings and the normalization of \(g_1\).

Here the symmetric bilinear pairing is defined by

$$
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2}.
$$

In the ordered basis

$$
(e_1,e_2,f_1,f_2,g_1),
$$

the expected Gram matrix is

```math
\mathrm{Gram}(\mathcal B)=
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}
```

The biquadratic implementation uses SymPy to verify these identities symbolically, avoiding floating-point approximations.

Symbolic verification confirms the algebraic identities of the constructed vectors. Existence of the required norm witnesses over the local field is justified by the mathematical norm conditions.

### 5. Software Organization

| File | Description |
|---|---|
| `Witt_Basis_Sum_Quadratics.py` | Constructs Witt bases for sums of two quadratic extensions. |
| `biquadratic.py` | Constructs Witt bases for biquadratic extensions. |
| `quartic.py` | Constructs symbolic Witt bases for three quartic cases. |
| `app.py` | Streamlit interface for all three families. |
| `web_runner.py` | Connects the quadratic construction to the web application. |
| `verification_panel.py` | Provides supporting verification functionality. |
| `tests/` | Contains automated tests for the mathematical constructions. |
| `tests/test_biquadratic.py` | Tests the biquadratic Witt basis constructions. |
| `tests/test_quartic.py` | Tests the quartic Witt basis constructions. |
| `.github/workflows/tests.yml` | Runs the Python test suite automatically using GitHub Actions. |
| `requirements.txt` | Lists the Python dependencies. |
| `README.md` | Documents the mathematical framework, implementation, and research objectives. |

The interactive application is available through the link at the beginning of this README.

## Research Development Plan

### Stage 1 — Explicit Witt Basis Constructions

- [x] Sums of two quadratic extensions — implemented and tested.
- [x] Biquadratic extensions — implemented with exact symbolic verification and automated tests.
- [x] Quartic field extensions with a unique quadratic subfield — implemented and tested.

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
- [x] Implement symbolic Witt basis constructions for the three quartic cases.
- [x] Verify their Gram matrices using SymPy.
- [x] Integrate the quartic constructions into the interactive application.
- [ ] Add worked examples and document the mathematical assumptions.



## Author
**Serine Bairakji**

PhD in Mathematics, University of Ottawa

Research interests: representation theory of p-adic groups, algebraic groups, elliptic tori, and Bruhat–Tits buildings.
