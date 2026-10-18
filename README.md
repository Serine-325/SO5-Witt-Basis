# Explicit Witt Bases for Elliptic Maximal Tori in Split SO(5)

This Python research project constructs and verifies **explicit Witt bases** for quadratic spaces associated with elliptic maximal tori in the split group $\mathrm{SO}(5)$. It is based on doctoral research in the representation theory of $p$-adic groups at the University of Ottawa.

**[Launch the interactive Witt Basis Constructor](https://so5-witt-basis.onrender.com/)** — explore parameters, inspect the basis vectors, and check their Gram matrices in a web interface. No local installation or GitHub account is needed to use the app.

**Implemented:** constructions for sums of two quadratic extensions, biquadratic extensions, and three quartic cases. **Verification:** symbolic Gram-matrix checks in SymPy, arithmetic and admissibility tests for supported local-field settings, and automated tests using pytest and GitHub Actions. **Reported test status (October 18, 2026): 188 passing tests.**

**Research boundary:** the Witt basis stage is implemented; the algorithms to locate corresponding points in the Bruhat–Tits building are planned work. The arithmetic routines currently cover specific local fields and examples, not all non-archimedean local fields.

## Quick start

```bash
git clone https://github.com/Serine-325/SO5-Witt-Basis.git
cd SO5-Witt-Basis
python -m pip install -r requirements.txt
python -m pytest -q
python -m streamlit run app.py
```

The test command runs the mathematical checks; the last command starts the Streamlit interface. Use a Python environment with the required packages installed. The repository's `requirements.txt` lists the application dependencies.

## Mathematical setting and research goal

Let $F$ be a non-archimedean local field of odd residue characteristic. We study elliptic maximal tori in the split special orthogonal group

$$
G=\mathrm{SO}(2\mathbb H\oplus\langle1\rangle)\simeq\mathrm{SO}(5).
$$

The long-term goal is to compute the points in the Bruhat–Tits building associated with these tori, directly from their defining parameters and quadratic forms. The first stage constructs explicit Witt bases, giving coordinates for the lattice calculations needed in that work.

The implementation treats three types of étale $F$-algebras:

1. **Two quadratic extensions:** $F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F$.
2. **Biquadratic extension:** $K\oplus F$, where $K/F$ is biquadratic.
3. **Quartic field extension:** $E\oplus F$, where $E/F$ is a quartic field extension with a unique quadratic subfield in the implemented cases.

The algorithms build Witt bases from torus parameters subject to the required admissibility and norm conditions. This project extends the author's doctoral research at the University of Ottawa.

## Implemented constructions

### 1. Two quadratic extensions

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

Under the relevant norm and admissibility conditions, the algorithm constructs an explicit Witt basis

$$
\mathcal B=(e_1,f_1,e_2,f_2,g_1).
$$

For the Gram-matrix display, the vectors are ordered as $(e_1,e_2,f_1,f_2,g_1)$.

The implementation distinguishes two cases:

- **Equal quadratic extensions:** $\alpha=\beta$, with explicit constructions using solutions of the corresponding norm equation.

- **Distinct quadratic extensions:** $\alpha\ne\beta$, using norm witnesses and the appropriate Hilbert-symbol conditions.

The original implementation has been tested on 36 admissible parameter configurations.

### 2. Biquadratic extensions

**Python file:** `biquadratic.py`

Let

$$
K=F(\sqrt{\epsilon},\sqrt{\varpi})
$$

be a biquadratic extension of $F$, and let

$$
E_1=F(\sqrt{\alpha})
$$

be one of its three quadratic intermediate fields.

Writing $\sigma$ for the nontrivial automorphism of $K/E_1$, consider the quadratic form

$$
q_{\mu,\lambda}(x,z)=\mathrm{Tr}_{K/F}\left(\frac{\mu}{4}x\sigma(x)\right)+\lambda z^2
$$

The construction uses the relative norm map

$$
N_{K/E_1}:K^\times\longrightarrow E_1^\times
$$

and the relevant norm-class criterion to determine suitable parameters for the split quadratic space.

The algorithm:

1. Allows the user to select a quadratic intermediate field and the associated square classes.

2. Distinguishes the cases where $-1$ is or is not a square in $F$.

3. Constructs explicit isotropic vectors and their normalized dual vectors.

4. Orders the hyperbolic pairs according to the adopted ramification convention.

5. Performs exact symbolic verification of the resulting Witt basis.

The implementation includes automated tests for the biquadratic construction.

### 3. Quartic extensions

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

The quartic tests are located in `tests/test_quartic.py`. The complete suite is run with the same command.

## Admissibility and verification

In supported quartic cases, the automated tests check the two conditions in the admissibility criterion of Theorem 3.6 of the underlying research:

$$
q_{\mu_4,\lambda}\cong 2\mathbb H\perp\langle1\rangle
\quad\Longleftrightarrow\quad
\begin{cases}
\lambda\equiv D\pmod{F^{\times2}},\\
(\mu_4,d)_{E_1}=\epsilon_H(q_1)(\lambda,-1)_F.
\end{cases}
$$

The Hilbert-symbol routines include computations over $\mathbb Q_p$ and $\mathbb Q_3(i)$, with $3$-adic approximations used to check the exceptional quartic examples. These computations support the specified examples; they do not constitute a general-purpose local-field package.

### Witt basis checks

The constructed vectors are required to satisfy

$$
q(e_i)=q(f_i)=0,\qquad B(e_i,f_j)=\delta_{ij},
$$

together with vanishing cross-pairings and the normalization of $g_1$.

Here the symmetric bilinear pairing is defined by

$$
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2}.
$$

In the ordered basis

$$
(e_1,e_2,f_1,f_2,g_1),
$$

the expected Gram matrix is

$$
\mathrm{Gram}(\mathcal B)=
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}
$$

The biquadratic implementation uses SymPy to verify these identities symbolically, avoiding floating-point approximations.

Symbolic verification confirms the algebraic identities of the constructed vectors. Existence of the required norm witnesses over the local field is justified by the mathematical norm conditions.

## Worked examples

### Pure quartic extension

Consider the pure quartic case with parameters

$$
\alpha=3,\qquad r=1,\qquad \mu_4=1,\qquad \lambda=-3.
$$

These parameters may be interpreted over the local field $\mathbb Q_5$. Here $-1$ is a square, so $-3$ and $3$ represent the same square class; the example uses the representative $\lambda=-3$.

The associated quadratic form is

$$
q(x_0,x_1,x_2,x_3,z)
=x_0^2+3x_1^2-6x_2x_3-3z^2.
$$

The Python implementation constructs the following Witt basis:

$$
\begin{aligned}
e_1&=(0,0,1,0,0),\\
e_2&=(0,1,0,0,1),\\
f_1&=(0,0,0,-1/3,0),\\
f_2&=(0,1/6,0,0,-1/6),\\
g_1&=(1,0,0,0,0).
\end{aligned}
$$

With respect to the normalized polar form

$$
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2},
$$

the Gram matrix in the ordered basis $(e_1,e_2,f_1,f_2,g_1)$ is

$$
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}.
$$

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

### Exceptional unramified quartic extension

We consider the local field $F=\mathbb Q_3$, for which $-1$ is not a square.

Let

$$
E_1=F(i),\qquad i^2=-1.
$$

Choose $t=1$ and let $s\in\mathbb Q_3$ satisfy

$$
s^4=-2,\qquad s\equiv1\pmod3.
$$

Such an element exists by Hensel's lemma, applied to $X^4+2$ at $X=1$.

Define

$$
\eta=s^2+i,\qquad K=E_1(\sqrt{\eta}).
$$

Then

$$
N_{E_1/F}(\eta)=s^4+1=-1.
$$

Consequently, $\eta$ is not a square in $E_1$, and $K/F$ is a cyclic quartic extension with unique quadratic subfield $E_1$.

For the exceptional unramified construction, the parameters are

$$
d=\eta,\qquad \mu_4=1,\qquad\lambda=-1.
$$

**Explicit Witt basis**

Using the coordinate order $(x_0,x_1,x_2,x_3,z)$, the Python implementation produces

$$
\begin{aligned}
e_1&=(1,1,0,0,0),\\
e_2&=(0,0,0,1/s,1),\\
f_1&=(1/2,-1/2,0,0,0),\\
f_2&=(0,0,0,1/(2s),-1/2),\\
g_1&=(0,0,s,-1/s,0).
\end{aligned}
$$

**Symbolic verification**

The Gram matrix is computed using the symmetric bilinear form

$$
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2}.
$$

Reducing its entries under the defining relation $s^4+2=0$ gives

$$
\mathrm{Gram}(e_1,e_2,f_1,f_2,g_1)=
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}.
$$

The construction can be reproduced using:

```python
import sympy as sp
from quartic import exceptional_quartic_basis, exceptional_quartic_gram
s = sp.Symbol("s", nonzero=True)
basis = exceptional_quartic_basis(s, 1)
gram = exceptional_quartic_gram(s, 1)
for name, vector in zip(("e1", "e2", "f1", "f2", "g1"), basis):
    print(name, list(vector))
print("Gram matrix:")
print(gram)
```

The entries of the Gram matrix are to be simplified using $s^4=-2$.

This example demonstrates the explicit Witt basis construction for an admissible exceptional unramified quartic torus. The symbolic verification establishes the required Gram identities, while the existence of the local-field parameter $s$ follows independently from Hensel's lemma.

### Exceptional twisted quartic extension

We continue working over the local field $F=\mathbb Q_3$, with

$$
E_1=F(i),\qquad i^2=-1.
$$

Let $s\in\mathbb Q_3$ satisfy $s^4=-2$ and $s\equiv1\pmod3$, and set

$$
t=1,\qquad \eta=s^2+i,\qquad \varpi=3.
$$

The existence of $s$ follows from Hensel's lemma, as established in the preceding example.

For the exceptional twisted construction, choose

$$
d=3\eta,\qquad \mu_4=\eta,\qquad \lambda=-1.
$$

Since

$$
N_{E_1/F}(d)=9(s^4+1)=-9,
$$

the norm has nonsquare class $-1$ in $F$. Thus $K=E_1(\sqrt{3\eta})$ is a cyclic quartic extension with unique quadratic subfield $E_1$.

**Explicit Witt basis**

In the coordinate order $(x_0,x_1,x_2,x_3,z)$, the construction produces the following vectors, simplified using $s^4=-2$:

$$
\begin{aligned}
e_1&=(0,0,-s^2-1,1-s^2,0),\\
e_2&=(-1/s,0,0,0,1),\\
f_1&=(0,0,(s^2-1)/6,-(s^2+1)/6,0),\\
f_2&=(-1/(2s),0,0,0,-1/2),\\
g_1&=(1/s,s,0,0,0).
\end{aligned}
$$

**Symbolic verification**

The symmetric bilinear form is defined by

$$
B(u,v)=\frac{q(u+v)-q(u)-q(v)}{2}.
$$

Reducing the entries under the relation $s^4+2=0$ yields

$$
\operatorname{Gram}(e_1,e_2,f_1,f_2,g_1)=
\begin{pmatrix}
0&0&1&0&0\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&1
\end{pmatrix}.
$$

The example can be reproduced using the existing Python implementation:

```python
import sympy as sp
from quartic import (
    exceptional_twisted_quartic_basis,
    exceptional_twisted_quartic_gram,
    TARGET_GRAM,
)
s = sp.Symbol("s", nonzero=True)
basis = exceptional_twisted_quartic_basis(s, 1, 3)
gram = exceptional_twisted_quartic_gram(s, 1, 3)
for name, vector in zip(("e1", "e2", "f1", "f2", "g1"), basis):
    print(name, list(vector))
difference = gram - TARGET_GRAM
def reduce_numerator(expression):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    assert sp.rem(denominator, s**4 + 2, s) != 0
    return sp.rem(numerator, s**4 + 2, s)
reduced_difference = difference.applyfunc(reduce_numerator)
assert reduced_difference == sp.zeros(5)
print("Gram matrix verification passed.")
```

This verifies the Witt basis identities exactly under the defining algebraic relation. The symbolic calculation does not approximate $3$-adic numbers; the local-field existence of the parameter $s$ is established separately by Hensel's lemma.

## Repository structure

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

The interactive application is available through the link at the beginning of this README. The last reported full test run had **188 passing tests** (October 18, 2026).

## Research development plan

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

- [x] Provide an interactive interface for all three implemented families.

- [ ] Add detailed worked examples for the two-quadratic and biquadratic constructions.

- [ ] Document the mathematical proofs and their relationship to the algorithms.

- [x] Implement symbolic Witt basis constructions for the three quartic cases.

- [x] Verify their Gram matrices using SymPy.

- [x] Integrate the quartic constructions into the interactive application.

- [x] Include worked quartic examples and document their mathematical assumptions.

## Author

**Serine Bairakji**

PhD in Mathematics, University of Ottawa

Research interests: representation theory of p-adic groups, algebraic groups, elliptic tori, and Bruhat–Tits buildings.
