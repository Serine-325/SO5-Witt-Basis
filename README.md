# Explicit Witt Bases and Bruhat–Tits Buildings for Elliptic Tori in Split SO(5)

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

We consider elliptic maximal tori arising from different étale algebra structures:

1. **Sums of two quadratic extensions:** $$ F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F$$.
2. **Biquadratic extensions:** $E\oplus F$, where $ E/F$ is biquadratic.
3. **Quartic field extensions:** \(E\oplus F\), where \(E/F\) is a quartic field extension.

For each case, the project aims to construct an explicit Witt basis and use it to study the corresponding point in the Bruhat–Tits building of the split group.

## Current Implementation

### Sum of Two Quadratic Extensions

**File:** `Witt_Basis_Sum_Quadratics.py`

For the quadratic form

\[
q(x,y,z)=\mu_1N_\alpha(x)+\mu_2N_\beta(y)+\lambda z^2,
\]

the program constructs a Witt basis

\[
\mathcal B^*=\{e_1,e_2,f_1,f_2,g_1\}
\]

satisfying the required isotropy, normalization, and orthogonality conditions.

The implementation covers equal quadratic extensions and distinct quadratic extensions, including the relevant Hilbert-symbol cases.

The current version has been tested on 36 admissible parameter configurations.

## Research Development Plan

### Stage 1 — Explicit Witt Basis Constructions
- Sum of two quadratic extensions — implemented and tested.
- Biquadratic extensions — planned.
- Quartic field extensions — planned.

### Stage 2 — Bruhat–Tits Building Computations
- Express the relevant lattice data in the constructed Witt bases.
- Determine the corresponding points in the Bruhat–Tits building.
- Implement the necessary valuation and coordinate calculations.
- Compare the computed results with the theoretical classification.

### Stage 3 — Verification and Documentation
- Verify the mathematical identities and computational outputs.
- Document the constructions and their proofs.
- Develop examples and tests for the implemented algorithms.

## Author

**Serine Bairakji**

PhD in Mathematics, University of Ottawa

Research interests: representation theory of p-adic groups, algebraic groups, elliptic tori, and Bruhat–Tits buildings.
