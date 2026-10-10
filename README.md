# Explicit Witt Basis Constructions for Elliptic Tori in SO(5)

## Overview

This repository contains Python implementations of explicit Witt basis constructions associated with elliptic maximal tori in the split special orthogonal group SO(5).

The project builds on my doctoral research in the representation theory of p-adic groups at the University of Ottawa.

## Mathematical Background

Let F be a non-archimedean local field of odd residue characteristic.

For the sum-of-two-quadratic-extensions case, consider the  algebra

\[
\mathcal A=F(\sqrt{\alpha})\oplus F(\sqrt{\beta})\oplus F.
\]

The associated quadratic form is

\[
q(x,y,z)=\mu_1N_\alpha(x)+\mu_2N_\beta(y)+\lambda z^2.
\]

The program constructs a Witt basis

\[
\mathcal B^*=\{e_1,e_2,f_1,f_2,g_1\}
\]

for admissible split quadratic forms, satisfying

\[
q(e_i)=q(f_i)=0,\quad B(e_i,f_j)=\delta_{ij},\quad q(g_1)=1,
\]

with the remaining pairings equal to zero.

## Current Implementation

**Witt_Basis_Sum_Quadratics.py**

Handles three constructions:

1. Equal quadratic extensions.
2. Distinct quadratic extensions with Hilbert symbol +1.
3. Distinct quadratic extensions with Hilbert symbol -1.

The implementation includes explicit norm-witness constructions and exceptional cases involving the unramified quadratic extension.

## Verification

The current implementation has been tested on 36 admissible parameter configurations, using algebraic checks of the Witt basis identities.

These tests provide computational verification of the implemented constructions; they do not replace the general mathematical proofs.

## Planned Extensions

- Explicit Witt basis constructions for biquadratic extensions.
- Explicit Witt basis constructions for quartic field extensions.
- Applications to the Bruhat–Tits building and elliptic tori of SO(5).

## Author

Serine Bairakji

PhD in Mathematics, University of Ottawa

Research interests: representation theory of p-adic groups, algebraic groups, and explicit constructions in SO(5).
