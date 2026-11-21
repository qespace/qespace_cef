---
id: CONV-0005
title: Basis: the ground J multiplet and the angular-momentum matrices
status: accepted
owner: qespace_cef
supersedes:
---

## Statement

The basis is the ground J multiplet of the free ion (LS coupling, Hund's
rules), 2J + 1 states, with the Landé `gJ`; there is no intermediate coupling.

The matrices of J_x, J_y and J_z are the root of every matrix in the package.
They are written in the eigenstates of J_z, ordered `|J>, |J−1>, ..., |−J>`,
so `J_z = diag(J, J−1, ..., −J)`, with the Condon–Shortley phases:

    <m+1|J₊|m> = √(J(J+1) − m(m+1)),  real and non-negative
    J₋ = J₊†,  J_x = (J₊ + J₋)/2,  J_y = (J₊ − J₋)/2i

Every other matrix is derived from these and states no convention of its
own: the Stevens operators (by the formulas of CONV-0001), the Hamiltonian,
and the eigenvectors and kets, whose components are coefficients on these
J_z eigenstates. The one phase the J matrices leave free, the overall phase
of time reversal, is fixed in CONV-0011.

## Rationale

For most rare-earth ions the spin–orbit splitting is much larger than the CEF,
so J is a good quantum number and the ground multiplet is enough. It fails
when an excited multiplet is close, as for Sm³⁺ and Eu³⁺, where the CEF mixes
multiplets (J-mixing).

The order and the phases are not an implementation detail. Writing the same
matrices in the order `|−J>, ..., |J>` maps J_x → J_x, J_y → −J_y, J_z → −J_z,
a rotation by π about x: 12 of the 27 even-rank Stevens matrices change sign
(cosine type with m odd, sine type with m even; checked on J = 4). The same
parameter values with arrays in the other order describe the ion rotated by
π. Stating the convention once, on J_x, J_y, J_z, lets it reach every derived
matrix by construction, and one test of the J matrices covers it.

## Consequences

- The package has one construction of J_x, J_y, J_z (and J₊, J₋), and every
  other matrix is built from it (AUDIT A5).
- Data exchanged with another code (matrices, eigenvectors) needs that code's
  order and phases, converted at the boundary.
- Intermediate coupling would need the full fⁿ basis and Wybourne parameters
  on spherical tensors. It may come later; a new CONV then supersedes this one.
- Tests, each named `test_conv_0005_*`:
  - `J_z` is diagonal with `J, J−1, ..., −J`, so `<J|J_z|J> = J` for the first
    state;
  - `J₊` has the Condon–Shortley elements on its first superdiagonal and
    nothing else;
  - `[J_x, J_y] = i J_z` and `J_x² + J_y² + J_z² = J(J+1)`;
  - the J matrices `StevensBase` uses are these.
