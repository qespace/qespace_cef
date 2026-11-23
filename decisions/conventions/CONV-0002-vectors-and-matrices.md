---
id: CONV-0002
title: Vectors and matrices of vectors
status: accepted
owner: qespace_cef
supersedes:
superseded_by:
---

## Statement

In formulas a vector is a column, `v`, and `vᵀ`, `v†` are rows. A matrix
element `M_ij` has i the row and j the column; `<i|M|f>` is `M_if`.

A matrix built from a set of vectors (a basis, eigenvectors, lattice
vectors) holds them as columns: its element (i, j) is component i of
vector j. For every such matrix:

- eigenvectors: `U_mn = <m|n>`, with m the basis of CONV-0005 and n the
  eigenstates, so `H U = U D` and `D = U† H U`; an operator in the
  eigenbasis is `O′ = U† O U`;
- a vector with components c in a basis A is `x = A c`;
- two bases with `a_i·b_j ∝ δ_ij` satisfy `Bᵀ A ∝ 1`.

The kets of CONV-0005 are orthonormal and U is unitary, so the dual of a
basis of kets is its conjugate transpose, `U⁻¹ = U†`. The frame of xyz
components follows CONV-0004.

## Rationale

Physics texts write vectors as columns and the eigenvalue problem as
`H U = U D`. Stacking the vectors as rows instead turns every formula into
its transpose (`D = U H U†`).

This differs from the crystallographic convention of the International
Tables, which is mixed: the direct basis is written as a row `(a, b, c)`
with the coordinates as a column, and the reciprocal basis as a column with
the Miller indices as a row. Its direct basis vectors are thus stacked as
columns and its reciprocal ones as rows, so the duality reads `B A ∝ 1`.
Here both bases are stacked as columns and the duality reads `Bᵀ A ∝ 1`,
like every other change of basis; reciprocal coordinates are a column too.
Data written in the crystallographic convention need one transpose, of the
reciprocal basis, where they are read.

A matrix such as U or a lattice matrix is both a set of vectors and a map,
a change of basis. This CONV takes the reading of the map, so the formulas
read as in the books. How arrays store matrices and stacks of independent
vectors is ADR-0007 (ADR-0004).

The choice concerns every change of basis here (operators in the
eigenbasis, CONV-0005, CONV-0011) and the lattices of the suite's
crystallography package, so it should become a suite convention; the 2π
of a reciprocal lattice is crystallography's.

## Consequences

- Every formula with a change of basis is written by this CONV.
- Input that stacks vectors as rows needs one transpose, where it is read.
- Tests: an operator taken to the eigenbasis and back is unchanged; `H` in
  its own eigenbasis is diagonal with the energies; for a basis and its
  dual, `Bᵀ A` is proportional to the identity.
