---
id: CONV-0002
title: Vectors and matrices of vectors
status: draft
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
`H U = U D`; NumPy's `eigh` returns the eigenvectors as columns. Stacking
the vectors as rows instead, as pymatgen and spglib's Python interface do
for lattices, turns every formula into its transpose (`D = U H U†`). The
International Tables mix the two, the basis as a row and the coordinates
as a column, and mixing within one formula puts transposes in the wrong
place: `B A = 1` with A in columns and B in rows.

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
- Converting from codes that stack vectors as rows needs one transpose,
  where the data are read.
- Tests: an operator taken to the eigenbasis and back is unchanged; `H` in
  its own eigenbasis is diagonal with the energies; for a basis and its
  dual, `Bᵀ A` is proportional to the identity.
