---
id: ADR-0007
title: Broadcasting and the layout of arrays
status: draft
date: 2026-11-23
supersedes: ADR-0003
superseded_by:
---

## Context

ADR-0003 made broadcasting the principle of v2 and said that a vector has
its component on the last axis. That holds for data, such as fields and Q
points, but not for eigenvectors: `eigh` returns U with
`U[..., m, n] = <m|n>` (CONV-0002), and ADR-0005 gives `Jxyz` the shape
`S + (3, 2J+1, 2J+1)`, the vector axis before the matrix axes. CONV-0012's
flat `(2J+1)²` entries are being replaced by matrices over the pairs (i, f)
(AUDIT A16), and ADR-0005 now holds the parameters by `(l, m)` rather than
in one `(K,)` axis. ADR-0003's measurement still holds: on Ho³⁺ with 8
parameters, 100 × 100 points, one batched diagonalisation is 15× faster
than a loop.

## Options

1. Keep ADR-0003 and handle eigenvectors, `Jxyz` and the pair matrices
   case by case.
2. One layout rule for data, matrices and their combinations.

## Decision

Option 2. Broadcasting stays as in ADR-0003, as in a NumPy gufunc: loop
axes lead and broadcast, core axes trail.

- **Loop shape.** Each parameter value of a `CrystalFieldParameters` has
  a shape `S_p`; the field `Bxyz` has `S_f + (3,)`. An `Ion` has the loop
  shape `S = broadcast(S_p, S_f)`, one ion per point. A single point is
  `()`. Temperature and Q in an observable broadcast against S; the caller
  chooses a grid or point-by-point pairing with `None` axes.
- **Data vectors.** A vector that is data (a field, a Q point, a
  position, a state vector) has its component on the last axis: `Bxyz`
  `(..., 3)`, a state `(..., 2J+1)`.
- **Matrices.** A matrix is stored as its matrix, element `[..., i, j]`,
  so code reads like the formula. A matrix of vectors holds them as
  columns (CONV-0002): `eigenvectors[..., m, n] = <m|n>`, the n-th
  eigenvector is `eigenvectors[..., :, n]`.
- **Hilbert-space axes trail.** Wherever a pair of state indices appears,
  its two axes are the last ones, and 3D component axes come right before
  them:
  - `hamiltonian`, `eigenvectors`: `S + (2J+1, 2J+1)`;
  - `Jxyz`: `S + (3, 2J+1, 2J+1)`;
  - `energies`: `S + (2J+1,)`;
  - `excitation_energies`: `S + (2J+1, 2J+1)`, element `[..., i, f]`
    the excitation from state i to state f (CONV-0012's ε_if);
  - scattering weights: `S + (3, 3, 2J+1, 2J+1)`, element
    `[..., α, β, i, f]` (CONV-0012's W_if^{αβ});
  - cross-section weights: `broadcast(S, Q) + (2J+1, 2J+1)`, with Q the
    loop shape of the Q points.
- **Products** between data and matrices, and over 3D components, are
  written with named indices (`einsum`), so no transpose is implied.

## Consequences

Every array of the v2 API has one place for each axis, and the
documentation of each states its shape. Eigenvectors are taken as `eigh`
returns them. A selection by the last index picks a final state: the
excitations from the ground state are `excitation_energies[..., 0, :]`.
A measurement log pairs `Bxyz (N, 3)` with `temperature (N,)` point by
point; a grid needs explicit new axes. A test checks batched results
against a loop. Degenerate levels only through basis-independent masks;
levels sorted per point are not branches where they cross. Memory grows
with S: compute only the matrix elements an observable needs, or work in
chunks.
