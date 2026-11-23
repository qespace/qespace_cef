---
id: ADR-0003
title: Broadcasting and the shape rule
status: superseded
date: 2026-11-20
supersedes:
superseded_by: ADR-0007
---

> Superseded by ADR-0007 on 2026-11-23: broadcasting is unchanged; it adds
> the layout of data vectors, matrices and Hilbert-space axes.

## Context

The Hamiltonian, `Σ_lm B_lm O_lm` plus a Zeeman term proportional to `J·B`, is
linear in the parameters `B_lm` and the field B, and `np.linalg.eigh`
diagonalises stacks, so a whole grid is one `einsum` and one call. Measured on Ho³⁺ with 8 parameters, 100 × 100 points: 15× faster than a
loop over `CEFion`. T and Q don't enter H; they act after the diagonalisation.

## Options

1. (a) Every input adds its own axes: T `(nT,)` gives `S + (nT, 3)`.
2. (b) numpy broadcasting: the inputs broadcast against each other.

## Decision

Broadcasting is the principle of v2, with rule (b), as in a numpy gufunc: core
dimensions trail, loop dimensions broadcast.

- A vector has its component on the last axis: field and Q `(..., 3)`.
- Parameters: `S_p + (K,)`, K the number of parameters `B_lm`.
- A Hamiltonian has the loop shape `S = broadcast(S_p, S_f)`: matrices
  `S + (n, n)`, energies `S + (n,)`, n = 2J + 1, one ion per Hamiltonian.
- T and Q in an observable broadcast against `S`.
- A single point is the shape `()`. The caller picks a grid or point-by-point
  pairing with `None` axes.

## Consequences

A measurement log pairs `field (N, 3)` with `T (N,)` point by point; a grid
needs explicit new axes. A test checks batched results against a loop. Degenerate
levels only through basis-independent masks; levels sorted per point are not
branches where they cross. Memory grows with `S`: compute only the matrix
elements an observable needs, or work in chunks.
