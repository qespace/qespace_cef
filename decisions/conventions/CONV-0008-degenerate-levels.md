---
id: CONV-0008
title: Degenerate levels
status: draft
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

Eigenvectors within a degenerate level are basis dependent: results and tests
use quantities that don't depend on that basis (projectors onto the level,
traces, sums over the level), never single eigenvectors.

`OPEN:` when two levels count as degenerate: a tolerance relative to the
energy scale of the Hamiltonian, from one config place.

## Rationale

If |a> and |b> have the same energy, so does every combination
`α|a> + β|b>`. The solver returns one orthonormal pair out of infinitely many,
chosen by rounding, so it can change between numpy versions or platforms.
Anything computed from one vector of the pair, e.g. `<a|J_z|a>` in a Kramers
doublet, can then take any value between −j and +j, where ±j are the
eigenvalues of J_z within the doublet. What is
fixed is the level as a whole: the projector `P = |a><a| + |b><b|`, and sums
such as `Σ_{i,j ∈ level} |<i|J_α|j>|²`. The g-tensor is an example: g depends
on the basis inside the doublet, `G = g gᵀ` does not. In a small field the
degeneracy lifts and the eigenvectors become well defined, but they jump as
the field direction changes.

## Consequences

- Observables sum over each degenerate level; tests compare projectors or sums,
  or rotate the basis inside the level and check that the result is unchanged.
- Level grouping uses the tolerance, never rounding to a number of decimals.
