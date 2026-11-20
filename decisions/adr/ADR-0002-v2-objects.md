---
id: ADR-0002
title: The v2 objects Ion, Blm and Hamiltonian
status: accepted
date: 2026-11-20
supersedes:
---

## Context

v1 has `Ion`, `CEFpars` and a stateful `CEFion`: `diagonalize()` changes the
meaning of its J operators, and a second call rotates them twice. The broadcasting
API needs parameter grids and solved states of any loop shape. A mock with bare
arrays, a `NamedTuple` of parameters and a general `Eigensystem` holding `gJ` and
J, was rejected: that solved state belongs to the CEF problem.

## Options

1. Keep the v1 objects and fix only the bugs.
2. Bare arrays and functions: `hamiltonian()`, `solve()` to an `Eigensystem`.
3. The v1 physics objects, renamed, with the arrays inside them.

## Decision

Option 3, used as `import qespace_cef as cef`:

- `cef.Ion` (unchanged name): the ion, J and `gJ`.
- `cef.Blm` (was `CEFpars`): the CEF parameters by name, with values that
  broadcast to a loop shape.
- `cef.Hamiltonian` (was `CEFion`): ion, parameters and field. Once solved it
  holds the energies, the eigenvectors and J in the eigenbasis, so observables
  don't rotate operators themselves.

`OPEN:` when a `Hamiltonian` is solved and whether it is immutable.
`OPEN:` observables as functions or methods.

## Consequences

Breaking, version 2.0. J in the eigenbasis costs memory, `S + (3, n, n)`:
139 MB per 10⁴ points at J = 8. Rotated from the kept free-ion J, it cannot be
rotated twice. The other v1 names are not settled by this.
