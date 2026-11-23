---
id: ADR-0005
title: The v2 objects FreeIon, CrystalFieldParameters and SingleIon
status: draft
date: 2026-11-23
supersedes: ADR-0002
---

## Context

ADR-0002 named the v2 objects `Ion`, `Blm` and `Hamiltonian`. The naming
audit (`NAMING.md`) found that `Ion` is the free ion, while the object one
studies is the ion in its crystal field; that `Blm` names one notation of
the parameters, while other notations (CONV-0001, Table 5) are to come; and
that `Hamiltonian` names one matrix of an object that holds the whole
solved problem. McPhase calls that object the single ion, Mantid
`CrystalField`, PyCrystalField `CFLevels`. The other parts of ADR-0002
still hold: the v1 physics objects, renamed, with the arrays inside them,
and no bare-array `Eigensystem`.

## Options

1. Keep ADR-0002's names.
2. `FreeIon`, `CrystalFieldParameters`, `SingleIon`, one container for every
   notation, the notation given by a keyword.
3. As 2, with one class per notation.
4. As 2, with the notation in the key's letter (`'B20'`, `'A20'`), as
   McPhase and PyCrystalField do.

## Decision

Option 2, used as `import qespace_cef as cef`:

- `cef.FreeIon(name)` (v1 `Ion`): the free ion, J, `gJ`, the Stevens
  factors and the form factor.
- `cef.CrystalFieldParameters(values, notation=..., units=...,
  point_group=...)` (v1 `CEFpars`): the parameters of the crystal field in one notation of
  CONV-0001. `values` maps `(l, m)` tuples to numbers or arrays that
  broadcast to a loop shape (ADR-0003); `notation` is McPhase's name of the
  notation. The keys never carry the notation. Not `CrystalField`: in some
  notations, Stevens among them, the parameters belong to an ion
  (CONV-0001).
- `cef.SingleIon(free_ion, parameters, field)` (v1 `CEFion`): the ion in
  its crystal field and a field. It converts the parameters to the notation
  of its Hamiltonian, with the factors of CONV-0001, and holds the solved
  state: the Hamiltonian matrix, the energies, the eigenvectors and `Jxyz`,
  J in the eigenbasis, so observables don't rotate operators themselves.

Names inside these objects follow `NAMING.md` until the naming rule is an
ADR.

`OPEN:` when a `SingleIon` is solved and whether it is immutable (Q29).
`OPEN:` observables as functions or methods (Q31).
`OPEN:` the name of the field argument (Q59).

## Consequences

Breaking, version 2.0: v1 scripts fail at the first call, not silently.
A notation is added by its conversion, not by a new class. A
`CrystalFieldParameters` is built without choosing an ion; the ion comes
with `SingleIon`. `Jxyz` costs memory, `S + (3, 2J+1, 2J+1)`: 139 MB per 10⁴
points at J = 8. Rotated from the kept free-ion J, it cannot be rotated
twice.
