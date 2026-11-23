---
id: ADR-0005
title: The v2 package qespace_cft and its objects IonData, CrystalFieldParameters and Ion
status: draft
date: 2026-11-23
supersedes: ADR-0002
superseded_by:
---

## Context

ADR-0002 named the v2 objects `Ion`, `Blm` and `Hamiltonian`. The naming
audit (`NAMING.md`) found that its `Ion` only holds data of the free ion;
that `Blm` names one notation of the parameters, while other notations
(CONV-0001, Table 5) are to come; and that `Hamiltonian` names one matrix
of the object that holds the whole problem. McPhase calls that object the
single ion, Mantid `CrystalField`, PyCrystalField `CFLevels`. The suite
splits into packages by theory: crystallography, linear spin-wave theory,
later the random phase approximation; this one holds crystal-field theory.
The other parts of ADR-0002 still hold: the v1 physics objects, renamed,
with the arrays inside them, and no bare-array `Eigensystem`.

## Options

1. Keep ADR-0002's names and the package name `qespace_cef`.
2. The package `qespace_cft`; `IonData`, `CrystalFieldParameters`, `Ion`;
   one container for every notation, the notation given by a keyword.
3. As 2, with one class per notation.
4. As 2, with the notation in the key's letter (`'B20'`, `'A20'`), as
   McPhase and PyCrystalField do.
5. As 2, with the main object named `SingleIon`, `CrystalFieldIon` or
   `CrystalField`.

## Decision

Option 2. The package is `qespace_cft`, crystal-field theory, used as
`import qespace_cft as cft`:

- `cft.IonData(name)` (v1 `Ion`): the data of the free ion: J, `gJ`, the
  Stevens factors, the form-factor coefficients, and what is derived from
  them, such as the form factor. It holds no operators.
- `cft.CrystalFieldParameters(values, notation=..., units=...,
  point_group=...)` (v1 `CEFpars`): the parameters of the crystal field in
  one notation of CONV-0001. `values` maps `(l, m)` tuples to numbers or
  arrays that broadcast to a loop shape (ADR-0003); `notation` is McPhase's
  name of the notation. The keys never carry the notation. Not
  `CrystalField`: in some notations, Stevens among them, the parameters
  belong to an ion (CONV-0001).
- `cft.Ion(ion_data, parameters, Bxyz)` (v1 `CEFion`): the ion in its
  crystal field and the field B of CONV-0006, its xyz components in the
  local frame (CONV-0004), and the entry point of the package. It converts
  the parameters to the notation of its Hamiltonian, with the factors of
  CONV-0001, and holds the solved state: the Hamiltonian matrix, the
  energies, the eigenvectors and `Jxyz`, J in the eigenbasis, so
  observables don't rotate operators themselves. A short name, because
  everything starts from it, also in packages that couple ions.

Public names are listed in `__all__` of the package and of every module and
imported explicitly; nothing uses `import *`. Names inside these objects
follow `NAMING.md` until the naming rule is an ADR.

`OPEN:` when an `Ion` is solved and whether it is immutable (Q29).
`OPEN:` observables as functions or methods (Q31).

## Consequences

Breaking, version 2.0, together with the new package name: every import
changes, and v1 scripts fail at the first call, not silently. `Ion` now
names the ion in its crystal field, not the free ion of v1. The repository
is renamed when the names are final. A notation is added by its
conversion, not by a new class. A `CrystalFieldParameters` is built
without choosing an ion; the ion comes with `Ion`. `Jxyz` costs memory,
`S + (3, 2J+1, 2J+1)`: 139 MB per 10⁴ points at J = 8. Rotated from the
kept free-ion J, it cannot be rotated twice.
