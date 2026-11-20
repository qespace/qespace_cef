---
id: CONV-0003
title: Units of energies, temperature, field, momentum and the observables
status: draft
owner: qespace_cef
enforced_by:
since:
supersedes:
---

## Statement

`OPEN:` the units of the API.

As implemented: energies and CEF parameters in meV (parameters can also be
given in K and cm⁻¹ and are converted); level energies shifted so the ground
state is at 0; temperature in K; magnetic field `B = μ0 H` in T; |Q| in Å⁻¹;
moments in μB per ion; the susceptibility in μB²/meV per ion with
`method='perturbation'` and in μB/T per ion with `method='magnetization'`.

`OPEN:` one unit for the susceptibility (μB/T per ion, or emu/mol), and
whether level energies are always shifted to the ground state.

## Rationale

meV, K, T and Å⁻¹ are the units of neutron spectroscopy and magnetometry, where
the results are compared. One unit per quantity avoids the factor errors of
mixed units: a heat capacity computed from energies in meV and temperatures
in K is off by about 11.6², the square of the meV-to-K conversion. Units other qespace packages rely on may come
from the suite.

## Consequences

- Every public function states the unit of every argument and return.
- Conversions happen at the boundary, with constants from one module that
  cites its CODATA set.
- A test checks each observable in its unit, e.g. the Curie constant at high T.
