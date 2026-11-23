---
id: CONV-0003
title: Units of energies, temperature, field, momentum and the observables
status: accepted
owner: qespace_cef
supersedes:
---

## Statement

As implemented: energies and CEF parameters in meV (parameters can also be
given in K and cm⁻¹ and are converted); temperature in K; magnetic field
`B = μ0 H` in T; |Q| in Å⁻¹; moments in μB per ion; the susceptibility
tensor in μB/T per ion (CONV-0009, accepted).

The zero of the energies is not a unit; it is CONV-0010.

## Rationale

meV, K, T and Å⁻¹ are the units of neutron spectroscopy and magnetometry, where
the results are compared. One unit per quantity avoids the factor errors of
mixed units: a heat capacity computed from energies in meV and temperatures
in K is off by about 11.6², the square of the meV-to-K conversion. McPhase
gives its CEF parameters in meV. Units other qespace packages rely on may
come from the suite.

## Consequences

- Every public function states the unit of every argument and return.
- Conversions happen at the boundary, with constants from one module that
  cites its CODATA set.
- A test checks each observable in its unit, e.g. the Curie constant at high T.
