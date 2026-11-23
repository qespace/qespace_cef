---
id: ADR-0006
title: Names of the public API
status: draft
date: 2026-11-23
supersedes:
superseded_by:
---

## Context

v1 had no naming rule (AUDIT Q2): camelCase next to snake_case, one
quantity under several names (temperature three times, the point group
four), and one name for several quantities (`J` the quantum number and the
operators, `E` the eigenvalues and the internal energy, `S` the entropy and
the scattering function). The naming audit (`NAMING.md`) listed every
public name of v1 and settled the v2 names with the maintainer. The suite
has no glossary yet; these names are its first version.

## Options

1. No rule; decide name by name.
2. PEP 8 only.
3. PEP 8 and the rules below, with the names of the v2 API.

## Decision

Option 3. The rules:

1. snake_case for modules, functions, arguments and attributes; CapWords
   for classes; no camelCase.
2. Class names are written out, without abbreviations.
3. A symbol of the formulas stays a symbol where physicists use it: `J`,
   `gJ`, `Jx`, `Jxyz`, `lam`, `u2`.
4. One name per quantity and one quantity per name, in the package and in
   its documentation.
5. A crystallographic vector carries its basis in its name: `Qxyz`, `Qhkl`,
   `Bxyz`; its length is `Qnorm`. A bare `Q` is not a name. Other vectors
   (moments, magnetisation) carry no suffix.
6. Thermodynamic quantities are written in words: `partition_function`,
   `internal_energy`, `free_energy`, `entropy`, `heat_capacity`.
7. No name stands for 2J + 1; shapes write `(2J+1, 2J+1)`. Other shape
   letters are placeholders local to one document.
8. `precision` is a number of decimal places; a convergence tolerance is
   `conv_tol`.
9. A physical constant is named after itself (`mu_B`, `k_B`); a unit
   conversion after its units (`invcm2meV`).
10. Every module lists its public names in `__all__`, and the package
    imports them explicitly; nothing uses `import *`. Private names start
    with `_`.
11. An observable takes the ion first, as `ion`, then `temperature` and the
    other inputs.

The names of the v2 API, in the modules of ADR-0005:

| Module | Names |
|---|---|
| `ion_data` | `IonData(name)`: `J`, `gJ`, `alpha`, `beta`, `gamma`, `form_factor(Qnorm)`, `form_factor_coefficients` |
| `parameters` | `CrystalFieldParameters(values, notation, units, point_group)` with `(l, m)` keys; `allowed(point_group)` |
| `ion` | `Ion(ion_data, parameters, Bxyz)`: `hamiltonian`, `energies`, `excitation_energies`, `eigenvectors`, `Jxyz`, `magnetic_moment`, `basis_labels`, `get_degeneracies()`, `summary(show_eigenvectors, precision)` |
| `operators` | `J_x`, `J_y`, `J_z`, `J_p`, `J_m`; `stevens_operator(J, l, m)`, which raises for an `(l, m)` it doesn't implement |
| `observables` | `boltzmann_population`, `scattering_weights`, `cross_section_crystal(ion, temperature, Qxyz, E_incident, u2)`, `cross_section_powder(ion, temperature, Qnorm, E_incident, u2)`, `magnetization`, `magnetization_exchange(…, lam, conv_tol, max_iterations)`, `susceptibility(ion, temperature, method='van_vleck' \| 'finite_difference')`, `g_tensor`, `thermodynamics` with the names of rule 6 |
| `constants` | `mu_B`, `k_B`, `invcm2meV`, `gamma_n`, `r0`, `gamma_r0_sq` |
| `config` | `rtol`, `E_floor` |
| `tables` | private lookups (point groups, ions) |

What each quantity is, its unit and sign, are the CONVs' (CONV-0003,
CONV-0006, CONV-0007, CONV-0012, CONV-0013); this ADR only names them.

## Consequences

Breaking, with ADR-0005, in version 2.0. `NAMING.md` keeps the map from
the v1 names. A new public name follows the rules, and a review checks
them. The suite can take the rules and the table as its first glossary.
