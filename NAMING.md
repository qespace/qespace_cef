# Naming audit

The names of v1 and the names they get in v2: step 1 and 2 of the naming
audit in [`PLAN.md`](PLAN.md). This file tracks; it is not a rule
(ADR-0004). The rule it leads to is a draft ADR-0006 (Q26), and it is the
first version of the suite glossary, which does not exist yet (Q56).

Taken from `main` at `a1fac60` by introspection (`dir`,
`inspect.signature`, instance `vars`, docstring sections) and from the
accepted CONVs and ADRs. **Decided** marks your answers of 2026-11-22;
everything else is a proposal. `symbolic/`, the examples and private names
are left out.

## Rules so far

1. snake_case for functions, arguments and attributes; CapWords for
   classes; no camelCase.
2. A physics symbol is kept where the formulas use it: `J`, `Jx`, `gJ`.
3. A crystallographic vector carries its basis in its name: `Qxyz`,
   `Qhkl`. Its length is `Qnorm`. A bare `Q` never appears. Other vectors
   (moments, magnetisation) don't (Q60).
4. Thermodynamic quantities are written in words: `entropy`,
   `heat_capacity`.
5. No letter stands for 2J + 1: shapes write `(2J+1, 2J+1)`. Other shape
   letters (`n`, `S`, ...) are placeholders local to one document (Q61).
6. `precision` is a number of decimal places; a convergence tolerance is
   `conv_tol`.
7. Class names are written out, no abbreviations: `CrystalFieldParameters`,
   not `CrystalFieldPars`.

## Objects and modules

| Concept | v1 | Where | v2 | Note |
|---|---|---|---|---|
| free ion | `Ion(ion_name)` | `ion` | **Decided:** `FreeIon(name)` | ADR-0002 said `Ion` |
| CEF parameters | `CEFpars(Bdict, units, pointGroup)` | `cefpars` | **Decided:** `CrystalFieldParameters(values, notation=..., units=..., point_group=...)` (Q64) | ADR-0002 said `Blm`; **decided:** keyword `notation`, as McPhase |
| ion in the CEF and a field | `CEFion(ion, cfp, Hfield, diagonalize)` | `cefion` | **Decided:** `SingleIon(free_ion, parameters, Bxyz)` (Q63) | ADR-0002 said `Hamiltonian`; McPhase's single ion, coupled in RPA later |
| Stevens operators | `StevensBase(J)`, `.O4m3()`, `['4m3']` | `cefmatrices` | one function, e.g. `stevens_operator(J, l, m)` | A6, A14 |
| modules | `ion`, `cefion`, `cefpars`, `cefmatrices`, `data_containers` | | e.g. `free_ion`, `single_ion`, `parameters`, `operators`, `tables` | the `cef` prefix repeats the package |
| leaked names | `np`, `C`, `ArrayLike`, the submodules | `qespace_cef.*` | `__all__` | A4 |

## Quantum numbers and operators

| Concept | v1 | Where | v2 | Note |
|---|---|---|---|---|
| total angular momentum J | `Ion.J`, `CEFion.Jval`, `StevensBase.J` | 3 places | **Decided:** `J`, the quantum number only | |
| 2J + 1 | `Ion.J2p1`, `StevensBase.J2p1` | | **Decided:** no name (rule 5) | |
| J(J+1) | `StevensBase.JJ` | | none | |
| J operators, functions | `J_x`, `J_y`, `J_z`, `J_plus`, `J_minus` | `cefmatrices` | **Decided:** `J_x`, `J_y`, `J_z`, `J_p`, `J_m` | CONV-0005 root |
| J operators, attributes | `StevensBase.Jx, Jy, Jz, Jp, Jm`; `CEFion.Jx, Jy, Jz` | | **Decided:** `Jx`, `Jy`, `Jz`, `Jp`, `Jm` | |
| vector of J operators | `CEFion.J`, shape (3, 2J+1, 2J+1) | `cefion` | **Decided:** `Jxyz` | `J` was also the quantum number |
| Landé factor | `Ion.gJ` | | `gJ` | |
| Stevens factors α, β, γ | `Ion.Alpha`, `.Beta`, `.Gamma` | | **Decided:** `alpha`, `beta`, `gamma` | γ kept in its context |
| form factor | `Ion.mff(Q)`, `Ion.mffCoefficients` | | `form_factor(Qnorm)`, `form_factor_coefficients` | |
| free-ion kets | `CEFion.freeionkets` | | `basis_labels` | |

## Parameters, symmetry, field

| Concept | v1 | Where | v2 | Note |
|---|---|---|---|---|
| parameter set | `Bdict` (arg), `Bpars` (attr), `cfp` (arg of `CEFion`) | `cefpars`, `cefion` | `CrystalFieldParameters`; argument `parameters` | three names |
| parameter keys | `'B20'`, `'B4m3'` (strings, parsed) | `CEFpars`, `StevensBase['4m3']` | **Decided:** `(l, m)` tuples, `(4, -3)`; the notation only in the keyword (Q65) | no `m` parsing (A14) |
| notation | Stevens only | | keyword `notation='stevens'`, `'wybourne'`, ... (CONV-0001, Table 5); converted with the `FreeIon` in `SingleIon` | |
| parameter units | `units` in {'meV', 'K', 'invcm'} | `CEFpars` | `units` | stored in meV |
| point group | `pointGroup` (arg), `pointGroupName` (arg), `symmetry` (attr), `point_group_synonyms` | `cefpars`, `data_containers` | `point_group`, on `CrystalFieldParameters` | four names |
| allowed parameters | `CEFpars.allowed_Bpars(pointGroupName)` | | e.g. `CrystalFieldParameters.allowed(point_group)` | |
| implemented operators | `StevensBase.implementedOps` | | `implemented_operators` | |
| magnetic field | `Hfield` | `CEFion`, `magnetization*` | **Decided:** `Bxyz`, the induction B in T (Q59) | CONV-0006 (draft): ℋ_Z = −μ·B, μ = −gJ μB J |

## Solved state

| Concept | v1 | Unit, shape | v2 | Note |
|---|---|---|---|---|
| Hamiltonian matrix | `CEFion.hamiltonian` | meV, (2J+1, 2J+1) | `hamiltonian` | |
| solve | `CEFion.diagonalize(Eimag_precision)` → `E`, `U` | `E` a diagonal (2J+1, 2J+1) matrix | **Decided:** no diagonal matrix; with Q29 nothing to call | |
| energies E_n | `CEFion.energies` | meV, (2J+1,) | `energies` | CONV-0010 |
| E_n − E_0 | `CEFion.excitation_energies` | meV, (2J+1,) | goes to the print | A16 |
| E_f − E_i | first return of the spectra | meV, ((2J+1)²,) | **Decided:** `excitation_energies`, (2J+1, 2J+1), as in CONV-0012 | A16 |
| eigenvectors | `CEFion.eigenvectors` | (2J+1, 2J+1), columns | `eigenvectors` | |
| moment per state | `CEFion.magnetic_moment` | μB, (2J+1, 3) | `magnetic_moment` | not crystallographic: no suffix (Q60) |
| levels | `degenerate_levels(energies)` → slices; `CEFion.get_degeneracies()` → [[E, g]] | | **Decided:** private helper `_degenerate_levels` on energies; `get_degeneracies` on `SingleIon` | renamed by you on the working tree |
| level print | `CEFion.summary(evecs, precision)` | str | `summary(show_eigenvectors, precision)` | rule 6 |

## Observables

| Concept | v1 | Unit, shape | v2 | Note |
|---|---|---|---|---|
| temperature | `temperature`, `temperatures`, `T` | K | **Decided:** `temperature` | three names |
| populations | `boltzman_population(energies, temperature)` | (2J+1,) | **Decided:** `boltzmann_population` | misspelt |
| scattering weights | `scattering_function(cefion, temperature)` | | **Decided:** `scattering_weights` (Q55) | |
| cross sections | `cross_section_crystal(…, Q, Ei, u2)`, `cross_section_powder(…, Qnorm, Ei, u2)` | barn sr⁻¹ per ion | `Qxyz`, `Qnorm`, `E_incident` (**decided**), `u2` | |
| magnetisation | `magnetization(cefion, temperature, Hfield)` | μB, tuple of 3 | **Decided:** returns an array; field `Bxyz` | |
| mean field | `magnetization_exchange(…, lam, precision, max_iterations)` | λ in T per μB | `lam`, **decided:** `conv_tol`, `max_iterations` | |
| susceptibility | `susceptibility(cefion, temperatures, method)` | μB/T per ion, (…, 3, 3) | `temperature` | |
| g-tensor | `gtensor(cefion)` | (3, 3) | `g_tensor` | C7 work |
| thermodynamics | `thermodynamics(cefion, T)` → `Z, E, F, S, Cv` | meV, k_B (P3) | **Decided:** in words: `partition_function`, `internal_energy`, `free_energy`, `entropy`, `heat_capacity` | |

## Constants and config

| v1 | Unit | v2 | Note |
|---|---|---|---|
| `uB_SI`, `uB` | J/T, meV/T | `mu_B` (meV/T); SI private | |
| `K2meV`, `meV2K` | meV/K, K/meV | `k_B` (meV/K) | |
| `invcm2meV` | meV per cm⁻¹ | `invcm2meV` | |
| `gamma_n`, `r0`, `gamma_r0_sq` | –, m, barn | keep | γ in its context |
| `config.rtol`, `config.E_floor` | –, meV | keep | |

## Collisions: one name, two quantities

| ID | Name | Quantities | Decided |
|---|---|---|---|
| C1 | excitation energies | E_n − E_0 (CONV-0010) and E_f − E_i (CONV-0012) | E_f − E_i |
| C2 | `E`, `E_i` | eigenvalues, internal energy, the matrix of `diagonalize`; initial state and incident energy (CONV-0007) | `energies`, `internal_energy`, no matrix, `E_incident` |
| C3 | `S` | entropy, scattering function, loop shape (ADR-0003) | the scattering function; entropy in words; shape letters are local (rule 5) |
| C4 | `n` | 2J + 1, a count of temperatures | a placeholder only (rule 5) |
| C5 | `J` | quantum number, operator vector | the quantum number; operators `Jxyz` |
| C6 | `precision` | decimals, convergence tolerance | decimals; `conv_tol` |
| C7 | γ, Γ | Stevens factor, irreps, `gamma_n` | each in its context |
| C8 | `hamiltonian`, `Hamiltonian` | matrix, v2 class | the class is `SingleIon`; `hamiltonian` stays the matrix |

## How other codes name them

| Code | Ion in the CEF | Parameters | Notation chosen by |
|---|---|---|---|
| McPhase | single ion (`so1ion`, single ion property file) | `B20 = …` lines | the key's letter (CONV-0001, Table 5) |
| PyCrystalField | `CFLevels`, `CFLevels.Bdict(ion, Bdict)` | dict `Bdict`, keys `'B20'` | the letter: `Adict` (`'A20'`) through `WybourneToStevens` |
| Mantid | `CrystalField('Ce', 'C2v', B20=…)`, ion, symmetry and solver in one | keywords `B20`, `IB22` | Stevens only |

## Decisions these need

- ADR-0005 (draft) supersedes ADR-0002: `FreeIon`, `CrystalFieldParameters`,
  `SingleIon`. What each notation means and its conversion factors stay in
  CONV-0001 (ADR-0004).
- The naming rule (Q26) as ADR-0006, from the rules above.

## Open

None.
