# qespace_cef v2 release and naming audit

**23 November 2026 · main at `a1fac603dadb` · package metadata 1.3.0 - by GPT-6.2 Astra**

**Release assessment.** Hold the v2 API freeze. This checkout has a working 1.x public interface and a substantial test suite. Six high-severity findings affect scientific results or symmetry selection. Two more affect the contract that downstream code would adopt. The most consequential naming changes concern which energy difference a name denotes, which axis selects an eigenstate, and which field an observable uses.

| Severity | Findings | Meaning in this review |
|---|---:|---|
| Critical | 0 | A failure invalidating the package broadly |
| High | 8 | A wrong scientific result on accepted input, or incompatible definitions of a central v2 API |
| Medium | 10 | A narrower correctness, migration, validation, or test gap with an identifiable consequence |
| Low | 4 | A bounded edge case, documentation, or release-traceability issue |

The review covers the current working tree, including the naming inventory. Before the report was created, `CLAUDE.md` differed from the local Git index and `NAMING.md` was untracked. The latter's SHA-256 begins `72ab73f7db6c24b1`. “Recorded” below means the inventory marks a maintainer choice as decided. Accepted ADRs and CONVs are identified separately. Numerical examples use synthetic parameters in meV unless otherwise indicated. [Naming inventory](NAMING.md#L3) · [Decision authority](decisions/adr/ADR-0004-rules-descriptions-tracking.md#L28)

## High severity

### H1. Zero-field susceptibility depends on the field of the supplied ion

Both methods promise the response at zero field and document that the input ion's field is ignored. The default perturbation method reads its already field-split energies and operators. The finite-difference method rebuilds the ion in small fields around zero. For free Ce³⁺ at 2 K, supplying an ion solved at 10 T along z gave **1.485173** for the default `χzz` and **0.719693** for the finite difference, in μB/T per ion. The analytic zero-field Curie value is **0.719693**. The perturbation formula also omits the connected-moment subtraction a finite-field response would require. Use a zero-field eigensystem for the advertised quantity and cover both methods with a finite-field input. [Field contract](qespace_cef/observables.py#L395) · [Perturbation branch](qespace_cef/observables.py#L426) · [CONV-0009](decisions/conventions/CONV-0009-susceptibility.md#L11) · [Existing method comparison](tests/test_observables.py#L222)

### H2. Thermodynamic entropy changes when the energy origin moves

`thermodynamics` combines internal energy in meV with `−T log Z` in kelvin. For Ce³⁺ with `B20=0.5, B40=0.01, B44=0.1` at 10 K, the raw entropy and heat capacity were **4.854638** and **0.00126422**. Direct Boltzmann sums give **1.292718** and **0.170246** in kB per ion. Adding only `B00=1` changed the returned entropy by **−1.060452**, although every population was identical up to rounding. This reproduces P3 with the current energy zero. The value 1.054 quoted in the existing audit and docstring is historical. Record the units under Q44, then check the offset invariance and derivatives in those units. [Calculation](qespace_cef/observables.py#L576) · [CONV-0010](decisions/conventions/CONV-0010-energy-zero.md#L22) · [P3](AUDIT.md#L30)

### H3. At a positive cryogenic temperature, finite thermodynamic observables become NaN

The routine exponentiates absolute eigenvalues before normalizing. With Ce³⁺, `B20=1` and `T=np.array([0.1])` K, it returned `(inf, nan, -inf, nan, nan)` for `(Z, U, F, S, Cv)`. The ground doublet is at −8 meV. Here the absolute partition function exceeds the float64 range, while the internal energy approaches −8 meV and the entropy approaches ln 2 in kB per ion. This failure survives a correction of the unit factors alone. Normalize using energy differences, retain the energy offset for U and F, and evaluate the variance about the thermal mean. Test the representable outputs even when Z overflows. [Exponential sums](qespace_cef/observables.py#L568) · [Existing partition-function test](tests/test_observables.py#L280)

### H4. The g-tensor can describe the wrong subspace and give wrong principal values

The ground-level guard accepts **any multiplicity ≥ 2**, then uses only states 0 and 1. Free Ce³⁺ has one sixfold level, yet `gtensor` returned diagonal entries `[1.916630, −1.916630, 4.285714]`. Optimized Python also removes the `assert` guard. For an actual axial Ce³⁺ ground doublet (`B20=1`), a unitary basis change moved the principal values from `[0.857143, 2.571429, 2.571429]` to `[0.857126, 2.431936, 2.703740]`. The latter reproduces P17. Enforce exactly two states with an ordinary runtime check, and settle C7's basis-independent return before publishing `g_tensor`. [Projection and guard](qespace_cef/observables.py#L506) · [Principal-value tests](tests/test_observables.py#L594) · [P4 and P17](AUDIT.md#L31)

### H5. The level summary can print an excited-state ket as a ground state

The constructor accepts odd-rank terms. `summary` assumes Kramers pairing from half-integer J and zero applied field alone, and constructs a partner by time reversal. For Ce³⁺ with `B20=1, B10=−6`, the exact ground states are `|+1/2>` and `|+3/2>` at −11 meV. The summary printed `|+1/2>` and **`|−1/2>`**, whose energy is −5 meV. This is a wrong state even for a diagonal Hamiltonian. Q12 can remain open while pairing is restricted to Hamiltonians that support it. A regression check can verify that every displayed state lies inside its reported eigenspace. [Allowed ranks](qespace_cef/cefpars.py#L83) · [Pair construction](qespace_cef/cefion.py#L123) · [Pairing condition](qespace_cef/cefion.py#L422)

### H6. The monoclinic lookup returns operators forbidden by its stated symmetry

`CEFpars.allowed_Bpars("C2")` includes `B43`, `B4m3`, `B63`, and `B6m3`. With the documented twofold axis along z, a rotation by π changes the sign of each order-three term. For J = 7/2, rotating `O43` gave `−O43` to **3.7 × 10⁻¹⁶** relative to its largest entry. A fit built from this list can therefore break its declared site symmetry. This is independent of the tracked absence of cubic parameter relations. The companion lookup also advertises all 32 point groups but contains **31**. Both `"-1"` and `"Ci"` raise `NameError` for the missing inversion group. Check all group identifiers and verify the listed operators against the generators of each stated setting. [Axis and lookup table](qespace_cef/cefpars.py#L93) · [Monoclinic entries](qespace_cef/cefpars.py#L131) · [Point-group records](qespace_cef/data_containers.py#L39) · [Current coverage](tests/test_cefpars.py#L36)

### H7. The recorded v2 object names and the accepted design disagree

| Role | Importable in this checkout | Accepted ADR-0002 | Recorded in NAMING.md |
|---|---|---|---|
| Free ion | `Ion` | `Ion` | `FreeIon` |
| Crystal-field parameters | `CEFpars` | `Blm` | `CrystalField` |
| Ion in the CEF and field | `CEFion` | `Hamiltonian` | `SingleIon` |

The successor ADR and proposed naming ADR-0005 are absent. ADR-0002 and the plan still leave solve timing, immutability, and functions versus methods open. The inventory assumes there will be nothing to call to solve the ion, and A16 assumes immutability. A downstream caller following either document can design for a different lifecycle. Record the maintainer's choices in a superseding ADR, including constructor keywords, field ownership, and the status of the old import paths. The working 1.x exports are the migration baseline. [Current exports](qespace_cef/__init__.py#L10) · [Accepted design](decisions/adr/ADR-0002-v2-objects.md#L25) · [Recorded names](NAMING.md#L31) · [Open questions](PLAN.md#L37) · [A16](AUDIT.md#L154)

### H8. The same excitation-energy name would acquire a different physical meaning and rank

Today `CEFion.excitation_energies[i]` means `E[i] − E[0]`. Spectral functions return flat signed `E[f] − E[i]` and weights, including elastic and neutron-energy-gain entries. The inventory assigns `excitation_energies[i, f]` to the pairwise quantity and moves the ground-relative values to the printout. For Ce³⁺ these are respectively six values, 36 flat entries, and a proposed 6 × 6 matrix. A consumer indexing `[i]` could silently get an entire row.

CONV-0010 specifies the ground-relative property and CONV-0012 explicitly specifies flat spectra. The plan assigns the new shapes to a CONV despite ADR-0004 assigning names and shapes to ADRs. Preserve the definitions of the three quantities, record the public spelling, axis order and broadcast shapes in the design decision, and give a migration example. A fixed pair of state indices identifies entries for one solve. Rotating a degenerate subspace changes individual weights. Only the appropriate sums over states are basis independent. [Current spectral layout](qespace_cef/observables.py#L79) · [Energy convention](decisions/conventions/CONV-0010-energy-zero.md#L18) · [Spectral convention](decisions/conventions/CONV-0012-spectral-weights.md#L24) · [Recorded choice](NAMING.md#L76) · [Plan](PLAN.md#L28) · [ADR-0004](decisions/adr/ADR-0004-rules-descriptions-tracking.md#L32)

## Medium severity

### M1. The migration inventory transposes stored eigenvectors

The inventory describes the v1 `eigenvectors` attribute as columns. Only the `U` returned by `diagonalize()` uses columns. The stored attribute is `U.T`, with one ket per **row**. For the Ce³⁺ parameters from finding H2 with `B4m2=0.03` and field `[0.1, 0.2, 0.3]` T, the largest component of `H v − E v` was **2.7 × 10⁻¹⁵ meV** using rows and **7.55 meV** using columns. The square array shape hides the error. Correct the v1 inventory and explicitly specify the v2 state axis using a complex, low-symmetry example. [Inventory row](NAMING.md#L77) · [Assignment](qespace_cef/cefion.py#L316) · [Row-eigenpair test](tests/test_cefion.py#L229)

### M2. Tuple keys leave the coefficient convention ambiguous

The recorded `(l, m)` key removes the notation letter. Proposed values such as `notation="stevens"` or `"wybourne"` identify families. CONV-0001's Table 5 lists multiple members of each family with different factors, units, or real/complex representations. An `Alm` coefficient also requires radial information that `Blm` already incorporates. Define the exact supported coefficient set for each keyword value, its units and its negative-m convention. A Stevens-`Blm`-only first release can state that scope explicitly. A migration example such as `"B4m3" → (4, -3)` should preserve the Hamiltonian matrix, including its sign. [Naming proposal](NAMING.md#L59) · [Table 5](decisions/conventions/CONV-0001-CEF-parameters.md#L67) · [Existing label-to-matrix check](tests/test_cefpars.py#L64)

### M3. Field, momentum and temperature need explicit frame and axis contracts

`SingleIon(..., field)`, `Bxyz`, and `Hxyz` all appear in the inventory. CONV-0004 places both Q and the field in the **local Cartesian frame of the CEF parameters**. An HKL vector needs reciprocal-lattice information and a local-frame transform. The `Qxyz` argument needs that local-frame definition in its docstring. Document the field in tesla, its owner, the supported Q representation, and which observable accepts a field override. Leave the sign decision to CONV-0006.

The scalar-temperature v1 contract also needs a concrete broadcast migration example. Passing six different temperatures to `scattering_function` for Ce³⁺ silently returned **one** set of 36 energies and 36 weight matrices. NumPy paired one temperature with each eigenstate inside the population formula. A v2 test should compare a temperature scan with six separate scalar calls and distinguish paired data from an outer grid. Specify whether Q-independent energies retain their smaller shape when the weights carry a Q axis. [Field and Q spellings](NAMING.md#L19) · [Unresolved field keyword](NAMING.md#L133) · [Local frame](decisions/conventions/CONV-0004-frames.md#L11) · [Population calculation](qespace_cef/observables.py#L36) · [Broadcasting decision](decisions/adr/ADR-0003-broadcasting.md#L23)

### M4. A solved ion can simultaneously describe two different Hamiltonians

Changing `h.cfp.Bpars["B20"]` from +1 to −1 after solving Ce³⁺ leaves the stored ground energy at **−8 meV**. Rebuilding from that same mapping gives **−10 meV**. On the original object at 2 K, the two susceptibility methods then give **0.061688** and **1.542200** μB/T per ion along z. A related path, `diagonalize=False`, still exposes a populated `magnetic_moment` array while `energies` is empty. Specify ownership of input arrays and parameters and when derived properties become valid. A lazy solution requires protection against mutation of its inputs and cached arrays. [Stored references](qespace_cef/cefion.py#L226) · [Moment property](qespace_cef/cefion.py#L331) · [Reconstruction in magnetization](qespace_cef/observables.py#L278) · [A2–A3](AUDIT.md#L45)

### M5. An identity shift can change both the reported degeneracy and susceptibility

With `B20=5e-11`, Ce³⁺ gives three doublets. Adding `B00=10` merges them into one sixfold level because the tolerance scales with the largest absolute eigenvalue. At 10⁻⁹ K, a deliberately sensitive boundary case, `χxx` changed from **1.616 × 10⁹** to **2.083 × 10⁹** μB/T per ion. The energy gaps are unchanged to floating-point accuracy. The grouping follows CONV-0008 as written, while CONV-0010 requires susceptibility to be independent of the energy zero. Reconcile the grouping scale and the evaluation of near-degenerate thermal weights before using grouped levels as v2 identifiers. [Grouping](qespace_cef/cefion.py#L50) · [Thermal weights](qespace_cef/observables.py#L428) · [CONV-0008](decisions/conventions/CONV-0008-degenerate-levels.md#L11) · [CONV-0010](decisions/conventions/CONV-0010-energy-zero.md#L22)

### M6. Several public inputs are silently reinterpreted

| Input | Observed consequence |
|---|---|
| `CEFion(..., Hfield=[1])` | Identical Hamiltonian to `Hfield=[1, 1, 1]`. The v1 signature promises three real components. |
| `CEFpars({"B20": 1, "O20": 1})` | Identical Hamiltonian to `B20=2`. An operator label is accepted as a second spelling for a parameter. |
| `degenerate_levels([0, nan, 100])` | Reports one three-state level spanning the entire input. |
| Crystal Q of `[nan, 0, 0]`, free Ce³⁺, `Ei=50` | All 36 weights are NaN although every transition is kinematically accessible. Powder Q validation rejects non-finite lengths. |
| Populations for energies `[0, 0, 1]` meV at 0 K | All three probabilities are NaN. Specify the supported temperature domain or implement its limiting behavior. |
| `J_z(0.75)`, `J_plus(0.75)` | A non-quantized J gives matrices with adjacent m values 1.5 apart. The largest error in `[J_z, J_+] = J_+` is 0.61237. |

Validate the field's component axis, parameter spelling and finite inputs at their public entry points. Hermiticity also uses two independent thresholds, default `allclose` at construction and decimal rounding on solve. A change of scale can alter which matrices pass. [Field and Hermiticity](qespace_cef/cefion.py#L249) · [Second Hermiticity check](qespace_cef/cefion.py#L300) · [Parameter spelling](qespace_cef/cefpars.py#L81) · [Level validation](qespace_cef/cefion.py#L50) · [Crystal Q](qespace_cef/observables.py#L187) · [Thermal probabilities](qespace_cef/observables.py#L36) · [Angular momentum inputs](qespace_cef/cefmatrices.py#L22)

### M7. The fitting example cannot run its advertised scan

The example imports `cef` but uses `qc.Ion`. It calls undefined module-level `solve` and `hamiltonian` after moving the implementation into a class. Even after those names are repaired, its three `(21,)` inputs broadcast to **21 paired points**, while the example claims **9,261 Cartesian combinations**. Its mock solver also subtracts the ground energy from `energies`, unlike the accepted absolute-energy convention. The file advertises a runnable command and the roadmap cites its timings. Restore a runnable example or mark it as a historical sketch outside the runnable examples. A check on the scan shape would catch the scientific error after the name errors are fixed. [Imports and run command](examples/fit_cef_mock.py#L14) · [Mock solve](examples/fit_cef_mock.py#L87) · [Model calls](examples/fit_cef_mock.py#L157) · [Scan inputs](examples/fit_cef_mock.py#L215) · [Roadmap](ROADMAP.md#L32)

### M8. Six strict expected failures use the wrong scale for the proposed thermodynamic units

If Q44 adopts meV per ion for U and F and kB per ion for S and Cv, the numerical identities are `S = −(dF/dT)/k_B` and `Cv = (dU/dT)/k_B`, with `k_B` in meV/K. Two parameterized tests compare the unconverted derivatives at three temperatures each. A physically correct P3 fix could therefore leave six tests xfailed. The separate identity `Cv = T dS/dT` is already dimensionally consistent under that choice. Align the oracles with the accepted units and include the invariant B00 shift from H2. [Derivative tests](tests/test_observables.py#L289) · [Proposed units](AUDIT.md#L146)

### M9. The Ce form factor uses a different charge state, and the Stevens factors lack provenance

All **140** stored `j0` and `j2` coefficients for the ten ions matched the bundled numerical table exactly. For Ce³⁺ the matched rows are explicitly **Ce²⁺**. The package documents this approximation. A nonzero-Q absolute cross section for Ce³⁺ inherits it. The tests check the near-unity limit and overall decrease for most ions. The Stevens factors also lack a source beside the values, which matters if v2 starts using them in notation conversion. Carry the charge-state provenance into the new API and validate any conversion against an independent reference before advertising it. [Ion data and warning](qespace_cef/ion.py#L40) · [Bundled Ce j0](qespace_cef/atoms-form-factors.txt#L78) · [Bundled Ce j2](qespace_cef/atoms-form-factors.txt#L187) · [Form-factor tests](tests/test_ion.py#L34)

### M10. The public export and migration inventory leaves several caller-visible choices open

The wildcard export makes `np`, `C`, `ArrayLike`, and `degenerate_levels` reachable at package root. A caller can already import the level helper there. Adding `__all__` controls star imports, but existing explicit paths also need a disposition. The inventory gives one v2 name, `k_B`, for two reciprocal constants. Multiplication by v1 `K2meV` becomes multiplication by `k_B`, while multiplication by `meV2K` becomes **division** by it. The five thermodynamic names could denote positional returns, named fields, or separate calls. Likewise `scattering_weights` leaves open whether it returns weights alone or the existing energy/weight pair. Specify import path, keyword, units, return structure, and loop/core axes for each public entry point. Preserve the distinction between a per-eigenstate moment inside a degenerate subspace and its basis-independent level sum. [Exports](qespace_cef/__init__.py#L15) · [Observable imports](qespace_cef/observables.py#L8) · [Migration table](NAMING.md#L84) · [Reciprocal definitions](qespace_cef/constants.py#L20) · [Per-state moment](qespace_cef/cefion.py#L344)

## Low severity

- **L1. A zero iteration limit raises an unrelated exception.** With `max_iterations=0`, the loop never assigns `change`, and formatting the failure raises `UnboundLocalError`. Define the valid range and check it before starting the iteration. [Mean-field iteration](qespace_cef/observables.py#L335)
- **L2. Zero-decimal summaries remove pure kets.** `summary(precision=0)` on axial Ce³⁺ prints empty `ev_0:` lines, including for a coefficient of exactly one. The strict `>` filter discards it at the threshold. Define the cutoff separately from output precision or include boundary coefficients. [Summary formatting](qespace_cef/cefion.py#L432)
- **L3. Convention coverage and descriptions have drifted.** The tests contain no test naming accepted CONV-0003 or CONV-0004 as required by the repository. Unit conversion has behavioral coverage already. A common-frame rotation check would exercise the latter. The architecture still directs readers to superseded ADR-0001, CONV-0008 describes an obsolete dependency on `numerical_zero` in CONV-0009, and the roadmap calls the tracked fitting example untracked. [Test rule](CLAUDE.md#L164) · [Unit tests](tests/test_constants.py#L27) · [Architecture](ARCHITECTURE.md#L144) · [Earlier susceptibility description](decisions/conventions/CONV-0008-degenerate-levels.md#L55) · [Current susceptibility convention](decisions/conventions/CONV-0009-susceptibility.md#L29) · [Roadmap](ROADMAP.md#L32)
- **L4. Version alone cannot identify the results of a development checkout.** The README installs from the moving `main` branch. Its metadata still says 1.3.0 while the unreleased changelog lists removed arguments, changed energy zeros, and a replacement neutron API. Record the commit with research results and give users a pinned release or identified prerelease when publishing the migration guide. [Install example](README.md#L14) · [Version](pyproject.toml#L3) · [Breaking changes](CHANGELOG.md#L14)

## Recorded naming choices

These spellings are marked as decided in the local inventory and can be carried into the proposed ADR. The last column gives an observable check for each migration. [Operator and ion names](NAMING.md#L44) · [Observable names](NAMING.md#L86)

| v1 spelling | Recorded v2 spelling | Check at the public entry point |
|---|---|---|
| `CEFion.Jval`, `CEFion.J` | `J`, `Jxyz` | `J` is the scalar quantum number. The three components of `Jxyz` are operators in a documented basis. |
| `J_plus`, `J_minus` | `J_p`, `J_m` | The function results agree with the named `Jp`, `Jm` attributes and obey the ladder-operator identities. |
| `Alpha`, `Beta`, `Gamma` | `alpha`, `beta`, `gamma` | The three values retain their ranks 2, 4 and 6 and their signs. |
| `boltzman_population`, `T` or `temperatures` | `boltzmann_population`, `temperature` | Normalize over eigenstates separately for each broadcast point. |
| `scattering_function`, `Ei` | `scattering_weights`, `E_incident` | Preserve the energy-transfer sign and the incident-energy cutoff when changing the spelling. |
| `precision` in mean field and in level printing | `conv_tol` for iteration, `precision` for decimal places | Changing display precision affects the text. Changing `conv_tol` affects the stopping criterion. |

The choice of package-wide or suite-wide scope (Q26), functions or methods (Q31), v1 coexistence (Q51), the field keyword (Q59), and the public representation of grouped levels are still recorded as open questions. [Plan](PLAN.md#L37) · [Open naming questions](NAMING.md#L133)

## Verification

| Check | Result |
|---|---|
| Unit tests in all five configured Python/NumPy pairs | **536 passed, 0 failed, 28 xfailed, 4 deselected per pair** |
| Package doctests in the same five pairs | **7 passed, 0 failed per pair** |
| Performance tests run separately on Python 3.14.2 / NumPy 2.5.5 | **4 passed** |
| README quick start on the same environment | Executed successfully, powder weights finite with shape `(3, 36)` and energies `(36,)` |
| Ruff 0.17.6, `check --no-cache .` | **7 diagnostics** in the fitting example (5 undefined-name occurrences, 1 unused import and 1 import-format diagnostic) |
| Ruff 0.17.6, `format --check --no-cache .` | **1 file** would be reformatted, the fitting example. 50 files already formatted |

The installed matrix was Python 3.12.12 with NumPy 1.26.0, 2.0.2 and 2.5.5, plus Python 3.13.11 and 3.14.2 with NumPy 2.5.5. Each run confirmed it imported this checkout. The 28 strict expected failures cover 16 reference-code sign differences, 9 thermodynamic identities, 2 g-tensor rotations and 1 unresolved cubic Γ3 display basis. The baseline tests include independent Stevens-operator comparisons, scattering sum rules, detailed balance, and powder/crystal projection checks. [Configured matrix](noxfile.py#L19) · [Sign comparisons](tests/test_cefmatrices.py#L160) · [Cubic basis case](tests/test_cefion.py#L374) · [Scattering checks](tests/test_observables.py#L341)

Tests used the existing environments with bytecode, pytest caches, disk capture and temporary-directory fixtures disabled. Fresh installation, a built wheel, Sphinx HTML, external literature and the experimental symbolic scripts remain unverified. The fitting example was checked statically and its broadcast dimensions were checked in memory. The optional optimizer was outside the executed checks.
