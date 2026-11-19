# Code audit of qespace_cef

One pass through the whole repository at commit `747c871` (version 1.3.0),
on 2026-11-19. It collects everything that needs addressing so the code does
not have to be read again; [`ARCHITECTURE.md`](ARCHITECTURE.md) describes
the code as it is. Line numbers refer to that commit. Remove an item, or
mark it done, when its change lands.

**Read:** every tracked file. The `symbolic/` scripts were only skimmed (they
stay experimental), and of `atoms-form-factors.txt` only the header and the Ce
lines.
**Ran:** `uv run pytest` (33 passed), plus 14 numerical checks against analytic
limits and exact formulas, in throwaway scripts outside the repository. Eight
of them show a bug (P1–P8). Six pass: the g-tensor of the `B20`-only case, the
cubic doublet and quartet, the high-T entropy limit, the direction of
magnetisation, and the powder and single-crystal sum rules within 1 %.

Severity: **bug** = wrong result, verified by running it; **bug (read)** =
wrong from reading, not run; **rule** = breaks a rule in CLAUDE.md;
**docs** = documentation does not match behaviour; **suggestion** = no rule
broken.

**Status, 2026-11-19.** The maintainer's decisions on this audit are in
section 7. Done: 6.1 (with P12), Q1, Q4, Q8, Q9 (`chore/cleanup`); X1–X3
(`build/numpy-only`); D1–D12, A13 (`docs/docstrings`); Q3, Q6, Q2 in text
(`style/quality`). Q7 in part (see its row). Every P item still in
the code carries an `# AUDIT Pn:` comment there and, except P15
(`constants`), a Warnings section in the docstring of the public object it
affects; the P items stay open until each fix is reviewed.

## 1. Physics and correctness

| ID | Severity | Where | Finding | Evidence | Proposal |
|---|---|---|---|---|---|
| P1 | bug | `observables.susceptibility`, `method='perturbation'` | The van Vleck term counts only `E_m > E_n` pairs (`jumps_with_positive_energy`). The docstring sums over every `m ≠ n`, and the `E_m < E_n` pairs contribute the same amount, so the van Vleck part is half what it should be | Ce, `B20=0.5, B40=0.01, B44=0.1`, H ∥ z. χ in μB/T: exact formula 0.13130 / 0.02981 / 0.00999 at 5 / 20 / 100 K, code (perturbation) 0.12704 / 0.02680 / 0.00847, code (magnetization) 0.13130 / 0.02981 / 0.00999 | Sum over both orderings, i.e. `(p_n − p_m)/(E_m − E_n)`; test it against `method='magnetization'` |
| P2 | bug | same, line 296 (`# TODO WTF is this projection?`) | The field direction enters as `Σ_α n_α \|J_α\|² / \|n\|` instead of `\|n·J\|²`. For n = (1,1,0) this gives √2 times the x value | Same ion: χ(110) / χ(100) = 1.414 with `perturbation`, 1.000 with `magnetization` (for this tetragonal case, rank-2 symmetry requires 1) | Project the operator: `Jn = Σ n_α J_α / \|n\|`, then `\|Jn_mn\|²` |
| P3 | bug | `observables.thermodynamics` | Energies in meV are mixed with temperatures in K: `F = −T ln Z` is in K, while `E` stays in meV, so `S = (E − F)/T` and `Cv = (E2 − E²)/T²` come out in mixed units | Entropy (k_B) at 10 / 20 / 40 K: exact 1.293 / 1.374 / 1.494, code 1.054 / 1.203 / 1.325. Cv (k_B): exact 0.170 / 0.108 / 0.239, code 0.00126 / 0.00080 / 0.00177 (low by a factor of 11.6², about 135) | Convert the energies to `k_B T` once and state the units in the return values; test against `−∂F/∂T` and against the high-T limit `ln(2J+1)` |
| P4 | bug, rule | `observables.gtensor` | The result depends on the arbitrary basis that the diagonalisation (`schur`, now `eigh`) picks inside the ground doublet. That goes against the CLAUDE.md section 3 row on degenerate levels (itself marked `OPEN:`): "compare projectors or sums over the subspace, never single eigenvectors" | Ce, `B20=1`: with no field, g = diag(2.571, 2.571, 0.857) with `schur`, which is correct, and diag(2.571, −2.571, 0.857) with `eigh`. A field of 1e-13 T along x (a different doublet basis) gives `[[-0.80, 0, 2.44], [0, 2.57, 0], [0.81, 0, 0.27]]` | Build g from basis-invariant quantities, e.g. `G = g gᵀ` from `Σ_ab <a\|J_α\|b><b\|J_β\|a>` over the doublet; the convention is a question (C7) |
| P5 | bug | `CEFion.diagonalize` | Calling it a second time rotates `Jx`, `Jy`, `Jz`, which are already in the eigenbasis, again | Ce, field (0.3, 0, 1) T: after a second `diagonalize()` the moments change by up to 2.58 μB | Keep the free-ion operators unchanged and store the rotated ones separately, or make diagonalisation run once (A2) |
| P6 | bug | `CEFion.diagonalize(sortWithE=False)` | `np.arange(self.Jval)` takes a float J and has the wrong length (2J+1 is needed) | `IndexError: arrays used as indices must be of integer (or boolean) type` | `np.arange(2*J+1)` |
| P7 | bug | `CEFion.get_degeneracies(precision)` | Line 183 overwrites the argument with `precision = 5`. `astype(int)` truncates instead of rounding, so levels at 0.999999 and 1.000001 come out as two levels | `get_degeneracies(precision=1)` returns 0.74073. With `eigh`, Pr³⁺ with `B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4` has its 28.584 meV doublet at 28.584000000000000 and 28.583999999999996 meV, which come out as two levels | Group levels by a relative tolerance from the config (Q3), not by decimals |
| P8 | bug | `Ion.__str__` | `%d` prints J = 5/2 as 2 | `str(Ion('Ce'))` gives `'Ce: J = 2, gJ = 0.86'` | Print the fraction |
| P9 | bug (read) | `observables.spectral_weight`, single crystal | `Σ_α (1 − Q̂_α²) \|J_α\|²` leaves out the cross terms `−Q̂_α Q̂_β J_α J_β*`. The correct form is `Σ_αβ (δ_αβ − Q̂_α Q̂_β) <i\|J_α\|f><f\|J_β\|i>` | Q ∥ (110): summed weight 5.6516 from the code, 5.6573 with the cross terms. That's small here; per transition it can be larger | Use the full projector; the commented-out "first implementation" (lines 132–135) went in that direction |
| P10 | bug (read) | `observables.susceptibility`, `method='magnetization'` | Returns `\|M\|/ε`, the magnitude, not the component along the field, and so loses the sign (and is wrong when M is not parallel to H). ε is `1e-5`, but the docstring says 10⁻⁸ | read | Return `M·n/ε` |
| P11 | bug (read) | `observables.magnetization_exchange` | No iteration limit, so it can loop forever when the mean field doesn't converge. M = 0 gives `nan` and silently ends the loop. `residual = 42` is a placeholder | read | Iteration cap, a relative tolerance from the config, raise on no convergence |
| P12 | bug (read) | `Ion._m_dynamic` | `gJ·J(J+1)` instead of `gJ·√(J(J+1))`; unused | read | Remove, or fix and test |
| P13 | data | `Ion` table, Ce | The `<j0>`/`<j2>` coefficients are Ce²⁺, because `composeTables` maps `ce` to `Ce2`, the only Ce entry in the ILL table. The other ions use the 3+ values | `atoms-form-factors.txt`: only `{F Ce2 ...}` lines | Record it as a known approximation next to the value, with its source (C9) |
| P14 | data, rule | `Ion` table | The values (J, gJ, α, β, γ, form factors) have no source next to them (CLAUDE.md 3.2.2); `composeTables` names "Czech guys" and "J.P. Brown" | read | Cite Stevens (1952) / Hutchings (1964) for α, β, γ and the ILL / ITC Vol. C tables (with DOI) for the form factors |
| P15 | rule | `constants.py` | No source or CODATA year. The values mix vintages (h, e from CODATA 1986; μB, k_B, N_A from 2010), `pi` is hard-coded to 11 digits, and `meV2K = 11.60452` and `K2meV = 0.086173` are separate numbers, not reciprocals. `R0`, `_gamn`, `C1`–`C3`, `NA`, `kB_SI` are unused | read | One CODATA set (2022), derived conversions, `np.pi`; drop what's unused (CLAUDE.md 3.2.1) |
| P16 | physics question | `cefmatrices`, `CEFpars` | Odd-rank operators (k = 1, 3, 5) are implemented and accepted as CEF parameters (`B52` in the `CEFpars` example, `B5m3` in `test_cefion.py`). Inside one J multiplet, parity removes odd k from the CEF Hamiltonian | read | Keep them as building blocks, but decide whether `CEFpars` accepts them (C3) |

## 2. Documentation against behaviour

| ID | Where | Finding |
|---|---|---|
| D1 | `CEFion` docstring | Says `ion: str` (it takes an `Ion`), `qespace_cef:ion`, `qespace_cef.reion.cfp` (old names), and "Examples: TODO" |
| D2 | `CEFion` docstring | Gives `H = ΣB O + g_J (H·J)`, with a + sign and no μB; the code has `−μB gJ J·H` |
| D3 | `CEFion.diagonalize` | Returns `E, U` from `schur`, unshifted and unsorted, not the stored energies; the docstring says `H = U E Uᵀ` (it is `U E U†`). "Raises ValueError if not real" can't happen for a Hermitian H |
| D4 | `spectral_weight` | The docstring gives `J⊥ = Σ(1 − Q_α/\|Q\|) J_α`, which is neither the code nor the physics. It promises sorted, flattened lists, but the function returns two (2J+1)² matrices; the sorting code is dead (lines 156–159). "List of Q vectors" (it takes one). A "RuntimeWarning for Q = 0" is never raised. The error message for a bad `scheme` speaks of "Q". The reference "Enderle book following Stephane Raymond article" has no DOI |
| D5 | `spectral_weight` | Debye–Waller: the docstring says `W = ½ B Q²`, the code computes `exp(−⅓ Q² B_iso)`. Neither matches the crystallographic `B_iso` (`exp(−B Q²/16π²)`). `Ei = 1e6` meV means that k_f/k_i ≈ 1 by default |
| D6 | `susceptibility` | Says the units are μB²; they are μB²/meV per ion (×`C.uB` gives μB/T). "Old implementation" is mentioned without saying which one. Typos throughout |
| D7 | `thermodynamics` | `Cv = (<E²> − <E>)/T²` is missing the square; no units given |
| D8 | `CEFpars` | The docstring says both that "symmetry handling was removed" and that it "simplifies creation from the point group". The example prints `symmetry=('2', 'C2', ...)` for the default C1. The `symmetry` attribute is typed `str` but is a NumPy record |
| D9 | `StevensBase` | The docstring shows `StevensBase(J)['O32']`; the key is `'32'`. The type hints say complex for `J_z` (it's float) |
| D10 | `point_group_synonyms` | Documented as returning a tuple; it returns a structured array, several rows for a lattice name |
| D11 | `Ion.mff` | No unit for Q (Å⁻¹ assumed); "List of Q" but scalars are used |
| D12 | `README.md` | The quick start and the API table are correct for 1.3.0 (checked by running them), but they say "magnetic form factor (`mff`)" without the dipole approximation or the Ce²⁺ caveat |

## 3. Architecture

| ID | Severity | Finding | Proposal |
|---|---|---|---|
| A1 | rule | `ARCHITECTURE.md` (now added), `DEV.md` and `CONVENTIONS.md` are required by CLAUDE.md (intro, section 3, section 7), and only the first exists | Write `DEV.md` (commands, release steps from the 1.3.0 tag) and `CONVENTIONS.md` once C1–C10 are decided |
| A2 | suggestion | `CEFion` is a stateful object: the constructor builds H, `diagonalize()` mutates `Jx/Jy/Jz` and changes their meaning (free-ion basis before, eigenbasis after; a list before, an array after), and `energies` is `[]` until diagonalised. This causes P5 and makes every observable depend on the call order | Make the result immutable: `CEFion` holds `ion`, `cfp`, field and H; a `solve()` (or the constructor) returns energies, eigenvectors and the operators in the eigenbasis as read-only arrays; observables are pure functions of that |
| A3 | suggestion | Observables in a field rebuild a whole `CEFion` from `ion` and `cfp` and ignore `cefion.Hfield`, so the field lives in two places | Pass the field to the solver only; the observables take a solved state |
| A4 | rule | `from .observables import *` exports `np`, `C`, `Tuple`, `CEFion` again, and the submodules. With no `__all__`, the public API is undefined (CLAUDE.md 6.2 "public entry point") | `__all__` in `__init__` and in every module |
| A5 | suggestion | The J matrices are built twice (`J_x`… functions and `StevensBase.__init__`); `symbolic/cefmatrices.py` is a third, sympy, copy | `StevensBase` uses the module functions |
| A6 | suggestion | `StevensBase` exposes 49 hand-written methods found through `cls.__dict__`; `implementedOps()` lists every name starting with `O`. It works and is tested, but adding k = 7… or checking the normalisation means editing 49 formulas | Keep it as the reference for now; later, generate `O_k^q` from one closed formula (e.g. via the Racah tensor operators) and keep the hand-written table only as a test fixture |
| A7 | suggestion | `CEFpars.from_symmetry` is broken: an instance method marked `# @classmethod`; `Bvalues` and `units` are ignored; it fills zeros; for `T`/`Th` it writes past the end. `_assignLattice` duplicates `point_group_synonyms`, and the constructor's symmetry branch is a `pass` | Decide what point-group support should be (C8), then rebuild it on one table with tests, or delete it |
| A8 | suggestion | `data_containers.py` is a mix: point groups (used), `composeTables` (a one-off generator), `sympy_expression_mapping` (unused), a commented-out `CEF_tables` class. The ion data sit in a dictionary literal rebuilt on every `Ion()` | Put the data tables (ions, point groups, allowed parameters) in one data module with their sources; move `composeTables` and `atoms-form-factors.txt` out of the package (e.g. `tools/`) |
| A9 | suggestion | Files that ship without a use: `thermodynamics-magnetism.png`, `qespace_cef/README.md`, `atoms-form-factors.txt` (kept for now, by your decision) | Move the reference to the README or the docs; keep the rest out of the wheel when the data question is settled |
| A10 | rule | Validation: `assert` in `gtensor` (removed by `python -O`), a bare `Exception` for a non-Hermitian H, `NameError` for an unknown point group, `KeyError` for an unknown ion; the Hermiticity check is absolute `allclose` (CLAUDE.md 5.5, 5.6) | `ValueError` at the public boundary, relative tolerances from one config place |
| A11 | rule | No central tolerance config: `numerical_zero = 1e-12` (meV, absolute), `Eimag_precision=10`, `precision=5`, a print threshold of `10⁻³`, `eps = 1e-5` T, `precision=1e-5` (CLAUDE.md 5.6) | One `config` with relative tolerances, read at call time |
| A12 | rule | A pure computation prints: `susceptibility` has `print(jumps)` at line 298 (CLAUDE.md 5.8) | Remove it |
| A13 | rule | Units missing in the API: `mff` Q, the `spectral_weight` returns, `lam` in `magnetization_exchange`, the `thermodynamics` returns, `boltzman_population` T (CLAUDE.md 3.3) | State them in every docstring |

## 4. Code quality (for the reformat)

| ID | Finding |
|---|---|
| Q1 | Tabs and spaces mixed: `cefmatrices.py` (tabs inside classes, spaces at module level and on lines 301–321), `data_containers.py` line 115 |
| Q2 | Naming: camelCase and snake_case mixed (`sortWithE`, `shiftToZero`, `pointGroup`, `Hfield`, `mffCoefficients`, `implementedOps`, `allowed_Bpars`), class attributes `Alpha/Beta/Gamma`, typos in public names (`boltzman_population`) and in text (`faiclitates`, `Ptobabitlity`, `ocntaining`, `tempreature`, `magneitzation`, `calulated`, `strenght`, `identifyer`) |
| Q3 | Type hints: `typing.Tuple`/`List`/`Union` instead of built-ins; a mutable default `Hfield=[0,0,0]`; wrong types (`J2p1: float` is int, `symmetry: str`) |
| Q4 | Docstrings mix NumPy style, Google style ("Parameters:"), Markdown (`## Class`) and RST `\|` blocks, and lack array shapes (CLAUDE.md 5.4) |
| Q5 | Unused imports: `conj`, `transpose`, `zeros` in `cefmatrices`; `np` and `Union` used only for small things; commented imports in `cefion` lines 10–11 and `observables` line 2 |
| Q6 | NumPy scalars leak into results and strings: `get_degeneracies` returns `np.float64`/`np.int64` (CLAUDE.md 5.7.3) |
| Q7 | `CEFion.__str__(evecs, precision)` takes arguments `str()` can't pass; `{energy:.{precision}}` gives significant digits, not decimals; the kets are created only for J % 1 in {0, 0.5}, so other values leave the attribute undefined. **Partly done:** energies are printed with decimals; `CEFion.summary(evecs, precision)` exists as a stub that raises `NotImplementedError`, to be implemented with A2; the kets remain (no `Ion` gives another J) |
| Q8 | Tests: `test_cefion.py` and `test_cefpars.py` contain only `__main__` scripts; `test_cefmatrices.py` keeps a verbatim copy of PyCrystalField's operator code (GPL-3, compatible, but without attribution or version) and an old implementation as reference; `test_implementation` asserts nothing; `tests/__init__.py` exists (the suite's DEV.md asks for none) |
| Q9 | No linter or formatter configured (CLAUDE.md 6.7 `OPEN:` ruff) |

## 5. Dependencies

| ID | Now | Needed by | Proposal |
|---|---|---|---|
| X1 | `matplotlib>=3.10.8` | nothing | Remove |
| X2 | `sympy>=1.14.0` | `symbolic/` only, not imported | Remove from runtime; an optional extra `symbolic` if the scripts should stay runnable |
| X3 | `scipy>=1.17.1` | `scipy.linalg.schur` in `CEFion.diagonalize`; the comment says "other functions do not produce perpendicular subspaces for degenerate eigenvalues" | Replace with `numpy.linalg.eigh`, which returns orthonormal eigenvectors for a Hermitian matrix, degenerate ones included. That removes scipy and leaves `numpy` alone, as CLAUDE.md 5.1 asks (and settles its `OPEN:` on scipy) |
| X4 | `numpy>=2.4.4` | everything | The floor is the latest release, not what's needed, and it excludes numpy 1, which CLAUDE.md 5.7 says must work. Set the real floor (e.g. `>=1.26`) and test numpy 1 and 2 (`OPEN:` nox) |
| X5 | `requires-python >=3.12` | – | Fine; keep it unless numpy 1 support needs 3.11 |

After X1–X4 the runtime needs only `numpy`.

## 6. Comment inventory

Every comment that carries information, grouped by what it asks for.

### 6.1 Old implementations and dead code: remove, or keep as test fixtures

| Where | What |
|---|---|
| `cefion.py` 10–11 | commented imports |
| `cefion.py` 91–92 | "Might take cc of everything ... `(H + H.T.conj())/2`" |
| `cefion.py` 97–101 | a commented-out warning for a non-symmetric H |
| `observables.py` 2, 186, 284, 292 | a commented import, an old `M = gJ·abs(...)`, `fill(nan)`, `jumps_with_negative_energy` ("obsolete but left for clarity") |
| `observables.py` 132–135 | "First implementation does not seem to work well" (the projector in `spectral_weight`, see P9) |
| `observables.py` 156–159 | dead sorting and the alternative `return` |
| `observables.py` 332–336 | "My guess" g-tensor formulas next to "PyCrystalField" |
| `data_containers.py` 4–5, 93–104, 107–178 | `sympy_expression_mapping`, the `CEF_tables` stub, `composeTables` and its `__main__` |
| `cefpars.py` 59–63, 66–102 | a symmetry `pass`, `from_symmetry` (A7) |
| `ion.py` 88–98 | `_LandeGFactor` ("maybe for the future"), `_m_in_vacuum`, `_m_dynamic` (P12), all unused |
| `__init__.py` 16 | `# from .cefmatrices import *` |
| `tests/test_cefmatrices.py` | the old `O20…O6m6` and `StevensOpPCF`: keep as reference fixtures, with their source |

### 6.2 Open questions and ideas: issues or decisions

| Where | Text (shortened) | Goes to |
|---|---|---|
| `cefmatrices.py` 4 | the matrices and conventions "need to be cross-checked against some reliable book" | C2, C4 |
| `cefion.py` 140 | `eigenvectors = U.T ### TRIPLE CHECK` | verified by reading: rows are eigenvectors, and `V* J Vᵀ` gives `<i\|J\|j>`. Replace the comment with a test |
| `cefion.py` 196 | "After refactoring the precision there might be cleanup needed" | P7, Q7 |
| `cefion.py` 32 | "Examples: TODO" | D1 |
| `observables.py` 296 | "TODO WTF is this projection?" | P2 |
| `cefpars.py` 80–86 | the cubic prefactors "from Martin's table have an additional factor 1/2, which ... does not reproduce the doublet–quartet scheme. Hopefully it's an implementation problem only" | C3, C8; the code (B44 = 5 B40) does give doublet + quartet (verified) |
| `ion.py` 89 | initialise from S, L instead of J | C5 (intermediate coupling) |
| `symbolic/README.md` | only the eigenvalues are fast; reduce parameters by setting small ones to zero; physical bounds such as "B20 > 0" | stays experimental |

### 6.3 Definitions that belong in conventions

These comments and docstrings define something every result depends on.
Each becomes a row in CLAUDE.md section 3 and a CONV, with a test.

| ID | Definition in the code | Where | CLAUDE.md row |
|---|---|---|---|
| C1 | Energies in meV, shifted so the ground state is at 0; temperatures in K; fields in T; Q in Å⁻¹ | `cefion` 146, `ion.mff` | Energy unit |
| C2 | Basis order `\|J>, ..., \|-J>` (the first component is m = +J); free-ion kets are labelled in that order | `cefmatrices` 6, 40; `cefion` 77 | Basis |
| C3 | Stevens `O_k^q` (McPhase, EasySpin), q < 0 as sine-type operators; which k are allowed as parameters (P16); the cubic relations B44 = 5 B40, B64 = −21 B60; the factor-½ remark | `cefmatrices`, `cefpars` | CEF parameters, normalisation |
| C4 | PyCrystalField agrees only up to a sign for some operators (the test accepts ±) | `tests/test_cefmatrices.py` | Normalisation, reference codes |
| C5 | Ground J multiplet only; `gJ` Landé; no intermediate coupling | `ion`, `cefion` | Basis |
| C6 | Zeeman `H = −μB gJ J·H`, moment `+gJ<J>`; the field in the frame of the J operators | `cefion` 84, 176 | Zeeman term, quantisation axis |
| C7 | g-tensor: definition from the ground doublet (PyCrystalField's convention, `2gJ<J>`), and what to do when the ground state is not a doublet | `gtensor` | new row |
| C8 | Point groups: Hermann–Mauguin ↔ Schoenflies table; the allowed parameters per group from the McPhase manual (node133); which frame (z along the main axis) | `data_containers`, `cefpars` | Quantisation axis |
| C9 | Form factor: dipole approximation `j0 + (2 − gJ)/gJ · j2`, the ILL coefficients, Ce²⁺ for Ce | `ion.mff` | Cross sections |
| C10 | Neutron cross section: `N(γr0)² = 1` (arbitrary units), the Debye–Waller form (D5), `k_f/k_i = √(1 − ΔE/Ei)`, positive ΔE = neutron energy loss, powder average `⅔ Σ_α \|J_α\|²` | `spectral_weight` | Cross sections |
| C11 | Susceptibility units (μB²/meV per ion, or μB/T, or emu/mol); the conversions in `constants` (`C1`–`C3`) | `susceptibility`, `constants` | new row |
| C12 | Degenerate levels: tolerance and grouping (P7); eigenvectors inside a level are arbitrary (P4) | `cefion` | Degenerate levels (the rule is written, but the row is marked `OPEN:`; the code doesn't follow it) |

## 7. Work packages

The maintainer's decisions on this audit (2026-11-19) replace the first
proposal: the physics fixes (P) wait until the maintainer has reviewed each
one, and the architecture items (A) are decided later. Each package is one
branch (CLAUDE.md 1.3).

1. **`chore/cleanup`** (merged): ruff with a light rule set and the format
   (Q9, Q1); dead code removed (6.1, with P12); tests of the current
   behaviour of every module, the reference implementations moved to
   `tests/reference/` with attribution (Q8, Q15); performance stubs.
2. **`build/numpy-only`** (merged): `eigh` for `schur`, compared with the
   old results on 7 ions and parameter sets; scipy and matplotlib dropped,
   sympy in the dependency group `symbolic` (X1–X3, Q14).
3. **`docs/docstrings`**: D1–D12 and A13 in NumPy style (Q4); markers and
   Warnings sections for the P items.
4. **`style/quality`**: Q3 type hints and the mutable default, Q6
   plain floats, Q7 `__str__` (energies with decimals, a `summary()` stub),
   Q2 typos in text only; renaming public names waits for A2.
5. Delayed: the P fixes, the A items except A13, the items of 6.2 and 6.3
   (C1–C12), X4 (the numpy floor, once tests show what it can be), and the
   questions Q11–Q13.

## 8. Questions for the maintainer

- **Q10.** Is the order of the work packages in section 7 right? In
  particular, should the numpy-only change come before or after the fixes?
  **Answer:** numpy-only first; the fixes wait for review (section 7).
- **Q11 (C1–C12).** Which conventions do you want to decide first? P4 (the
  g-tensor) and P9 (the cross section) wait for C7 and C10.
- **Q12 (P16).** Should `CEFpars` accept odd-rank parameters?
- **Q13 (A2).** Is an immutable solved-state object (section 3, A2) the
  direction you want for the API, or should `CEFion` stay as it is with only
  the bugs fixed?
- **Q14 (X2).** Should the symbolic scripts keep `sympy` as an optional
  extra, or have no declared dependency at all?
  **Answer:** a dependency group `symbolic`.
- **Q15 (Q8).** Keep the verbatim PyCrystalField code in the tests (with
  attribution), or replace it with stored reference matrices?
  **Answer:** keep it, in `tests/reference/`, with attribution.
