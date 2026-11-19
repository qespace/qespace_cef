# Code audit of qespace_cef

The open findings of the audit of 2026-11-19 (commit `747c871`, version
1.3.0) and of the review of the work packages that followed, as of `main` at
`610c884` (2026-11-20). What has been fixed is gone from this file: it is in
[`CHANGELOG.md`](CHANGELOG.md) (Unreleased) and in the history of the
branches `chore/cleanup`, `build/numpy-only`, `docs/docstrings` and
`style/quality`. [`ARCHITECTURE.md`](ARCHITECTURE.md) describes the code as
it is. Remove an item when its change lands; IDs are never reused.

Places are named by module and function. Every P item still in the code has
an `# AUDIT Pn:` comment where it happens and, except P15, a Warnings section
in the docstring of the public object it affects.

Severity: **bug** = wrong result, verified by running it; **bug (read)** =
wrong from reading, not run; **rule** = breaks a rule in CLAUDE.md;
**docs** = documentation does not match behaviour; **suggestion** = no rule
broken. The evidence was rerun on `610c884` where it says so; the rest is
from `747c871` and does not depend on what changed since.

## 1. Physics and correctness

| ID | Severity | Where | Finding | Evidence | Proposal |
|---|---|---|---|---|---|
| P1 | bug | `observables.susceptibility`, `method='perturbation'` | The van Vleck term counts only the pairs with `E_m > E_n`. The `E_m < E_n` pairs contribute the same amount, so the van Vleck part is half what it should be | Ce, `B20=0.5, B40=0.01, B44=0.1`, H ∥ z. χ in μB/T: exact formula 0.13130 / 0.02981 / 0.00999 at 5 / 20 / 100 K, code (perturbation, × `C.uB`) 0.12704 / 0.02680 / 0.00847, code (magnetization) 0.13130 / 0.02981 / 0.00999 | Sum over both orderings, i.e. `(p_n − p_m)/(E_m − E_n)`; test it against `method='magnetization'` |
| P2 | bug | same | The field direction enters as `Σ_α n_α \|J_α\|² / \|n\|` instead of `\|n̂·J\|²`. For n = (1,1,0) this gives √2 times the x value | Same ion, 20 K: χ(110) / χ(100) = 1.414 with `perturbation`, 1.000 with `magnetization` (for this tetragonal case, rank-2 symmetry requires 1). Rerun on `610c884` | Project the operator: `Jn = Σ n_α J_α / \|n\|`, then `\|Jn_mn\|²` |
| P3 | bug | `observables.thermodynamics` | Energies in meV are mixed with temperatures in K: `F = −T ln Z` is in K, while `E` stays in meV, so `S = (E − F)/T` and `Cv` come out in mixed units | Entropy (k_B) at 10 / 20 / 40 K: exact 1.293 / 1.374 / 1.494, code 1.054 / 1.203 / 1.325. Cv (k_B): exact 0.170 / 0.108 / 0.239, code 0.00126 / 0.00080 / 0.00177 (low by a factor of 11.6², about 135) | Convert the energies to `k_B T` once and state the units of the returns (C1, C11); test against `−∂F/∂T` and the high-T limit `ln(2J+1)` |
| P4 | bug, rule | `observables.gtensor` | Even with a correct formula, g itself depends on the arbitrary basis inside the ground doublet: a change of basis turns g into g R with R a rotation, and only `G = g gᵀ` (the principal values and axes) is basis independent. The function returns g in whatever basis `eigh` picks. That goes against the CLAUDE.md section 3 row on degenerate levels (itself marked `OPEN:`) | Ce, `B20=1`, no field: diag(2.571, 2.571, 0.857) with `schur`, diag(2.571, −2.571, 0.857) with `eigh` | Return G, or the principal values and axes; together with P17 and C7 |
| P5 | bug | `CEFion.diagonalize` | Calling it a second time rotates `Jx`, `Jy`, `Jz`, which are already in the eigenbasis, again | Ce, `B20=0.5, B40=0.01, B44=0.1`, field (0.3, 0, 1) T: after a second `diagonalize()` the moments change by up to 2.94 μB. Rerun on `610c884` | Rotate from the free-ion operators every time (keep them), or make diagonalisation run once (A2) |
| P6 | bug | `CEFion.diagonalize(sortWithE=False)` | `np.arange(self.Jval)` takes a float J and has the wrong length (2J+1 is needed) | `IndexError: arrays used as indices must be of integer (or boolean) type`. Rerun on `610c884` | `np.arange(2*J+1)` |
| P7 | bug | `CEFion.get_degeneracies(precision)` | The argument is overwritten with `precision = 5`, and `astype(int)` truncates instead of rounding, so two degenerate states whose energies differ by rounding noise can come out as two levels | `get_degeneracies(precision=1)` still gives 5 decimals (0.78397). Pr³⁺ with `B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4`: the 28.584 meV doublet comes out of `eigh` as 28.584000000000000 and 28.583999999999996 meV, two levels. Rerun on `610c884` | Group levels by a relative tolerance from one config place (A11), not by decimals (Q24) |
| P8 | bug | `Ion.__str__` | `%d` prints J = 5/2 as 2 | `str(Ion('Ce'))` gives `'Ce: J = 2, gJ = 0.86'`. Rerun on `610c884` | Print the fraction |
| P9 | bug (read) | `observables.spectral_weight`, `'single-crystal'` | `Σ_α (1 − Q̂_α²) \|J_α\|²` leaves out the cross terms `−Q̂_α Q̂_β J_α J_β*` | Q ∥ (110): summed weight 5.6516 from the code, 5.6573 with the cross terms. That's small here; per transition it can be larger | Use the full projector `Σ_αβ (δ_αβ − Q̂_α Q̂_β) <i\|J_α\|f><f\|J_β\|i>` (C10) |
| P10 | bug | `observables.susceptibility`, `method='magnetization'` | Returns `\|M\|/ε`, the length, not the component along the field: wrong when M is not parallel to H, and it loses the sign | Ce, `B20=0.5, B40=0.01, B44=0.1`, 20 K, field along (1,0,1): the code gives 0.05824 μB/T, `M·n̂/ε` is 0.05329 μB/T. Run on `610c884` | Return `M·n̂/ε` |
| P11 | bug (read) | `observables.magnetization_exchange` | No iteration limit, so it can loop forever when the mean field doesn't converge. M = 0 gives `nan` and silently ends the loop. `residual = 42` is a placeholder | read | Iteration cap, a relative tolerance from the config, raise on no convergence |
| P13 | data | `Ion` table, Ce | The `<j0>`/`<j2>` coefficients for Ce are those of Ce²⁺, the only Ce entry in the ILL table (`atoms-form-factors.txt`); the generator that built the table (`composeTables`, in the history at `747c871`) mapped `ce` to `Ce2`. The other ions use the 3+ values | `atoms-form-factors.txt`: only `{F Ce2 ...}` lines | Record it as a known approximation next to the value, with its source (C9) |
| P14 | data, rule | `Ion` table | The values (J, gJ, α, β, γ, form factors) have no source next to them (CLAUDE.md 3.2.2); the generator named "Czech guys" and "J.P. Brown" | read | Cite Stevens (1952) / Hutchings (1964) for α, β, γ and the ILL / ITC Vol. C tables (with DOI) for the form factors |
| P15 | rule | `constants.py` | No source or CODATA year. The values mix vintages (h, e from CODATA 1986; μB, k_B, N_A from 2010), `pi` is hard-coded to 11 digits, and `meV2K = 11.60452` and `K2meV = 0.086173` are separate numbers, not reciprocals. `R0`, `_gamn`, `C1`–`C3`, `NA`, `kB_SI` are unused | read | One CODATA set (2022), derived conversions, `np.pi`; drop what's unused (CLAUDE.md 3.2.1) |
| P16 | physics question | `cefmatrices`, `CEFpars` | Odd-rank operators (k = 1, 3, 5) are implemented and accepted as CEF parameters. Inside one J multiplet, parity removes odd k from the CEF Hamiltonian | `CEFpars(dict(B52=0.1, B3m1=0.2))` is accepted. Rerun on `610c884` | Keep them as building blocks, but decide whether `CEFpars` accepts them (C3, Q12) |
| P17 | bug | `observables.gtensor` | The z column uses `\|<0\|J_α\|0>\|`. The `abs()` makes even `G = g gᵀ` depend on the basis inside the doublet, so the principal g values are wrong in a general basis. (From the maintainer's notes: "I feel implementation of `gtensor` is wrong") | Ce, 5 random bases of the ground doublet: the spread of g gᵀ is 2.93 (`B20=1`) and 2.15 (low symmetry) with the code's formula, 4e-15 with `g_αβ = gJ Tr(P J_α σ_β)`. For `B20=1` the principal values of the latter are 0.857, 2.571, 2.571 (correct). Run on `610c884` | `g_αβ = gJ Tr(P J_α σ_β)` over the doublet, reported as G or principal values (P4); the convention with C7; non-Kramers doublets need their own treatment |

## 2. Documentation against behaviour

| ID | Where | Finding | Proposal |
|---|---|---|---|
| D13 | `StevensBase` docstring, Notes, "Observations" | Observations 1 and 2 say O(n, ±(n−1)) = ½{O(n, ±n), Jz}, which cannot hold: the two sides change m by n−1 and by n. The form that holds, and that the code comments use (`O32 = 0.5 (O22 Jz + Jz O22)`), is O(n, ±(n−1)) = ½{O(n−1, ±(n−1)), Jz}: checked for J = 1/2 … 8, n = 2 … 6, both signs, 160 cases, on `610c884`. Observation 4 (q < 0: purely imaginary, not symmetric) holds; 3 is too general to check as written. (From the maintainer's notes: "validate observations ... maybe can go to tests?") | Correct 1 and 2; turn 1, 2 and 4 into tests in `tests/test_cefmatrices.py`: an independent check of the hand-written table, also for A6 |

## 3. Architecture

| ID | Severity | Finding | Proposal |
|---|---|---|---|
| A1 | rule | `DEV.md` and `CONVENTIONS.md` are required by CLAUDE.md (intro, section 3, section 7) and don't exist | Write `DEV.md` (commands, release steps from the 1.3.0 tag) and `CONVENTIONS.md` once C1–C12 are decided |
| A2 | suggestion | `CEFion` is a stateful object: the constructor builds H, `diagonalize()` mutates `Jx/Jy/Jz` and changes their meaning (free-ion basis before, eigenbasis after; a list before, an array after), and `energies` is `[]` until diagonalised. This causes P5 and makes every observable depend on the call order | Make the result immutable: `CEFion` holds `ion`, `cfp`, field and H; a `solve()` (or the constructor) returns energies, eigenvectors and the operators in the eigenbasis as read-only arrays; observables are pure functions of that (Q13) |
| A3 | suggestion | Observables in a field rebuild a whole `CEFion` from `ion` and `cfp` and ignore `cefion.Hfield`, so the field lives in two places | Pass the field to the solver only; the observables take a solved state |
| A4 | rule | `from .observables import *` exports `np`, `C`, `ArrayLike`, `CEFion` again, and the submodules. With no `__all__`, the public API is undefined (CLAUDE.md 6.2 "public entry point") | `__all__` in `__init__` and in every module |
| A5 | suggestion | The J matrices are built twice (`J_x`… functions and `StevensBase.__init__`); `symbolic/cefmatrices.py` is a third, sympy, copy | `StevensBase` uses the module functions |
| A6 | suggestion | `StevensBase` has 49 hand-written methods. It works and is tested, but adding k = 7… or checking the normalisation means editing 49 formulas | Keep it as the reference for now; later, generate `O_k^q` from one closed formula (e.g. via the Racah tensor operators) and keep the hand-written table only as a test fixture |
| A7 | suggestion | `CEFpars` stores the point group but doesn't use it; `allowed_Bpars` is a lookup table without the relations between parameters (cubic B44 = 5 B40) | Decide what point-group support should be (C8), then build it on the one table, with tests |
| A8 | suggestion | The ion data sit in a dictionary literal in `Ion.__init__`, rebuilt on every `Ion()`; `data_containers.py` holds only the point groups | Put the data tables (ions, point groups, allowed parameters) in one data module with their sources; move `atoms-form-factors.txt` out of the package (e.g. `tools/`) |
| A9 | suggestion | Files that ship without a use: `thermodynamics-magnetism.png`, `qespace_cef/README.md`, `atoms-form-factors.txt` (kept for now, by your decision) | Move the reference to the README or the docs; keep the rest out of the wheel when the data question is settled |
| A10 | rule | Validation: `assert` in `gtensor` (removed by `python -O`), a bare `Exception` for a non-Hermitian H in `CEFion`, `NameError` for an unknown point group, `KeyError` for an unknown ion; the Hermiticity check in the constructor is an absolute `allclose` (CLAUDE.md 5.5, 5.6) | `ValueError` at the public boundary, relative tolerances from one config place |
| A11 | rule | No central tolerance config: `numerical_zero = 1e-12` (meV, absolute), `Eimag_precision=10`, `precision=5`, a print threshold of `10⁻³`, `eps = 1e-5` T, `precision=1e-5` (CLAUDE.md 5.6) | One `config` with relative tolerances, read at call time |
| A12 | rule | `susceptibility(method='perturbation')` prints the matrix of transition energies (CLAUDE.md 5.8) | Remove the print (Q25) |
| A14 | suggestion | Stevens operators can be had two ways, `StevensBase(J).O4m3()` and `StevensBase(J)['4m3']`, and `implementedOps` finds them through `cls.__dict__`: two access paths for one thing, with string keys to parse. (From the maintainer's notes: "Should we settle on just one of them?") | Settle on one with A2 and A6, e.g. `stevens(J, k, q)` or `ops[k, q]` with integers; the methods go when A6 generates the operators. Breaking (v2) |

## 4. Code quality

| ID | Finding | Proposal |
|---|---|---|
| Q2 | Naming: camelCase and snake_case mixed (`sortWithE`, `shiftToZero`, `pointGroup`, `Hfield`, `mffCoefficients`, `implementedOps`, `allowed_Bpars`), class attributes `Alpha/Beta/Gamma`, and the typo in the public name `boltzman_population` | Rename with A2, in v2 |
| Q7 | `CEFion.summary(evecs, precision)` is a stub that raises `NotImplementedError`; `__str__` still takes `evecs` and `precision`, which `str()` can't pass | Implement `summary` with A2 and give `__str__` no arguments |

## 5. Dependencies

The runtime needs only `numpy`.

| ID | Now | Proposal |
|---|---|---|
| X4 | `numpy>=2.4.4`: the floor is the latest release, not what's needed, and it excludes numpy 1, which CLAUDE.md 5.7 says must work | Find the real floor (e.g. `>=1.26`) by running the tests on numpy 1 and 2 (Q23), then set it |

## 6. Comments and conventions

### 6.1 Comments that ask for something

| Where | Text (shortened) | Goes to |
|---|---|---|
| `cefmatrices` module docstring | the matrices and conventions "need to be cross-checked against some reliable book" | C2, C4 |
| `CEFion.__str__` | "After refactoring the precision there might be cleanup needed" | P7, Q7 |

### 6.2 Definitions that belong in conventions

These define something every result depends on. Each becomes a row in
CLAUDE.md section 3 and a CONV, with a test. The ones another qespace package
also relies on (units, the cross-section conventions) are suite conventions
and go to the workspace first.

| ID | Definition in the code | Where | CLAUDE.md row |
|---|---|---|---|
| C1 | Energies in meV, shifted so the ground state is at 0; temperatures in K; fields in T; Q in Å⁻¹ | `CEFion.diagonalize`, `Ion.mff` | Energy unit |
| C2 | Basis order `\|J>, ..., \|-J>` (the first component is m = +J); free-ion kets are labelled in that order | `cefmatrices`, `CEFion.freeionkets` | Basis |
| C3 | Stevens `O_k^q` (McPhase, EasySpin), q < 0 as sine-type operators; which k are allowed as parameters (P16); the cubic relations B44 = 5 B40, B64 = −21 B60. A removed comment said "Martin's table" has an extra factor ½ that does not reproduce the doublet–quartet scheme; the code's B44 = 5 B40 does (verified) | `cefmatrices`, `CEFpars` | CEF parameters, normalisation |
| C4 | PyCrystalField agrees only up to a sign for some operators (the test accepts ±) | `tests/test_cefmatrices.py` | Normalisation, reference codes |
| C5 | Ground J multiplet only; `gJ` Landé; no intermediate coupling (a removed comment wanted to initialise from S, L) | `Ion`, `CEFion` | Basis |
| C6 | Zeeman `H = −μB gJ J·H`, moment `+gJ<J>`; the field in the frame of the J operators | `CEFion` Hamiltonian, `magnetic_moment` | Zeeman term, quantisation axis |
| C7 | g-tensor: its definition from the ground doublet (P4, P17), what is returned (g, G or principal values and axes), and what to do when the ground state is not a Kramers doublet | `gtensor` | new row |
| C8 | Point groups: Hermann–Mauguin ↔ Schoenflies table; the allowed parameters per group from the McPhase manual (node133); which frame (z along the main axis) | `data_containers`, `CEFpars` | Quantisation axis |
| C9 | Form factor: dipole approximation `j0 + (2 − gJ)/gJ · j2`, the ILL coefficients, Ce²⁺ for Ce | `Ion.mff` | Cross sections |
| C10 | Neutron cross section: `N(γr0)² = 1` (arbitrary units), the Debye–Waller form (`exp(−B_iso Q²/3)` in the code, not the crystallographic one), `k_f/k_i = √(1 − ΔE/Ei)`, positive ΔE = neutron energy loss, powder average `⅔ Σ_α \|J_α\|²`, the single-crystal projector (P9). The orientation of the transition matrices: `jumps[i, j] = E_j − E_i` and `Sqw[i, j]` for the transition from initial i to final j, rows = initial states (from the maintainer's notes: "this should be a matter of convention, documented somewhere") | `spectral_weight` | Cross sections |
| C11 | Susceptibility units (μB²/meV per ion from `perturbation`, μB/T from `magnetization`, or emu/mol); the conversions in `constants` (`C1`–`C3`) | `susceptibility`, `constants` | new row |
| C12 | Degenerate levels: tolerance and grouping (P7); eigenvectors inside a level are arbitrary (P4) | `CEFion` | Degenerate levels (the rule is written, but the row is marked `OPEN:`; the code doesn't follow it) |

## 7. Next work packages

Done so far: `chore/cleanup`, `build/numpy-only`, `docs/docstrings`,
`style/quality`. The order of what follows: tests first, so every fix is
driven by a failing test (your answer to Q19); then the fixes that need no
convention; then the conventions; then the fixes that depend on them; the
API redesign last, as version 2.0. Each package is one branch (CLAUDE.md
1.3) and stops for your review after it, or after about 5 commits. Each P
fix is its own commit with a physics explanation, so you can review them one
by one.

1. **`test/known-problems`.** A strict `xfail` test for each P item with an
   expected result that does not come from the code under test: P1 (the two
   susceptibility methods agree, in μB/T), P2 (χ(110) = χ(100) for a
   tetragonal case), P3 (`S = −∂F/∂T`, `Cv = T ∂S/∂T`, `S → ln(2J+1)` at high
   T, in k_B), P5 (a second `diagonalize()` changes nothing), P6
   (`sortWithE=False` runs), P7 (the Pr cubic doublet is one level), P8
   (`str` shows 5/2), P9 (single-crystal weights against the full projector,
   computed in the test), P10 (field along (1,0,1)), P11 (M = 0 gives 0, not
   NaN), P17 (g gᵀ unchanged under a rotation of the doublet basis). Tests
   that pass now: the sum rules of `spectral_weight`, and the D13 identities
   (with the docstring corrected). P4 gets its test with C7. About 4–5
   commits. Needs Q22.
2. **`build/numpy-floor`** (X4): run the tests on numpy 1 and 2, set the
   floor that passes, record the command in `DEV.md`. 1–2 commits. Needs
   Q23.
3. **`fix/susceptibility`**: P1, P2, P10, and A12 if Q25 says so. Each fix
   turns its `xfail` green and carries `!` and `BREAKING CHANGE`, as the
   results change. 3–4 commits.
4. **`fix/cefion`**: P5 (rotate from the kept free-ion operators, without
   A2), P6, P7 (Q24), P8, P11. About 5 commits.
5. **`docs/conventions`**: I prepare a decision sheet for C1–C12 (options,
   what each changes, what the code does now); you decide (Q11); I write
   CLAUDE.md section 3, `CONVENTIONS.md` and `DEV.md` (A1). Suite-wide
   conventions go to the workspace as CONV drafts first.
6. **`fix/conventions`**: the fixes that need a convention: P3 (C1, C11), P4
   and P17 (C7), P9 (C10), P13 and P14 (C9, sources), P15 (CODATA), P16 (C3,
   Q12). Probably two branches.
7. **`refactor/api-v2`**, version 2.0: A2, A3, A4, A5, A10, A11, A14, Q2,
   Q7. Planned first (CLAUDE.md 2.3), after Q13.
8. Later: A6 (generated operators), A7 and A8 (point groups, data tables),
   A9 (package data).

## 8. Questions for the maintainer

The IDs continue after the questions asked in the review sessions (Q15–Q21),
so they don't clash.

- **Q11 (C1–C12).** Which conventions do you want to decide first? P4, P17
  (the g-tensor) and P9 (the cross section) wait for C7 and C10.
- **Q12 (P16).** Should `CEFpars` accept odd-rank parameters?
- **Q13 (A2).** Is an immutable solved-state object (section 3, A2) the
  direction you want for the API, or should `CEFion` stay as it is with only
  the bugs fixed?
- **Q22 (WP 1).** Which independent references should the tests use
  (CLAUDE.md 6.5, `OPEN:`)? Analytic limits and identities only, or also
  PyCrystalField as a numerical reference in a test dependency group, which
  is a new dependency?
- **Q23 (X4).** Check numpy 1 and 2 with `nox` (a new dev dependency,
  CLAUDE.md 6.7 `OPEN:`), or with a plain `uv run --with numpy==...` command
  recorded in `DEV.md`?
- **Q24 (P7).** Fix `get_degeneracies` by rounding now, or wait for the
  central relative tolerance (A11) and group the levels by it?
- **Q25 (A12).** Remove the print in `susceptibility` together with P1 (same
  function), although the A items are delayed?
