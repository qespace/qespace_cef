# Code audit of qespace_cef

Open findings from reviews of the code: the audit of 2026-11-19 (commit
`747c871`, version 1.3.0) and the reviews since, as of `main` at `423827b`
(2026-11-20). This file tracks; it is not a rule (ADR-0001). It cites the
rules (`CLAUDE.md`, the ADRs and CONVs in [`decisions/`](decisions/)) by ID.
The step being worked on is in [`PLAN.md`](PLAN.md), feature ideas are in
[`ROADMAP.md`](ROADMAP.md), what was fixed is in
[`CHANGELOG.md`](CHANGELOG.md). Remove an item when its change lands; IDs
are never reused.

Places are named by module and function. Every P item still in the code has
an `# AUDIT Pn:` comment where it happens and a Warnings section in the
docstring of the public object it affects. P3 and P17
also have strict `xfail` tests with expected values that don't come from
the code (`test/known-problems`); a fix turns its tests green and removes
the markers. P4 gets its test with the g-tensor convention (C7); P13, P14 and
P16 are data and conventions.

Severity: **bug** = wrong result, verified by running it; **bug (read)** =
wrong from reading, not run; **rule** = breaks a rule (CLAUDE.md, an ADR or a CONV);
**docs** = documentation does not match behaviour; **suggestion** = no rule
broken. The evidence was rerun on `610c884` where it says so; the rest is
from `747c871` and does not depend on what changed since.

## 1. Physics and correctness

| ID | Severity | Where | Finding | Evidence | Proposal |
|---|---|---|---|---|---|
| P3 | bug | `observables.thermodynamics` | Energies in meV are mixed with temperatures in K: `F = −T ln Z` is in K, while `E` stays in meV, so `S = (E − F)/T` and `Cv` come out in mixed units | Entropy (k_B) at 10 / 20 / 40 K: exact 1.293 / 1.374 / 1.494, code 1.054 / 1.203 / 1.325. Cv (k_B): exact 0.170 / 0.108 / 0.239, code 0.00126 / 0.00080 / 0.00177 (low by a factor of 11.6², about 135) | Convert the energies to `k_B T` once and state the units of the returns (CONV-0003); test against `−∂F/∂T` and the high-T limit `ln(2J+1)` |
| P4 | bug, rule | `observables.gtensor` | Even with a correct formula, g itself depends on the arbitrary basis inside the ground doublet: a change of basis turns g into g R with R a rotation, and only `G = g gᵀ` (the principal values and axes) is basis independent. The function returns g in whatever basis `eigh` picks. That goes against the draft on degenerate levels (CONV-0008) | Ce, `B20=1`, no field: diag(2.571, 2.571, 0.857) with `schur`, diag(2.571, −2.571, 0.857) with `eigh` | Return G, or the principal values and axes; together with P17 and C7 (CONV-0008) |
| P13 | data | `Ion` table, Ce | The `<j0>`/`<j2>` coefficients for Ce are those of Ce²⁺, the only Ce entry in the ILL table (`atoms-form-factors.txt`); the generator that built the table (`composeTables`, in the history at `747c871`) mapped `ce` to `Ce2`. The other ions use the 3+ values | `atoms-form-factors.txt`: only `{F Ce2 ...}` lines | Record it as a known approximation next to the value, with its source (CONV-0007) |
| P14 | data, rule | `Ion` table | The values (J, gJ, α, β, γ, form factors) have no source next to them (CLAUDE.md 3.2.2); the generator named "Czech guys" and "J.P. Brown" | read | Cite Stevens (1952) / Hutchings (1964) for α, β, γ and the ILL / ITC Vol. C tables (with DOI) for the form factors |
| P16 | physics question | `cefmatrices`, `CEFpars` | Odd-rank operators (l = 1, 3, 5) are implemented and accepted as CEF parameters. Inside one J multiplet, parity removes odd l from the CEF Hamiltonian | `CEFpars(dict(B52=0.1, B3m1=0.2))` is accepted. Rerun on `610c884` | Keep them as building blocks, but decide whether `CEFpars` accepts them (Q12) |
| P17 | bug | `observables.gtensor` | The z column uses `\|<0\|J_α\|0>\|`. The `abs()` makes even `G = g gᵀ` depend on the basis inside the doublet, so the principal g values are wrong in a general basis. (From the maintainer's notes: "I feel implementation of `gtensor` is wrong") | Ce, 5 random bases of the ground doublet: the spread of g gᵀ is 2.93 (`B20=1`) and 2.15 (low symmetry) with the code's formula, 4e-15 with `g_αβ = gJ Tr(P J_α σ_β)`. For `B20=1` the principal values of the latter are 0.857, 2.571, 2.571 (correct). Run on `610c884` | `g_αβ = gJ Tr(P J_α σ_β)` over the doublet, reported as G or principal values (P4); the convention with C7; non-Kramers doublets need their own treatment |
| P21 | convention | `CEFion` (Zeeman term), `CEFion.magnetic_moment`, `magnetization`, `magnetization_exchange`, `susceptibility(method='magnetization')` | The code uses the moment-parallel convention, `−gJ μB J·B` and the moment `+gJ μB <J>`; the draft CONV-0006 takes the electron's sign, `+gJ μB J·B` and `−gJ μB <J>` | read (ARCHITECTURE.md, Physics conventions) | Flip both signs in one change, with CONV-0006's tests; energies, M along B, χ and the cross sections stay, ⟨J⟩ and the eigenvectors in a field become their time-reversed partners. Probably in `refactor/api-v2` |

## 2. Documentation against behaviour

None open.

## 3. Architecture

| ID | Severity | Finding | Proposal |
|---|---|---|---|
| A2 | suggestion | `CEFion` is a stateful object: the constructor builds H, `diagonalize()` mutates `Jx/Jy/Jz` and changes their meaning (free-ion basis before, eigenbasis after; a list before, an array after), and `energies` is `[]` until diagonalised. This made a second `diagonalize()` rotate J twice (P5, fixed) and makes every observable depend on the call order | ADR-0002 (`Hamiltonian`). Open there: when it is solved and whether it is immutable (Q29), observables as functions or methods (Q31) |
| A3 | suggestion | Observables in a field rebuild a whole `CEFion` from `ion` and `cfp` and ignore `cefion.Hfield`, so the field lives in two places | Pass the field to the solver only; the observables take a solved state |
| A4 | rule | `from .observables import *` exports `np`, `C`, `ArrayLike`, `CEFion` again, and the submodules. With no `__all__`, the public API is undefined (CLAUDE.md 6.2 "public entry point") | `__all__` in `__init__` and in every module |
| A5 | rule (CONV-0005) | `symbolic/cefmatrices.py` is a second, sympy, construction of the J matrices besides the module functions `J_x`, `J_y`, `J_z`, `J_plus`, `J_minus`. It is experimental, not imported and not tested. (`StevensBase` takes its J matrices from the module functions since `refactor/j-matrices`, with all its matrices bit-identical.) | Waits with A9 |
| A6 | suggestion | `StevensBase` has 49 hand-written methods. It works and is tested, but adding l = 7… or checking the normalisation means editing 49 formulas | Keep it as the reference for now; later, generate `O_lm` from one closed formula (e.g. via the Racah tensor operators) and keep the hand-written table only as a test fixture |
| A7 | suggestion | `CEFpars` stores the point group but doesn't use it; `allowed_Bpars` is a lookup table without the relations between parameters (cubic B44 = 5 B40) | Decide what point-group support should be (CONV-0004, `ROADMAP.md`), then build it on the one table, with tests |
| A8 | suggestion | The ion data sit in a dictionary literal in `Ion.__init__`, rebuilt on every `Ion()`; `data_containers.py` holds only the point groups | Put the data tables (ions, point groups, allowed parameters) in one data module with their sources; move `atoms-form-factors.txt` out of the package (e.g. `tools/`) |
| A9 | suggestion | Files that ship without a use: `thermodynamics-magnetism.png`, `qespace_cef/README.md`, `atoms-form-factors.txt` (kept for now, by your decision) | Move the reference to the README or the docs; keep the rest out of the wheel when the data question is settled |
| A10 | rule | Validation: `assert` in `gtensor` (removed by `python -O`), a bare `Exception` for a non-Hermitian H in `CEFion`, `NameError` for an unknown point group, `KeyError` for an unknown ion; the Hermiticity check in the constructor is an absolute `allclose` (CLAUDE.md 5.5, 5.6) | `ValueError` at the public boundary, relative tolerances from one config place |
| A11 | rule | `qespace_cef.config` holds only the level tolerance of CONV-0008 (`rtol`, `E_floor`; `rtol · J` also for the J_z degeneracy of CONV-0011). Still outside it: `Eimag_precision=10` in `diagonalize`, the coefficient threshold `10⁻precision` of `CEFion.summary`, `eps = 1e-5` T, `precision=1e-5` (CLAUDE.md 5.6) | Move them into `config` as relative tolerances, read at call time |
| A14 | suggestion | Stevens operators can be had two ways, `StevensBase(J).O4m3()` and `StevensBase(J)['4m3']`, and `implementedOps` finds them through `cls.__dict__`: two access paths for one thing, with string keys to parse. (From the maintainer's notes: "Should we settle on just one of them?") | Settle on one with A2 and A6, e.g. `stevens(J, l, m)` or `ops[l, m]` with integers; the methods go when A6 generates the operators. Breaking (v2) |
| A15 | suggestion | Every call solves one point: a scan over parameters or fields is a Python loop over `CEFion` objects, 15× slower than one batched diagonalisation (ADR-0003, measured on `afa1dd5`) | ADR-0002 and ADR-0003, in `refactor/api-v2`; the feature in `ROADMAP.md` |

## 4. Code quality

| ID | Finding | Proposal |
|---|---|---|
| Q2 | Naming is not uniform, and no rule says what it should be. camelCase: `pointGroup`, `pointGroupName`, `mffCoefficients`, `implementedOps`. Words fused with symbols or capitalised: `Hfield`, `Hfield_direction`, `Bpars`, `Bdict`, `allowed_Bpars`, `Jval`, `J2p1`, `Eimag_precision`, `Alpha/Beta/Gamma`, the classes `CEFion`, `CEFpars`. The typo in `boltzman_population`. Physics symbols (`J`, `Q`, `T`, `Ei`, `B_iso`, `gJ`, `Jx`, `O4m3`) follow the formulas. The code and the docstrings index the operators and parameters with k, q; CONV-0001 and the documentation use l, m. Scanned at run time on `8a39456` | The classes are renamed in ADR-0002 (`CEFpars` → `Blm`, `CEFion` → `Hamiltonian`, `import qespace_cef as cef`). For the rest, a written rule (Q26): snake_case for words, physics symbols as in the formulas, listed; then rename with A2, in v2, possibly with deprecated aliases in a last 1.x |

## 5. Dependencies

None open. The runtime needs only `numpy>=1.26`.

## 6. Comments and conventions

### 6.1 Comments that ask for something

None open.

### 6.2 Conventions the findings wait for

The definitions found in the code (C1–C12 of the audit) are now CONVs in
`decisions/conventions/`, except C7. Each accepted CONV gets a test that
names its ID (CLAUDE.md 6.6).

| CONV | Status | Open | Findings waiting |
|---|---|---|---|
| CONV-0001 parameters and operators | accepted | Q12, Q33 | P16 |
| CONV-0003 units | accepted | its Statement still says `OPEN:` the units of the API (Q42) | P3 |
| CONV-0004 frames | accepted | – (one local frame; the crystal frame later) | A7 |
| CONV-0005 basis | accepted | – (the J matrices are the root: order `\|J>, ..., \|−J>`, Condon–Shortley phases; everything else follows) | A5 |
| CONV-0006 Zeeman term | accepted | – (the electron's sign, μ = −gJ μB J, ℋ_Z = −μ·B) | P21 |
| CONV-0007 cross sections | accepted | – | P13 |
| CONV-0012 spectral weights | accepted | – | – |
| CONV-0008 degenerate levels | accepted | – (the tolerance `rtol · max(max\|E\|, E_floor)`, defaults 1e-10 and 1 meV, from the config) | P4 |
| CONV-0011 canonical basis | accepted | – (J_z inside the level, J_x then J_y, Kramers partners by time reversal; integer-J levels unpaired until the irreps replace it) | the eigenvectors `summary` prints |
| CONV-0013 susceptibility (supersedes CONV-0009) | accepted | – | P21 (the finite-difference test with the moment of CONV-0006) |
| CONV-0010 zero of energy | accepted | – | P3 (the zero of U and F) |
| C7 g-tensor | no draft | its definition from the ground doublet, what is returned, non-Kramers doublets | P4, P17 |

## 7. Work packages for the open findings

Each is one branch (CLAUDE.md 1.3); the one being worked on is in
`PLAN.md`. Each P fix is its own commit with a physics explanation, turns
its `xfail` green and carries `!` and `BREAKING CHANGE` when results change.

1. **`fix/conventions`**: the fixes that wait for a CONV (section 6.2) or
   need sources: P3, P4, P13, P14, P16, P17. Probably two branches.
2. **`refactor/api-v2`**, version 2.0, on ADR-0002 and ADR-0003: A2–A4, A10,
   A11, A14, A15, Q2 (with Q26). Planned first (CLAUDE.md 2.3),
   after Q29 and Q31.
3. Later: A6 and A7 (in `ROADMAP.md`), A8, A9.

## 8. Questions for the maintainer

Answered questions are removed; their answers are in the ADRs, the CONVs and
`PLAN.md`. The IDs continue after those of the review sessions.

- **Q11 (section 6.2).** Which draft do you want to decide next? P4 and P17
  wait for C7 (no draft yet) and CONV-0008, P3 for the
  units of the thermodynamic quantities (Q44).
- **Q12 (P16, CONV-0001).** Should `CEFpars` accept odd-rank parameters
  (l = 1, 3, 5)? CONV-0001 tabulates only even l.
- **Q26 (Q2).** The naming rule for everything but the classes (ADR-0002):
  snake_case for words and physics symbols as in the formulas (which ones?),
  and is it a rule for this package or for the whole suite (the workspace's
  `DEV.md` and glossary)?
- **Q29 (ADR-0002).** When is a `Hamiltonian` solved? (a) It is immutable and
  solves itself on first access to `energies` (cached); a new `Blm` means a
  new `Hamiltonian`, so a fit loop never sees a stale solution (A2). (b) An
  explicit `H.solve()` that fills the object in place. I recommend (a).
- **Q31 (ADR-0002).** Observables as functions, `cef.magnetization(H, T)` as
  in v1, or as methods, `H.magnetization(T)`? I recommend functions: new
  observables then don't change the class.
- **Q32 (ADR-0001).** IDs: this package numbers its own `ADR-xxxx` and
  `CONV-xxxx`, and a suite ADR copied here keeps its suite ID (CLAUDE.md 3.4),
  so two different ADR-0001 can meet. Prefix ours (e.g. `CEF-ADR-0001`), or
  keep suite copies in their own folder?
- **Q33 (CONV-0001).** CONV-0001 labels the sine-type parameters `B4m3`,
  the package's current name. McPhase itself writes `B43S` (manual node123:
  "negative q parameters such as B_2^{-2}, are specified as: B22S"). Keep
  `B4m3` and map McPhase files on input, or switch to `B43S`?
- **Q34 (CONV-0001).** The CONV template has an `enforced_by` line naming the
  test; CONV-0001 has none, and an accepted CONV is never edited. Leave it,
  with the tests naming the CONV (`test_conv_0001_*`), or supersede CONV-0001
  with a copy that adds the line?
- **Q42 (CONV-0003).** The accepted CONV-0003 still says "`OPEN:` the
  units of the API" and lists the units "As implemented" (meV, K, T, Å⁻¹,
  μB per ion, χ in μB/T per ion). Is that list the decision? An accepted
  CONV is never edited, so the clean form is a new CONV that supersedes it
  without the `OPEN:` line.
- **Q44 (P3).** The units of the thermodynamic quantities, which no CONV
  states: U and F in meV per ion from the zero of CONV-0010, S and C_v in
  k_B per ion? With that, P3 can be fixed.

## 9. Needs later attention

- **A16 (CONV-0012).** The excitations are returned flat, n² entries without
  level information, as CONV-0012 says; that came from a misunderstanding
  (maintainer, 2026-11-22). Intended: `energies` are the levels of the
  Hamiltonian, (n,); the transition energies `E_f − E_i` keep the pairs,
  (n, n), and so do the weights, (n, n, 3, 3) and (n, n). As a solved
  Hamiltonian is immutable (Q29), the pair (i, f) traces back to the levels;
  the shape costs nothing, as a reshape is a view. Needs a CONV superseding
  CONV-0012, then `scattering_function`, the cross sections and the
  `test_conv_0012_*` tests follow.
