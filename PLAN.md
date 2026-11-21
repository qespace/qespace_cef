# Plan: the next implementation step

The step being worked on now, in the form of CLAUDE.md 2.3. It borrows from
[`AUDIT.md`](AUDIT.md) or comes from a discussion, follows the rules
(`CLAUDE.md` and `decisions/`), and is replaced when the step is done. It is
tracking, not a rule (ADR-0001). Ideas for later are in
[`ROADMAP.md`](ROADMAP.md).

## Now: `feat/energy-zero`, CONV-0010

**Goal.** `CEFion` follows CONV-0010 (accepted):

- `energies` are the eigenvalues of `hamiltonian`, not shifted. With
  `B00 = 0` their zero is the centre of gravity of the J multiplet.
- `excitation_energies = energies − energies[0]` sits next to them and starts
  at 0.
- `diagonalize` loses `shiftToZero`.
- `CEFion.__str__` goes. The `summary` stub stays as the placeholder for a
  print of the level scheme (AUDIT Q7, `ROADMAP.md`).

**Steps.** Branch `feat/energy-zero`, 4 commits.

1. `feat(cefion)!: keep the energies unshifted, add excitation_energies
   (CONV-0010)`.
   - `diagonalize` stores the eigenvalues as `energies`, and `energies −
     energies[0]` as `excitation_energies` (`[]` before diagonalisation,
     like `energies`). `shiftToZero` goes.
   - The class docstring, its doctest and the README line describe both.
   - Tests that encoded the shift move to the new convention, with the same
     assertions:
     - `test_axial_levels` checks `energies` against its own docstring
       formula, `B20 (3m² − J(J+1))` = −8, −2, 10 meV, and
       `excitation_energies` against 0, 6, 18;
     - `test_energies_sorted_from_zero` becomes ascending `energies` plus
       `excitation_energies[0] == 0`;
     - `test_eigenvectors_belong_to_energies` loses its `E0` correction;
     - `test_degeneracies_are_plain_python_numbers` expects −8, −2, 10.
   - Results change: `energies`, `get_degeneracies`, and U and F of
     `thermodynamics` (P3, still open) shift by `E_0`. Everything else is
     shift invariant: populations, χ, M, the spectra, S and C_v.
2. `test(cefion): check the zero of energy against CONV-0010`.
   - `test_conv_0010_*`, from its Consequences:
     - `energies` equal `eigvalsh(hamiltonian)` and average to 0, with and
       without a field;
     - free Ho³⁺ at 2 T along z gives `E_0 = −gJ μB J B`;
     - `excitation_energies` is `energies − energies[0]`;
     - Hellmann–Feynman, `∂E_n/∂B40 = <n|O_40|n>`, against a finite
       difference on low-symmetry Ce³⁺ in a field.
3. `refactor(cefion)!: remove CEFion.__str__`.
   - `__str__` goes, with `test_str_prints_energies_with_decimals`, which
     tests only the removed method.
   - The `summary` stub stays, and its message points to the planned print
     of the level scheme instead of `__str__`.
4. `docs: describe the energies of CONV-0010`.
   - `ARCHITECTURE.md`: the data flow, the Energy row of the conventions
     and the test tree.
   - `AUDIT.md`: Q7 and section 6.2 (CONV-0010 accepted).
   - `ROADMAP.md`: the print of the level scheme.
   - `CHANGELOG.md`.

**Files.** `qespace_cef/cefion.py`, `tests/test_cefion.py`, `README.md`,
`ARCHITECTURE.md`, `AUDIT.md`, `ROADMAP.md`, `CHANGELOG.md`.

**Checked in advance** (scratch run on `c94cb46`, with a pytest plugin that
forces `shiftToZero=False`, the package unchanged):
- 9 tests and the `CEFion` doctest fail, and each of them encodes the shift:
  - `test_axial_levels`;
  - `test_energies_sorted_from_zero` (3 ions);
  - `test_eigenvectors_belong_to_energies` (3 ions);
  - `test_degeneracies_are_plain_python_numbers`;
  - `test_str_prints_energies_with_decimals`.
- All other 347 tests pass. The 30 strict xfails stay xfails, P3 and P7
  among them.
- The four checks of step 2 hold: mean 0 to 10⁻¹²; E_0 = −1.15768 meV
  exactly; Hellmann–Feynman to 1.4·10⁻⁷ on a largest value of 60.5.

**Risks.**
- **`get_degeneracies` (P7)** truncates `energies` to 5 decimals. With
  negative energies, truncation goes towards zero, so a level at −7.9999999…
  and one at −8.0000000… can land in different groups. The P7 xfail stays.
  The fix of P7 waits for A11 and CONV-0008.
- **`thermodynamics`:** U and F change by `E_0` while P3 is still open. The
  P3 xfails don't depend on the zero and stay xfails.
- **`str(cefion)`** prints Python's default representation after step 3. Any
  script that printed a `CEFion` changes.

**Tests.** `uv run pytest`: the rewritten tests and the `test_conv_0010_*`
tests pass, and the xfails stay at 30. Then `ruff check`,
`ruff format --check` and `uv run nox`.

## In parallel: decisions only you can take

- The `OPEN:` parts of the drafts CONV-0004 to CONV-0008 and CONV-0011
  (AUDIT section 6.2, Q11), and a draft for the g-tensor (C7).
- Q29 and Q31 (ADR-0002), which block `refactor/api-v2`.

## Next

Most remaining fixes wait for a convention (AUDIT section 6.2, Q11):
- P3 for CONV-0003 (its zero of U and F is CONV-0010, accepted);
- P4 and P17 for C7 and CONV-0008;
- P9 and P13 for CONV-0007;
- P7 for A11 and CONV-0008.

The others need no CONV:
- P14 and P15 need sources: the ion table and one CODATA set;
- P16 needs your answer to Q12.

`refactor/api-v2` waits for Q29 and Q31.
