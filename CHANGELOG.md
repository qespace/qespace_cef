# Changelog

`qespace_cef` continues `crysfipy`; its history and version numbers carry on.

## Unreleased

- New `qespace_cef.config`, the one place for tolerances, read at call
  time: `rtol = 1e-10` and `E_floor = 1` meV (CONV-0008, AUDIT A11). New
  `qespace_cef.cefion.degenerate_levels(energies)` groups ascending
  energies into the levels of CONV-0008, as slices.
- **Breaking:** `qespace_cef.constants` takes its values from CODATA 2022
  (e, h, c, k exact; μB = 9.2740100657(29)·10⁻²⁴ J/T), cited in the module,
  and derives the conversions from them: `K2meV = 1e3 k/e`,
  `meV2K = 1/K2meV`, `uB = 1e3 μB/e`, `invcm2meV = 1e3 · 100 h c/e`
  (AUDIT P15). Before, they mixed CODATA 1986 and 2010 and were typed
  separately, so `K2meV` and `meV2K` were not reciprocal. Results move,
  relative, largest in the cases checked, by the unit of the parameters:
  meV: energies 0 (zero field), M 2.4·10⁻⁷, χ 3.2·10⁻⁷; cm⁻¹: energies
  3.8·10⁻⁷, M 1.9·10⁻⁷, χ 5.8·10⁻⁷; K (`K2meV` had 5 digits): energies
  3.9·10⁻⁶, M 8.8·10⁻⁷, χ 1.5·10⁻⁶. The
  unused names `pi`, `R0`, `kB_SI`, `NA`, `C1`, `C2`, `C3` and the private
  `_m`, `_gn`, `_hq`, `_gamn`, `_r0` are removed.
- `susceptibility(method='perturbation')` evaluates the van Vleck weight
  `(p_n − p_m)/(E_m − E_n)` as `p_n (−expm1(−(E_m − E_n)/k_B T))/(E_m −
  E_n)` on the pairs `E_m > E_n`, mirrored to the others, instead of a
  plain difference that cancelled for two levels just above the degeneracy
  rule (AUDIT P20, CONV-0008). Ce³⁺ at 5 K with the ground doublet split by
  6·10⁻¹¹ meV: χ_zz was off by 6.9·10⁻⁷ relative, now by 2·10⁻¹⁵. Other
  results change only by rounding.
- The J matrices are built in one place, `J_x`, `J_y`, `J_z`, `J_plus`
  and `J_minus`, and `StevensBase` takes its `Jx`, `Jy`, `Jz`, `Jp`, `Jm`
  from them instead of its own copy of the formulas (AUDIT A5). They
  follow CONV-0005 (order `|J>, ..., |−J>`, Condon–Shortley phases) and
  are tested against it. No change to the results: every operator matrix
  is bit-identical.
- **Breaking:** `CEFion.__str__` is removed (CONV-0010); `str(cefion)`
  gives Python's default text. `CEFion.summary(evecs, precision)` prints
  the level scheme instead: the levels above the ground state with their
  degeneracies and eigenvectors, and E_0 in the header. For this print
  only, states within 10⁻¹⁰ meV count as one level.
- **Breaking:** `CEFion.energies` are the eigenvalues of the Hamiltonian,
  not shifted (CONV-0010): with `B00 = 0` they average to zero, the centre
  of gravity of the J multiplet. The new `CEFion.excitation_energies`
  (`energies - energies[0]`) start at 0, as `energies` did. `diagonalize`
  has no `shiftToZero` argument any more. `get_degeneracies` and U and F of
  `thermodynamics` move with `energies`; populations, χ, M, the spectra, S
  and C_v do not depend on the zero.
- `magnetization_exchange` stops after `max_iterations` steps (new keyword,
  default 1000) and raises `RuntimeError` if the mean field has not
  converged; it could loop forever. A zero moment converges at once
  instead of giving 0/0 (AUDIT P11).
- `str(Ion)` shows a half-integer J as a fraction, `'Ce: J = 5/2, gJ =
  0.86'`; it printed 2 (AUDIT P8).
- **Breaking:** `CEFion.diagonalize` has no `sortWithE` argument any more.
  `numpy.linalg.eigh` already returns the levels ascending, so sorting
  changed nothing, and `sortWithE=False` raised `IndexError` (AUDIT P6).
- `CEFion.diagonalize` can be called again: it rotates the free-ion J
  operators, not the already rotated ones (AUDIT P5). The first call is
  unchanged.
- **Breaking:** `susceptibility(cefion, temperatures, method)` returns the
  tensor `chi_ab`, shape `S + (3, 3)` for temperatures of shape `S`, in
  μB/T per ion from both methods (CONV-0009). `Hfield_direction` is
  removed; `n @ chi @ n` gives the value along `n`. `'perturbation'` was
  in μB²/meV and took the direction wrongly (AUDIT P2); `'magnetization'`
  returned |M|/ε instead of the component along the field (P10).
- `susceptibility(method='perturbation')` no longer prints the matrix of
  transition energies (AUDIT A12).
- **Breaking:** `susceptibility(method='perturbation')` counts the van
  Vleck term over every pair of levels in both orders (CONV-0009); it
  was half as large (AUDIT P1). Ce³⁺ with `B20=0.5, B40=0.01, B44=0.1`
  along z: +3 %, +11 %, +18 % at 5, 20, 100 K.
- Decisions are recorded in `decisions/`: design choices as ADRs, physics
  conventions as CONVs. The CEF parameters and the Stevens operators follow
  McPhase (CONV-0001), checked by tests against its table, signs included.
  No change to the code or its results.
- The numpy floor is 1.26 (it was 2.4.4); numpy 1 and 2 are tested on
  Python 3.12–3.14 with `nox`. `DEV.md` lists the commands.
- **Breaking:** the only runtime dependency is `numpy`. `scipy` and
  `matplotlib` are dropped; `sympy` moves to the dependency group
  `symbolic`, for the experimental scripts.
- **Breaking:** `CEFion.diagonalize` uses `numpy.linalg.eigh` instead of
  `scipy.linalg.schur`. Energies and basis-independent results are
  unchanged; inside degenerate levels the eigenvectors, and what depends on
  them (`gtensor`, single moments and transitions, `str()`), differ.
  `diagonalize` raises `ValueError` for a non-Hermitian Hamiltonian.
- **Breaking:** `CEFpars.from_symmetry`, which never worked, is removed,
  with other dead code.
- Code formatted with `ruff`; tests for `CEFion`, `CEFpars` and `Ion`;
  performance stubs (`pytest -m perf -s`).
- Docstrings in NumPy style with units and shapes; known problems from
  `AUDIT.md` marked in the code and in Warnings sections.
- `str(CEFion)` prints energies with `precision` decimals (`9.968 meV`,
  was `9.97 meV`) and "2-fold degenerate". `CEFion.summary()` is a stub
  for later. `get_degeneracies`, `magnetization` and `CEFpars.Bpars` give
  Python floats and ints instead of NumPy scalars.

## 1.3.0

- Rename the package from `crysfipy` to `qespace_cef`, part of the qespace
  suite.
- Build with `uv_build`; licence GPL-3.0-or-later shipped as `LICENSE`.
- The version is set in `pyproject.toml` only; `__version__` reads it from
  the package metadata.
- Tests of the Stevens matrices run and assert (J = 1/2 ... 8).
- Documentation sources pruned; build with `uv run python docs/build.py`.
- Symbolic scripts in `qespace_cef/symbolic/` marked experimental.

## 1.2

- Implement additional Stevens operators and reimplement the way to
  initialise crystal-field parameters.

## 1.1

- Implement magnetisation.

## 1.0

- Hello World!
