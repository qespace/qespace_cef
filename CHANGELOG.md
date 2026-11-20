# Changelog

`qespace_cef` continues `crysfipy`; its history and version numbers carry on.

## Unreleased

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
