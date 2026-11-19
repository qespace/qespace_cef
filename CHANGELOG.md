# Changelog

`qespace_cef` continues `crysfipy`; its history and version numbers carry on.

## Unreleased

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
