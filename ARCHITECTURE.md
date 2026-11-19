# Architecture of qespace_cef

How the package is built today (version 1.3.0). Known problems and the plan
to fix them are in [`AUDIT.md`](AUDIT.md); this file describes the code as it
is and is updated with every change that moves a module, a name or a
convention (CLAUDE.md, section 7).

## Data flow

```
Ion('Ce')                 J, gJ, Stevens factors, form-factor coefficients
   |
CEFpars({'B20': ...})     CEF parameters B_k^q, converted to meV
   |
CEFion(ion, cfp, Hfield)  H = sum B_k^q O_k^q  - muB gJ J.H   (2J+1 x 2J+1)
   |  diagonalize()       energies (meV, ground state = 0), eigenvectors,
   |                      J operators rotated into the eigenbasis
   v
observables               populations, S(Q, w), M(H, T), chi(T), g-tensor,
                          Z, U, F, S, Cv
```

Every observable takes a diagonalised `CEFion`. Observables in a field
(`magnetization`, `susceptibility(method='magnetization')`) build a new
`CEFion` with that field from the same `Ion` and `CEFpars`.

## Module map

| Module | Contents | Depends on |
|---|---|---|
| `__init__.py` | Public API: `Ion`, `CEFpars`, `CEFion` and `from .observables import *`; `__version__` read from the package metadata | all below |
| `constants.py` | Physical constants and unit conversions (meV, K, T, cm⁻¹, emu) | – |
| `cefmatrices.py` | Angular-momentum matrices `J_x`, `J_y`, `J_z`, `J_plus`, `J_minus`; class `StevensBase(J)`, which gives the Stevens operators `O_k^q` for k = 0 to 6, every q, as methods `O20()`, `O4m3()` and as `StevensBase(J)['20']` | numpy |
| `ion.py` | `Ion(name)`: rare-earth ions Ce, Pr, Nd, Sm, Tb, Dy, Ho, Er, Tm, Yb (3+). J, `gJ`, Stevens factors α, β, γ, `<j0>` and `<j2>` coefficients, and `mff(Q)`, the dipole magnetic form factor | numpy |
| `cefpars.py` | `CEFpars(Bdict, units, pointGroup)`: dictionary of `B_k^q` in meV, checked against the operators `StevensBase` implements; `allowed_Bpars(pointGroup)`, a table of the non-zero parameters for each point group | `cefmatrices`, `constants`, `data_containers` |
| `data_containers.py` | `point_group_synonyms` (Hermann–Mauguin, Schoenflies, lattice); `composeTables`, a developer tool that generated the ion table in `ion.py` from `atoms-form-factors.txt` | numpy |
| `cefion.py` | `CEFion`: builds the Hamiltonian, checks that it is Hermitian, diagonalises it with `scipy.linalg.schur`, sorts the levels, rotates `Jx`, `Jy`, `Jz` into the eigenbasis; `magnetic_moment`, `get_degeneracies`, `__str__` | `ion`, `cefpars`, `cefmatrices`, `constants`, scipy |
| `observables.py` | `boltzman_population`, `spectral_weight`, `magnetization`, `magnetization_exchange`, `susceptibility`, `gtensor`, `thermodynamics` | `cefion`, `constants` |
| `symbolic/` | Experimental `sympy` scripts and their own `cefmatrices.py`. Not imported, not tested, not API | sympy |
| `atoms-form-factors.txt` | J. Brown's form-factor table (ILL), the input of `composeTables` | – |
| `thermodynamics-magnetism.png`, `README.md` (in the package) | A figure and the book reference. Neither is used by the code | – |

Layering inside the package, from the bottom up: `constants` and
`data_containers` → `cefmatrices` and `ion` → `cefpars` → `cefion` →
`observables`. Nothing imports upward.

## Conventions as implemented

These are what the code does today, not decisions: every row marked `OPEN:`
in CLAUDE.md section 3 is still open, and no row here settles it.

| Topic | Implemented as | Where |
|---|---|---|
| Energy | meV internally; `CEFpars` accepts `meV`, `K`, `invcm`. After diagonalisation the ground state is set to 0 | `cefpars`, `cefion.diagonalize` |
| Temperature, field | K; T (a 3-vector in the frame of the J operators) | `observables`, `CEFion` |
| Momentum transfer | \|Q\| in Å⁻¹ (assumed: the code computes s = Q/4π, as in the ILL tables) | `Ion.mff` |
| Basis | Ground J multiplet, `\|J>, \|J-1>, ..., \|-J>` (the first vector component is m = +J) | `cefmatrices` |
| Operators | Stevens `O_k^q`; q < 0 written `m`, e.g. `O4m3`. Taken from McPhase and EasySpin; the tests compare them with PyCrystalField up to a sign | `cefmatrices`, `tests/test_cefmatrices.py` |
| Parameters | `B_k^q` in the Stevens convention, `H = Σ B_k^q O_k^q`; named `B20`, `B4m3`, ... | `cefpars`, `cefion` |
| Zeeman term | `H_Z = -μB gJ J·H`; the moment is `+gJ <J>` | `cefion` |
| Degenerate levels | The eigenvectors inside a degenerate level are whatever `schur` returns | `cefion.diagonalize` |
| Tolerance | Absolute: `constants.numerical_zero = 1e-12`, plus parameters for the decimal places in `diagonalize` and `get_degeneracies` | `constants`, `cefion` |

## Dependencies

Runtime dependencies as declared: `numpy`, `scipy`, `sympy`, `matplotlib`. What
the core actually imports: `numpy`, plus `scipy` for one call
(`scipy.linalg.schur`). `sympy` is used only by `symbolic/`, and `matplotlib`
by nothing. Dev: `pytest`, `sphinx`, `furo`. Build: `uv_build`.

## Tests

```
tests/
  __init__.py
  test_cefmatrices.py   Stevens matrices vs the old implementation and
                        PyCrystalField, Hermiticity; J = 1/2 ... 8
  test_cefion.py        script under __main__ only, no tests
  test_cefpars.py       script under __main__ only, no tests
```

`cefion`, `cefpars`, `ion` and `observables` have no tests.

## Documentation

`docs/` is a Sphinx project: `index.rst` and one page per public object
(`source/*.rst`, autodoc), built with `uv run python docs/build.py` into
`docs/_build/html`. The docstrings mix the NumPy and Google styles; napoleon
reads both.
