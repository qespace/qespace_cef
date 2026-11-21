# Architecture of qespace_cef

How the package is built today (version 1.3.0). Known problems are in
[`AUDIT.md`](AUDIT.md), the next step in [`PLAN.md`](PLAN.md); this file
describes the code as it is and is updated with every change that moves a
module, a name or a convention (CLAUDE.md, section 7).

## Data flow

```
Ion('Ce')                 J, gJ, Stevens factors, form-factor coefficients
   |
CEFpars({'B20': ...})     CEF parameters B_lm, converted to meV
   |
CEFion(ion, cfp, Hfield)  H = sum B_lm O_lm  - muB gJ J.H     (2J+1 x 2J+1)
   |  diagonalize()       energies (meV, eigenvalues of H, zero at the centre
   |                      of gravity), excitation energies, eigenvectors,
   |                      J operators rotated into the eigenbasis
   v
observables               populations, S(Q, w), M(H, T), chi(T), g-tensor,
                          Z, U, F, S, Cv
```

Every observable takes a diagonalised `CEFion`. Observables in a field
(`magnetization`, `susceptibility(method='magnetization')`) build a new
`CEFion` with that field from the same `Ion` and `CEFpars`;
`susceptibility(method='magnetization')` builds three per temperature, one
per field axis.

## Module map

| Module | Contents | Depends on |
|---|---|---|
| `__init__.py` | Public API: `Ion`, `CEFpars`, `CEFion` and `from .observables import *`; `__version__` read from the package metadata | all below |
| `constants.py` | Physical constants from CODATA 2022 (e, h, c, k exact; μB), cited in the module docstring, and the unit conversions derived from them: `K2meV` and `meV2K = 1 / K2meV` (meV/K, K/meV), `uB` (meV/T), `invcm2meV` (meV per cm⁻¹); `numerical_zero`, the degeneracy rule of CONV-0009 | – |
| `cefmatrices.py` | Angular-momentum matrices `J_x`, `J_y`, `J_z`, `J_plus`, `J_minus`, the one construction of J (CONV-0005); class `StevensBase(J)`, which takes its `Jx`, `Jy`, `Jz`, `Jp`, `Jm` from these functions and gives the Stevens operators `O_lm` for l = 0 to 6, every m, as methods `O20()`, `O4m3()` and as `StevensBase(J)['20']` | numpy |
| `ion.py` | `Ion(name)`: rare-earth ions Ce, Pr, Nd, Sm, Tb, Dy, Ho, Er, Tm, Yb (3+). J, `gJ`, Stevens factors α, β, γ, `<j0>` and `<j2>` coefficients, and `mff(Q)`, the dipole magnetic form factor | numpy |
| `cefpars.py` | `CEFpars(Bdict, units, pointGroup)`: dictionary of `B_lm` in meV, checked against the operators `StevensBase` implements; `allowed_Bpars(pointGroup)`, a table of the non-zero parameters for each point group | `cefmatrices`, `constants`, `data_containers` |
| `data_containers.py` | `point_group_synonyms` (Hermann–Mauguin, Schoenflies, lattice) | numpy |
| `cefion.py` | `CEFion`: builds the Hamiltonian, checks that it is Hermitian, diagonalises it with `numpy.linalg.eigh` (levels ascending), rotates `Jx`, `Jy`, `Jz` into the eigenbasis; `energies`, `excitation_energies`, `magnetic_moment`, `get_degeneracies`, `summary` (the level scheme as text) | `ion`, `cefpars`, `cefmatrices`, `constants`, numpy |
| `observables.py` | `boltzman_population`, `spectral_weight`, `magnetization`, `magnetization_exchange`, `susceptibility`, `gtensor`, `thermodynamics` | `cefion`, `constants` |
| `symbolic/` | Experimental `sympy` scripts and their own `cefmatrices.py`. Not imported, not tested, not API | sympy (dependency group `symbolic`) |
| `atoms-form-factors.txt` | J. Brown's form-factor table (ILL), from which the ion table in `ion.py` was generated (the generator, `composeTables`, is in the history at 747c871) | – |
| `thermodynamics-magnetism.png`, `README.md` (in the package) | A figure and the book reference. Neither is used by the code | – |

Layering inside the package, from the bottom up: `constants` and
`data_containers` → `cefmatrices` and `ion` → `cefpars` → `cefion` →
`observables`. Nothing imports upward.

## Conventions as implemented

These are what the code does today. The conventions themselves are the CONVs
in `decisions/conventions/`; no row here sets one.

| Topic | Implemented as | Where |
|---|---|---|
| Energy | meV internally; `CEFpars` accepts `meV`, `K`, `invcm`. `energies` are the eigenvalues of H, not shifted: with `B00 = 0` they average to zero, the centre of gravity of the J multiplet. `excitation_energies` are measured from the ground state (CONV-0010, tested) | `cefpars`, `cefion.diagonalize` |
| Temperature, field | K; T (a 3-vector in the frame of the J operators) | `observables`, `CEFion` |
| Momentum transfer | \|Q\| in Å⁻¹ (assumed: the code computes s = Q/4π, as in the ILL tables) | `Ion.mff` |
| Basis | Ground J multiplet, `\|J>, \|J-1>, ..., \|-J>` (the first vector component is m = +J), Condon–Shortley phases; the J matrices are built once, by the module functions, and every other matrix from them (CONV-0005, tested) | `cefmatrices` |
| Operators | Stevens `O_lm`; m < 0 written `m`, e.g. `O4m3`, the sine-type operators. As in McPhase (CONV-0001), tested against its table; PyCrystalField has the opposite sign for l + m odd | `cefmatrices`, `tests/test_cefmatrices.py` |
| Parameters | `B_lm` in the Stevens convention, `H = Σ B_lm O_lm`; named `B20`, `B4m3`, ... | `cefpars`, `cefion` |
| Zeeman term | `H_Z = -μB gJ J·H`; the moment is `+gJ <J>` | `cefion` |
| Susceptibility | The tensor `χ_αβ` in μB/T per ion, shape `S + (3, 3)` for temperatures of shape `S`, from Van Vleck's formula or a finite difference (CONV-0009, tested) | `observables.susceptibility` |
| Degenerate levels | The eigenvectors inside a degenerate level are whatever `eigh` returns (orthonormal, otherwise arbitrary) | `cefion.diagonalize` |
| Tolerance | Absolute: `constants.numerical_zero = 1e-12` (meV), also the degeneracy rule of CONV-0009, plus parameters for the decimal places in `diagonalize` and `get_degeneracies` | `constants`, `cefion` |

## Dependencies

Runtime: `numpy>=1.26` only, tested with numpy 1.26.0, 2.0.2 and the newest
on Python 3.12–3.14 (`noxfile.py`, see `DEV.md`). Dependency groups (not
published): `dev` with `pytest`, `ruff`, `nox`, `sphinx`, `furo`; `symbolic`
with `sympy`, for the scripts in `symbolic/`. Build: `uv_build`.

## Tests

```
tests/
  reference/
    mcphase.py          Stevens operators from the formulas and table in
                        CONV-0001 (even l), with its own J matrices
    old_stevens.py      first implementation of the Stevens operators
    pycrystalfield.py   PyCrystalField's StevensOp, with attribution
  test_cefmatrices.py   the J matrices (CONV-0005: J_z, the Condon–Shortley
                        J+, the commutator, J², the matrices StevensBase
                        uses); Stevens matrices vs the three references
                        (CONV-0001 for even l, signs included), Hermiticity,
                        the StevensBase observations; J = 1/2 ... 8; strict
                        xfail where PyCrystalField's sign differs
  test_cefion.py        analytic limits: free ion, B20 alone, cubic, Kramers,
                        Zeeman, ascending energies, orthonormality, H = U E
                        U^dagger, eigenvectors paired with energies, the
                        zero of energy (CONV-0010), the summary print,
                        repeated diagonalisation; strict
                        xfail for P7
  test_cefpars.py       unit conversion, point-group tables, the parameter
                        labels of CONV-0001
  test_constants.py     the conversions against CODATA 2022 (values written
                        in the test), K2meV and meV2K reciprocal, uB
                        against NIST's value in eV/T
  test_ion.py           J and gJ from Hund's rules, form factor, str()
  test_observables.py   susceptibility tensor against the exact formula and
                        CONV-0009, thermodynamic identities, sum rules, g
                        values, the mean-field iteration; strict xfail for
                        P3, P9, P17
  test_performance.py   timing stubs, marker `perf`, deselected by default
```

`pyproject.toml` puts `tests/` on the path, so the tests import the
references as `reference.*`. A known problem from AUDIT.md is a strict `xfail`
that names its ID; it fails (XPASS) once the problem is fixed, and the fix
removes the marker. A test of a CONV names its ID (`test_conv_0001_*`); a
reference code that disagrees with an accepted CONV is a strict `xfail`
naming the CONV. The expected values come from formulas in the tests, not
from the code under test.

## Documentation

`docs/` is a Sphinx project: `index.rst` and one page per public object
(`source/*.rst`, autodoc), built with `uv run python docs/build.py` into
`docs/_build/html`. The docstrings are in NumPy style (napoleon), with units
and array shapes. Their examples are doctests, run with
`uv run pytest --doctest-modules qespace_cef --ignore=qespace_cef/symbolic`
(not part of the default `pytest` run). Known problems are marked in the
code as `# AUDIT Pn:` and in Warnings sections of the docstrings; the
findings themselves are in [`AUDIT.md`](AUDIT.md).

`docs/` documents the implementation; the decisions are kept apart in
`decisions/` (ADR-0001): design choices in `decisions/adr/` (`ADR-xxxx`),
physics conventions in `decisions/conventions/` (`CONV-xxxx`), one file per
decision, each from the `_template.md` in its folder. They are not part of
the Sphinx build. `AUDIT.md`, `PLAN.md` and `ROADMAP.md` track the work.
