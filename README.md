# qespace_cef

Crystal electric field (CEF) calculations for a magnetic rare-earth ion: build
the CEF Hamiltonian from Stevens operators, diagonalise it, and compute the
observables that follow from it (energy levels, neutron cross sections,
magnetisation, susceptibility, g-tensor, thermodynamics).

The theory follows A. T. Boothroyd, *Principles of Neutron Scattering from
Condensed Matter* (Oxford University Press, 2020),
DOI: [10.1093/oso/9780198862314.001.0001](https://doi.org/10.1093/oso/9780198862314.001.0001).

## Installation

With [`uv`](https://docs.astral.sh/uv/), from GitHub:

```bash
uv add git+https://github.com/qespace/qespace_cef@main
```

or from a local checkout, editable:

```bash
uv add --editable path/to/qespace_cef
```

Requires Python ≥ 3.12 and `numpy` ≥ 1.26; numpy 1 and 2 are both tested.

## Quick start

```python
import numpy as np
import qespace_cef as qc

ion = qc.Ion("Ce")  # Ce3+, J = 5/2
cfp = qc.CEFpars(dict(B20=0.5, B40=0.01, B44=0.1), units="meV")
cef = qc.CEFion(ion, cfp, Hfield=[0, 0, 0])  # field in T; diagonalised on creation

cef.energies  # eigenvalues of H in meV, not shifted (CONV-0010)
cef.excitation_energies  # the same, measured from the ground state
cef.get_degeneracies()  # [[energy, degeneracy], ...]: three Kramers doublets here

qc.gtensor(cef)  # g-tensor of the ground doublet
qc.magnetization(cef, 2.0, [0, 0, 1])  # T = 2 K, H = 1 T along z
Z, U, F, S, Cv = qc.thermodynamics(cef, np.linspace(1, 300, 300))  # T in K
energies, w = qc.cross_section_powder(cef, 10.0, [0.5, 1.0, 2.0], Ei=50.0)
# 10 K, |Q| in 1/A, Ei in meV: energies (36,) in meV, weights (3, 36) in barn/sr per ion
```

## What is in the package

| Name | What it is |
|---|---|
| `Ion` | Free rare-earth ion: `J`, Landé `gJ`, Stevens factors α, β, γ, magnetic form factor in the dipole approximation (`mff`; for Ce the coefficients are those of Ce²⁺) |
| `CEFpars` | Set of CEF parameters `Blm` (e.g. `B20`, `B44`, `B6m1` for m = −1); input in `meV`, `K` or `invcm`, stored in meV |
| `CEFion` | Ion in the CEF and an external field: Hamiltonian, `J` operators, energies, eigenvectors, degeneracies |
| `boltzman_population` | Thermal populations of the levels |
| `scattering_function` | The magnetic scattering function as excitations: energies and scattering weights `p_i <i\|J_α\|f><f\|J_β\|i>`, no neutron factors (CONV-0007, CONV-0012) |
| `cross_section_crystal`, `cross_section_powder` | Neutron cross section per ion as excitations: energies and weights in barn sr⁻¹, at a temperature, Q (a vector, or \|Q\| for a powder) and incident energy `Ei` (CONV-0007) |
| `magnetization`, `magnetization_exchange` | Magnetic moment in a field, optionally with a mean-field exchange term |
| `susceptibility` | Magnetic susceptibility tensor `χ_αβ` in μB/T per ion versus temperature (Van Vleck, CONV-0009) |
| `gtensor` | g-tensor of the ground state |
| `thermodynamics` | Partition function, internal energy, free energy, entropy, heat capacity |

`qespace_cef/symbolic/` holds experimental scripts that diagonalise the CEF
Hamiltonian symbolically with `sympy`; they are not part of the public API.
`sympy` is not installed with the package; in a checkout, get it with
`uv sync --group symbolic`.

## Conventions

The conventions are recorded in [`decisions/conventions/`](decisions/conventions/),
one file per convention. The CEF parameters and the Stevens operators follow
McPhase: $\mathcal{H} = \sum_{lm} B_{lm} \hat{O}_{lm}$, with the parameters
named `B20`, `B44`, `B4m3` (m = −3), and so on. The code works with energies in
meV, temperatures in K and magnetic fields in T.

## Documentation

The API documentation is built with Sphinx from the docstrings; the
sources are in `docs/`. Build it into `docs/_build/html` with:

```bash
uv run python docs/build.py
```

## Development

```bash
uv sync
```

```bash
uv run pytest
```

## Authors and licence

Michal Stekiel and Petr Čermák, based on earlier code by Petr Čermák, Jan
Zubáč and Karel Pajskr (2014–2018), distributed under the GNU General Public
License, version 3 or later.
