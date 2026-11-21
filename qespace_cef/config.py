"""The one config place of qespace_cef: numerical tolerances (AUDIT A11).

The functions of the package read these values when they are called, not when
they are imported, so an assignment changes every later call, e.g.
``qespace_cef.config.rtol = 1e-8``.
"""

rtol = 1e-10
"""Relative tolerance of CONV-0008, dimensionless: two energies are one level
when they differ by at most ``rtol * max(max|E|, E_floor)``."""

E_floor = 1.0
"""Energy scale of CONV-0008 in meV: the tolerance never falls below
``rtol * E_floor``, e.g. for a Hamiltonian H = 0."""
