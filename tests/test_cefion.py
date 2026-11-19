"""Tests of ``qespace_cef.CEFion`` against analytic limits.

Only behaviour that AUDIT.md does not question is tested here; the known
problems (P1-P16) get their own tests when they are fixed.
"""

import numpy as np
import pytest

import qespace_cef as qc
from qespace_cef import constants as C

KRAMERS_IONS = ["Ce", "Nd", "Sm", "Dy", "Er", "Yb"]  # odd number of 4f electrons
ALL_IONS = ["Ce", "Pr", "Nd", "Sm", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"]
LOW_SYMMETRY = dict(B20=0.3, B22=-0.1, B40=0.01, B4m3=0.02, B44=0.03, B60=1e-4, B6m6=2e-4)


def _cefion(name, pars, field=(0, 0, 0)):
    return qc.CEFion(qc.Ion(name), qc.CEFpars(pars), Hfield=list(field))


@pytest.mark.parametrize("name", ALL_IONS)
def test_free_ion_is_degenerate(name):
    """Without CEF and field all 2J+1 states are at zero energy."""
    cef = _cefion(name, {})
    assert len(cef.energies) == int(2 * cef.ion.J + 1)
    assert np.allclose(cef.energies, 0)


def test_axial_levels():
    """B20 alone: E(m) = B20 (3 m^2 - J(J+1)), here Ce (J = 5/2), B20 = 1 meV."""
    cef = _cefion("Ce", dict(B20=1.0))
    expected = np.array([0, 0, 6, 6, 18, 18], dtype=float)  # m = ±1/2, ±3/2, ±5/2
    assert np.allclose(cef.energies, expected, atol=1e-12 * expected.max())


def test_cubic_doublet_quartet():
    """Cubic CEF on Ce (B44 = 5 B40) splits J = 5/2 into a doublet and a quartet."""
    cef = _cefion("Ce", dict(B40=0.1, B44=0.5))
    e = cef.energies
    assert np.allclose(e[:2], e[0]) and np.allclose(e[2:], e[2])
    assert e[2] - e[0] > 1e-3 * e[2]


@pytest.mark.parametrize("name", KRAMERS_IONS)
def test_kramers_degeneracy(name):
    """Without a field every level of a Kramers ion is (at least) doubly degenerate."""
    e = _cefion(name, LOW_SYMMETRY).energies
    assert np.allclose(e[0::2], e[1::2], atol=1e-9 * e.max())


@pytest.mark.parametrize("name", ALL_IONS)
def test_zeeman_levels_equally_spaced(name):
    """Field only: the levels are spaced by gJ muB B (spacing only, no sign)."""
    B = 2.0  # T
    cef = _cefion(name, {}, field=(0, 0, B))
    spacing = np.diff(cef.energies)
    assert np.allclose(spacing, cef.ion.gJ * C.uB * B, rtol=1e-10)


@pytest.mark.parametrize("name", ["Ce", "Tb", "Er"])
def test_energies_sorted_from_zero(name):
    """Energies are sorted ascending and the ground state is at zero."""
    e = _cefion(name, LOW_SYMMETRY, field=(0.1, 0.2, 0.3)).energies
    assert e[0] == 0 and np.all(np.diff(e) >= 0)


@pytest.mark.parametrize("name", ["Ce", "Tb", "Er"])
def test_eigenvectors_orthonormal(name):
    """The eigenvectors (rows of ``eigenvectors``) are orthonormal."""
    v = _cefion(name, LOW_SYMMETRY, field=(0.1, 0.2, 0.3)).eigenvectors
    assert np.allclose(v.conj() @ v.T, np.eye(len(v)), atol=1e-12)


@pytest.mark.parametrize("name", ALL_IONS)
def test_saturated_moment_in_strong_field(name):
    """Free ion in a field: the ground state carries a moment of size gJ J."""
    cef = _cefion(name, {}, field=(0, 0, 1.0))
    assert np.isclose(np.linalg.norm(cef.magnetic_moment[0]), cef.ion.gJ * cef.ion.J)
