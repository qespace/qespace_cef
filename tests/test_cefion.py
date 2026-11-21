"""Tests of ``qespace_cef.CEFion`` against analytic limits.

A known problem from AUDIT.md is a strict ``xfail`` naming its ID: it turns
into a failure (XPASS) when the problem is fixed, so the marker must then be
removed.
"""

import numpy as np
import pytest

import qespace_cef as qc
from qespace_cef import constants as C
from qespace_cef.cefmatrices import StevensBase

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
    expected = np.array([-8, -8, -2, -2, 10, 10], dtype=float)  # m = ±1/2, ±3/2, ±5/2
    assert np.allclose(cef.energies, expected, atol=1e-12 * expected.max())
    excitations = expected - expected[0]
    assert np.allclose(cef.excitation_energies, excitations, atol=1e-12 * expected.max())


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
def test_energies_ascending_excitations_from_zero(name):
    """Energies are ascending; the excitation energies start at zero."""
    cef = _cefion(name, LOW_SYMMETRY, field=(0.1, 0.2, 0.3))
    assert np.all(np.diff(cef.energies) >= 0)
    assert cef.excitation_energies[0] == 0 and np.all(np.diff(cef.excitation_energies) >= 0)


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


@pytest.mark.parametrize("name", ["Ce", "Tb"])
def test_diagonalize_returns_decomposition(name):
    """``diagonalize`` returns E diagonal and U unitary with H = U E U^dagger."""
    cef = qc.CEFion(
        qc.Ion(name), qc.CEFpars(LOW_SYMMETRY), Hfield=[0.1, 0.2, 0.3], diagonalize=False
    )
    E, U = cef.diagonalize()
    scale = np.max(np.abs(cef.hamiltonian))
    assert np.allclose(E, np.diag(np.diag(E)))
    assert np.allclose(U.conj().T @ U, np.eye(len(U)), atol=1e-12)
    assert np.allclose(U @ E @ U.conj().T, cef.hamiltonian, atol=1e-12 * scale)


def test_diagonalize_rejects_non_hermitian_hamiltonian():
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)), diagonalize=False)
    cef.hamiltonian = cef.hamiltonian + 1j * np.eye(len(cef.hamiltonian))
    with pytest.raises(ValueError):
        cef.diagonalize()


def test_degeneracies_are_plain_python_numbers():
    """``get_degeneracies`` returns Python float and int, not NumPy scalars (CLAUDE.md 5.7.3)."""
    levels = _cefion("Ce", dict(B20=1.0)).get_degeneracies()
    assert levels == [[-8.0, 2], [-2.0, 2], [10.0, 2]]
    assert all(type(e) is float and type(d) is int for e, d in levels)


def test_summary_lists_levels_from_the_ground_state():
    """``summary`` prints excitation energies with degeneracies, and E_0 in its header."""
    text = _cefion("Ce", dict(B20=1.0)).summary()
    assert "(E_0 = -8.000 meV)" in text
    for level, energy in enumerate(["0.000", "6.000", "18.000"]):
        assert f"E({level}) =\t{energy} meV\t 2-fold degenerate" in text
    assert text.count("ev_") == 6  # one line per state


def test_summary_options():
    """``precision`` sets the decimals; ``evecs=False`` leaves out the eigenvectors."""
    text = _cefion("Ce", dict(B20=1.0)).summary(evecs=False, precision=1)
    assert "E(2) =\t18.0 meV" in text and "ev_" not in text


def test_summary_groups_levels_for_printing():
    """States within 1e-10 meV print as one level; a larger splitting does not.

    Cubic Pr3+: eigh returns the 28.584 meV doublet as two values that differ
    by rounding noise (AUDIT P7); the print groups them. Ce3+ in 1e-6 T splits
    each doublet by about 5e-8 meV, which the print keeps apart.
    """
    pr = _cefion("Pr", CUBIC_PR).summary(evecs=False)
    assert [line.split("\t")[-1] for line in pr.splitlines() if line.startswith("E(")] == [
        " 3-fold degenerate",
        " 1-fold degenerate",
        " 3-fold degenerate",
        " 2-fold degenerate",
    ]
    ce = _cefion("Ce", dict(B20=1.0), field=(0, 0, 1e-6)).summary(evecs=False)
    assert ce.count(" 1-fold degenerate") == 6


def test_summary_before_diagonalisation():
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)), diagonalize=False)
    assert cef.summary() == "Undiagonalized Hamiltonian"


# --- CONV-0010: the zero of energy ----------------------------------------
# Expected values from formulas in the tests (CONV-0010, Consequences).

# meV/T, CODATA 2022 (NIST): mu_B = 9.2740100657(29)e-24 J/T, e = 1.602176634e-19 C exact.
MUB_CODATA = 1e3 * 9.2740100657e-24 / 1.602176634e-19


@pytest.mark.parametrize("name", ["Ce", "Ho"])
@pytest.mark.parametrize("field", [(0, 0, 0), (0.3, 0.0, 1.0)])
def test_conv_0010_energies_are_eigenvalues_centred_on_zero(name, field):
    """``energies`` are the eigenvalues of ``hamiltonian``; with B00 = 0 they average to 0."""
    cef = _cefion(name, LOW_SYMMETRY, field=field)
    scale = np.abs(cef.energies).max()
    assert np.allclose(
        cef.energies, np.linalg.eigvalsh(cef.hamiltonian), rtol=0, atol=1e-12 * scale
    )
    assert abs(np.mean(cef.energies)) <= 1e-12 * scale


def test_conv_0010_B00_shifts_every_level():
    """B00 multiplies the identity, so it moves every energy by B00."""
    shift = 0.7  # meV, a test value
    plain = _cefion("Ce", LOW_SYMMETRY).energies
    shifted = _cefion("Ce", dict(LOW_SYMMETRY, B00=shift)).energies
    assert np.allclose(shifted - plain, shift, rtol=0, atol=1e-12 * np.abs(plain).max())


@pytest.mark.parametrize("B", [0.5, 2.0])  # T
def test_conv_0010_free_ion_ground_state_in_field(B):
    """Free Ho3+ in a field along z: E_0 = -gJ muB J B, measured from the free-ion energy."""
    cef = _cefion("Ho", {}, field=(0, 0, B))
    expected = -cef.ion.gJ * MUB_CODATA * cef.ion.J * B
    assert cef.energies[0] == pytest.approx(expected, rel=1e-6)


def test_conv_0010_excitation_energies():
    """``excitation_energies`` is ``energies - energies[0]`` and starts at 0."""
    cef = _cefion("Er", LOW_SYMMETRY, field=(0.1, 0.2, 0.3))
    assert np.array_equal(cef.excitation_energies, cef.energies - cef.energies[0])
    assert cef.excitation_energies[0] == 0


def test_conv_0010_hellmann_feynman():
    """dE_n/dB40 = <n|O_40|n> for the unshifted energies, against a finite difference.

    Low-symmetry Ce3+ in a field, so no level is degenerate. With energies
    shifted to the ground state the derivative would gain -<0|O_40|0>.
    """
    field = (0.3, 0.0, 1.0)
    cef = _cefion("Ce", LOW_SYMMETRY, field=field)
    O40 = StevensBase(cef.Jval)["40"]
    expected = np.einsum("ni,ij,nj->n", cef.eigenvectors.conj(), O40, cef.eigenvectors).real
    h = 1e-6  # meV
    up = _cefion("Ce", dict(LOW_SYMMETRY, B40=LOW_SYMMETRY["B40"] + h), field=field).energies
    down = _cefion("Ce", dict(LOW_SYMMETRY, B40=LOW_SYMMETRY["B40"] - h), field=field).energies
    assert np.min(np.diff(cef.energies)) > 1e-3  # no degenerate level
    np.testing.assert_allclose(
        (up - down) / (2 * h), expected, rtol=0, atol=1e-5 * np.abs(expected).max()
    )


# --- known problems (AUDIT.md), strict xfail until fixed ------------------


def test_second_diagonalize_changes_nothing():
    """Calling ``diagonalize`` again leaves energies and moments unchanged."""
    cef = _cefion("Ce", LOW_SYMMETRY, field=(0.3, 0, 1))
    energies, moments = cef.energies.copy(), cef.magnetic_moment.copy()
    cef.diagonalize()
    assert np.allclose(cef.energies, energies)
    assert np.allclose(cef.magnetic_moment, moments, atol=1e-10)


@pytest.mark.parametrize("name", ["Ce", "Tb", "Er"])
def test_eigenvectors_belong_to_energies(name):
    """Row i of ``eigenvectors`` is the eigenvector of ``energies[i]``: H v_i = E_i v_i."""
    cef = _cefion(name, LOW_SYMMETRY, field=(0.1, 0.2, 0.3))
    H = cef.hamiltonian
    scale = np.max(np.abs(H))
    for energy, v in zip(cef.energies, cef.eigenvectors):
        assert np.allclose(H @ v, energy * v, atol=1e-10 * scale)


CUBIC_PR = dict(B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4)  # B44 = 5 B40, B64 = -21 B60


@pytest.mark.xfail(strict=True, reason="AUDIT P7: energies truncated, a doublet split in two")
def test_degeneracies_of_cubic_pr():
    """Pr3+ in a cubic field: levels grouped within 1e-9 of the largest energy."""
    e = _cefion("Pr", CUBIC_PR).energies
    new_level = np.diff(e) > 1e-9 * e.max()
    expected = np.diff(np.flatnonzero(np.r_[True, new_level, True])).tolist()
    assert [d for _, d in _cefion("Pr", CUBIC_PR).get_degeneracies()] == expected
