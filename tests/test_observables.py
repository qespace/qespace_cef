"""Tests of ``qespace_cef.observables``.

The expected values come from formulas written out here, not from the code
under test. A known problem from AUDIT.md is a strict ``xfail`` naming its
ID: it turns into a failure (XPASS) when the problem is fixed, so the marker
must then be removed.
"""

import numpy as np
import pytest

import qespace_cef as qc

# CODATA 2018: k_B = 8.617333262e-5 eV/K, mu_B = 5.7883818060e-5 eV/T.
KB = 0.08617333262  # meV/K
MUB = 0.057883818060  # meV/T

TETRAGONAL = dict(B20=0.5, B40=0.01, B44=0.1)  # Ce3+, as in AUDIT P1, P2, P10
T = np.array([5.0, 20.0, 100.0])  # K


def _ce(pars=TETRAGONAL):
    return qc.CEFion(qc.Ion("Ce"), qc.CEFpars(pars))


def _chi_exact(cef, temperature, direction):
    """Zero-field susceptibility along ``direction``, in muB^2/meV per ion.

    chi = gJ^2 [ sum_{E_m = E_n} p_n |<n|J_n|m>|^2 / k_B T
               + sum_{E_m != E_n} (p_n - p_m) / (E_m - E_n) |<n|J_n|m>|^2 ],
    with J_n the component of J along the unit vector of ``direction``, from
    the energies and the J matrices in the eigenbasis.
    """
    n = np.asarray(direction, float) / np.linalg.norm(direction)
    Jn = np.einsum("a,aij->ij", n, cef.J)
    E = cef.energies
    p = np.exp(-(E - E.min()) / (KB * temperature))
    p /= p.sum()
    dE = E[None, :] - E[:, None]  # E_m - E_n, n rows
    degenerate = np.abs(dE) <= 1e-9 * max(np.abs(E).max(), 1.0)
    weight = np.where(
        degenerate,
        p[:, None] / (KB * temperature),
        (p[:, None] - p[None, :]) / np.where(degenerate, 1.0, dE),
    )
    return cef.ion.gJ**2 * float(np.sum(weight * np.abs(Jn) ** 2))


def test_magnetization_returns_python_floats():
    """``magnetization`` returns a tuple of three Python floats (CLAUDE.md 5.7.3)."""
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)))
    M = qc.magnetization(cef, 2.0, [0.0, 0.0, 1.0])
    assert len(M) == 3 and all(type(m) is float for m in M)


# --- susceptibility -------------------------------------------------------


@pytest.mark.parametrize("temperature", T)
def test_susceptibility_magnetization_method_along_z(temperature):
    """``method='magnetization'`` along z equals the exact formula, in muB/T.

    Along a symmetry axis M is parallel to H, so P10 does not show here; the
    test checks the exact formula as much as the code.
    """
    chi = qc.susceptibility(_ce(), [temperature], [0, 0, 1], method="magnetization")[0]
    assert chi == pytest.approx(_chi_exact(_ce(), temperature, [0, 0, 1]) * MUB, rel=1e-5)


@pytest.mark.xfail(strict=True, reason="AUDIT P1: the van Vleck term is halved")
@pytest.mark.parametrize("temperature", T)
def test_susceptibility_perturbation_along_z(temperature):
    """``method='perturbation'`` along z equals the exact formula, in muB^2/meV."""
    chi = qc.susceptibility(_ce(), [temperature], [0, 0, 1])[0]
    assert chi == pytest.approx(_chi_exact(_ce(), temperature, [0, 0, 1]), rel=1e-6)


@pytest.mark.xfail(strict=True, reason="AUDIT P2: the field direction is not projected")
def test_susceptibility_perturbation_tetragonal_in_plane():
    """Tetragonal symmetry: chi is the same along (1,0,0) and (1,1,0)."""
    chi_100 = qc.susceptibility(_ce(), T, [1, 0, 0])
    chi_110 = qc.susceptibility(_ce(), T, [1, 1, 0])
    assert np.allclose(chi_110, chi_100, rtol=1e-10)


@pytest.mark.xfail(strict=True, reason="AUDIT P10: returns |M|/eps, not M.n/eps")
def test_susceptibility_magnetization_method_off_axis():
    """Off a symmetry axis, ``'magnetization'`` gives n.chi.n, in muB/T."""
    chi = qc.susceptibility(_ce(), [20.0], [1, 0, 1], method="magnetization")[0]
    assert chi == pytest.approx(_chi_exact(_ce(), 20.0, [1, 0, 1]) * MUB, rel=1e-5)


# --- thermodynamics -------------------------------------------------------
# Identities that hold in any consistent set of units, so they do not choose
# the units of the returns (C1, C11).


def _thermo(temperatures):
    return qc.thermodynamics(_ce(), np.asarray(temperatures, float))


def _derivative(i, temperature):
    """Central difference d(return value i)/dT at ``temperature``."""
    h = 1e-4 * temperature
    hi, lo = _thermo([temperature + h])[i][0], _thermo([temperature - h])[i][0]
    return (hi - lo) / (2 * h)


def test_partition_function_and_energy():
    """Z = sum exp(-E/kT) and <E> = sum E exp(-E/kT) / Z, energies in meV."""
    E = _ce().energies
    Z, U = _thermo(T)[:2]
    w = np.exp(-E[None, :] / (KB * T[:, None]))
    assert np.allclose(Z, w.sum(1), rtol=1e-5)
    assert np.allclose(U, (w * E).sum(1) / w.sum(1), rtol=1e-5)


@pytest.mark.xfail(strict=True, reason="AUDIT P3: meV and K mixed in S")
@pytest.mark.parametrize("temperature", T)
def test_entropy_is_minus_dF_dT(temperature):
    """S = -dF/dT."""
    S = _thermo([temperature])[3][0]
    assert S == pytest.approx(-_derivative(2, temperature), rel=1e-4)


@pytest.mark.xfail(strict=True, reason="AUDIT P3: meV and K mixed in Cv")
@pytest.mark.parametrize("temperature", T)
def test_heat_capacity_is_dE_dT(temperature):
    """Cv = dE/dT."""
    Cv = _thermo([temperature])[4][0]
    assert Cv == pytest.approx(_derivative(1, temperature), rel=1e-4)


@pytest.mark.xfail(strict=True, reason="AUDIT P3: meV and K mixed in S and Cv")
@pytest.mark.parametrize("temperature", T)
def test_heat_capacity_is_T_dS_dT(temperature):
    """Cv = T dS/dT."""
    Cv = _thermo([temperature])[4][0]
    assert Cv == pytest.approx(temperature * _derivative(3, temperature), rel=1e-4)
