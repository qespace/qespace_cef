"""Tests of ``qespace_cef.observables``.

The expected values come from formulas written out here, not from the code
under test. A known problem from AUDIT.md is a strict ``xfail`` naming its
ID: it turns into a failure (XPASS) when the problem is fixed, so the marker
must then be removed.
"""

import copy
import warnings

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


# --- spectral weight ------------------------------------------------------

LOW_SYMMETRY = dict(B20=0.3, B22=-0.1, B40=0.01, B4m3=0.02, B44=0.03)
TEMPERATURE = 10.0  # K


def _weights_without_kfki(cef, Q, scheme):
    """Sqw / (kf/ki), with kf/ki = sqrt(1 - (E_j - E_i)/Ei) for the default Ei."""
    jumps, Sqw = qc.spectral_weight(cef, TEMPERATURE, Q, scheme)
    return Sqw / np.sqrt(1 - jumps / 1e6)


def _populations(cef):
    p = np.exp(-(cef.energies - cef.energies.min()) / (KB * TEMPERATURE))
    return p / p.sum()


def test_powder_sum_rule():
    """sum_f |<i|J_perp|f>|^2 = 2/3 J(J+1) for every initial state i (powder)."""
    cef = _ce(LOW_SYMMETRY)
    Q = [1.0, 0, 0]
    rows = _weights_without_kfki(cef, Q, "powder").sum(axis=1)
    norm = cef.ion.mff(1.0) ** 2 * _populations(cef)
    J = cef.ion.J
    assert np.allclose(rows / norm, 2 / 3 * J * (J + 1), rtol=1e-5)


def test_single_crystal_trace_sum_rule():
    """sum_if |<i|J_perp|f>|^2 = 2/3 J(J+1)(2J+1), the trace over i (single crystal)."""
    cef = _ce(LOW_SYMMETRY)
    Q = np.array([0.3, 0.5, 0.7])
    W = _weights_without_kfki(cef, Q, "single-crystal")
    norm = cef.ion.mff(np.linalg.norm(Q)) ** 2 * _populations(cef)
    J = cef.ion.J
    assert np.sum(W / norm[:, None]) == pytest.approx(2 / 3 * J * (J + 1) * (2 * J + 1), rel=1e-5)


@pytest.mark.xfail(strict=True, reason="AUDIT P9: the cross terms of the projector are missing")
@pytest.mark.parametrize("Q", [[1.0, 1.0, 0.0], [0.3, 0.5, 0.7]])
def test_single_crystal_weights_full_projector(Q):
    """Sqw_if = f^2 p_i sum_ab (delta_ab - Qa Qb / Q^2) <i|J_a|f><f|J_b|i> (kf/ki divided out)."""
    cef = _ce(LOW_SYMMETRY)
    q = np.asarray(Q) / np.linalg.norm(Q)
    projector = np.eye(3) - np.outer(q, q)
    M = np.einsum("ab,aif,bfi->if", projector, cef.J, cef.J).real
    expected = cef.ion.mff(np.linalg.norm(Q)) ** 2 * _populations(cef)[:, None] * M
    W = _weights_without_kfki(cef, Q, "single-crystal")
    assert np.allclose(W, expected, rtol=1e-5, atol=1e-8 * expected.max())


# --- g-tensor -------------------------------------------------------------


def _principal_g(cef):
    g = qc.gtensor(cef)
    return np.sort(np.sqrt(np.linalg.eigvalsh(g @ g.T)))


def _rotate_ground_doublet(cef, axis, angle):
    """A copy of ``cef`` whose ground-doublet basis is rotated by exp(-i angle n.sigma/2)."""
    n = np.asarray(axis, float) / np.linalg.norm(axis)
    sigma = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    u = np.cos(angle / 2) * np.eye(2) - 1j * np.sin(angle / 2) * sum(
        c * s for c, s in zip(n, sigma)
    )
    V = np.eye(len(cef.energies), dtype=complex)
    V[:2, :2] = u
    rotated = copy.copy(cef)
    rotated.J = np.einsum("ij,ajk,kl->ail", V.conj().T, cef.J, V)
    return rotated


AXIAL_CE = dict(B20=1.0)  # ground doublet |+-1/2>: g_par = gJ, g_perp = 3 gJ


def test_gtensor_principal_values_axial():
    """Ce3+ with B20 > 0: principal g values gJ, 3 gJ, 3 gJ in the basis eigh picks."""
    cef = _ce(AXIAL_CE)
    gJ = cef.ion.gJ
    assert np.allclose(_principal_g(cef), [gJ, 3 * gJ, 3 * gJ], rtol=1e-10)


@pytest.mark.xfail(strict=True, reason="AUDIT P17: abs() makes g g^T depend on the doublet basis")
@pytest.mark.parametrize(("axis", "angle"), [([1, 0, 0], 1.0), ([1, 1, 1], 2.0)])
def test_gtensor_principal_values_any_doublet_basis(axis, angle):
    """The principal g values do not depend on the basis inside the ground doublet."""
    cef = _ce(AXIAL_CE)
    gJ = cef.ion.gJ
    rotated = _rotate_ground_doublet(cef, axis, angle)
    assert np.allclose(_principal_g(rotated), [gJ, 3 * gJ, 3 * gJ], rtol=1e-10)


# --- magnetization with a molecular field ---------------------------------


class _Runaway(Exception):
    """Raised by the stand-in below when the iteration has not stopped."""


@pytest.mark.xfail(strict=True, reason="AUDIT P11: no iteration limit")
def test_magnetization_exchange_stops_without_convergence(monkeypatch):
    """A mean field that never converges ends in an error, not in an endless loop.

    ``magnetization`` is replaced by a stand-in that alternates between two
    moments, so the iteration cannot converge; it raises ``_Runaway`` after
    10000 calls, which is what happens while there is no iteration limit.
    """
    calls = []

    def alternating(cefion, temperature, Hfield):
        calls.append(1)
        if len(calls) > 10_000:
            raise _Runaway
        return (0.0, 0.0, 1.0 if len(calls) % 2 else -1.0)

    monkeypatch.setattr(qc.observables, "magnetization", alternating)
    with pytest.raises(Exception) as error:
        qc.magnetization_exchange(_ce(), 5.0, [0, 0, 1.0], lam=1.0)
    assert not isinstance(error.value, _Runaway)


@pytest.mark.xfail(strict=True, raises=RuntimeWarning, reason="AUDIT P11: 0/0 when M = 0")
def test_magnetization_exchange_zero_field():
    """In zero field the moment is zero, without a division 0/0 on the way."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        M = qc.magnetization_exchange(_ce(), 5.0, [0.0, 0.0, 0.0], lam=1.0)
    assert M == (0.0, 0.0, 0.0)
