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
from qespace_cef.cefion import degenerate_levels

# CODATA 2022 (NIST, physics.nist.gov/cuu/Constants): e = 1.602176634e-19 C and
# k = 1.380649e-23 J/K exact, mu_B = 9.2740100657(29)e-24 J/T.
KB = 1e3 * 1.380649e-23 / 1.602176634e-19  # meV/K
MUB = 1e3 * 9.2740100657e-24 / 1.602176634e-19  # meV/T

TETRAGONAL = dict(B20=0.5, B40=0.01, B44=0.1)  # Ce3+, as in AUDIT P1, P2, P10
T = np.array([5.0, 20.0, 100.0])  # K


def _ce(pars=TETRAGONAL):
    return qc.CEFion(qc.Ion("Ce"), qc.CEFpars(pars))


def _chi_exact(cef, temperature, direction, tolerance=None):
    """Zero-field susceptibility along ``direction``, in muB^2/meV per ion.

    chi = gJ^2 [ sum_{E_m = E_n} p_n |<n|J_n|m>|^2 / k_B T
               + sum_{E_m != E_n} (p_n - p_m) / (E_m - E_n) |<n|J_n|m>|^2 ],
    with J_n the component of J along the unit vector of ``direction``, from
    the energies and the J matrices in the eigenbasis. E_m = E_n when they
    differ by at most ``tolerance`` (meV), by default 1e-9 max(max|E|, 1 meV).
    """
    n = np.asarray(direction, float) / np.linalg.norm(direction)
    Jn = np.einsum("a,aij->ij", n, cef.J)
    E = cef.energies
    p = np.exp(-(E - E.min()) / (KB * temperature))
    p /= p.sum()
    dE = E[None, :] - E[:, None]  # E_m - E_n, n rows
    if tolerance is None:
        tolerance = 1e-9 * max(np.abs(E).max(), 1.0)
    degenerate = np.abs(dE) <= tolerance
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


def _along(chi, direction):
    """n.chi.n for the unit vector of ``direction``; chi of shape (..., 3, 3)."""
    n = np.asarray(direction, float) / np.linalg.norm(direction)
    return np.einsum("a,...ab,b->...", n, chi, n)


@pytest.mark.parametrize("temperature", T)
def test_susceptibility_magnetization_method_along_z(temperature):
    """``method='magnetization'``: chi_zz equals the exact formula, in muB/T."""
    chi = qc.susceptibility(_ce(), [temperature], method="magnetization")[0]
    assert chi[2, 2] == pytest.approx(_chi_exact(_ce(), temperature, [0, 0, 1]) * MUB, rel=1e-5)


@pytest.mark.parametrize("temperature", T)
def test_susceptibility_perturbation_along_z(temperature):
    """``method='perturbation'``: chi_zz equals the exact formula, in muB/T."""
    chi = qc.susceptibility(_ce(), [temperature])[0]
    assert chi[2, 2] == pytest.approx(_chi_exact(_ce(), temperature, [0, 0, 1]) * MUB, rel=1e-6)


@pytest.mark.parametrize("direction", [[1, 0, 0], [1, 1, 0], [1, 0, 1]])
def test_susceptibility_perturbation_off_axis(direction):
    """``method='perturbation'``: n.chi.n equals the exact formula off the axes (P2)."""
    chi = qc.susceptibility(_ce(), [20.0])[0]
    assert _along(chi, direction) == pytest.approx(
        _chi_exact(_ce(), 20.0, direction) * MUB, rel=1e-6
    )


def test_susceptibility_perturbation_tetragonal_in_plane():
    """Tetragonal symmetry: chi is the same along (1,0,0) and (1,1,0)."""
    chi = qc.susceptibility(_ce(), T)
    assert np.allclose(_along(chi, [1, 1, 0]), _along(chi, [1, 0, 0]), rtol=1e-10)


def test_susceptibility_magnetization_method_off_axis():
    """Off a symmetry axis, ``'magnetization'`` gives n.chi.n, in muB/T (P10)."""
    chi = qc.susceptibility(_ce(), [20.0], method="magnetization")[0]
    assert _along(chi, [1, 0, 1]) == pytest.approx(
        _chi_exact(_ce(), 20.0, [1, 0, 1]) * MUB, rel=1e-5
    )


@pytest.mark.parametrize("method", ["perturbation", "magnetization"])
@pytest.mark.parametrize(("temperatures", "shape"), [(5.0, (3, 3)), ([[5.0, 20.0]], (1, 2, 3, 3))])
def test_susceptibility_shape(method, temperatures, shape):
    """Temperatures of shape S give a tensor of shape S + (3, 3) (CONV-0009)."""
    assert qc.susceptibility(_ce(), temperatures, method=method).shape == shape


def test_susceptibility_prints_nothing(capsys):
    """A computation writes nothing to stdout (CLAUDE.md 5.8)."""
    qc.susceptibility(_ce(), T)
    assert capsys.readouterr().out == ""


def _level_tolerance(energies):
    """The tolerance of CONV-0008 with its defaults rtol = 1e-10, E_floor = 1 meV, in meV."""
    return 1e-10 * max(np.abs(energies).max(), 1.0)


# The field along x splits the ground doublet by about 6e-12, 6e-11, 6e-9 and
# 6e-8 meV; the tolerance is 6.4e-10 meV.
@pytest.mark.parametrize(
    ("field", "below"), [(1e-10, True), (1e-9, True), (1e-7, False), (1e-6, False)]
)  # T
def test_susceptibility_continuous_through_degeneracy(field, below):
    """chi is continuous as a splitting shrinks through the tolerance (CONV-0008, P20).

    A field along x splits the ground Kramers doublet of the tetragonal Ce3+
    ion. Below the tolerance of CONV-0008 the pair is one level and takes the
    Curie weight p_n / k_B T; above it, the van Vleck weight
    (p_n - p_m) / (E_m - E_n), whose limit is the Curie weight. The change of
    chi_zz with the field is of second order, about (muB B / k_B T)^2 < 1e-13
    here, so on both sides chi_zz must equal its zero-field value.
    """
    temperature = 5.0  # K
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(TETRAGONAL), Hfield=[field, 0.0, 0.0])
    split = cef.energies[1] - cef.energies[0]
    assert (split <= _level_tolerance(cef.energies)) == below  # the side the case is about
    assert split < 1e-6 * KB * temperature
    chi_zz = qc.susceptibility(cef, temperature)[2, 2]
    assert chi_zz == pytest.approx(qc.susceptibility(_ce(), temperature)[2, 2], rel=1e-12)


def test_conv_0008_susceptibility_levels_from_config(monkeypatch):
    """The Curie and van Vleck pairs follow the levels of CONV-0008, read at call time.

    With rtol = 0.5 the tolerance, 0.5 max|E| = 3.2 meV, merges the two lower
    doublets of the tetragonal Ce3+ ion (0.78 meV apart) into one level, so
    their pairs take the Curie weight. J_x couples the two doublets (J_z does
    not), so chi_xx sees it.
    """
    temperature = 5.0  # K
    cef = _ce()
    tolerance = 0.5 * np.abs(cef.energies).max()  # meV
    assert np.diff(cef.energies)[1] < tolerance < np.diff(cef.energies)[3]
    monkeypatch.setattr(qc.config, "rtol", 0.5)
    default = _chi_exact(cef, temperature, [1, 0, 0]) * MUB
    expected = _chi_exact(cef, temperature, [1, 0, 0], tolerance=tolerance) * MUB
    assert abs(expected / default - 1) > 1e-2  # the merge changes chi_xx
    chi = qc.susceptibility(cef, temperature)
    assert chi[0, 0] == pytest.approx(expected, rel=1e-6)


def test_susceptibility_low_temperature_finite():
    """At 0.05 K, where (E_m - E_n) / k_B T reaches 10^4, chi is finite and warns of nothing.

    Pairs with E_m < E_n would overflow in exp((E_n - E_m) / k_B T), so the
    formula must never evaluate it (P20).
    """
    temperature = 0.05  # K
    cef = _ce()
    assert (cef.energies[-1] - cef.energies[0]) / (KB * temperature) > 1e3  # exp(+x) overflows
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        chi = qc.susceptibility(cef, temperature)
    assert np.all(np.isfinite(chi))
    assert chi[2, 2] == pytest.approx(_chi_exact(cef, temperature, [0, 0, 1]) * MUB, rel=1e-6)


# --- CONV-0009: the susceptibility tensor ---------------------------------
# Expected values from formulas in the tests (CONV-0009, Consequences).

CHI_LOW_SYMMETRY = dict(B20=0.5, B22=0.2, B40=0.01, B43=0.05, B4m2=0.03)  # Ce3+, meV, test values


def _curie(ion, temperature):
    """Free-ion Curie law gJ^2 muB J(J+1) / (3 k_B T), in muB/T per ion."""
    return ion.gJ**2 * MUB * ion.J * (ion.J + 1) / (3 * KB * temperature)


@pytest.mark.parametrize(("method", "rel"), [("perturbation", 1e-6), ("magnetization", 1e-5)])
def test_conv_0009_free_ion_curie_law(method, rel):
    """H_CEF = 0: chi_ab = delta_ab gJ^2 muB J(J+1) / (3 k_B T)."""
    ion = qc.Ion("Ho")
    cef = qc.CEFion(ion, qc.CEFpars(dict(B20=0.0)))
    for temperature in [2.0, 50.0]:
        chi = qc.susceptibility(cef, temperature, method=method)
        expected = _curie(ion, temperature) * np.eye(3)
        np.testing.assert_allclose(chi, expected, rtol=0, atol=rel * expected.max())


def test_conv_0009_high_temperature_limit():
    """Tr chi / 3 approaches the free-ion Curie law as T grows."""
    ion = qc.Ion("Ce")
    cef = qc.CEFion(ion, qc.CEFpars(CHI_LOW_SYMMETRY))
    T_high = np.array([1e2, 1e3, 1e4])
    deviation = np.abs(
        np.trace(qc.susceptibility(cef, T_high), axis1=-2, axis2=-1) / 3 / _curie(ion, T_high) - 1
    )
    assert deviation[-1] < 1e-4 and np.all(np.diff(deviation) < 0)


@pytest.mark.parametrize("temperature", [2.0, 20.0, 200.0])
def test_conv_0009_perturbation_equals_finite_difference(temperature):
    """The formula and the finite difference M(eps e_b)/eps agree, off-diagonal terms included."""
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(CHI_LOW_SYMMETRY))
    chi = qc.susceptibility(cef, temperature)
    chi_fd = qc.susceptibility(cef, temperature, method="magnetization")
    off_diagonal = chi[~np.eye(3, dtype=bool)]
    assert np.abs(off_diagonal).max() > 1e-3 * np.abs(chi).max()  # the test sees them
    np.testing.assert_allclose(chi_fd, chi, rtol=0, atol=1e-5 * np.abs(chi).max())


def test_conv_0009_symmetric_positive_semidefinite():
    """chi is real, symmetric and positive semi-definite."""
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(CHI_LOW_SYMMETRY))
    for chi in qc.susceptibility(cef, [1.0, 10.0, 100.0]):
        assert np.allclose(chi, chi.T, rtol=0, atol=1e-14 * np.abs(chi).max())
        assert np.linalg.eigvalsh(chi).min() >= -1e-14 * np.abs(chi).max()


def test_conv_0009_invariant_inside_degenerate_levels():
    """chi is unchanged when the basis inside each Kramers doublet is rotated."""
    cef = _ce()
    chi = qc.susceptibility(cef, T)
    rng = np.random.default_rng(1)
    n = len(cef.energies)
    U = np.zeros((n, n), complex)
    for k in range(0, n, 2):  # Ce3+: three Kramers doublets
        Q, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        U[k : k + 2, k : k + 2] = Q
    cef.J = np.einsum("ki,akl,lj->aij", U.conj(), np.asarray(cef.J), U)
    np.testing.assert_allclose(
        qc.susceptibility(cef, T), chi, rtol=0, atol=1e-12 * np.abs(chi).max()
    )


def test_conv_0009_tetragonal_tensor():
    """Tetragonal field, z the fourfold axis: chi_xx = chi_yy and no off-diagonal terms."""
    for chi in qc.susceptibility(_ce(), T):
        scale = np.abs(chi).max()
        assert chi[0, 0] == pytest.approx(chi[1, 1], rel=1e-12)
        assert np.abs(chi[~np.eye(3, dtype=bool)]).max() <= 1e-14 * scale


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


# --- CONV-0012: the scattering function as excitations --------------------
# Tests compare sums over the entries at one energy, never single entries.

SPECTRUM_CASES = [  # (ion, parameters in meV, field in T), test values
    ("Ce", dict(B20=0.3, B22=-0.1, B40=0.01, B4m3=0.02, B44=0.03), (0, 0, 0)),
    ("Pr", dict(B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4), (0, 0, 0)),  # cubic, degenerate
    ("Er", dict(B20=-0.02, B22=0.01, B40=1e-4, B43=3e-4), (0.1, 0.2, 0.3)),
]


def _cef(name, pars, field):
    return qc.CEFion(qc.Ion(name), qc.CEFpars(pars), Hfield=list(field))


def _distinct(energies, tolerance):
    """The distinct values of ``energies``, merged within ``tolerance``."""
    values = np.sort(energies)
    keep = np.r_[True, np.diff(values) > tolerance]
    return values[keep]


def _summed_at(energies, weights, energy, tolerance):
    """The sum of the weights of the entries within ``tolerance`` of ``energy``."""
    return weights[np.abs(energies - energy) <= tolerance].sum(axis=0)


@pytest.mark.parametrize(("name", "pars", "field"), SPECTRUM_CASES)
@pytest.mark.parametrize("temperature", [1.0, 30.0, 300.0])
def test_conv_0012_sum_rule(name, pars, field, temperature):
    """sum_k Tr W_k = sum_i p_i <i|J^2|i> = J(J+1)."""
    cef = _cef(name, pars, field)
    energies, W = qc.scattering_function(cef, temperature)
    n = len(cef.energies)
    assert energies.shape == (n * n,) and W.shape == (n * n, 3, 3)
    J = cef.ion.J
    assert np.trace(W, axis1=1, axis2=2).sum().real == pytest.approx(J * (J + 1), rel=1e-12)


@pytest.mark.parametrize(("name", "pars", "field"), SPECTRUM_CASES)
def test_conv_0012_detailed_balance(name, pars, field):
    """Summed over each energy: S^ab(-e) = exp(-e / k_B T) S^ba(e)."""
    temperature = 20.0  # K
    cef = _cef(name, pars, field)
    energies, W = qc.scattering_function(cef, temperature)
    tolerance = 1e-9 * np.abs(cef.energies).max()
    scale = np.abs(W).max()
    for e in _distinct(energies, tolerance):
        if e <= tolerance:
            continue
        gain = _summed_at(energies, W, -e, tolerance)
        loss = _summed_at(energies, W, e, tolerance)
        assert np.allclose(
            gain, np.exp(-e / (KB * temperature)) * loss.T, rtol=1e-9, atol=1e-13 * scale
        )


@pytest.mark.parametrize(("name", "pars", "field"), SPECTRUM_CASES)
def test_conv_0012_summed_weights_invariant_inside_levels(name, pars, field):
    """A unitary rotation inside each level leaves the summed weights at every energy unchanged."""
    temperature = 20.0  # K
    cef = _cef(name, pars, field)
    energies, W = qc.scattering_function(cef, temperature)
    rng = np.random.default_rng(12)
    n = len(cef.energies)
    U = np.zeros((n, n), complex)
    for level in degenerate_levels(cef.energies):
        d = level.stop - level.start
        U[level, level], _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    cef.J = np.einsum("ki,akl,lj->aij", U.conj(), np.asarray(cef.J), U)
    energies_r, W_r = qc.scattering_function(cef, temperature)
    tolerance = 1e-9 * np.abs(cef.energies).max()
    scale = np.abs(W).max()
    for e in _distinct(energies, tolerance):
        assert np.allclose(
            _summed_at(energies_r, W_r, e, tolerance),
            _summed_at(energies, W, e, tolerance),
            rtol=0,
            atol=1e-12 * scale,
        )


@pytest.mark.parametrize(("name", "pars", "field"), SPECTRUM_CASES)
def test_conv_0012_weights_hermitian_positive_semidefinite(name, pars, field):
    """Each W_k is a Hermitian, positive semidefinite 3x3 matrix."""
    _, W = qc.scattering_function(_cef(name, pars, field), 20.0)
    scale = np.abs(W).max()
    assert np.allclose(W, np.conj(np.swapaxes(W, 1, 2)), rtol=0, atol=1e-14 * scale)
    assert np.linalg.eigvalsh(W).min() >= -1e-14 * scale


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


@pytest.mark.parametrize(("field", "doublet"), [(1e-10, True), (1e-6, False)])  # T
def test_conv_0008_gtensor_ground_doublet(field, doublet):
    """The ground doublet is a level of CONV-0008: a split below the tolerance still counts.

    A field along x splits the ground doublet of the tetragonal Ce3+ ion by
    6e-12 meV at 1e-10 T (below the tolerance, 6.4e-10 meV) and by 6e-8 meV at
    1e-6 T (above it).
    """
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(TETRAGONAL), Hfield=[field, 0.0, 0.0])
    split = cef.energies[1] - cef.energies[0]
    assert (split <= _level_tolerance(cef.energies)) == doublet
    if doublet:
        g = qc.gtensor(cef)
        assert g.shape == (3, 3) and np.all(np.isfinite(g)) and np.abs(g).max() > 0
    else:
        with pytest.raises(AssertionError):
            qc.gtensor(cef)


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
    with pytest.raises(RuntimeError) as error:
        qc.magnetization_exchange(_ce(), 5.0, [0, 0, 1.0], lam=1.0)
    assert not isinstance(error.value, _Runaway)


def test_magnetization_exchange_zero_field():
    """In zero field the moment is zero, without a division 0/0 on the way."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        M = qc.magnetization_exchange(_ce(), 5.0, [0.0, 0.0, 0.0], lam=1.0)
    assert M == (0.0, 0.0, 0.0)
