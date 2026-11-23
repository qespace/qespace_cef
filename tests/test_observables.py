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
from qespace_cef.cefion import _degenerate_levels

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
    """The tolerance of CONV-0015 with its defaults rtol = 1e-10, E_floor = 1 meV, in meV."""
    return 1e-10 * max(np.max(energies) - np.min(energies), 1.0)


def test_conv_0015_susceptibility_independent_of_B00():
    """chi does not change when B00 shifts all energies (CONV-0014, CONV-0015, M5).

    Ce3+ with B20 = 5e-11 meV at 1e-9 K, a case sensitive to the grouping:
    with a tolerance from max|E|, B00 = 10 meV merged the three doublets, and
    chi_xx moved from 1.616e9 to 2.083e9 muB/T (29 %).

    What remains is rounding: shifted by 10 meV, the eigenvalues carry an
    error of about eps * 10 meV on splits of 3e-10 meV, and the van Vleck
    weight goes as 1 / split, so chi agrees to eps * max|E| / min split.
    """
    temperature = 1e-9  # K
    cefs = [_ce(dict(B20=5e-11, B00=B00)) for B00 in (0.0, 10.0)]
    chi = [qc.susceptibility(cef, temperature) for cef in cefs]
    E = cefs[1].energies
    split = np.diff(E)[np.diff(E) > 1e-13].min()  # meV, between the doublets
    bound = np.finfo(float).eps * np.abs(E).max() / split  # about 7e-6
    assert np.abs(chi[1] - chi[0]).max() <= bound * np.abs(chi[0]).max()


TETRAGONAL_PR = dict(B20=0.3, B40=0.01, B44=0.05)  # Pr3+, meV, test values: two doublets


# B22 splits the two non-Kramers doublets of tetragonal Pr3+ by about 9e-12 and
# 3e-11 meV (B22 = 1e-12 meV), 9e-11 and 3e-10, 9e-9 and 3e-8, 9e-8 and 3e-7 meV
# (1e-8 meV); the tolerance is about 2e-9 meV.
@pytest.mark.parametrize(
    ("B22", "below"), [(1e-12, True), (1e-11, True), (1e-9, False), (1e-8, False)]
)  # meV
def test_susceptibility_continuous_through_degeneracy(B22, below):
    """chi is continuous as a splitting shrinks through the tolerance (CONV-0015, P20).

    A tiny B22 splits the two doublets of the tetragonal Pr3+ ion, in zero
    field (chi is a zero-field quantity, CONV-0013). Below the tolerance each
    doublet is one level and takes the Curie weight p_n / k_B T; above it, the
    van Vleck weight (p_n - p_m) / (E_m - E_n), whose limit is the Curie
    weight. B22 -> -B22 is a rotation by 45 degrees about z, so chi_zz and
    Tr chi are even in B22 and change by O(B22^2): on both sides they must
    equal their values at B22 = 0.
    """
    temperature = 5.0  # K
    cef = qc.CEFion(qc.Ion("Pr"), qc.CEFpars(dict(TETRAGONAL_PR, B22=B22)))
    gaps = np.diff(cef.energies)
    splits = gaps[gaps < 1e-3]  # meV: the two split doublets
    assert len(splits) == 2
    assert np.all((splits <= _level_tolerance(cef.energies)) == below)  # the side
    assert splits.max() < 1e-6 * KB * temperature
    chi = qc.susceptibility(cef, temperature)
    chi0 = qc.susceptibility(qc.CEFion(qc.Ion("Pr"), qc.CEFpars(TETRAGONAL_PR)), temperature)
    assert chi[2, 2] == pytest.approx(chi0[2, 2], rel=1e-12)
    assert np.trace(chi) == pytest.approx(np.trace(chi0), rel=1e-12)


@pytest.mark.parametrize("method", ["perturbation", "magnetization"])
def test_conv_0013_susceptibility_is_zero_field(method):
    """chi is the zero-field response whatever field the given ion was solved in (H1).

    Free Ce3+ solved at 10 T along z: chi_zz at 2 K is the zero-field Curie
    value gJ^2 muB J(J+1) / (3 k_B T) = 0.7197 muB/T. The field-split
    eigensystem gave 1.4852 by perturbation.
    """
    ion = qc.Ion("Ce")
    cef = qc.CEFion(ion, qc.CEFpars(dict(B20=0.0)), Hfield=[0.0, 0.0, 10.0])
    chi = qc.susceptibility(cef, 2.0, method=method)
    expected = _curie(ion, 2.0) * np.eye(3)
    rel = 1e-10 if method == "perturbation" else 1e-5
    np.testing.assert_allclose(chi, expected, rtol=0, atol=rel * expected.max())


def test_conv_0015_susceptibility_levels_from_config(monkeypatch):
    """The Curie and van Vleck pairs follow the levels of CONV-0015, read at call time.

    With rtol = 0.5 the tolerance, 0.5 max|E| = 3.2 meV, merges the two lower
    doublets of the tetragonal Ce3+ ion (0.78 meV apart) into one level, so
    their pairs take the Curie weight. J_x couples the two doublets (J_z does
    not), so chi_xx sees it.
    """
    temperature = 5.0  # K
    cef = _ce()
    tolerance = 0.5 * (cef.energies.max() - cef.energies.min())  # meV
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
# U and F in meV per ion, S and C in k_B per ion (CONV-0003). The identities
# S = -(dF/dT) / k_B and C = (dU/dT) / k_B carry k_B (meV/K) for those units.


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


@pytest.mark.parametrize("temperature", T)
def test_entropy_is_minus_dF_dT(temperature):
    """S = -(dF/dT) / k_B, S in k_B and F in meV per ion."""
    S = _thermo([temperature])[3][0]
    assert S == pytest.approx(-_derivative(2, temperature) / KB, rel=1e-4)


@pytest.mark.parametrize("temperature", T)
def test_heat_capacity_is_dE_dT(temperature):
    """Cv = (dU/dT) / k_B, Cv in k_B and U in meV per ion."""
    Cv = _thermo([temperature])[4][0]
    assert Cv == pytest.approx(_derivative(1, temperature) / KB, rel=1e-4)


@pytest.mark.parametrize("temperature", T)
def test_heat_capacity_is_T_dS_dT(temperature):
    """Cv = T dS/dT."""
    Cv = _thermo([temperature])[4][0]
    assert Cv == pytest.approx(temperature * _derivative(3, temperature), rel=1e-4)


def _boltzmann_entropy_and_heat_capacity(energies, temperature):
    """S = -sum p ln p and C = <(E - <E>)^2> / (k_B T)^2, in k_B, from the populations."""
    x = (energies - energies.min()) / (KB * temperature)
    p = np.exp(-x) / np.exp(-x).sum()
    S = -np.sum(p[p > 0] * np.log(p[p > 0]))
    C = np.sum(p * (x - np.sum(p * x)) ** 2)
    return S, C


@pytest.mark.parametrize("B00", [0.0, 1.0])  # meV
def test_entropy_and_heat_capacity_independent_of_B00(B00):
    """S and C in k_B per ion do not depend on the zero of energy (CONV-0014, H2).

    Ce3+, B20=0.5, B40=0.01, B44=0.1 meV at 10 K: S = 1.2927, C = 0.1702 k_B
    (AUDIT_GPT H2) from the Boltzmann sums.
    """
    temperature = 10.0  # K
    cef = _ce(dict(TETRAGONAL, B00=B00))
    _, U, F, S, C = qc.thermodynamics(cef, np.array([temperature]))
    S_ref, C_ref = _boltzmann_entropy_and_heat_capacity(cef.energies, temperature)
    assert S[0] == pytest.approx(S_ref, rel=1e-12)
    assert C[0] == pytest.approx(C_ref, rel=1e-12)
    assert S[0] == pytest.approx(1.292718, abs=1e-6) and C[0] == pytest.approx(0.170246, abs=1e-6)
    _, U0, F0, _, _ = qc.thermodynamics(_ce(), np.array([temperature]))
    assert U[0] - U0[0] == pytest.approx(B00, abs=1e-12) and F[0] - F0[0] == pytest.approx(
        B00, abs=1e-12
    )


def test_thermodynamics_finite_at_low_temperature():
    """Ce3+ with B20=1 meV at 0.1 K: U -> -8 meV, F -> -8 meV - k_B T ln 2, S -> ln 2, C -> 0 (H3).

    The ground doublet lies at -8 meV and the next level 6 meV above it, so
    the excited levels carry exp(-6 meV / k_B T) ~ 1e-302 of the weight.
    """
    temperature = 0.1  # K
    cef = _ce(dict(B20=1.0))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _, U, F, S, C = qc.thermodynamics(cef, np.array([temperature]))
    assert U[0] == pytest.approx(-8.0, abs=1e-12)
    assert F[0] == pytest.approx(-8.0 - KB * temperature * np.log(2), abs=1e-12)
    assert S[0] == pytest.approx(np.log(2), rel=1e-12)
    assert C[0] == pytest.approx(0.0, abs=1e-12)


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
    for level in _degenerate_levels(cef.energies):
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


# --- CONV-0007: the neutron cross section ---------------------------------
# Expected values from formulas written here: (gamma r0)^2 from CODATA 2022,
# the projector, k_f/k_i and the matrix elements of J in the eigenbasis.

LOW_SYMMETRY = dict(B20=0.3, B22=-0.1, B40=0.01, B4m3=0.02, B44=0.03)
TEMPERATURE = 10.0  # K
EI = 1e3  # meV, test value
# (gamma r0)^2, CODATA 2022: mu_n/mu_N = -1.91304276(45), r_e = 2.8179403205(13)e-15 m.
GAMMA_R0_SQ = (1.91304276 * 2.8179403205e-15) ** 2 / 1e-28  # barn


def _populations(cef):
    p = np.exp(-(cef.energies - cef.energies.min()) / (KB * TEMPERATURE))
    return p / p.sum()


def _neutron_factor(cef, Qnorm, u2=0.0):
    """(gamma r0)^2 (gJ/2)^2 F(|Q|)^2 exp(-|Q|^2 u2 / 3), in barn."""
    return (
        GAMMA_R0_SQ * (cef.ion.gJ / 2) ** 2 * cef.ion.mff(Qnorm) ** 2 * np.exp(-(Qnorm**2) * u2 / 3)
    )


def _bare_weights(cef, energies, weights, Qnorm, Ei=EI):
    """weights / (neutron factor k_f/k_i), as an (n, n) array of pairs (i, f)."""
    n = len(cef.energies)
    kf_ki = np.sqrt(1 - energies / Ei)
    return (weights / (_neutron_factor(cef, Qnorm) * kf_ki)).reshape(n, n)


def _projected_matrix_elements(cef, P):
    """p_i sum_ab P_ab <i|J_a|f><f|J_b|i>, as (n, n), from the J matrices."""
    M = np.einsum("ab,aif,bfi->if", P, cef.J, cef.J).real
    return _populations(cef)[:, None] * M


def test_powder_sum_rule():
    """sum_f |<i|J_perp|f>|^2 = 2/3 J(J+1) for every initial state i (powder)."""
    cef = _ce(LOW_SYMMETRY)
    energies, w = qc.cross_section_powder(cef, TEMPERATURE, 1.0, Ei=EI)
    rows = _bare_weights(cef, energies, w, 1.0).sum(axis=1)
    J = cef.ion.J
    assert np.allclose(rows / _populations(cef), 2 / 3 * J * (J + 1), rtol=1e-10)


def test_single_crystal_trace_sum_rule():
    """sum_if |<i|J_perp|f>|^2 = 2/3 J(J+1)(2J+1), the trace over i (single crystal)."""
    cef = _ce(LOW_SYMMETRY)
    Q = np.array([0.3, 0.5, 0.7])
    energies, w = qc.cross_section_crystal(cef, TEMPERATURE, Q, Ei=EI)
    W = _bare_weights(cef, energies, w, np.linalg.norm(Q))
    J = cef.ion.J
    assert np.sum(W / _populations(cef)[:, None]) == pytest.approx(
        2 / 3 * J * (J + 1) * (2 * J + 1), rel=1e-10
    )


@pytest.mark.parametrize("Q", [[1.0, 1.0, 0.0], [0.3, 0.5, 0.7]])
def test_single_crystal_weights_full_projector(Q):
    """The weights summed at each energy equal, summed at the same energy,
    neutron factor (kf/ki) p_i sum_ab (delta_ab - Qa Qb / Q^2) <i|J_a|f><f|J_b|i> (P9).

    Sums over equal energies, never single entries (CONV-0012); energies
    grouped as in the CONV-0012 tests.
    """
    cef = _ce(LOW_SYMMETRY)
    Qnorm = np.linalg.norm(Q)
    q = np.asarray(Q) / Qnorm
    E = cef.energies
    expected_energies = (E[None, :] - E[:, None]).ravel()  # E_f - E_i, pairs (i, f)
    expected = (
        _neutron_factor(cef, Qnorm)
        * np.sqrt(1 - expected_energies / EI)
        * _projected_matrix_elements(cef, np.eye(3) - np.outer(q, q)).ravel()
    )
    energies, w = qc.cross_section_crystal(cef, TEMPERATURE, Q, Ei=EI)
    tolerance = 1e-9 * np.abs(E).max()
    scale = expected.max()
    for e in _distinct(expected_energies, tolerance):
        assert _summed_at(energies, w, e, tolerance) == pytest.approx(
            _summed_at(expected_energies, expected, e, tolerance), rel=1e-10, abs=1e-12 * scale
        )


@pytest.mark.parametrize("name", ["Nd", "Tb", "Ho"])  # F(0) = 1 to rounding in the table
def test_conv_0007_absolute_sum_rule(name):
    """Powder, Q = 0, u2 = 0: sum_k w_k / (kf/ki)_k = (gamma r0)^2 (gJ/2)^2 (2/3) J(J+1)."""
    cef = _cef(name, LOW_SYMMETRY, (0, 0, 0))
    energies, w = qc.cross_section_powder(cef, TEMPERATURE, 0.0, Ei=EI)
    J, gJ = cef.ion.J, cef.ion.gJ
    expected = GAMMA_R0_SQ * (gJ / 2) ** 2 * 2 / 3 * J * (J + 1)  # barn / sr
    assert np.sum(w / np.sqrt(1 - energies / EI)) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("kind", ["powder", "crystal"])
def test_conv_0007_nan_exactly_above_Ei(kind):
    """The weights are NaN exactly where the excitation energy exceeds Ei, without a warning."""
    cef = _ce(LOW_SYMMETRY)
    E = cef.energies
    Ei = 0.5 * (E.max() - E.min())  # meV, between the excitations
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if kind == "powder":
            energies, w = qc.cross_section_powder(cef, TEMPERATURE, 1.0, Ei=Ei)
        else:
            energies, w = qc.cross_section_crystal(cef, TEMPERATURE, [0.3, 0.5, 0.7], Ei=Ei)
    assert np.any(energies > Ei) and np.any(energies < Ei)
    assert np.array_equal(np.isnan(w), energies > Ei)


def test_conv_0007_powder_is_crystal_averaged_over_directions():
    """The crystal weights averaged over three orthogonal Q directions equal the powder weights.

    Any orthonormal triple q_1, q_2, q_3 has sum_j q_j q_j^T = 1, so the mean
    projector is (2/3) delta, the average over all directions.
    """
    cef = _ce(LOW_SYMMETRY)
    Qnorm = 1.3  # 1/A
    axes, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))
    energies, crystal = qc.cross_section_crystal(cef, TEMPERATURE, Qnorm * axes.T, Ei=EI)
    _, powder = qc.cross_section_powder(cef, TEMPERATURE, Qnorm, Ei=EI)
    assert crystal.shape == (3, len(energies))
    finite = ~np.isnan(powder)
    assert np.allclose(crystal.mean(axis=0)[finite], powder[finite], rtol=1e-12, atol=1e-15)


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_conv_0007_projector_along_an_axis(axis):
    """Q along an axis: the full projector is diagonal, sum over the two other components."""
    cef = _ce(LOW_SYMMETRY)
    Q = np.zeros(3)
    Q[axis] = 1.2  # 1/A
    P = np.diag([0.0 if a == axis else 1.0 for a in range(3)])
    expected = _projected_matrix_elements(cef, P)
    energies, w = qc.cross_section_crystal(cef, TEMPERATURE, Q, Ei=EI)
    W = _bare_weights(cef, energies, w, 1.2)
    assert np.allclose(W, expected, rtol=1e-10, atol=1e-12 * expected.max())


def test_conv_0007_debye_waller_and_loop_shape():
    """u2 scales the weights by exp(-|Q|^2 u2 / 3); Q of shape S + (3,) gives S + (n^2,)."""
    cef = _ce(LOW_SYMMETRY)
    Q = np.array([[[0.3, 0.5, 0.7]], [[1.0, 0.0, 2.0]]])  # S = (2, 1)
    u2 = 0.02  # A^2, test value
    energies, w0 = qc.cross_section_crystal(cef, TEMPERATURE, Q, Ei=EI)
    _, w = qc.cross_section_crystal(cef, TEMPERATURE, Q, Ei=EI, u2=u2)
    n2 = len(energies)
    assert w.shape == (2, 1, n2)
    factor = np.exp(-np.sum(Q**2, axis=-1) * u2 / 3)[..., None]
    assert np.allclose(w, w0 * factor, rtol=1e-14, atol=0)


@pytest.mark.parametrize(
    ("func", "Q", "Ei", "u2"),
    [
        ("cross_section_crystal", [1.0, 0.0], EI, 0.0),
        ("cross_section_crystal", [0.0, 0.0, 0.0], EI, 0.0),
        ("cross_section_crystal", [1.0, 0.0, 0.0], 0.0, 0.0),
        ("cross_section_crystal", [1.0, 0.0, 0.0], EI, np.nan),
        ("cross_section_powder", -1.0, EI, 0.0),
        ("cross_section_powder", 1.0, -5.0, 0.0),
        ("cross_section_powder", 1.0, EI, -0.01),
    ],
)
def test_cross_sections_reject_bad_input(func, Q, Ei, u2):
    with pytest.raises(ValueError):
        getattr(qc, func)(_ce(LOW_SYMMETRY), TEMPERATURE, Q, Ei=Ei, u2=u2)


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
def test_conv_0015_gtensor_ground_doublet(field, doublet):
    """The ground doublet is a level of CONV-0015: a split below the tolerance still counts.

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
        with pytest.raises(ValueError):
            qc.gtensor(cef)


@pytest.mark.parametrize(
    ("name", "pars"),
    [("Ce", {}), ("Nd", dict(B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4))],
)  # meV: free Ce3+ (one sixfold level), cubic Nd3+ (a quartet among its levels)
def test_gtensor_needs_exactly_one_ground_doublet(name, pars):
    """A ground level that is not exactly twofold raises ValueError (H4, A10).

    The free ion has one level of 2J + 1 states; a ground quartet has four.
    states 0 and 1 then span only part of the level, and g would describe a
    wrong subspace.
    """
    cef = qc.CEFion(qc.Ion(name), qc.CEFpars(pars))
    ground = _degenerate_levels(cef.energies)[0]
    assert ground.stop - ground.start == {"Ce": 6, "Nd": 4}[name]  # the case the test is about
    with pytest.raises(ValueError):
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


# --- validation of public inputs (M3, M6, L1) -----------------------------


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf])
def test_populations_need_a_positive_temperature(temperature):
    """Populations at T <= 0 were NaN (M6)."""
    with pytest.raises(ValueError):
        qc.boltzman_population(np.array([0.0, 0.0, 1.0]), temperature)


@pytest.mark.parametrize("temperatures", [[0.0], [5.0, -1.0], [np.nan]])
def test_thermodynamics_needs_positive_temperatures(temperatures):
    with pytest.raises(ValueError):
        qc.thermodynamics(_ce(), np.array(temperatures))


def test_scattering_function_needs_a_scalar_temperature():
    """Six temperatures gave one spectrum, paired with the six states (M3)."""
    with pytest.raises(ValueError):
        qc.scattering_function(_ce(), np.linspace(1.0, 50.0, 6))


@pytest.mark.parametrize("Q", [[np.nan, 0.0, 0.0], [1.0, np.inf, 0.0]])
def test_cross_section_crystal_needs_a_finite_Q(Q):
    """A NaN component made every weight NaN (M6)."""
    with pytest.raises(ValueError):
        qc.cross_section_crystal(_ce(), 10.0, Q, Ei=50.0)


@pytest.mark.parametrize("max_iterations", [0, -1])
def test_magnetization_exchange_needs_one_iteration(max_iterations):
    """max_iterations = 0 raised UnboundLocalError (L1)."""
    with pytest.raises(ValueError):
        qc.magnetization_exchange(
            _ce(), 5.0, [0.0, 0.0, 1.0], lam=1.0, max_iterations=max_iterations
        )
