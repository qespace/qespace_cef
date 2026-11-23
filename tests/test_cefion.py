"""Tests of ``qespace_cef.CEFion`` against analytic limits.

A known problem from AUDIT.md is a strict ``xfail`` naming its ID: it turns
into a failure (XPASS) when the problem is fixed, so the marker must then be
removed.
"""

import re
import warnings

import numpy as np
import pytest

import qespace_cef as qc
from qespace_cef import constants as C
from qespace_cef.cefion import _degenerate_levels
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
    """States in one level of CONV-0008 print as one level; a larger splitting does not.

    Cubic Pr3+: eigh returns the 28.584 meV doublet as two values that differ
    by rounding noise; the print groups them. Ce3+ in 1e-6 T splits each
    doublet by about 5e-8 meV, above the tolerance of 1e-9 meV, which the
    print keeps apart.
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


def test_degeneracies_of_cubic_pr():
    """Pr3+ in a cubic field: levels grouped within 1e-9 of the largest energy."""
    e = _cefion("Pr", CUBIC_PR).energies
    new_level = np.diff(e) > 1e-9 * e.max()
    expected = np.diff(np.flatnonzero(np.r_[True, new_level, True])).tolist()
    assert [d for _, d in _cefion("Pr", CUBIC_PR).get_degeneracies()] == expected


@pytest.mark.parametrize(("name", "pars"), [("Pr", CUBIC_PR), ("Er", LOW_SYMMETRY)])
def test_degeneracies_report_mean_energy(name, pars):
    """Each level is reported with the mean of its energies (meV, CONV-0010), not truncated."""
    e = _cefion(name, pars).energies
    new_level = np.diff(e) > 1e-9 * np.abs(e).max()
    groups = np.split(e, np.flatnonzero(new_level) + 1)
    levels = _cefion(name, pars).get_degeneracies()
    assert [d for _, d in levels] == [len(g) for g in groups]
    for (energy, _), group in zip(levels, groups):
        assert energy == pytest.approx(group.mean(), rel=1e-12, abs=1e-12 * np.abs(e).max())


# --- CONV-0008: degenerate levels -----------------------------------------
# Expected values from the physics (Kramers' theorem, the free ion), not from
# the code under test.

CUBIC_HO = dict(B40=1e-3, B44=5e-3, B60=1e-5, B64=-21e-5)  # meV, test values, cubic
CEF_SCALE = {2: 0.5, 4: 0.01, 6: 1e-4}  # meV, typical size of B_lm per rank l, test values


def _random_low_symmetry(rng):
    """Every B_lm of rank 2, 4, 6, random: no symmetry left but time reversal."""
    pars = {}
    for rank in (2, 4, 6):
        for m in range(-rank, rank + 1):
            label = f"B{rank}{m}" if m >= 0 else f"B{rank}m{-m}"
            pars[label] = CEF_SCALE[rank] * rng.normal()
    return pars


def _degeneracies(energies):
    return [s.stop - s.start for s in _degenerate_levels(energies)]


def test_conv_0008_kramers_doublets_are_one_level():
    """In 200 low-symmetry Hamiltonians of Kramers ions every level is one exact doublet."""
    rng = np.random.default_rng(8)
    for k in range(200):
        cef = _cefion(KRAMERS_IONS[k % len(KRAMERS_IONS)], _random_low_symmetry(rng))
        assert _degeneracies(cef.energies) == [2] * (len(cef.energies) // 2)


@pytest.mark.parametrize("name", ALL_IONS)
def test_conv_0008_free_ion_is_one_level(name):
    """H = 0: all 2J+1 states form one level (the tolerance is rtol E_floor)."""
    cef = _cefion(name, {})
    assert _degeneracies(cef.energies) == [int(2 * cef.ion.J + 1)]


@pytest.mark.parametrize(
    ("name", "pars"),
    [("Ce", LOW_SYMMETRY), ("Er", LOW_SYMMETRY), ("Pr", CUBIC_PR), ("Ho", CUBIC_HO)],
)
def test_conv_0008_grouping_stable_under_perturbation(name, pars):
    """A random Hermitian perturbation of size eps max|E| leaves the levels unchanged."""
    H = _cefion(name, pars).hamiltonian
    levels = _degenerate_levels(np.linalg.eigvalsh(H))
    scale = np.finfo(float).eps * np.abs(np.linalg.eigvalsh(H)).max()
    rng = np.random.default_rng(80)
    for _ in range(20):
        A = rng.normal(size=H.shape) + 1j * rng.normal(size=H.shape)
        dH = (A + A.conj().T) / 2
        dH *= scale / np.abs(np.linalg.eigvalsh(dH)).max()  # spectral norm eps max|E|
        assert _degenerate_levels(np.linalg.eigvalsh(H + dH)) == levels


def test_conv_0008_tolerance_read_from_config_at_call_time(monkeypatch):
    """rtol and E_floor come from qespace_cef.config when the function is called."""
    E = np.array([-1e-3, -1e-3 + 5e-11, 1e-3])  # meV; split 5e-11 meV
    assert _degeneracies(E) == [2, 1]  # tolerance 1e-10 * max(1e-3, 1) meV
    monkeypatch.setattr(qc.config, "rtol", 1e-11)
    assert _degeneracies(E) == [1, 1, 1]  # 1e-11 meV
    monkeypatch.setattr(qc.config, "rtol", 1e-10)
    monkeypatch.setattr(qc.config, "E_floor", 1e-6)
    assert _degeneracies(E) == [1, 1, 1]  # 1e-10 * 1e-3 meV


@pytest.mark.parametrize("energies", [[1.0, 0.0], [[0.0, 1.0]]])
def test_degenerate_levels_rejects_unsorted_or_not_1d(energies):
    with pytest.raises(ValueError):
        _degenerate_levels(energies)


# --- CONV-0011: the canonical basis printed by summary --------------------
# The basis is read back from the text summary prints, its public output.

CANONICAL_CASES = [  # (ion, parameters in meV, field in T), test values
    ("Ce", dict(B20=0.5, B40=0.01, B44=0.1), (0, 0, 0)),  # Kramers doublets
    ("Ce", LOW_SYMMETRY, (0, 0, 0)),
    ("Nd", CUBIC_PR, (0, 0, 0)),  # cubic: Kramers doublets and quartets
    ("Pr", dict(B20=0.3, B40=0.01, B44=0.05), (0, 0, 0)),  # integer J: two doublets
    ("Pr", CUBIC_PR, (0, 0, 0)),  # integer J: triplets
    ("Er", LOW_SYMMETRY, (0.1, 0.2, 0.3)),  # field: no degeneracy
]


def _printed_basis(cef, precision=12):
    """The eigenvectors summary prints, one array per level, rows in the free-ion basis."""
    levels = []
    for block in cef.summary(precision=precision).strip().split("\n\n")[1:]:
        rows = []
        for line in block.splitlines()[1:]:
            v = np.zeros(len(cef.freeionkets), complex)
            for coefficient, ket in re.findall(r"\(([^()]*)\)(\|[^>]*>)", line):
                v[cef.freeionkets.index(ket)] = complex(coefficient)
            rows.append(v)
        levels.append(np.array(rows))
    return levels


def _jz(J):
    """J_z in the basis |J>, ..., |-J> (CONV-0005)."""
    return np.diag(np.arange(J, -J - 1, -1))


def _time_reversal(v, J):
    """Theta|m> = (-1)^(J-m) |-m>, then complex conjugation; v in the basis |J>, ..., |-J>."""
    m = np.arange(J, -J - 1, -1)
    out = np.zeros_like(v)
    for i, mi in enumerate(m):
        out[list(m).index(-mi)] = (-1) ** round(J - mi) * np.conj(v[i])
    return out


GAMMA3_UNFIXED = pytest.mark.xfail(
    strict=True,
    reason="CONV-0011 rule 1: cubic Pr3+ has a non-magnetic Gamma3 doublet in which J_x, J_y "
    "and J_z all vanish (|<i|J|j>| ~ 2e-15), so no J component fixes its basis",
)


@pytest.mark.parametrize(
    ("name", "pars", "field"),
    [
        pytest.param(*case, marks=GAMMA3_UNFIXED)
        if case[1] is CUBIC_PR and case[0] == "Pr"
        else case
        for case in CANONICAL_CASES
    ],
)
def test_conv_0011_basis_unchanged_by_rotation_inside_levels(name, pars, field):
    """A random unitary rotation inside each level does not change the printed basis."""
    cef = _cefion(name, pars, field)
    before = _printed_basis(cef)
    rng = np.random.default_rng(11)
    for level in _degenerate_levels(cef.energies):
        d = level.stop - level.start
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
        cef.eigenvectors[level] = Q @ cef.eigenvectors[level]
    for a, b in zip(before, _printed_basis(cef)):
        assert np.allclose(a, b, rtol=0, atol=1e-9)


@pytest.mark.parametrize(("name", "pars", "field"), CANONICAL_CASES)
def test_conv_0011_jz_ascending_and_phases(name, pars, field):
    """<J_z> ascends inside each level; the states that get their own phase have their
    largest component real and positive (CONV-0011, rules 1-3)."""
    cef = _cefion(name, pars, field)
    J = cef.ion.J
    kramers = J % 1 == 0.5 and not np.any(field)
    for states in _printed_basis(cef):
        jz = np.einsum("ik,kl,il->i", states.conj(), _jz(J), states).real
        assert np.all(np.diff(jz) >= -1e-9)
        own_phase = states[jz < 0] if kramers else states
        for v in own_phase:
            largest = v[np.argmax(np.abs(v))]
            assert abs(largest.imag) <= 1e-9 and largest.real > 0


@pytest.mark.parametrize(("name", "pars", "field"), CANONICAL_CASES[:3])
def test_conv_0011_kramers_partners_by_time_reversal(name, pars, field):
    """Kramers ion in zero field: the partner of each state with <J_z> < 0 is Theta of it."""
    cef = _cefion(name, pars, field)
    J = cef.ion.J
    for states in _printed_basis(cef):
        d = len(states)
        assert d % 2 == 0
        for i in range(d // 2):
            assert np.real(states[i].conj() @ _jz(J) @ states[i]) < 0
            assert np.allclose(_time_reversal(states[i], J), states[d - 1 - i], rtol=0, atol=1e-9)


# --- summary: time reversal and the coefficient cutoff (H5, L2) -----------


@pytest.mark.parametrize(
    ("name", "pars", "field"),
    [
        ("Ce", dict(B20=1.0, B10=-6.0), (0, 0, 0)),  # odd rank: not time-even
        ("Ce", dict(B20=0.5, B40=0.01, B44=0.1), (0, 0, 0)),
        ("Er", LOW_SYMMETRY, (0.1, 0.2, 0.3)),
    ],
)
def test_summary_states_lie_in_their_level(name, pars, field):
    """Every printed state lies inside the eigenspace of its level (H5).

    Ce3+ with B20=1, B10=-6 meV: the ground level is |+1/2>, |+3/2> at
    -11 meV. B10 J_z is odd under time reversal, so pairing a state with its
    time reverse put |-1/2> (at -5 meV) into the ground level.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # odd-rank parameters warn
        cef = _cefion(name, pars, field)
    for level, states in zip(_degenerate_levels(cef.energies), _printed_basis(cef)):
        V = cef.eigenvectors[level]  # rows span the level
        residual = states - (states @ V.conj().T) @ V
        assert np.abs(residual).max() < 1e-9


def test_summary_precision_zero_keeps_pure_kets():
    """``precision=0`` keeps coefficients at the cutoff 10^0 = 1: pure kets are printed (L2)."""
    text = _cefion("Ce", dict(B20=1.0)).summary(precision=0)
    lines = [line for line in text.splitlines() if line.startswith("ev_")]
    assert len(lines) == 6 and all("|" in line for line in lines)
