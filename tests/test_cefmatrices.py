"""Tests of the J matrices and the Stevens operator matrices in ``qespace_cef.cefmatrices``."""

import math

import numpy as np
import pytest
from reference import mcphase
from reference import old_stevens as old
from reference.pycrystalfield import StevensOpPCF

from qespace_cef.cefmatrices import J_minus, J_plus, J_x, J_y, J_z, StevensBase

# Total angular momenta J = 1/2 ... 8, covering the rare-earth ions.
J_VALUES = [float(J) for J in np.arange(1, 17) / 2]

# Every operator O_k^q with k = 0 ... 6 and q = -k ... k.
ALL_OPERATORS = [f"O{k}{q}" if q >= 0 else f"O{k}m{-q}" for k in range(7) for q in range(-k, k + 1)]

# Operators compared with PyCrystalField, and with the old implementation where
# one exists: name -> (old implementation or None, PyCrystalField (k, q)).
REFERENCES = {
    "O11": (J_x, (1, 1)),
    "O10": (J_z, (1, 0)),
    "O1m1": (J_y, (1, -1)),
    "O22": (old.O22, (2, 2)),
    "O21": (None, (2, 1)),
    "O20": (old.O20, (2, 0)),
    "O2m1": (None, (2, -1)),
    "O2m2": (old.O2m2, (2, -2)),
    "O33": (None, (3, 3)),
    "O32": (None, (3, 2)),
    "O31": (None, (3, 1)),
    "O30": (None, (3, 0)),
    "O3m1": (None, (3, -1)),
    "O3m2": (None, (3, -2)),
    "O3m3": (None, (3, -3)),
    "O44": (old.O44, (4, 4)),
    "O43": (old.O43, (4, 3)),
    "O42": (old.O42, (4, 2)),
    "O41": (None, (4, 1)),
    "O40": (old.O40, (4, 0)),
    "O4m1": (None, (4, -1)),
    "O4m2": (old.O4m2, (4, -2)),
    "O4m3": (old.O4m3, (4, -3)),
    "O4m4": (old.O4m4, (4, -4)),
    "O66": (old.O66, (6, 6)),
    "O65": (None, (6, 5)),
    "O64": (old.O64, (6, 4)),
    "O63": (old.O63, (6, 3)),
    "O62": (old.O62, (6, 2)),
    "O61": (None, (6, 1)),
    "O60": (old.O60, (6, 0)),
    "O6m1": (None, (6, -1)),
    "O6m2": (old.O6m2, (6, -2)),
    "O6m3": (old.O6m3, (6, -3)),
    "O6m4": (old.O6m4, (6, -4)),
    "O6m5": (None, (6, -5)),
    "O6m6": (old.O6m6, (6, -6)),
}


def _close(a, b):
    """True if ``a`` equals ``b`` to 1e-10 of the largest element of ``a``."""
    return np.allclose(a, b, atol=np.max(np.abs(a)) * 1e-10)


def _m_values(J):
    """The m of the basis states, in the order of CONV-0005: J, J-1, ..., -J."""
    return [J - i for i in range(int(2 * J + 1))]


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0005_jz_diagonal(J):
    """CONV-0005: J_z = diag(J, J-1, ..., -J), so <J|J_z|J> = J for the first state."""
    Jz = J_z(J)
    assert Jz.shape == (int(2 * J + 1),) * 2
    assert Jz[0, 0] == J
    assert _close(Jz, np.diag(_m_values(J)))
    assert np.count_nonzero(Jz - np.diag(np.diag(Jz))) == 0


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0005_jplus_condon_shortley(J):
    """CONV-0005: <m+1|J+|m> = sqrt(J(J+1) - m(m+1)), real and non-negative, alone on the first superdiagonal.

    Row i is the state m_i = J - i, so <m+1|J+|m> sits in row i, column i+1,
    with m = m_{i+1}.
    """
    ms = _m_values(J)
    n = len(ms)
    expected = np.zeros((n, n))
    for i in range(n - 1):
        m = ms[i + 1]
        expected[i, i + 1] = math.sqrt(J * (J + 1) - m * (m + 1))
    Jp = J_plus(J)
    assert np.isrealobj(Jp)
    assert np.all(np.diag(Jp, 1) >= 0)
    assert _close(Jp, expected)
    off_superdiagonal = ~np.eye(n, k=1, dtype=bool)
    assert np.count_nonzero(Jp[off_superdiagonal]) == 0


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0005_jminus_jx_jy(J):
    """CONV-0005: J- = J+^dagger, J_x = (J+ + J-)/2, J_y = (J+ - J-)/2i."""
    Jp = J_plus(J)
    assert _close(J_minus(J), Jp.conj().T)
    assert _close(J_x(J), (Jp + Jp.conj().T) / 2)
    assert _close(J_y(J), (Jp - Jp.conj().T) / 2j)


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0005_commutator_and_j_squared(J):
    """CONV-0005: [J_x, J_y] = i J_z and J_x^2 + J_y^2 + J_z^2 = J(J+1)."""
    Jx, Jy, Jz = J_x(J), J_y(J), J_z(J)
    assert _close(Jx @ Jy - Jy @ Jx, 1j * Jz)
    assert _close(Jx @ Jx + Jy @ Jy + Jz @ Jz, J * (J + 1) * np.eye(int(2 * J + 1)))


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0005_stevens_base_uses_these(J):
    """CONV-0005: the J matrices StevensBase uses are the module's J_x, J_y, J_z, J_plus, J_minus."""
    ops = StevensBase(J)
    pairs = {"Jx": J_x, "Jy": J_y, "Jz": J_z, "Jp": J_plus, "Jm": J_minus}
    different = [a for a, f in pairs.items() if not np.array_equal(getattr(ops, a), f(J))]
    assert different == []


def test_implemented_operators():
    """Every O_k^q with k <= 6 is implemented, and nothing else."""
    assert sorted(StevensBase.implementedOps()) == sorted(ALL_OPERATORS)


@pytest.mark.parametrize("J", J_VALUES)
def test_conv_0001_operators(J):
    """CONV-0001: every even-rank O_lm equals the McPhase definition, signs included.

    The reference (tests/reference/mcphase.py) is built from the formulas and
    the polynomials copied into CONV-0001, for l = 2, 4, 6 and m = -l ... l.
    """
    ops = StevensBase(J)
    mismatches = [lm for lm in mcphase.LM if not _close(ops[_name(*lm)], mcphase.stevens(J, *lm))]
    assert mismatches == []


@pytest.mark.parametrize("J", J_VALUES)
def test_matrices_match_old_implementation(J):
    """Stevens matrices equal the first implementation (tests/reference)."""
    ops = StevensBase(J)
    mismatches = [
        name
        for name, (reference, _) in REFERENCES.items()
        if reference is not None and not _close(ops[name[1:]], reference(J))
    ]
    assert mismatches == []


# PyCrystalField's operators with the opposite sign, all with l + m odd. The
# package follows McPhase (CONV-0001), so the difference is recorded, not
# accepted up to a sign: exactly -1 (a passing test) and a strict xfail.
PCF_OPPOSITE_SIGN = {
    "O10", "O21", "O2m1", "O30", "O32", "O3m2", "O41", "O4m1",
    "O43", "O4m3", "O61", "O6m1", "O63", "O6m3", "O65", "O6m5",
}  # fmt: skip


@pytest.mark.parametrize("J", J_VALUES)
def test_matrices_match_pycrystalfield(J):
    """Stevens matrices equal PyCrystalField's, signs included, except PCF_OPPOSITE_SIGN."""
    ops = StevensBase(J)
    mismatches = [
        name
        for name, (_, (k, q)) in REFERENCES.items()
        if name not in PCF_OPPOSITE_SIGN and not _close(ops[name[1:]], StevensOpPCF(J, k, q))
    ]
    assert mismatches == []


@pytest.mark.parametrize("J", J_VALUES)
def test_pycrystalfield_opposite_sign(J):
    """The operators of PCF_OPPOSITE_SIGN are exactly -1 times PyCrystalField's."""
    ops = StevensBase(J)
    mismatches = [
        name
        for name in sorted(PCF_OPPOSITE_SIGN)
        if not _close(ops[name[1:]], -StevensOpPCF(J, *REFERENCES[name][1]))
    ]
    assert mismatches == []


@pytest.mark.xfail(
    strict=True,
    reason="PyCrystalField's O_lm with l + m odd have the opposite sign to McPhase, "
    "which the package follows (CONV-0001)",
)
@pytest.mark.parametrize("name", sorted(PCF_OPPOSITE_SIGN))
def test_pycrystalfield_sign_conv_0001(name):
    """CONV-0001: PyCrystalField's O_lm equals the package's, signs included."""
    k, q = REFERENCES[name][1]
    mismatches = [
        J for J in J_VALUES if not _close(StevensBase(J)[name[1:]], StevensOpPCF(J, k, q))
    ]
    assert mismatches == []


@pytest.mark.parametrize("J", J_VALUES)
def test_matrices_hermitian(J):
    """Every Stevens matrix is Hermitian."""
    ops = StevensBase(J)
    not_hermitian = [n for n in ALL_OPERATORS if not _close(ops[n[1:]], ops[n[1:]].T.conj())]
    assert not_hermitian == []


@pytest.mark.parametrize("J", J_VALUES)
def test_getitem_and_method_agree(J):
    """``StevensBase(J)['42']`` and ``StevensBase(J).O42()`` give the same matrix."""
    ops = StevensBase(J)
    assert all(np.array_equal(ops[n[1:]], getattr(ops, n)()) for n in ALL_OPERATORS)


def _name(k, q):
    return f"{k}{q}" if q >= 0 else f"{k}m{-q}"


@pytest.mark.parametrize("J", J_VALUES)
@pytest.mark.parametrize("sign", [1, -1])
def test_observations_1_and_2(J, sign):
    """O(n, +-(n-1)) = 1/2 {O(n-1, +-(n-1)), Jz} for n = 2 ... 6 (StevensBase notes)."""
    ops = StevensBase(J)
    Jz = ops.Jz
    mismatches = []
    for n in range(2, 7):
        lower = ops[_name(n - 1, sign * (n - 1))]
        if not _close(0.5 * (lower @ Jz + Jz @ lower), ops[_name(n, sign * (n - 1))]):
            mismatches.append(n)
    assert mismatches == []


@pytest.mark.parametrize("J", J_VALUES)
def test_observation_4_negative_q(J):
    """Every non-zero O_k^q with q < 0 is purely imaginary and not symmetric (StevensBase notes)."""
    ops = StevensBase(J)
    wrong = []
    for k in range(1, 7):
        for q in range(-k, 0):
            op = ops[_name(k, q)]
            if np.abs(op).max() > 0 and (not np.allclose(op.real, 0) or np.allclose(op, op.T)):
                wrong.append(_name(k, q))
    assert wrong == []


@pytest.mark.parametrize("J", [0.75, -0.5, 1.2, float("nan")])
@pytest.mark.parametrize("func", ["J_z", "J_y", "J_x", "J_plus", "J_minus"])
def test_j_matrices_reject_j_not_a_multiple_of_one_half(func, J):
    """J must be a non-negative multiple of 1/2 (M6): J_z(0.75) had m values 1.5 apart."""
    from qespace_cef import cefmatrices

    with pytest.raises(ValueError):
        getattr(cefmatrices, func)(J)


@pytest.mark.parametrize("J", [0, 0.5, 4, 7.5])
def test_j_matrices_accept_multiples_of_one_half(J):
    from qespace_cef.cefmatrices import J_z

    assert J_z(J).shape == (int(2 * J + 1),) * 2
