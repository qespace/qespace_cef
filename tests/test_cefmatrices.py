"""Tests of the Stevens operator matrices in ``qespace_cef.cefmatrices``."""

import numpy as np
import pytest
from reference import old_stevens as old
from reference.pycrystalfield import StevensOpPCF

from qespace_cef.cefmatrices import J_x, J_y, J_z, StevensBase

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


def test_implemented_operators():
    """Every O_k^q with k <= 6 is implemented, and nothing else."""
    assert sorted(StevensBase.implementedOps()) == sorted(ALL_OPERATORS)


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


@pytest.mark.parametrize("J", J_VALUES)
def test_matrices_match_pycrystalfield(J):
    """Stevens matrices equal PyCrystalField's, up to an overall sign."""
    ops = StevensBase(J)
    mismatches = []
    for name, (_, (k, q)) in REFERENCES.items():
        new, pcf = ops[name[1:]], StevensOpPCF(J, k, q)
        if not (_close(new, pcf) or _close(new, -pcf)):
            mismatches.append(name)
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
