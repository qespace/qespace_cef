"""Tests of ``qespace_cef.CEFpars``."""

import numpy as np
import pytest

import qespace_cef as qc

# k_B / e = 8.617333262e-5 eV/K and h c / e = 1.239841984e-4 eV cm (CODATA 2018),
# so 1 K = 0.08617333 meV and 1 cm^-1 = 0.1239842 meV.
K_IN_MEV = 0.08617333262
INVCM_IN_MEV = 0.1239841984


@pytest.mark.parametrize(
    ("units", "meV_per_unit"), [("meV", 1.0), ("K", K_IN_MEV), ("invcm", INVCM_IN_MEV)]
)
def test_units_converted_to_meV(units, meV_per_unit):
    """Parameters given in K or cm^-1 are stored in meV."""
    cfp = qc.CEFpars(dict(B20=2.0, B4m3=-0.5), units=units)
    assert cfp.Bpars["B20"] == pytest.approx(2.0 * meV_per_unit, rel=1e-5)
    assert cfp.Bpars["B4m3"] == pytest.approx(-0.5 * meV_per_unit, rel=1e-5)


def test_unknown_unit_raises():
    with pytest.raises(ValueError):
        qc.CEFpars(dict(B20=1.0), units="eV")


def test_unknown_operator_raises():
    with pytest.raises(KeyError):
        qc.CEFpars(dict(B80=1.0))


def test_allowed_parameters_tetragonal():
    """4mm (C4v): B20, B40, B44, B60, B64, as in the McPhase manual table."""
    assert qc.CEFpars.allowed_Bpars("4mm") == ["B20", "B40", "B44", "B60", "B64"]


@pytest.mark.parametrize(("hm", "schoenflies"), [("4mm", "C4v"), ("m-3m", "Oh"), ("-3m", "D3d")])
def test_allowed_parameters_same_for_both_notations(hm, schoenflies):
    assert qc.CEFpars.allowed_Bpars(hm) == qc.CEFpars.allowed_Bpars(schoenflies)


def test_triclinic_allows_every_operator():
    allowed = qc.CEFpars.allowed_Bpars("1")
    assert len(allowed) == 49 and "B00" in allowed and "B6m6" in allowed


def test_parameters_stored_as_python_floats():
    """NumPy scalars in the input are stored as Python floats (AUDIT Q6)."""
    cfp = qc.CEFpars(dict(B20=np.float32(0.5), B40=np.int64(2)), units="K")
    assert all(type(v) is float for v in cfp.Bpars.values())
