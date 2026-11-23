"""Tests of ``qespace_cef.CEFpars``."""

import warnings

import numpy as np
import pytest
from reference import mcphase

import qespace_cef as qc

# CODATA 2022 (NIST, physics.nist.gov/cuu/Constants), all exact: e = 1.602176634e-19 C,
# k = 1.380649e-23 J/K, h = 6.62607015e-34 J s, c = 299792458 m/s.
# 1 K = 1e3 k / e = 0.08617333 meV and 1 cm^-1 = 1e3 * 100 h c / e = 0.1239842 meV.
K_IN_MEV = 1e3 * 1.380649e-23 / 1.602176634e-19
INVCM_IN_MEV = 1e3 * 100 * 6.62607015e-34 * 299792458.0 / 1.602176634e-19


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
    """NumPy scalars in the input are stored as Python floats (CLAUDE.md 5.7.3)."""
    cfp = qc.CEFpars(dict(B20=np.float32(0.5), B40=np.int64(2)), units="K")
    assert all(type(v) is float for v in cfp.Bpars.values())


def _label(rank, m):
    """The CONV-0001 label of B_lm: B20, B43, B4m3, ..."""
    return f"B{rank}{m}" if m >= 0 else f"B{rank}m{-m}"


@pytest.mark.parametrize("ion", ["Yb", "Pr", "Ho"])  # J = 7/2, 4, 8
@pytest.mark.parametrize("lm", mcphase.LM, ids=[_label(*lm) for lm in mcphase.LM])
def test_conv_0001_labels(ion, lm):
    """CONV-0001: the parameter ``Blm`` multiplies O_lm, ``B4m3`` the sine-type O_4,-3.

    The Hamiltonian matrix is compared, not its energies: a rotation about z by
    pi/m turns O_lm into -O_lm and leaves the spectrum unchanged, so the
    energies cannot see the sign of a single B_lm with m != 0.
    """
    value = 0.3  # meV, an arbitrary test value
    cef_ion = qc.CEFion(qc.Ion(ion), qc.CEFpars({_label(*lm): value}))
    expected = value * mcphase.stevens(cef_ion.Jval, *lm)
    np.testing.assert_allclose(
        cef_ion.hamiltonian, expected, rtol=0, atol=1e-10 * np.abs(expected).max()
    )


@pytest.mark.parametrize("label", ["B10", "B3m2", "B55"])
def test_odd_rank_parameters_warn(label):
    """Odd-rank parameters are accepted with a UserWarning (Q12, P16): parity removes them."""
    with pytest.warns(UserWarning, match="odd rank"):
        cfp = qc.CEFpars({label: 0.1})
    assert cfp.Bpars == {label: pytest.approx(0.1)}


def test_even_rank_parameters_do_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        qc.CEFpars(dict(B00=1.0, B20=0.1, B4m3=0.01, B66=1e-4))


@pytest.mark.parametrize("label", ["O20", "b20", "20", "X44"])
def test_parameter_labels_start_with_B(label):
    """Only B... labels name parameters (M6): O20 was a second spelling of B20."""
    with pytest.raises(ValueError):
        qc.CEFpars({label: 1.0})
