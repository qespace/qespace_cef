"""Tests of ``qespace_cef.constants`` against CODATA 2022 (AUDIT P15).

The reference values are written here from CODATA 2022, read from NIST
(https://physics.nist.gov/cuu/Constants/index.html, "Source: 2022 CODATA"),
not imported from the module under test. P. J. Mohr, D. B. Newell,
B. N. Taylor, and E. Tiesinga, Rev. Mod. Phys. 97, 025002 (2025),
doi:10.1103/RevModPhys.97.025002.
"""

import numpy as np
import pytest

from qespace_cef import constants as C

# CODATA 2022, SI units.
E_CHARGE = 1.602176634e-19  # C, exact (SI 2019)
PLANCK = 6.62607015e-34  # J s, exact (SI 2019)
LIGHT = 299792458.0  # m/s, exact (SI 2019)
BOLTZMANN = 1.380649e-23  # J/K, exact (SI 2019)
BOHR_MAGNETON = 9.2740100657e-24  # J/T, uncertainty 0.0000000029e-24

# NIST's own derived value: mu_B = 5.7883817982(18)e-5 eV/T.
BOHR_MAGNETON_EV = 5.7883817982e-5  # eV/T
BOHR_MAGNETON_EV_UNCERTAINTY = 0.0000000018e-5  # eV/T


def test_K2meV_and_meV2K_are_reciprocal():
    """``meV2K`` is ``1 / K2meV``, so converting there and back is the identity."""
    assert abs(C.K2meV * C.meV2K - 1.0) <= np.finfo(float).eps


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("K2meV", 1e3 * BOLTZMANN / E_CHARGE),  # meV/K
        ("meV2K", E_CHARGE / (1e3 * BOLTZMANN)),  # K/meV
        ("uB", 1e3 * BOHR_MAGNETON / E_CHARGE),  # meV/T
        ("invcm2meV", 1e3 * 100 * PLANCK * LIGHT / E_CHARGE),  # meV per cm^-1
    ],
)
def test_conversions_from_codata_2022(name, expected):
    """The conversions follow from the CODATA 2022 values, to rounding."""
    assert getattr(C, name) == pytest.approx(expected, rel=1e-15, abs=0)


def test_uB_matches_nist_in_eV():
    """``uB`` (meV/T) agrees with NIST's mu_B in eV/T within its uncertainty."""
    assert abs(C.uB * 1e-3 - BOHR_MAGNETON_EV) <= BOHR_MAGNETON_EV_UNCERTAINTY


# CODATA 2022, neutron scattering (CONV-0007).
NEUTRON_MOMENT = -1.91304276  # mu_n / mu_N, uncertainty 0.00000045
ELECTRON_RADIUS = 2.8179403205e-15  # m, classical electron radius, uncertainty 0.0000000013e-15
BARN = 1e-28  # m^2


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("gamma_n", abs(NEUTRON_MOMENT)),
        ("r0", ELECTRON_RADIUS),  # m
        ("gamma_r0_sq", (NEUTRON_MOMENT * ELECTRON_RADIUS) ** 2 / BARN),  # barn
    ],
)
def test_neutron_constants_from_codata_2022(name, expected):
    """gamma, r0 and (gamma r0)^2 in barn follow from CODATA 2022, to rounding."""
    assert getattr(C, name) == pytest.approx(expected, rel=1e-15, abs=0)
