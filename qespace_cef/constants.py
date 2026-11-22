"""Physical constants and the unit conversions derived from them.

All values are CODATA 2022 (CONV-0003): P. J. Mohr, D. B. Newell,
B. N. Taylor, and E. Tiesinga, "CODATA recommended values of the fundamental
physical constants: 2022", Rev. Mod. Phys. 97, 025002 (2025),
doi:10.1103/RevModPhys.97.025002; companion article J. Phys. Chem. Ref.
Data 54, 033105 (2025). Values at
https://physics.nist.gov/cuu/Constants/index.html, read on 2026-11-21.

The conversions are computed from these values, not typed separately.
"""

# CODATA 2022, SI units.
_e = 1.602176634e-19  # C, elementary charge, exact (SI 2019)
_h = 6.62607015e-34  # J s, Planck constant, exact (SI 2019)
_c = 299792458.0  # m/s, speed of light in vacuum, exact (SI 2019)
_k = 1.380649e-23  # J/K, Boltzmann constant, exact (SI 2019)
uB_SI = 9.2740100657e-24  # J/T, Bohr magneton, standard uncertainty 0.0000000029e-24

# Conversions used by the package.
K2meV = 1e3 * _k / _e  # meV/K, k_B
meV2K = 1 / K2meV  # K/meV, 1 / k_B
uB = 1e3 * uB_SI / _e  # meV/T, Bohr magneton
invcm2meV = 1e3 * 100 * _h * _c / _e  # meV per cm^-1, h c

# Neutron scattering (CONV-0007), CODATA 2022.
gamma_n = 1.91304276  # |mu_n / mu_N|, neutron magnetic moment in nuclear magnetons, -1.91304276(45)
r0 = 2.8179403205e-15  # m, classical electron radius, standard uncertainty 0.0000000013e-15
gamma_r0_sq = (gamma_n * r0) ** 2 / 1e-28  # barn, (gamma r0)^2, about 0.2906
