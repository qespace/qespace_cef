"""Observables computed from a :class:`~qespace_cef.CEFion`.

Units: energies in meV, temperatures in K, magnetic fields in T, momentum
transfer in 1/Å, moments in :math:`\\mu_B` per ion. Functions that take a
``CEFion`` expect it diagonalised (the default when it is created).
"""

from typing import Tuple

import numpy as np

from . import constants as C
from .cefion import CEFion


def boltzman_population(energies: np.ndarray, temperature: float) -> np.ndarray:
    r"""Boltzmann population of each eigenstate.

    .. math:: p_i = \frac{1}{Z} e^{-E_i/k_B T}, \qquad Z = \sum_i e^{-E_i/k_B T}

    The sum runs over the given states, so pass every eigenvalue of the
    Hilbert space, each degenerate state separately; ``Z`` is computed here.

    Parameters
    ----------
    energies : ndarray of float, shape (n,)
        Eigenenergies in meV.
    temperature : float
        Temperature in K.

    Returns
    -------
    ndarray of float, shape (n,)
        Occupation probability of each state; they sum to 1.
    """

    p = np.exp(-np.array(energies - min(energies)) * C.meV2K / temperature)
    Z = sum(p)
    return p / Z


def spectral_weight(
    cefion: CEFion,
    temperature: float,
    Q: Tuple[float, float, float],
    scheme: str,
    Ei: float = 1e6,
    B_iso: float = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Neutron-scattering energies and spectral weights of every transition.

    For each pair of eigenstates, initial ``i`` and final ``j``, the code
    computes

    .. math::

        \Delta E_{ij} = E_j - E_i, \qquad
        S_{ij} = f^2(Q) \, e^{-B_{iso} Q^2 / 3} \,
                 \sqrt{1 - \Delta E_{ij} / E_\mathrm{inc}} \; p_i \; M_{ij},

    with :math:`f(Q)` the form factor (:meth:`Ion.mff`), :math:`p_i` the
    Boltzmann population (:func:`boltzman_population`), :math:`E_\mathrm{inc}` the
    incident neutron energy ``Ei``, and the matrix elements

    - ``'powder'``: :math:`M_{ij} = \frac{2}{3} \sum_\alpha |\langle i|\hat{J}_\alpha|j\rangle|^2`;
    - ``'single-crystal'``: :math:`M_{ij} = \sum_\alpha (1 - \hat{Q}_\alpha^2) |\langle i|\hat{J}_\alpha|j\rangle|^2`,
      with :math:`\hat{Q} = \vec{Q}/|Q|` in the frame of the J operators.

    The prefactor :math:`N (\gamma r_0)^2` is set to 1, so ``Sqw`` is in
    arbitrary units. Positive :math:`\Delta E` is neutron energy loss. The
    transitions are neither sorted nor merged over degenerate levels.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field.
    temperature : float
        Temperature in K.
    Q : array_like of float, shape (3,)
        Momentum transfer in 1/Å, in the frame of the J operators. For
        ``'powder'`` only its length is used.
    scheme : {'powder', 'single-crystal'}
        How the matrix elements are averaged, see above.
    Ei : float, optional
        Incident neutron energy in meV. The default, 1e6 meV, makes
        :math:`k_f/k_i \approx 1`.
    B_iso : float, optional
        Isotropic displacement parameter in Å², entering as
        :math:`e^{-B_{iso} Q^2/3}`. Default 0 (no Debye-Waller factor).

    Returns
    -------
    jumps : ndarray, shape (2J+1, 2J+1)
        ``jumps[i, j]`` is the transition energy :math:`E_j - E_i` in meV.
    Sqw : ndarray, shape (2J+1, 2J+1)
        ``Sqw[i, j]`` is the spectral weight of the transition from ``i`` to
        ``j``, in arbitrary units. NaN, with a NumPy ``RuntimeWarning``,
        where :math:`E_j - E_i > E_\mathrm{inc}`, and for ``'single-crystal'`` at
        Q = 0.

    Raises
    ------
    ValueError
        If ``scheme`` is unknown, or ``Q`` does not have 3 components for
        ``'single-crystal'``.

    Warnings
    --------
    Known problems, not fixed yet (AUDIT.md):

    - P9: the single-crystal matrix elements leave out the cross terms
      :math:`-\hat{Q}_\alpha \hat{Q}_\beta \langle i|\hat{J}_\alpha|j\rangle \langle j|\hat{J}_\beta|i\rangle`
      of the full projector. For Q along (110) in the audit's example the
      summed weight was 5.6516 instead of 5.6573; single transitions can be
      off by more.
    - D5: :math:`e^{-B_{iso} Q^2/3}` is not the crystallographic
      Debye-Waller factor :math:`e^{-B Q^2/16\pi^2}`; the convention is open
      (AUDIT C10).
    - The weights of single transitions inside a degenerate level depend on
      the arbitrary eigenbasis; sums over the level do not.

    Notes
    -----
    The original code cites "the Enderle book following Stephane Raymond
    article" for the formula; the full reference is missing.
    """

    # Magnetic form factor
    f2m = cefion.ion.mff(np.linalg.norm(Q)) ** 2

    # Debye-Waller factor
    eDW = np.exp(-1 / 3 * np.linalg.norm(Q) ** 2 * B_iso)

    # Tricky way to create a 2D array of energies associated with transitions between levels
    jumps = cefion.energies - cefion.energies[np.newaxis].T

    # Calculate the |<\Gamma_f|J_\perp|\Gamma_i>|^2 matrix
    if scheme == "single-crystal":
        if len(Q) != 3:
            raise ValueError("Dimension of the `Q` vector is not 3")

        J2 = np.square(np.abs(cefion.J))  # J is already diagonalized and in the eigenvector basis
        # AUDIT P9: only the diagonal of the projector; the cross terms are missing.
        projection = 1 - (Q / np.linalg.norm(Q)) ** 2

        J2_perp = np.einsum("i,ijk", projection, J2)
    elif scheme == "powder":
        J2_perp = 2.0 / 3 * np.einsum("ijk->jk", np.square(np.abs(cefion.J)))
    else:
        raise ValueError(f"Unknown scheme {scheme!r}; use 'powder' or 'single-crystal'.")

    # kf/ki factor, which is actually a matrix
    kfki = np.sqrt(1 - jumps / Ei)

    # Occupation
    transition_probs = boltzman_population(cefion.energies, temperature)

    # Multiply the factors, vectors and matrices properly to get spectral weight.
    Sqw = f2m * eDW * kfki * J2_perp * transition_probs[:, np.newaxis]

    return (jumps, Sqw)


def magnetization(
    cefion: CEFion, temperature: float, Hfield: Tuple[float, float, float]
) -> Tuple[float, float, float]:
    r"""Thermal average of the magnetic moment of one ion in a field.

    Builds and diagonalises a new :class:`~qespace_cef.CEFion` from
    ``cefion.ion`` and ``cefion.cfp`` with the field ``Hfield``, then

    .. math:: M_\alpha = g_J \sum_n p_n \, \mathrm{Re} \langle n|\hat{J}_\alpha|n\rangle

    with :math:`p_n` the Boltzmann populations in that field.

    Parameters
    ----------
    cefion : CEFion
        Ion in the crystal field. Its own ``Hfield`` is not used.
    temperature : float
        Temperature in K.
    Hfield : array_like of float, shape (3,)
        Total magnetic field in T, in the frame of the J operators.

    Returns
    -------
    tuple of 3 float
        :math:`(M_x, M_y, M_z)` in :math:`\mu_B` per ion.
    """

    # Solve the Hamiltonian of a copy of the given ion
    cefion_inH = CEFion(ion=cefion.ion, cfp=cefion.cfp, Hfield=Hfield, diagonalize=True)

    # The diagonalized Hamiltonians' operators are already transformed into the sorted eigenvector base
    p = boltzman_population(cefion_inH.energies, temperature)
    M = np.einsum("ij,i", cefion_inH.magnetic_moment, p)

    return tuple(M)


def magnetization_exchange(
    cefion: CEFion,
    temperature: float,
    Hfield: Tuple[float, float, float],
    lam: float,
    precision: float = 1e-5,
) -> Tuple[float, float, float]:
    r"""Magnetisation with a molecular (mean) field, solved self-consistently.

    Calls :func:`magnetization` with the effective field

    .. math:: \vec{H}_\mathrm{eff} = \vec{H} + \lambda \vec{M}

    and repeats with the new :math:`\vec{M}` until
    :math:`|\vec{M}_{k+1} - \vec{M}_k| / |\vec{M}_{k+1}| \le` ``precision``.

    Parameters
    ----------
    cefion : CEFion
        Ion in the crystal field. Its own ``Hfield`` is not used.
    temperature : float
        Temperature in K.
    Hfield : array_like of float, shape (3,)
        External magnetic field in T, in the frame of the J operators.
    lam : float
        Molecular-field constant :math:`\lambda` in T per :math:`\mu_B`.
    precision : float, optional
        Relative change of :math:`|\vec{M}|` at which the iteration stops.
        Default 1e-5.

    Returns
    -------
    tuple of 3 float
        :math:`(M_x, M_y, M_z)` in :math:`\mu_B` per ion.

    Warnings
    --------
    P11 (AUDIT.md), not fixed yet: there is no iteration limit, so the loop
    can run forever when the mean field does not converge. M = 0 gives NaN
    and ends the loop without an error.
    """

    # Solve the Hamiltonian of a copy of the given ion
    M = [magnetization(cefion, temperature, Hfield)]
    residual = 42

    # AUDIT P11: no iteration limit; M = 0 gives NaN and ends the loop.
    while residual > precision:
        Hexchange = lam * np.array(M[-1])
        M_next = magnetization(cefion, temperature, Hfield + Hexchange)
        M.append(M_next)
        dM = [m1 - m2 for m1, m2 in zip(M[-1], M[-2])]
        residual = np.linalg.norm(dM) / np.linalg.norm(M[-1])

    return M[-1]


def susceptibility(
    cefion: CEFion,
    temperatures: np.ndarray,
    Hfield_direction: Tuple[float, float, float],
    method: str = "perturbation",
) -> np.ndarray:
    r"""Single-ion magnetic susceptibility along a direction, versus temperature.

    The two methods return **different units**.

    ``'perturbation'`` (default) evaluates, in zero field,

    .. math::

        \chi = g_J^2 \sum_n p_n \Big[
            \sum_{m:\,E_m = E_n} \frac{|J^{\hat{n}}_{nm}|^2}{k_B T}
          + \sum_{m:\,E_m > E_n} \frac{1 - e^{-(E_m - E_n)/k_B T}}{E_m - E_n} |J^{\hat{n}}_{nm}|^2
        \Big],
        \qquad
        |J^{\hat{n}}_{nm}|^2 = \frac{\sum_\alpha n_\alpha |\langle n|\hat{J}_\alpha|m\rangle|^2}{|\vec{n}|},

    with :math:`\vec{n}` = ``Hfield_direction``, :math:`p_n` the Boltzmann
    populations and levels equal within ``constants.numerical_zero`` (meV).
    The result is in :math:`\mu_B^2`/meV per ion; multiply by
    ``constants.uB`` (meV/T) for :math:`\mu_B`/T.

    ``'magnetization'`` takes the finite difference
    :math:`\chi = |\vec{M}(\epsilon \hat{n})| / \epsilon` with
    :math:`\epsilon = 10^{-5}` T and :func:`magnetization`, in
    :math:`\mu_B`/T per ion.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field; its ``Hfield`` is not used.
    temperatures : array_like of float, shape (n,)
        Temperatures in K.
    Hfield_direction : array_like of float, shape (3,)
        Field direction in the frame of the J operators; its length does
        not matter.
    method : {'perturbation', 'magnetization'}, optional
        See above.

    Returns
    -------
    ndarray, shape (n,)
        Susceptibility at each temperature, in :math:`\mu_B^2`/meV
        (``'perturbation'``) or :math:`\mu_B`/T (``'magnetization'``) per ion.

    Raises
    ------
    ValueError
        If ``method`` is unknown.

    Warnings
    --------
    Known problems, not fixed yet (AUDIT.md):

    - P1 (``'perturbation'``): the van Vleck sum takes only the pairs with
      :math:`E_m > E_n`; the pairs with :math:`E_m < E_n` contribute the
      same, so that term is half what it should be. Ce3+ with
      ``B20=0.5, B40=0.01, B44=0.1``, field along z, 5 K: 0.12704 instead
      of 0.13130 :math:`\mu_B`/T (after multiplying by ``constants.uB``).
    - P2 (``'perturbation'``): the direction enters as
      :math:`\sum_\alpha n_\alpha |J_\alpha|^2 / |\vec{n}|` instead of
      :math:`|\hat{n} \cdot \hat{\vec{J}}|^2`; along (1, 1, 0) this gives
      :math:`\sqrt{2}` times the value along x.
    - P10 (``'magnetization'``): returns :math:`|\vec{M}|/\epsilon`, not the
      component along the field, so the sign is lost and the result is
      wrong when :math:`\vec{M}` is not parallel to the field.

    Notes
    -----
    ``'perturbation'`` prints the matrix of transition energies.
    """

    if method == "magnetization":
        susceptibility = np.zeros(len(temperatures))
        eps = 1e-5

        for it, temperature in enumerate(temperatures):
            Hfield = eps * np.array(Hfield_direction) / np.linalg.norm(Hfield_direction)
            M = magnetization(cefion, temperature, Hfield)
            # AUDIT P10: |M|/eps, not the component along the field.
            susceptibility[it] = np.linalg.norm(M) / eps
    elif method == "perturbation":
        susceptibility = np.empty(len(temperatures))

        # Tricky way to create a 2D array of energies associated with transitions between levels
        jumps = cefion.energies - cefion.energies[np.newaxis].T

        # Clean up
        jumps_with_zero_energy = np.where(np.abs(jumps) < C.numerical_zero)
        # AUDIT P1: only the E_m > E_n pairs; the van Vleck term comes out halved.
        jumps_with_positive_energy = np.where(jumps > 0)

        # Calculate the J^2 matrix
        J2 = np.square(np.abs(cefion.J))
        # AUDIT P2: sum_a n_a |J_a|^2 / |n| instead of |n.J|^2 / |n|^2.
        J2_directed = np.einsum("i,ijk", Hfield_direction, J2) / np.linalg.norm(Hfield_direction)

        print(jumps)
        for it, temperature in enumerate(temperatures):
            # Define the transition matrix taking into account the diagonal and negative energy jumps
            Tmx = np.zeros(np.shape(jumps))
            Tmx[jumps_with_zero_energy] = 1 * C.meV2K / temperature
            Tmx[jumps_with_positive_energy] = (
                -np.expm1(-jumps[jumps_with_positive_energy] * C.meV2K / temperature)
                / jumps[jumps_with_positive_energy]
            )

            # Include the Boltzmann factor
            transition_probs = boltzman_population(cefion.energies, temperature)
            Tmx = Tmx * transition_probs[:, np.newaxis]

            susceptibility[it] = (cefion.ion.gJ) ** 2 * np.sum(Tmx * J2_directed)
    else:
        raise ValueError(f"Unknown method {method!r}; use 'perturbation' or 'magnetization'.")

    return susceptibility


def gtensor(cefion: CEFion) -> np.ndarray:
    r"""g-tensor of the ground doublet, following PyCrystalField's construction.

    With :math:`|0\rangle, |1\rangle` the first two eigenstates, column by
    column:

    .. math::

        g_{\alpha x} = 2 g_J \, \mathrm{Re} \langle 0|\hat{J}_\alpha|1\rangle, \quad
        g_{\alpha y} = 2 g_J \, \mathrm{Im} \langle 0|\hat{J}_\alpha|1\rangle, \quad
        g_{\alpha z} = 2 g_J \, |\langle 0|\hat{J}_\alpha|0\rangle|.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion whose two lowest states are degenerate.

    Returns
    -------
    ndarray, shape (3, 3)
        Dimensionless; row :math:`\alpha` = x, y, z of the J operator.

    Raises
    ------
    AssertionError
        If the two lowest energies differ by ``constants.numerical_zero``
        (meV) or more.

    Warnings
    --------
    P4 (AUDIT.md), not fixed yet: the result depends on the arbitrary basis
    that :meth:`CEFion.diagonalize` picks inside the ground doublet, against
    CLAUDE.md, "Degenerate levels". For Ce3+ with ``B20=1`` in zero field it
    is currently ``diag(2.571, -2.571, 0.857)``; a different doublet basis
    gives a different matrix. How g should be defined is open (AUDIT C7).
    """

    assert np.abs(cefion.energies[0] - cefion.energies[1]) < C.numerical_zero, (
        "Ground state is not a doublet, g-tensor is not defined."
    )
    Jeff = cefion.J  # J operator projected to the ground state doublet

    # PyCrystalField. AUDIT P4: depends on the basis picked inside the doublet.
    gx = Jeff[:, 0, 1].real
    gy = Jeff[:, 0, 1].imag
    gz = np.abs(Jeff[:, 0, 0])

    return np.array([gx, gy, gz]).T * cefion.ion.gJ * 2


def thermodynamics(
    cefion: CEFion, T: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Single-ion thermodynamics from the partition function.

    With :math:`E_n` the eigenenergies (meV) and :math:`\beta = 1/k_B T`,
    the code computes

    .. math::

        Z = \sum_n e^{-\beta E_n}, \quad
        E = \langle E \rangle = \frac{1}{Z} \sum_n E_n e^{-\beta E_n}, \quad
        F = -T \ln Z, \quad
        S = \frac{E - F}{T}, \quad
        C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{T^2}.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field.
    T : ndarray of float, shape (n,)
        Temperatures in K.

    Returns
    -------
    Z : ndarray, shape (n,)
        Partition function, dimensionless.
    E : ndarray, shape (n,)
        Internal energy in meV per ion.
    F : ndarray, shape (n,)
        Free energy :math:`-T \ln Z`, in K per ion.
    S : ndarray, shape (n,)
        :math:`(E - F)/T`: meV/K and dimensionless terms mixed.
    Cv : ndarray, shape (n,)
        Heat capacity in meV²/K² per ion.

    Warnings
    --------
    P3 (AUDIT.md), not fixed yet: energies in meV are mixed with
    temperatures in K, so ``S`` and ``Cv`` are not in units of
    :math:`k_B`. In the audit's example, at 10 K: S = 1.054 instead of
    1.293 :math:`k_B`, and Cv low by a factor of about 135 (11.6²).
    """
    Z = np.zeros(len(T))
    E = np.zeros(len(T))
    E2 = np.zeros(len(T))

    for En in cefion.energies:
        Z += np.exp(-En * C.meV2K / T)
        E += En * np.exp(-En * C.meV2K / T)
        E2 += En * En * np.exp(-En * C.meV2K / T)

    E /= Z
    E2 /= Z

    # AUDIT P3: F in K, E in meV; S and Cv mix the two.
    F = -T * np.log(Z)
    S = (E - F) / T
    Cv = (E2 - E**2) / T**2

    return Z, E, F, S, Cv
