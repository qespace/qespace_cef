"""Observables computed from a :class:`~qespace_cef.CEFion`.

Units: energies in meV, temperatures in K, magnetic fields in T, momentum
transfer in 1/Å, moments in :math:`\\mu_B` per ion. Functions that take a
``CEFion`` expect it diagonalised (the default when it is created).
"""

import numpy as np
from numpy.typing import ArrayLike

from . import constants as C
from .cefion import CEFion, degenerate_levels


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


def scattering_function(cefion: CEFion, temperature: float) -> tuple[np.ndarray, np.ndarray]:
    r"""The magnetic scattering function of one ion as excitations (CONV-0007, CONV-0012).

    .. math::

        S^{\alpha\beta}(\hbar\omega) = \sum_k W_k^{\alpha\beta}\,
            \delta(\hbar\omega - \varepsilon_k), \qquad
        \varepsilon_{if} = E_f - E_i, \quad
        W_{if}^{\alpha\beta} = p_i \langle i|\hat{J}_\alpha|f\rangle
            \langle f|\hat{J}_\beta|i\rangle,

    one entry per pair of eigenstates (i, f), elastic ones included, with
    :math:`p_i` the Boltzmann populations. Positive :math:`\varepsilon` is
    neutron energy loss. No neutron factor (form factor, :math:`k_f/k_i`,
    :math:`(\gamma r_0)^2`) is applied; see :func:`cross_section_crystal`
    and :func:`cross_section_powder`.

    The entries form an unordered set (CONV-0012): they carry no level
    information, may share an energy and may have zero weight, and only
    sums over the entries at one energy (CONV-0008) are meaningful. (Here
    entry ``k`` is the pair ``(i, f) = divmod(k, n)``.)

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field.
    temperature : float
        Temperature in K.

    Returns
    -------
    energies : ndarray of float, shape (n**2,)
        Excitation energies :math:`\varepsilon_k` in meV, n = 2J + 1.
    weights : ndarray of complex, shape (n**2, 3, 3)
        Scattering weights :math:`W_k^{\alpha\beta}`, dimensionless, α and β
        in the frame of the J operators (CONV-0004). Each is Hermitian and
        positive semidefinite, and :math:`\sum_k \mathrm{Tr}\,W_k = J(J+1)`.
    """
    E = cefion.energies
    J = np.asarray(cefion.J)  # (3, n, n): J[a, i, f] = <i|J_a|f>
    n = len(E)
    p = boltzman_population(E, temperature)
    energies = (E[np.newaxis, :] - E[:, np.newaxis]).reshape(n * n)
    weights = np.einsum("i,aif,bfi->ifab", p, J, J).reshape(n * n, 3, 3)
    return energies, weights


def _cross_section_unprojected(
    cefion: CEFion,
    temperature: float,
    Qnorm: np.ndarray,
    P: np.ndarray,
    Ei: float,
    u2: float,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Cross-section weights of CONV-0007 for a given projector P.

    The shared kernel of :func:`cross_section_crystal` and
    :func:`cross_section_powder`, which build P and check Q:

    .. math::

        w_k = (\gamma r_0)^2 \, (k_f/k_i)_k \, [\tfrac{1}{2} g_J F(|Q|)]^2
              e^{-|Q|^2 \langle u^2 \rangle / 3}
              \sum_{\alpha\beta} P_{\alpha\beta} W_k^{\alpha\beta},
        \qquad (k_f/k_i)_k = \sqrt{1 - \varepsilon_k / E_i},

    NaN where :math:`\varepsilon_k > E_i`. ``Qnorm`` has the loop shape S
    (Å⁻¹), ``P`` the shape S + (3, 3); the weights get S + (n²,), in barn
    sr⁻¹ per ion.
    """
    if not np.isfinite(Ei) or Ei <= 0:
        raise ValueError(f"Ei must be a positive incident energy in meV, got {Ei!r}.")
    energies, W = scattering_function(cefion, temperature)
    ion_factor = (cefion.ion.gJ / 2) ** 2 * np.asarray(cefion.ion.mff(Qnorm)) ** 2
    debye_waller = np.exp(-(Qnorm**2) * u2 / 3)
    ratio = 1 - energies / Ei
    kf_ki = np.sqrt(np.where(ratio >= 0, ratio, np.nan))  # NaN where the neutron cannot lose it
    projected = np.einsum("...ab,kab->...k", P, W).real  # P real symmetric, W Hermitian
    prefactor = C.gamma_r0_sq * ion_factor * debye_waller
    return energies, prefactor[..., np.newaxis] * kf_ki * projected


def cross_section_crystal(
    cefion: CEFion,
    temperature: float,
    Q: ArrayLike,
    Ei: float,
    u2: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Magnetic neutron cross section of one ion in a single crystal (CONV-0007).

    The cross section as excitations (CONV-0012),
    :math:`d^2\sigma/d\Omega dE = \sum_k w_k \, \delta(\hbar\omega - \varepsilon_k)`,
    with the scattering weights :math:`W_k` of :func:`scattering_function`
    and

    .. math::

        w_k = (\gamma r_0)^2 \, (k_f/k_i)_k \, [\tfrac{1}{2} g_J F(|Q|)]^2
              e^{-|Q|^2 \langle u^2 \rangle / 3}
              \sum_{\alpha\beta} (\delta_{\alpha\beta} - \hat{Q}_\alpha \hat{Q}_\beta)
              W_k^{\alpha\beta},

    :math:`(k_f/k_i)_k = \sqrt{1 - \varepsilon_k / E_i}` at fixed incident
    energy, :math:`F` the dipole form factor (:meth:`Ion.mff`), and
    :math:`(\gamma r_0)^2` from ``constants.gamma_r0_sq`` (CODATA 2022).
    Positive :math:`\varepsilon_k` is neutron energy loss. The sum of the
    weights in an energy window is dσ/dΩ of that window.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field.
    temperature : float
        Temperature in K.
    Q : array_like of float, shape S + (3,)
        Momentum transfer in Å⁻¹, in the local frame of the J operators
        (CONV-0004); any loop shape S. Not zero.
    Ei : float
        Incident neutron energy in meV.
    u2 : float, optional
        Mean-square displacement :math:`\langle u^2 \rangle` in Å², in the
        Debye–Waller factor :math:`e^{-|Q|^2 \langle u^2 \rangle/3}`.
        Default 0.

    Returns
    -------
    energies : ndarray of float, shape (n**2,)
        Excitation energies in meV, n = 2J + 1, as in
        :func:`scattering_function`.
    weights : ndarray of float, shape S + (n**2,)
        Cross-section weights in barn sr⁻¹ per ion. NaN, without a
        warning, where the excitation energy exceeds ``Ei``.

    Raises
    ------
    ValueError
        If the last axis of ``Q`` is not 3, any Q is zero, or ``Ei`` is not
        positive.
    """
    Q = np.asarray(Q, dtype=float)
    if Q.shape[-1:] != (3,):
        raise ValueError(f"Q must have shape (..., 3), got {Q.shape}.")
    Qnorm = np.linalg.norm(Q, axis=-1)
    if np.any(Qnorm == 0):
        raise ValueError("Q must not be zero: the direction of Q sets the projector.")
    q = Q / Qnorm[..., np.newaxis]
    P = np.eye(3) - q[..., :, np.newaxis] * q[..., np.newaxis, :]
    return _cross_section_unprojected(cefion, temperature, Qnorm, P, Ei, u2)


def cross_section_powder(
    cefion: CEFion,
    temperature: float,
    Qnorm: ArrayLike,
    Ei: float,
    u2: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Magnetic neutron cross section of one ion in a powder (CONV-0007).

    As :func:`cross_section_crystal`, with the projector averaged over the
    directions of Q, :math:`\langle \delta_{\alpha\beta} - \hat{Q}_\alpha
    \hat{Q}_\beta \rangle = \tfrac{2}{3} \delta_{\alpha\beta}`:

    .. math::

        w_k = (\gamma r_0)^2 \, (k_f/k_i)_k \, [\tfrac{1}{2} g_J F(|Q|)]^2
              e^{-|Q|^2 \langle u^2 \rangle / 3}
              \tfrac{2}{3} \sum_\alpha W_k^{\alpha\alpha}.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field.
    temperature : float
        Temperature in K.
    Qnorm : array_like of float, shape S
        Length of the momentum transfer in Å⁻¹, any loop shape S.
    Ei : float
        Incident neutron energy in meV.
    u2 : float, optional
        Mean-square displacement in Å². Default 0.

    Returns
    -------
    energies : ndarray of float, shape (n**2,)
        Excitation energies in meV, n = 2J + 1.
    weights : ndarray of float, shape S + (n**2,)
        Cross-section weights in barn sr⁻¹ per ion. NaN, without a
        warning, where the excitation energy exceeds ``Ei``.

    Raises
    ------
    ValueError
        If any ``Qnorm`` is negative or not finite, or ``Ei`` is not
        positive.
    """
    Qnorm = np.asarray(Qnorm, dtype=float)
    if np.any(~np.isfinite(Qnorm) | (Qnorm < 0)):
        raise ValueError("Qnorm must be finite and not negative.")
    P = np.broadcast_to(2 / 3 * np.eye(3), Qnorm.shape + (3, 3))
    return _cross_section_unprojected(cefion, temperature, Qnorm, P, Ei, u2)


def magnetization(
    cefion: CEFion, temperature: float, Hfield: ArrayLike
) -> tuple[float, float, float]:
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

    # The diagonalized Hamiltonian's operators are already transformed into the sorted eigenvector base
    p = boltzman_population(cefion_inH.energies, temperature)
    M = np.einsum("ij,i", cefion_inH.magnetic_moment, p)

    return tuple(float(m) for m in M)


def magnetization_exchange(
    cefion: CEFion,
    temperature: float,
    Hfield: ArrayLike,
    lam: float,
    precision: float = 1e-5,
    max_iterations: int = 1000,
) -> tuple[float, float, float]:
    r"""Magnetisation with a molecular (mean) field, solved self-consistently.

    Calls :func:`magnetization` with the effective field

    .. math:: \vec{H}_\mathrm{eff} = \vec{H} + \lambda \vec{M}

    and repeats with the new :math:`\vec{M}` until
    :math:`|\vec{M}_{k+1} - \vec{M}_k| \le` ``precision`` :math:`\cdot
    \max(|\vec{M}_{k+1}|, |\vec{M}_k|)`. A zero moment, as in zero field,
    converges at once.

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
    max_iterations : int, optional
        Number of mean-field steps after which the iteration gives up.
        Default 1000.

    Returns
    -------
    tuple of 3 float
        :math:`(M_x, M_y, M_z)` in :math:`\mu_B` per ion.

    Raises
    ------
    RuntimeError
        If the moment has not converged after ``max_iterations`` steps.
    """

    Hfield = np.asarray(Hfield, dtype=float)
    M = magnetization(cefion, temperature, Hfield)
    for _ in range(max_iterations):
        M_next = magnetization(cefion, temperature, Hfield + lam * np.asarray(M))
        change = np.linalg.norm(np.subtract(M_next, M))
        if change <= precision * max(np.linalg.norm(M_next), np.linalg.norm(M)):
            return M_next
        M = M_next

    raise RuntimeError(
        f"The mean field did not converge in {max_iterations} steps "
        f"(last change of |M|: {float(change):.3g} muB)."
    )


def susceptibility(
    cefion: CEFion,
    temperatures: ArrayLike,
    method: str = "perturbation",
) -> np.ndarray:
    r"""Single-ion magnetic susceptibility tensor in zero field, versus temperature.

    The susceptibility of CONV-0009,
    :math:`\chi_{\alpha\beta} = \partial \langle m_\alpha \rangle / \partial B_\beta`
    at :math:`B = 0`, with :math:`m = g_J \mu_B J`, in :math:`\mu_B`/T per
    ion. Both methods return the same quantity.

    ``'perturbation'`` (default) evaluates Van Vleck's formula,

    .. math::

        \chi_{\alpha\beta} = g_J^2 \mu_B \sum_{n,m} W_{nm}\,
            \mathrm{Re}\,\langle n|\hat{J}_\alpha|m\rangle\langle m|\hat{J}_\beta|n\rangle,
        \qquad
        W_{nm} = \begin{cases}
            p_n / k_B T & E_m = E_n \
            (p_n - p_m) / (E_m - E_n) & E_m \ne E_n
        \end{cases}

    with :math:`p_n` the Boltzmann populations, both sums over all states,
    energies in meV, :math:`\mu_B` in meV/T, and :math:`E_m = E_n` when the
    two states belong to one level of CONV-0008
    (:func:`qespace_cef.cefion.degenerate_levels`). The van Vleck
    weight is evaluated as
    :math:`p_n (1 - e^{-(E_m - E_n)/k_B T}) / (E_m - E_n)` on the pairs
    :math:`E_m > E_n`, with ``expm1``, and mirrored, :math:`W_{mn} = W_{nm}`,
    so it stays accurate as two levels approach each other (CONV-0008).

    ``'magnetization'`` takes the finite difference
    :math:`\chi_{\alpha\beta} = M_\alpha(\epsilon \hat{e}_\beta) / \epsilon`
    with :math:`\epsilon = 10^{-5}` T and :func:`magnetization`, one field
    direction per axis.

    The susceptibility along a unit vector :math:`\hat{n}` is
    ``n @ chi @ n``; the powder average is ``np.trace(chi, axis1=-2,
    axis2=-1) / 3``.

    Parameters
    ----------
    cefion : CEFion
        Diagonalised ion in the crystal field; its ``Hfield`` is not used.
    temperatures : array_like of float, shape S
        Temperatures in K, any shape (a scalar gives ``S = ()``).
    method : {'perturbation', 'magnetization'}, optional
        See above.

    Returns
    -------
    ndarray of float, shape S + (3, 3)
        :math:`\chi_{\alpha\beta}` in :math:`\mu_B`/T per ion, at each
        temperature; α is the moment component and β the field component,
        both in the frame of the J operators. Real and symmetric.

    Raises
    ------
    ValueError
        If ``method`` is unknown.
    """

    temperatures = np.asarray(temperatures, dtype=float)
    # One row per temperature, in the order of temperatures.ravel(); reshaped
    # to S + (3, 3) only on return, so no array is filled through a view.
    chi = np.empty((temperatures.size, 3, 3))

    if method == "magnetization":
        eps = 1e-5  # T
        for it, temperature in enumerate(temperatures.ravel()):
            for beta in range(3):
                Hfield = np.zeros(3)
                Hfield[beta] = eps
                chi[it, :, beta] = np.array(magnetization(cefion, temperature, Hfield)) / eps
    elif method == "perturbation":
        J = np.asarray(cefion.J)  # (3, n, n), in the eigenbasis
        # jumps[n, m] = E_m - E_n; degenerate pairs are the pairs of states in
        # one level of CONV-0008 (CONV-0009 defers to it).
        jumps = cefion.energies[np.newaxis, :] - cefion.energies[:, np.newaxis]
        level = np.empty(len(cefion.energies), dtype=np.int64)
        for k, states in enumerate(degenerate_levels(cefion.energies)):
            level[states] = k
        degenerate = level[:, np.newaxis] == level[np.newaxis, :]
        # The pairs (n, m) with E_m > E_n that are not degenerate; their
        # mirror images (m, n) are all the other pairs that are not.
        up_n, up_m = np.nonzero((jumps > 0) & ~degenerate)
        up_jumps = jumps[up_n, up_m]

        for it, temperature in enumerate(temperatures.ravel()):
            # W[n, m]: Curie p_n / k_B T for degenerate pairs, van Vleck
            # (p_n - p_m) / (E_m - E_n) for every other pair, in both orders.
            # The van Vleck weight is p_n (1 - exp(-(E_m - E_n) / k_B T)) /
            # (E_m - E_n), with expm1 on the pairs E_m > E_n only: it does
            # not cancel for a pair near the degeneracy rule (CONV-0008),
            # and its argument is never positive, so it cannot overflow.
            # W is symmetric, so the pairs E_m < E_n take W[m, n] = W[n, m].
            p = boltzman_population(cefion.energies, temperature)
            W = np.where(degenerate, p[:, np.newaxis] * C.meV2K / temperature, 0.0)
            w_up = p[up_n] * -np.expm1(-up_jumps * C.meV2K / temperature) / up_jumps
            W[up_n, up_m] = w_up
            W[up_m, up_n] = w_up
            # <n|J_a|m><m|J_b|n> = J_a[n, m] conj(J_b[n, m]) for Hermitian J_b
            chi[it] = np.einsum("nm,anm,bnm->ab", W, J, J.conj()).real
        chi *= cefion.ion.gJ**2 * C.uB
    else:
        raise ValueError(f"Unknown method {method!r}; use 'perturbation' or 'magnetization'.")

    return chi.reshape(temperatures.shape + (3, 3))


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
        If the two lowest energies are not one level of CONV-0008
        (:func:`qespace_cef.cefion.degenerate_levels`).

    Warnings
    --------
    Known problems, not fixed yet (AUDIT.md):

    - P4: the result depends on the arbitrary basis that
      :meth:`CEFion.diagonalize` picks inside the ground doublet, against
      CLAUDE.md, "Degenerate levels". For Ce3+ with ``B20=1`` in zero field it
      is currently ``diag(2.571, -2.571, 0.857)``; a different doublet basis
      gives a different matrix.
    - P17: the absolute value in the z column makes even :math:`g g^T`
      basis dependent, so the principal g values are wrong in a general
      doublet basis.

    How g should be defined is open (AUDIT C7).
    """

    ground = degenerate_levels(cefion.energies)[0]  # CONV-0008
    assert ground.stop - ground.start >= 2, (
        "Ground state is not a doublet, g-tensor is not defined."
    )
    Jeff = cefion.J  # J operator projected to the ground state doublet

    # PyCrystalField. AUDIT P4: depends on the basis picked inside the doublet.
    gx = Jeff[:, 0, 1].real
    gy = Jeff[:, 0, 1].imag
    # AUDIT P17: abs() makes even g g^T depend on the doublet basis.
    gz = np.abs(Jeff[:, 0, 0])

    return np.array([gx, gy, gz]).T * cefion.ion.gJ * 2


def thermodynamics(
    cefion: CEFion, T: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
