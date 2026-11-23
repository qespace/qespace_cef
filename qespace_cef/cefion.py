import numpy as np
from numpy.typing import ArrayLike

from . import config
from . import constants as C
from .cefmatrices import J_x, J_y, J_z, StevensBase
from .cefpars import CEFpars
from .ion import Ion


def _degenerate_levels(energies: ArrayLike) -> list[slice]:
    r"""The levels of CONV-0015, as slices of the ascending energies.

    Two consecutive energies belong to the same level when

    .. math:: E_{i+1} - E_i \le r_\mathrm{tol} \, \max(\max_k E_k - \min_k E_k, E_\mathrm{floor}),

    with ``rtol`` and ``E_floor`` read from :mod:`qespace_cef.config` at
    call time (defaults 10⁻¹⁰ and 1 meV). A level ends wherever two
    consecutive energies differ by more than that. The scale is the spread
    of the energies, so a shift of all of them (``B00``) does not change
    the grouping.

    Parameters
    ----------
    energies : array_like of float, shape (n,)
        Energies in meV, ascending, with the zero of CONV-0010, e.g.
        ``CEFion.energies``.

    Returns
    -------
    list of slice
        One slice per level, in ascending energy: ``energies[s]`` are the
        energies of a level, and ``s.stop - s.start`` is its degeneracy.
        Empty for no energies.

    Raises
    ------
    ValueError
        If ``energies`` is not one-dimensional, not finite or not ascending.

    Examples
    --------
    Ce3+ with ``B20`` alone has three Kramers doublets:

    >>> import qespace_cef as qc
    >>> from qespace_cef.cefion import _degenerate_levels
    >>> cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)))
    >>> [s.stop - s.start for s in _degenerate_levels(cef.energies)]
    [2, 2, 2]
    """
    E = np.asarray(energies, dtype=float)
    if E.ndim != 1:
        raise ValueError(f"energies must be one-dimensional, got shape {E.shape}.")
    if E.size == 0:
        return []
    if not np.all(np.isfinite(E)):
        raise ValueError("energies must be finite.")
    gaps = np.diff(E)
    if np.any(gaps < 0):
        raise ValueError("energies must be ascending.")
    # The spread, not max|E|: a shift of all energies (B00) leaves it unchanged.
    tolerance = config.rtol * max(float(E[-1] - E[0]), config.E_floor)
    bounds = [0, *(np.flatnonzero(gaps > tolerance) + 1).tolist(), E.size]
    return [slice(start, stop) for start, stop in zip(bounds[:-1], bounds[1:])]


def _time_reversal(v: np.ndarray) -> np.ndarray:
    """Θ|m> = (-1)^(J-m) |-m>, then complex conjugation (CONV-0011).

    ``v`` is a state in the basis |J>, ..., |-J>; component ``i`` has
    ``J - m = i``, and ``|-m>`` is component ``2J - i``.
    """
    return ((-1.0) ** np.arange(len(v)) * v.conj())[::-1]


def _fix_phase(v: np.ndarray) -> np.ndarray:
    """The state times the phase that makes its largest component real and positive.

    Components equal to the largest within ``config.rtol`` count as largest,
    and the first of them in basis order is used, so rounding cannot switch
    between two of the same size.
    """
    size = np.abs(v)
    k = int(np.flatnonzero(size >= size.max() * (1 - config.rtol))[0])
    return v * (size[k] / v[k])


def _diagonalise_inside(states: np.ndarray, operators: list, tolerance: float) -> np.ndarray:
    """Diagonalise ``operators[0]`` in the span of the rows of ``states``, ascending.

    Where its eigenvalues agree within ``tolerance``, the next operator is
    diagonalised inside that subspace (CONV-0011, rule 1). Returns the new
    orthonormal rows.
    """
    if len(states) == 1 or not operators:
        return states
    values, U = np.linalg.eigh(states.conj() @ operators[0] @ states.T)
    states = U.T @ states  # row mu = sum_k U[k, mu] states[k]
    bounds = [0, *(np.flatnonzero(np.diff(values) > tolerance) + 1).tolist(), len(values)]
    return np.vstack(
        [
            _diagonalise_inside(states[start:stop], operators[1:], tolerance)
            for start, stop in zip(bounds[:-1], bounds[1:])
        ]
    )


def _canonical_basis(states: np.ndarray, J: float, kramers: bool) -> np.ndarray:
    """The basis of CONV-0011 for one level, from any orthonormal basis of it.

    Parameters
    ----------
    states : ndarray of complex, shape (d, 2J+1)
        Orthonormal rows spanning the level, in the basis |J>, ..., |-J>.
    J : float
        Total angular momentum quantum number.
    kramers : bool
        True for a half-integer J in zero field: the states are then paired
        by time reversal (rule 2); otherwise each gets its own phase (rule 3).

    Returns
    -------
    ndarray of complex, shape (d, 2J+1)
        The canonical states, ordered by <J_z> ascending (then <J_x>, <J_y>
        where J_z is degenerate within ``config.rtol * J``).
    """
    states = _diagonalise_inside(states, [J_z(J), J_x(J), J_y(J)], config.rtol * J)
    if not kramers:
        return np.array([_fix_phase(v) for v in states])
    # The order is antisymmetric under time reversal, which flips J: the
    # partner of state i (<J_z> < 0) is state d - 1 - i.
    d = len(states)
    for i in range(d // 2):
        states[i] = _fix_phase(states[i])
        states[d - 1 - i] = _time_reversal(states[i])
    return states


class CEFion:
    r"""
    Rare-earth ion in a crystal electric field and a magnetic field.

    Builds the Hamiltonian of the ground J multiplet,

    .. math:: \hat{H} = \sum_{k,q} B_k^q \hat{O}_k^q - \mu_B g_J \, \vec{H} \cdot \hat{\vec{J}},

    and, by default, diagonalises it (:meth:`diagonalize`). The basis is
    :math:`|J\rangle, |J-1\rangle, \dots, |-J\rangle`, so the first vector
    component belongs to :math:`m = +J`. Energies are in meV, with the zero
    of CONV-0010: the eigenvalues of the Hamiltonian, not shifted, so with
    ``B00 = 0`` they average to zero (the centre of gravity of the J
    multiplet).

    Parameters
    ----------
    ion : Ion
        The free ion: J, :math:`g_J` and the form-factor coefficients.
    cfp : CEFpars
        Crystal-field parameters :math:`B_k^q` in meV (Stevens operators).
    Hfield : array_like of float, shape (3,), optional
        Magnetic field in T, in the frame of the J operators (z is the
        quantisation axis). Default ``(0, 0, 0)``.
    diagonalize : bool, optional
        If True (default), call :meth:`diagonalize` at once.

    Attributes
    ----------
    ion : Ion
        The ion passed in.
    Jval : float
        Total angular momentum quantum number J (e.g. 2.5). Not the operator
        ``J``.
    Hfield : ndarray, shape (3,)
        Magnetic field in T.
    cfp : CEFpars
        The crystal-field parameters passed in.
    hamiltonian : ndarray, shape (2J+1, 2J+1)
        Hamiltonian matrix in meV, in the free-ion basis.
    Jx, Jy, Jz : ndarray, shape (2J+1, 2J+1)
        Angular-momentum matrices (units of :math:`\hbar`). In the free-ion
        basis before :meth:`diagonalize`; afterwards in the eigenbasis,
        ``Jx[i, j]`` :math:`= \langle i|\hat{J}_x|j\rangle`.
    J : list or ndarray
        ``[Jx, Jy, Jz]``: a list before :meth:`diagonalize`, an ndarray of
        shape (3, 2J+1, 2J+1) after.
    energies : ndarray, shape (2J+1,)
        Eigenvalues of ``hamiltonian`` in meV, ascending, not shifted
        (CONV-0010). An empty list before diagonalisation.
    excitation_energies : ndarray, shape (2J+1,)
        ``energies - energies[0]`` in meV: the energy of each state above the
        ground state, starting at 0, as in a level scheme or a neutron
        spectrum. An empty list before diagonalisation.
    eigenvectors : ndarray, shape (2J+1, 2J+1)
        Row ``i`` is the eigenvector of ``energies[i]`` in the free-ion basis.
        An empty list before diagonalisation.
    freeionkets : list of str
        Labels of the free-ion basis states, e.g. ``'|5/2>'``, in basis order.

    Raises
    ------
    ValueError
        If ``Hfield`` is not three finite components.
    Exception
        If the Hamiltonian is not Hermitian (``numpy.allclose``).

    Warnings
    --------
    Inside a degenerate level the eigenvectors are an arbitrary orthonormal
    basis (CLAUDE.md, "Degenerate levels"), and so is anything computed from
    single eigenvectors of that level; sums over the level are not.

    Examples
    --------
    Ce3+ (J = 5/2) in an axial field splits into three Kramers doublets,
    at :math:`B_2^0 (3m^2 - J(J+1))` for :math:`|m| = 1/2, 3/2, 5/2`:

    >>> import qespace_cef as qc
    >>> cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)))
    >>> cef.energies.round(6).tolist()
    [-8.0, -8.0, -2.0, -2.0, 10.0, 10.0]
    >>> cef.excitation_energies.round(6).tolist()
    [0.0, 0.0, 6.0, 6.0, 18.0, 18.0]
    """

    def __init__(
        self,
        ion: Ion,
        cfp: CEFpars,
        Hfield: ArrayLike = (0.0, 0.0, 0.0),
        diagonalize: bool = True,
    ):
        self.ion = ion
        self.Jval = self.ion.J
        Jval = self.Jval

        self.Hfield = np.array(Hfield, dtype=float)
        if self.Hfield.shape != (3,) or not np.all(np.isfinite(self.Hfield)):
            raise ValueError(f"Hfield must be 3 finite components in T, got {Hfield!r}.")
        self.cfp = cfp

        # Prepare the rest of the fields based on the main input parameters

        # Main Ji matrices
        self.Jx = J_x(Jval)
        self.Jy = J_y(Jval)
        self.Jz = J_z(Jval)
        self.J = [self.Jx, self.Jy, self.Jz]

        # Prepare a list of kets that form the basis
        if Jval % 1 == 0:
            self.freeionkets = [f"|{int(x)}>" for x in np.linspace(Jval, -Jval, int(2 * Jval + 1))]
        elif Jval % 1 == 0.5:
            self.freeionkets = [
                f"|{int(2 * x):d}/2>" for x in np.linspace(Jval, -Jval, int(2 * Jval + 1))
            ]

        # THE HAMILTONIAN
        H = -C.uB * self.ion.gJ * np.einsum("ijk,i", self.J, self.Hfield)

        # Store Stevens operators in the dictionary containing pointers to functions
        StevensOp = StevensBase(self.Jval)
        for Bname, Bvalue in self.cfp.Bpars.items():
            H += Bvalue * StevensOp[Bname[1:]]

        self.hamiltonian = H
        if not np.allclose(H, H.T.conj()):
            raise Exception("ERROR, Hamiltonian is not Hermitian")

        # Uncomputed, empty fields
        self.eigenvectors = []
        self.energies = []
        self.excitation_energies = []

        # Diagonalize the Hamiltonian
        if diagonalize:
            self.diagonalize()

    def diagonalize(self, Eimag_precision: int = 10):
        r"""Diagonalise the Hamiltonian and rotate the J operators into its eigenbasis.

        Sets ``energies``, ``excitation_energies``, ``eigenvectors``, ``Jx``,
        ``Jy``, ``Jz`` and ``J`` (see the class attributes). After the call,
        state ``i`` of every matrix is the eigenstate with energy
        ``energies[i]``. The states are in the order of ``numpy.linalg.eigh``:
        energies ascending.

        Parameters
        ----------
        Eimag_precision : int, optional
            Number of decimals to which the Hamiltonian must equal its
            Hermitian conjugate (in meV). Default 10.

        Returns
        -------
        E : ndarray, shape (2J+1, 2J+1)
            Diagonal matrix of the eigenvalues in meV, ascending, the same as
            ``energies``.
        U : ndarray, shape (2J+1, 2J+1)
            Unitary matrix whose columns are the eigenvectors, with
            :math:`H = U E U^\dagger`.

        Raises
        ------
        ValueError
            If the Hamiltonian is not Hermitian to ``Eimag_precision`` decimals.
        """

        # A Hermitian Hamiltonian has real energies; eigh assumes Hermiticity, so check it first.
        H = self.hamiltonian
        if np.any(np.around(H - H.conj().T, Eimag_precision) != 0):
            raise ValueError(
                "The Hamiltonian is not Hermitian, its energies would be complex. "
                "Check `Eimag_precision` option."
            )

        # eigh returns the eigenvalues ascending and orthonormal eigenvectors,
        # also within degenerate subspaces, as columns of U, with
        # H = U diag(eigenvalues) U^dagger. No further sorting is needed.
        eigenvalues, U = np.linalg.eigh(H)
        E = np.diag(eigenvalues)

        self.energies = eigenvalues  # not shifted (CONV-0010)
        self.excitation_energies = eigenvalues - eigenvalues[0]
        self.eigenvectors = U.T  # row i belongs to energies[i]

        # Rotate the free-ion operators, built anew, into the eigenbasis
        # (`self.eigenvectors[i]` is the vector of the i-th energy). Not the
        # current Jx, Jy, Jz: after a first call they are already rotated.
        V = self.eigenvectors
        self.Jx = V.conj() @ J_x(self.Jval) @ V.T
        self.Jy = V.conj() @ J_y(self.Jval) @ V.T
        self.Jz = V.conj() @ J_z(self.Jval) @ V.T

        # Handy for calculations is the overall J operator, shape = (3, 2J+1, 2J+1)
        self.J = np.array([self.Jx, self.Jy, self.Jz])

        return E, U

    @property
    def magnetic_moment(self):
        r"""Magnetic moment of each eigenstate, in :math:`\mu_B`.

        :math:`\vec{m}_i = g_J \, \mathrm{Re} \langle i|\hat{\vec{J}}|i\rangle`,
        from the J operators in the eigenbasis, so valid only after
        :meth:`diagonalize`.

        Returns
        -------
        ndarray, shape (2J+1, 3)
            Row ``i`` is :math:`(m_x, m_y, m_z)` of eigenstate ``i``.

        Warnings
        --------
        Inside a degenerate level the moments of single states depend on the
        arbitrary basis that :meth:`diagonalize` picks; their sum over the
        level does not (CLAUDE.md, "Degenerate levels").
        """
        # After diagonalization the J operator is in the eigenstate basis,
        # so we just trace it.
        mu = self.ion.gJ * np.einsum("ijj->ji", self.J)
        return np.real(mu + mu.conj()) / 2

    def get_degeneracies(self) -> list[list]:
        """Energy levels and their degeneracies, grouped by CONV-0015.

        The ``energies`` are grouped into levels by
        :func:`_degenerate_levels`, with the tolerance from
        :mod:`qespace_cef.config` at call time.

        Returns
        -------
        list of [float, int]
            ``[energy, degeneracy]`` for each level, ascending in energy:
            the mean of the level's ``energies`` in meV (CONV-0010) and the
            number of states in it.
        """
        return [
            [float(np.mean(self.energies[level])), level.stop - level.start]
            for level in _degenerate_levels(self.energies)
        ]

    def summary(self, evecs: bool = True, precision: int = 3) -> str:
        """The level scheme as text: levels, degeneracies and eigenvectors.

        The levels are the ``excitation_energies``, measured from the ground
        state, whose energy ``energies[0]`` is given in the first line
        (CONV-0010). The states are grouped into levels by
        :func:`_degenerate_levels` (CONV-0015). Under each level its states
        are printed in the canonical basis of CONV-0011, written in the
        free-ion kets, leaving out coefficients smaller than
        ``10**-precision`` (one at the cutoff is kept):

        1. J_z is diagonalised inside the level, states ordered by <J_z>
           ascending; where J_z is degenerate there (within
           ``config.rtol * J``), J_x and then J_y are diagonalised inside
           that subspace.
        2. Half-integer J, zero field and no parameter of odd rank (a
           Hamiltonian even under time reversal): each state with <J_z> < 0
           gets the phase that makes its largest component real and
           positive, and its partner is its time reverse,
           Θ|m> = (-1)^(J-m) |-m> followed by complex conjugation.
        3. Otherwise each state gets that phase on its own.

        Every printed state lies inside its level.

        ``eigenvectors`` are not changed.

        Parameters
        ----------
        evecs : bool, optional
            If True (default), list the eigenvectors under each level.
        precision : int, optional
            Decimals of the energies (meV) and of the eigenvector
            coefficients. Default 3.

        Returns
        -------
        str
            One block per level, e.g. ``'E(1) =\\t6.000 meV\\t 2-fold
            degenerate'`` followed by its eigenvectors. Before
            :meth:`diagonalize`: ``'Undiagonalized Hamiltonian'``.

        Warnings
        --------
        Where J_z, J_x and J_y are all degenerate inside a level (no J matrix
        element splits it, e.g. a non-magnetic doublet), rule 1 does not fix
        the basis, and the printed states there remain an arbitrary basis of
        the level.
        """
        if not len(self.energies):
            return "Undiagonalized Hamiltonian"

        # Time reversal pairs the states only of a time-even Hamiltonian: no
        # field and no odd-rank parameter (O_lm with odd l is odd under it).
        odd_rank = any(
            int(name[1]) % 2 == 1 and value != 0 for name, value in self.cfp.Bpars.items()
        )
        kramers = self.Jval % 1 == 0.5 and not np.any(self.Hfield) and not odd_rank
        # Coefficients at the cutoff are kept, within config.rtol, so that a
        # pure ket survives precision = 0.
        cutoff = 10.0**-precision * (1 - config.rtol)

        lines = [
            f"Energy levels{' and eigenvectors' if evecs else ''}, energies above the "
            f"ground state (E_0 = {self.energies[0]:.{precision}f} meV):"
        ]
        for level, states in enumerate(_degenerate_levels(self.energies)):
            energy = self.excitation_energies[states.start]
            degeneracy = states.stop - states.start
            block = [f"E({level:d}) =\t{energy:.{precision}f} meV\t{degeneracy:2d}-fold degenerate"]
            if evecs:
                canonical = _canonical_basis(self.eigenvectors[states], self.Jval, kramers)
                for ev in canonical:
                    components = [
                        f"({c:.{precision}f}){ket}"
                        for c, ket in zip(ev, self.freeionkets)
                        if np.abs(c) >= cutoff
                    ]
                    block.append(f"ev_{level}: " + " + ".join(components))
            lines.append("\n".join(block))
        return "\n\n".join(lines) + "\n"
