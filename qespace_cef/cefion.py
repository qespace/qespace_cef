import numpy as np
from numpy.typing import ArrayLike

from . import constants as C
from .cefmatrices import J_x, J_y, J_z, StevensBase
from .cefpars import CEFpars
from .ion import Ion


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

        self.Hfield = np.array(Hfield)
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

    def get_degeneracies(self, precision: int = 5):
        """Energy levels and their degeneracies.

        Eigenenergies are grouped by their value truncated to 5 decimals
        (in meV).

        Parameters
        ----------
        precision : int, optional
            Meant as the number of decimals; currently ignored, 5 is used.

        Returns
        -------
        list of [float, int]
            ``[energy, degeneracy]`` for each level, ascending in energy;
            energy in meV.

        Warnings
        --------
        P7 (AUDIT.md), not fixed yet: ``precision`` is overwritten with 5,
        and the energies are truncated, not rounded. Two degenerate states
        whose energies differ by rounding noise can fall on both sides of a
        1e-5 meV step and come out as two levels, e.g. the 28.584 meV doublet
        of Pr3+ with ``B40=0.01, B44=0.05, B60=1e-4, B64=-21e-4``.
        """
        # AUDIT P7: overwrites the argument; astype(int) truncates instead of rounding.
        precision = 5
        energies = (self.energies * 10**precision).astype(int)  # Convert to int for unique
        unique_energies, counts = np.unique(energies, return_counts=True)
        return [[float(e), int(c)] for e, c in zip(unique_energies / 10**precision, counts)]

    def summary(self, evecs: bool = True, precision: int = 3) -> str:
        """Energy levels and eigenvectors as text, with options. Not implemented yet.

        Placeholder for a print of the level scheme: the levels, their
        degeneracies and eigenvectors, with the ``evecs`` and ``precision``
        options (AUDIT Q7, ``ROADMAP.md``). ``CEFion`` has no ``str()`` of
        its own (CONV-0010); what the print shows is decided with it.

        Raises
        ------
        NotImplementedError
            Always, for now.
        """
        # AUDIT Q7: stub, to be implemented with A2.
        raise NotImplementedError(
            "CEFion.summary, a print of the level scheme, is planned (AUDIT Q7, ROADMAP.md)."
        )
