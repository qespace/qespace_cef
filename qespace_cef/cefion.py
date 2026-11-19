from typing import Tuple

import numpy as np
from scipy.linalg import schur

from . import constants as C
from .cefmatrices import J_x, J_y, J_z, StevensBase
from .cefpars import CEFpars
from .ion import Ion


class CEFion:
    r"""
    Object representing a rare-earth ion in CF potential. It is internally calculated in the meV units.

    Parameters
    ----------
    ion: str
        Name of the ion. They are tabulated in :obj:`ion` with their parameters.
    Hfield: 1D array of floats
        External magnetic field in *T* units.
    cfp : ``qespace_cef.CEFpars``
        Crystal field parameters
    diagonalize : bool, optional
        If true (default) then it automatically diagonalizes Hamiltonian, calculates energy levels and sorts all matrices so that the first eigenvector corresponds to the lowest energy level and the last one to the highest.

    Examples:

    TODO

    Attributes
    ----------
    ion : ``qespace_cef:ion``
        The ``ion`` object that represents an isolated rare-earth ion.
    Jval : int/2
        The J value corresponding to the L+S quantum numbers. Not to be confused with the ``J`` operator.
    Hfield : array_like
        Vector in real space corresponding to the external magnetic field in T units.
    cfp : ``qespace_cef.reion.cfp``
        Crystal field parameters
    hamiltonian : ndarray
        Hamiltonian operator. :math:`\hat{\mathcal{H}} = \sum_{ij} B_i^j \hat{O}_i^j + g_J (H_x \hat{J}_x + H_y \hat{J}_y + H_z \hat{J}_z)`
    Jx, Jy, Jz : ndarray
        Matrices representing the prinicpal quantum operators :math:`\hat{J}_\\alpha`` with :math:`\\alpha=[x,y,z]`.
    J : ndarray
        Total angular momentum vector operator :math:`\hat{J}=[\hat{J}_x, \hat{J}_y, \hat{J}_z]`. It's multidimensionality takes some time to get used to.
    energies : ndarray
        Eigenenergies of the Hamiltonian.
    eigenvectors : ndarray
        Eigenvectors of the Hamiltonian.
    degeneracies : list
        List containing entries like [energy, degeneracy], which describes how degenerated is the given energy level.
    freeionkets : list
        List of kets corresponding to the J basis of the free ion problem.
    """

    def __init__(
        self,
        ion: Ion,
        cfp: CEFpars,
        Hfield: Tuple[float, float, float] = [0, 0, 0],
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

        # Prepare a list os kets that form the basis
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
            raise Exception("ERROR, Hamiltonian is not hermitian")

        # Uncomputed, empty fields
        self.eigenvectors = []
        self.energies = []

        # Diagonalize the Hamiltonian
        if diagonalize:
            self.diagonalize()

    def diagonalize(
        self, sortWithE: bool = True, shiftToZero: bool = True, Eimag_precision: int = 10
    ):
        """Diagonalize the Hamiltonian, and change to the sorted eigenvector base.
        The default sorting is done according to eigenenergies, so that the first vector [1,0,...,0]
        is the lowest eigenstate, and the last one [0,...,0,1] is the highest one.
        Changing the base greatly faiclitates further calculations based on the evaluation of matrix elements.

        Updates
        -------
        ``Jx``,``Jy``,``Jz``,``J``,``energies``,``eigenvectors``,``degeneracies``

        Returns
        -------
            E, U
                Products of diagonalization, H = U E U^T.

        Raises
        ------
            ValueError
                If the calculated eigenenergies are not real.

        """

        # Diagonalize the Hamiltonian
        # Other functions than scipy.linalg.schur do not produce perpendicular subspaces for degenerated eigenvalues
        E, U = schur(self.hamiltonian)
        eigenvalues = np.diag(E)
        # With schur definition H = U D Ut we take columns as eigenvectors
        eigenvectors = U.T  ### TRIPLE CHECK

        # Check if energies are largely imaginary
        if np.any(np.abs(np.imag(np.around(eigenvalues, Eimag_precision))) > 0):
            raise ValueError(
                f"Final energies are complex! {eigenvalues} Check `Eimag_precision` option."
            )

        energies = np.real(eigenvalues) - int(shiftToZero) * min(
            np.real(eigenvalues)
        )  # shift to zero level and make real

        # Sorting
        if sortWithE:
            sortedIndices = energies.argsort()
        else:
            sortedIndices = np.arange(self.Jval)

        # Assign main values
        self.energies = energies[sortedIndices]
        self.eigenvectors = eigenvectors[sortedIndices]

        # Change the basis of principal operators to the eigenstate basis with specified sorting scheme
        # However obtained, keep notation that `self.eigenvectors[i]` is a vector of ith eigenvalue.
        self.Jx = self.eigenvectors.conj() @ self.Jx @ self.eigenvectors.T
        self.Jy = self.eigenvectors.conj() @ self.Jy @ self.eigenvectors.T
        self.Jz = self.eigenvectors.conj() @ self.Jz @ self.eigenvectors.T

        # Handy for calculations is the overall J operator, shape = (3, 2J+1, 2J+1)
        self.J = np.array([self.Jx, self.Jy, self.Jz])

        return E, U

    @property
    def magnetic_moment(self):
        """Expectation value of magnetic moment operator.

        :math:`m = gJ * [<Jx>, <Jy>, <Jz>]`"""
        # After diagonalization the J operator is in the eigenstate basis,
        # so we just trace it.
        mu = self.ion.gJ * np.einsum("ijj->ji", self.J)
        return np.real(mu + mu.conj()) / 2

    def get_degeneracies(self, precision: int = 5):
        """List of energy levels and their degeneracies.
        Follows the given precision, default is 5 decimal places."""
        precision = 5
        energies = (self.energies * 10**precision).astype(int)  # Convert to int for unique
        unique_energies, counts = np.unique(energies, return_counts=True)
        return [[e, c] for e, c in zip(unique_energies / 10**precision, counts)]

    def __str__(self, evecs: bool = True, precision: int = 3):
        """
        Nice printout of calculated parameters

        Parameters:
            precision : int, optional
                How many significant digits should be shown for printout of energy and eigenvector coefficients.
        """
        # After refactoring the precision there might be cleanup needed

        if not len(self.energies):
            return "Undiagonalized Hamiltonian"

        ret = ""
        ret += "Energy levels and corresponding eigenfunctions:\n"

        n = 0
        levels = []
        for level, x in enumerate(self.get_degeneracies(precision=precision)):
            level_str = ""
            energy = x[0]
            degeneracy = int(x[1])
            level_str += (
                f"E({level:d}) =\t{energy:.{precision}} meV\t{degeneracy:2d}fold-degenerated\n"
            )

            # List degenerated eigenvectors
            if evecs:
                for ev in self.eigenvectors[n : n + degeneracy]:
                    ev_components = []
                    for c, ket in zip(ev, self.freeionkets):
                        if np.abs(c) > 1 / np.power(10, precision):  # Arbitrary tolerance
                            ev_components.append(f"({c:.{precision}f}){ket}")

                    level_str += f"ev_{level}: " + " + ".join(ev_components) + "\n"

            levels.append(level_str)
            n += degeneracy

        ret += "\n".join(levels)
        return ret
