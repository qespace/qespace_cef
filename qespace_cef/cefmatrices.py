r"""Angular-momentum matrices and Stevens operators of a J multiplet.

All matrices are in the basis :math:`|J\rangle, |J-1\rangle, \dots, |-J\rangle`
(the first component belongs to the highest m; for J = 3/2,
:math:`|3/2\rangle = [1, 0, 0, 0]`), in units of :math:`\hbar = 1`, with the
Condon–Shortley phases (CONV-0005). The matrices of :math:`\hat{J}_x`,
:math:`\hat{J}_y`, :math:`\hat{J}_z`, :math:`\hat{J}_\pm` built by the
functions below are the root of every matrix in the package;
:class:`StevensBase` takes them from these functions. The operators are the Stevens operators
:math:`\hat{O}_k^q` as tabulated by McPhase (CONV-0001).

Additions by Michal Stekiel.
"""

import numpy as np
from numpy import diag, dot, eye, linspace, sqrt
from numpy.linalg import matrix_power as mp
from numpy.typing import NDArray


# J matrices
def _check_J(J: float) -> None:
    """Raise ValueError unless J is a finite, non-negative multiple of 1/2."""
    if not (np.isfinite(J) and J >= 0 and 2 * J == round(2 * J)):
        raise ValueError(f"J must be a non-negative multiple of 1/2, got {J!r}.")


def J_z(J: float) -> NDArray[np.float64]:
    r"""Matrix of :math:`\hat{J}_z`, shape (2J+1, 2J+1), real diagonal.

    Raises ValueError unless J is a non-negative multiple of 1/2.
    """
    _check_J(J)
    return diag(linspace(J, -J, int(2 * J + 1)))


def J_y(J: float) -> NDArray[np.complex128]:
    r"""Matrix of :math:`\hat{J}_y = (\hat{J}_+ - \hat{J}_-)/2i`, shape (2J+1, 2J+1)."""
    return -0.5j * (J_plus(J) - J_minus(J))


def J_x(J: float) -> NDArray[np.float64]:
    r"""Matrix of :math:`\hat{J}_x = (\hat{J}_+ + \hat{J}_-)/2`, shape (2J+1, 2J+1), real."""
    return 0.5 * (J_plus(J) + J_minus(J))


def J_plus(J: float) -> NDArray[np.float64]:
    r"""Matrix of the raising operator :math:`\hat{J}_+`, shape (2J+1, 2J+1), real.

    The only non-zero elements are on the first superdiagonal,
    :math:`\langle m+1|\hat{J}_+|m\rangle = \sqrt{J(J+1) - m(m+1)}` for
    m = J-1, ..., -J, the m of the lower state (CONV-0005). Raises
    ValueError unless J is a non-negative multiple of 1/2; ``J_x``, ``J_y``
    and ``J_minus`` are built from it.
    """
    _check_J(J)
    m = linspace(J - 1, -J, int(2 * J))
    return diag(sqrt(J * (J + 1) - m * (m + 1)), 1)


def J_minus(J: float) -> NDArray[np.float64]:
    r"""Matrix of the lowering operator :math:`\hat{J}_- = \hat{J}_+^\dagger`, shape (2J+1, 2J+1)."""
    return J_plus(J).conj().transpose()


class StevensBase:
    r"""Stevens operators :math:`\hat{O}_k^q` of one J multiplet.

    Every operator with k = 0 ... 6 and q = -k ... k is a method, named
    ``O{k}{q}`` for q >= 0 and ``O{k}m{|q|}`` for q < 0 (``O20``, ``O4m3``).
    Index the object with the name without the ``O`` to get the same
    matrix: ``StevensBase(J)['4m3']`` equals ``StevensBase(J).O4m3()``.
    Each call returns a matrix of shape (2J+1, 2J+1), dimensionless.

    Parameters
    ----------
    J : float
        Total angular momentum quantum number, integer or half-integer.

    Attributes
    ----------
    J : float
        The J passed in.
    JJ : float
        J(J+1).
    J2p1 : int
        Dimension 2J+1.
    Jx, Jy, Jz, Jp, Jm : ndarray, shape (2J+1, 2J+1)
        Angular-momentum matrices from :func:`J_x`, :func:`J_y`,
        :func:`J_z`, :func:`J_plus`, :func:`J_minus` (CONV-0005); ``Jy`` is
        complex, the others real.

    Examples
    --------
    >>> from qespace_cef.cefmatrices import StevensBase
    >>> ops = StevensBase(1)
    >>> ops["20"].tolist()
    [[1.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, 1.0]]
    >>> ops.O20().tolist() == ops["20"].tolist()
    True

    Notes
    -----
    Observations (1, 2 and 4 are checked in ``tests/test_cefmatrices.py``):

    1. O(n, n-1) = 0.5 (O(n-1,n-1) Jz + Jz O(n-1,n-1)).
    2. O(n, -(n-1)) = 0.5 (O(n-1,-(n-1)) Jz + Jz O(n-1,-(n-1))).
    3. Generally, O(n,k) is 0.5 (Okk Y + Y Okk), where Y is a polynomial in
       Jz of order n-k.
    4. All Stevens matrices are Hermitian; the ones with q < 0 are not
       symmetric and purely imaginary.

    Sources:

    1. https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node132.html
    2. https://easyspin.org/easyspin/documentation/stevensoperators.html and
       references therein.

    Performance:

    1. Faster than a dictionary as in PyCrystalField, as forming an operator
       goes directly to its method.
    2. The defining matrices (Jx, Jy, Jz, Jp, Jm) are built once per object,
       so build one object per ion and take all operators from it.
    """

    J: float
    JJ: float
    J2p1: int
    Jx: NDArray[np.float64]
    Jy: NDArray[np.complex128]
    Jz: NDArray[np.float64]

    def __init__(self, J: float):
        self.J = J
        self.JJ = J * (J + 1.0)
        self.J2p1 = int(2 * J + 1)

        self.Jz = J_z(J)
        self.Jp = J_plus(J)
        self.Jm = J_minus(J)
        self.Jx = J_x(J)
        self.Jy = J_y(J)

    @classmethod
    def implementedOps(cls):
        """Names of the implemented operators, e.g. ``['O00', 'O11', ...]`` (49 of them)."""
        return [Oname for Oname in cls.__dict__.keys() if Oname[0] == "O"]

    def __getitem__(self, nm: str) -> NDArray[np.complex128]:
        """Matrix of the operator named ``'O' + nm``, e.g. ``ops['4m3']``."""
        return self.__class__.__dict__[f"O{nm}"](self)

    ##################################################################
    def O00(self):
        return eye(self.J2p1)

    ##################################################################
    def O11(self):
        return self.Jx

    def O10(self):
        return self.Jz

    def O1m1(self):
        return self.Jy

    ##################################################################
    # Helper functions for higher order terms
    def _JpJm(self, k):
        return mp(self.Jp, k) + mp(self.Jm, k)

    def _JpmJm(self, k):
        return mp(self.Jp, k) - mp(self.Jm, k)

    ##################################################################
    def O2m2(self):
        return dot(self.Jx, self.Jy) + dot(self.Jy, self.Jx)

    def O2m1(self):
        return 0.5 * (dot(self.Jy, self.Jz) + dot(self.Jz, self.Jy))

    def O20(self):
        return 3 * dot(self.Jz, self.Jz) - self.JJ * eye(self.J2p1)

    def O21(self):
        return 0.5 * (dot(self.Jx, self.Jz) + dot(self.Jz, self.Jx))

    def O22(self):
        # Faster than mp(Jx, 2)
        return dot(self.Jx, self.Jx) - dot(self.Jy, self.Jy)

    ##################################################################
    def O3m3(self):
        return -0.5j * self._JpmJm(3)

    def O3m2(self):
        # O3m2 = 0.5* (O2m2*Jz + Jz*O2m2)
        O2m2 = self.O2m2()
        return 0.5 * (dot(O2m2, self.Jz) + dot(self.Jz, O2m2))

    def O3m1(self):
        # O3m1 = 0.5* (Jy*Y + Y*Jy)
        # Y = 5Jz**2 -J(J+1) -0.5
        Y = 5 * dot(self.Jz, self.Jz) - (self.JJ + 0.5) * eye(self.J2p1)
        return 0.5 * (dot(self.Jy, Y) + dot(Y, self.Jy))

    def O30(self):
        return 5 * mp(self.Jz, 3) - (3 * self.JJ - 1) * self.Jz

    def O31(self):
        # O31 = 0.5* (Jx*Y + Y*Jx)
        # Y = 5Jz**2 -J(J+1) -0.5
        Y = 5 * dot(self.Jz, self.Jz) - (self.JJ + 0.5) * eye(self.J2p1)
        return 0.5 * (dot(self.Jx, Y) + dot(Y, self.Jx))

    def O32(self):
        # O32 = 0.5* (O22*Jz + Jz*O22)
        O22 = self.O22()
        return 0.5 * (dot(O22, self.Jz) + dot(self.Jz, O22))

    def O33(self):
        return 0.5 * self._JpJm(3)

    ##################################################################
    def O4m4(self):
        return -0.5j * self._JpmJm(4)

    def O4m3(self):
        O3m3 = self.O3m3()
        return 0.5 * (dot(O3m3, self.Jz) + dot(self.Jz, O3m3))

    def O4m2(self):
        # O4m2 = 0.5* (O2m2*Y + Y*O2m2)
        # Y = 7Jz**2 - J(J+1) -5
        O2m2 = self.O2m2()
        Y = 7 * dot(self.Jz, self.Jz) - (self.JJ + 5) * eye(self.J2p1)
        return 0.5 * (dot(O2m2, Y) + dot(Y, O2m2))

    def O4m1(self):
        # O4m1 = 0.5* (Jy*Y + Y*Jy)
        # Y = 7Jz**3 - (3*J(J+1) + 1) Jz
        Y = 7 * mp(self.Jz, 3) - (3 * self.JJ + 1) * self.Jz
        return 0.5 * (dot(self.Jy, Y) + dot(Y, self.Jy))

    def O40(self):
        return (
            35 * mp(self.Jz, 4)
            - (30 * self.JJ - 25) * mp(self.Jz, 2)
            + (3 * self.JJ**2 - 6 * self.JJ) * eye(self.J2p1)
        )

    def O41(self):
        # O41 = 0.5* (Jx*Y + Y*Jx)
        # Y = 7Jz**3 - (3*J(J+1) + 1) Jz
        Y = 7 * mp(self.Jz, 3) - (3 * self.JJ + 1) * self.Jz
        return 0.5 * (dot(self.Jx, Y) + dot(Y, self.Jx))

    def O42(self):
        # O42 = 0.5* (O22*Y + Y*O22)
        # Y = 7Jz**2 - J(J+1) -5
        O22 = self.O22()
        Y = 7 * dot(self.Jz, self.Jz) - (self.JJ + 5) * eye(self.J2p1)
        return 0.5 * (dot(O22, Y) + dot(Y, O22))

    def O43(self):
        O33 = self.O33()
        return 0.5 * (dot(O33, self.Jz) + dot(self.Jz, O33))

    def O44(self):
        return 0.5 * self._JpJm(4)

    ##################################################################
    def O55(self):
        return 0.5 * self._JpJm(5)

    def O5m5(self):
        return -0.5j * self._JpmJm(5)

    def O54(self):
        O44 = self.O44()
        return 0.5 * (dot(O44, self.Jz) + dot(self.Jz, O44))

    def O5m4(self):
        O4m4 = self.O4m4()
        return 0.5 * (dot(O4m4, self.Jz) + dot(self.Jz, O4m4))

    def O53(self):
        # O53 = 0.5* (O33*Y + Y*O33)
        # Y = 9*Jz**2 - J(J+1) - 33/2
        O33 = self.O33()
        Y = 9 * dot(self.Jz, self.Jz) - (self.JJ + 33.0 / 2) * eye(self.J2p1)
        return 0.5 * (dot(O33, Y) + dot(Y, O33))

    def O5m3(self):
        # O5m3 = 0.5* (O3m3*Y + Y*O3m3)
        # Y = 9*Jz**2 - J(J+1) - 33/2
        O3m3 = self.O3m3()
        Y = 9 * dot(self.Jz, self.Jz) - (self.JJ + 33.0 / 2) * eye(self.J2p1)
        return 0.5 * (dot(O3m3, Y) + dot(Y, O3m3))

    def O52(self):
        # O52 = 0.5* (O22*Y + Y*O22)
        # Y = 3*Jz**3 - (J(J+1) + 6)*Jz
        O22 = self.O22()
        Y = 3 * mp(self.Jz, 3) - (self.JJ + 6) * self.Jz
        return 0.5 * (dot(O22, Y) + dot(Y, O22))

    def O5m2(self):
        # O5m2 = 0.5* (O2m2*Y + Y*O2m2)
        # Y = 3*Jz**3 - (J(J+1) + 6)*Jz
        O2m2 = self.O2m2()
        Y = 3 * mp(self.Jz, 3) - (self.JJ + 6) * self.Jz
        return 0.5 * (dot(O2m2, Y) + dot(Y, O2m2))

    def O51(self):
        # O51 = 0.5* (Jx*Y + Y*Jx)
        # Y = 21*Jz**4 - 14*JJ*Jz**2 + JJ**2 - JJ + 3/2
        Y = (
            21 * mp(self.Jz, 4)
            - 14 * self.JJ * mp(self.Jz, 2)
            + (self.JJ**2 - self.JJ + 1.5) * eye(self.J2p1)
        )
        return 0.5 * (dot(self.Jx, Y) + dot(Y, self.Jx))

    def O5m1(self):
        # O51 = 0.5* (Jy*Y + Y*Jy)
        # Y = 21*Jz**4 - 14*JJ*Jz**2 + JJ**2 - JJ + 3/2
        Y = (
            21 * mp(self.Jz, 4)
            - 14 * self.JJ * mp(self.Jz, 2)
            + (self.JJ**2 - self.JJ + 1.5) * eye(self.J2p1)
        )
        return 0.5 * (dot(self.Jy, Y) + dot(Y, self.Jy))

    def O50(self):
        return (
            63 * mp(self.Jz, 5)
            - (70 * self.JJ - 105) * mp(self.Jz, 3)
            + (15 * self.JJ**2 - 50 * self.JJ + 12) * self.Jz
        )

    ##################################################################
    def O66(self):
        return 0.5 * self._JpJm(6)

    def O6m6(self):
        return -0.5j * self._JpmJm(6)

    def O65(self):
        O55 = self.O55()
        return 0.5 * (dot(O55, self.Jz) + dot(self.Jz, O55))

    def O6m5(self):
        O5m5 = self.O5m5()
        return 0.5 * (dot(O5m5, self.Jz) + dot(self.Jz, O5m5))

    def O64(self):
        O44 = self.O44()
        Y = 11 * dot(self.Jz, self.Jz) - (self.JJ + 38) * eye(self.J2p1)
        return 0.5 * (dot(O44, Y) + dot(Y, O44))

    def O6m4(self):
        O4m4 = self.O4m4()
        Y = 11 * dot(self.Jz, self.Jz) - (self.JJ + 38) * eye(self.J2p1)
        return 0.5 * (dot(O4m4, Y) + dot(Y, O4m4))

    def O63(self):
        O33 = self.O33()
        Y = 11 * mp(self.Jz, 3) - (3 * self.JJ + 59) * self.Jz
        return 0.5 * (dot(O33, Y) + dot(Y, O33))

    def O6m3(self):
        O3m3 = self.O3m3()
        Y = 11 * mp(self.Jz, 3) - (3 * self.JJ + 59) * self.Jz
        return 0.5 * (dot(O3m3, Y) + dot(Y, O3m3))

    def O62(self):
        O22 = self.O22()
        Y = (
            33 * mp(self.Jz, 4)
            - (18 * self.JJ + 123) * dot(self.Jz, self.Jz)
            + (self.JJ**2 + 10 * self.JJ + 102) * eye(self.J2p1)
        )
        return 0.5 * (dot(O22, Y) + dot(Y, O22))

    def O6m2(self):
        O2m2 = self.O2m2()
        Y = (
            33 * mp(self.Jz, 4)
            - (18 * self.JJ + 123) * dot(self.Jz, self.Jz)
            + (self.JJ**2 + 10 * self.JJ + 102) * eye(self.J2p1)
        )
        return 0.5 * (dot(O2m2, Y) + dot(Y, O2m2))

    def O61(self):
        Y = (
            33 * mp(self.Jz, 5)
            - (30 * self.JJ - 15) * mp(self.Jz, 3)
            + (5 * self.JJ**2 - 10 * self.JJ + 12) * self.Jz
        )
        return 0.5 * (dot(self.Jx, Y) + dot(Y, self.Jx))

    def O6m1(self):
        Y = (
            33 * mp(self.Jz, 5)
            - (30 * self.JJ - 15) * mp(self.Jz, 3)
            + (5 * self.JJ**2 - 10 * self.JJ + 12) * self.Jz
        )
        return 0.5 * (dot(self.Jy, Y) + dot(Y, self.Jy))

    def O60(self):
        return (
            231 * mp(self.Jz, 6)
            - (315 * self.JJ - 735) * mp(self.Jz, 4)
            + (105 * self.JJ**2 - 525 * self.JJ + 294) * dot(self.Jz, self.Jz)
            - (5 * self.JJ**3 - 40 * self.JJ**2 + 60 * self.JJ) * eye(self.J2p1)
        )
