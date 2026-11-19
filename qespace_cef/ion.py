import numpy as np
from numpy.typing import ArrayLike


class Ion:
    r"""
    Free trivalent rare-earth ion: ground J multiplet and form-factor data.

    Parameters
    ----------
    ion_name : str
        Element symbol, case-insensitive: Ce, Pr, Nd, Sm, Tb, Dy, Ho, Er, Tm
        or Yb.

    Attributes
    ----------
    name : str
        Element symbol in lower case.
    J : float
        Total angular momentum of the ground multiplet (Hund's rules).
    J2p1 : int
        Degeneracy 2J+1 of the multiplet.
    gJ : float
        Landé factor, dimensionless.
    Alpha, Beta, Gamma : float
        Stevens factors :math:`\alpha_J, \beta_J, \gamma_J` of the rank 2, 4
        and 6 operators, dimensionless. Not used elsewhere in the package.
    mffCoefficients : dict
        ``{'j0': [A, a, B, b, C, c, D], 'j2': [A, a, B, b, C, c, D]}``, the
        coefficients of :math:`\langle j_0 \rangle` and
        :math:`\langle j_2 \rangle` used by :meth:`mff`.

    Raises
    ------
    KeyError
        If the ion is not in the table.

    Warnings
    --------
    Known problems, not fixed yet (AUDIT.md):

    - P13: the form-factor coefficients for Ce are those of Ce2+, the only
      Ce entry in the ILL table; the other ions use 3+ values.
    - P14: the values in the table (J, gJ, Stevens factors, form-factor
      coefficients) have no source next to them.
    - P8: ``str()`` prints J with ``%d``, so J = 5/2 shows as 2.

    Examples
    --------
    >>> import qespace_cef as qc
    >>> ce = qc.Ion("Ce")
    >>> ce.J, round(ce.gJ, 4)
    (2.5, 0.8571)
    """

    def __init__(self, ion_name: str):
        self.name = ion_name.lower()

        # AUDIT P14: no source for these values.
        atomicDatabase = {
            # AUDIT P13: the form-factor coefficients are Ce2+, not Ce3+.
            "ce": {
                "J": 2.5,
                "gJ": 0.8571428571428571,
                "alpha": -5.7142857143e-02,
                "beta": 6.3492063492e-03,
                "gamma": 0.0000000000e00,
                "js": {
                    "j0": [0.2953, 17.6846, 0.2923, 6.7329, 0.4313, 5.3827, -0.0194],
                    "j2": [0.9809, 18.063, 1.8413, 7.7688, 0.9905, 2.8452, 0.012],
                },
            },
            "pr": {
                "J": 4.0,
                "gJ": 0.8,
                "alpha": -2.1010101010e-02,
                "beta": -7.3461891644e-04,
                "gamma": 6.0994000388e-05,
                "js": {
                    "j0": [0.0504, 24.9989, 0.2572, 12.0377, 0.7142, 5.0039, -0.0219],
                    "j2": [0.8734, 18.9876, 1.5594, 6.0872, 0.8142, 2.415, 0.0111],
                },
            },
            "nd": {
                "J": 4.5,
                "gJ": 0.7272727272727273,
                "alpha": 6.4279155188e-03,
                "beta": -2.9110772912e-04,
                "gamma": -3.7987959158e-05,
                "js": {
                    "j0": [0.054, 25.0293, 0.3101, 12.102, 0.6575, 4.7223, -0.0216],
                    "j2": [0.6751, 18.3421, 1.6272, 7.26, 0.9644, 2.6016, 0.015],
                },
            },
            "sm": {
                "J": 2.5,
                "gJ": 0.2857142857142857,
                "alpha": 4.1269841270e-02,
                "beta": 2.5012025012e-03,
                "gamma": 0.0000000000e00,
                "js": {
                    "j0": [0.0288, 25.2068, 0.2973, 11.8311, 0.6954, 4.2117, -0.0213],
                    "j2": [0.4707, 18.4301, 1.4261, 7.0336, 0.9574, 2.4387, 0.0182],
                },
            },
            "tb": {
                "J": 6.0,
                "gJ": 1.5,
                "alpha": -1.0101010101e-02,
                "beta": 1.2243648607e-04,
                "gamma": 1.1212132424e-06,
                "js": {
                    "j0": [0.0177, 25.5095, 0.2921, 10.5769, 0.7133, 3.5122, -0.0231],
                    "j2": [0.2892, 18.4973, 1.1678, 6.7972, 0.9437, 2.2573, 0.0232],
                },
            },
            "dy": {
                "J": 7.5,
                "gJ": 1.3333333333333333,
                "alpha": -6.3492063492e-03,
                "beta": -5.9200059200e-05,
                "gamma": 1.0349660699e-06,
                "js": {
                    "j0": [0.1157, 15.0732, 0.327, 6.7991, 0.5821, 3.0202, -0.0249],
                    "j2": [0.2523, 18.5172, 1.0914, 6.7362, 0.9345, 2.2082, 0.025],
                },
            },
            "ho": {
                "J": 8.0,
                "gJ": 1.25,
                "alpha": -2.2222222222e-03,
                "beta": -3.3300033300e-05,
                "gamma": -9.9515968263e-08,
                "js": {
                    "j0": [0.0566, 18.3176, 0.3365, 7.688, 0.6317, 2.9427, -0.0248],
                    "j2": [0.2188, 18.5157, 1.024, 6.707, 0.9251, 2.1614, 0.0268],
                },
            },
            "er": {
                "J": 7.5,
                "gJ": 1.2,
                "alpha": 2.5396825397e-03,
                "beta": 4.4400044400e-05,
                "gamma": 2.6909117818e-05,
                "js": {
                    "j0": [0.0586, 17.9802, 0.354, 7.0964, 0.6126, 2.7482, -0.0251],
                    "j2": [0.171, 18.5337, 0.9879, 6.6246, 0.9044, 2.1004, 0.0278],
                },
            },
            "tm": {
                "J": 6.0,
                "gJ": 1.1666666666666667,
                "alpha": 1.0101010101e-02,
                "beta": 1.6324864810e-04,
                "gamma": -5.6060662121e-06,
                "js": {
                    "j0": [0.0581, 15.0922, 0.2787, 7.8015, 0.6854, 2.7931, -0.0224],
                    "j2": [0.176, 18.5417, 0.9105, 6.5787, 0.897, 2.0622, 0.0294],
                },
            },
            "yb": {
                "J": 3.5,
                "gJ": 1.1428571428571428,
                "alpha": 3.1746031746e-02,
                "beta": -1.7316017316e-03,
                "gamma": 1.4800014800e-04,
                "js": {
                    "j0": [0.0416, 16.0949, 0.2849, 7.8341, 0.6961, 2.6725, -0.0229],
                    "j2": [0.157, 18.5553, 0.8484, 6.5403, 0.888, 2.0367, 0.0318],
                },
            },
        }

        if self.name not in atomicDatabase.keys():
            raise KeyError(f"{ion_name} - ion with that name not implemented")

        elementData = atomicDatabase[self.name]

        self.J: float = elementData["J"]
        self.J2p1: float = int(2 * self.J + 1)
        self.gJ: float = elementData["gJ"]
        self.Alpha: float = elementData["alpha"]
        self.Beta: float = elementData["beta"]
        self.Gamma: float = elementData["gamma"]
        self.mffCoefficients = elementData["js"]

    def __str__(self):
        """Name, J and gJ, e.g. ``'Ho: J = 8, gJ = 1.25'``.

        Warnings
        --------
        P8 (AUDIT.md), not fixed yet: J is printed with ``%d``, so J = 5/2
        shows as 2.
        """
        # AUDIT P8: %d truncates a half-integer J.
        return "%s: J = %d, gJ = %.2f" % (self.name.title(), self.J, self.gJ)

    def mff(self, Q: ArrayLike) -> float | np.ndarray:
        r"""Magnetic form factor in the dipole approximation.

        .. math::

            f(Q) = \langle j_0 \rangle + \frac{2 - g_J}{g_J} \langle j_2 \rangle,

        with :math:`s = Q / 4\pi` and the coefficients in ``mffCoefficients``:

        .. math::

            \langle j_0 \rangle = A e^{-a s^2} + B e^{-b s^2} + C e^{-c s^2} + D, \qquad
            \langle j_2 \rangle = s^2 \left( A e^{-a s^2} + B e^{-b s^2} + C e^{-c s^2} + D \right).

        The dipole approximation holds at small Q.

        Parameters
        ----------
        Q : float or array_like of float
            Length of the momentum transfer in 1/Å (not a vector).

        Returns
        -------
        float or ndarray
            Form factor, dimensionless, with the shape of ``Q``;
            :math:`f(0) \approx 1`.

        Warnings
        --------
        P13 (AUDIT.md), not fixed yet: for Ce the coefficients are those of
        Ce2+.
        """
        s = np.divide(Q, 4 * np.pi)
        coefs0 = self.mffCoefficients["j0"]
        coefs2 = self.mffCoefficients["j2"]

        j0 = (
            coefs0[0] * np.exp(-coefs0[1] * s**2)
            + coefs0[2] * np.exp(-coefs0[3] * s**2)
            + coefs0[4] * np.exp(-coefs0[5] * s**2)
            + coefs0[6]
        )

        j2 = s**2 * (
            coefs2[0] * np.exp(-coefs2[1] * s**2)
            + coefs2[2] * np.exp(-coefs2[3] * s**2)
            + coefs2[4] * np.exp(-coefs2[5] * s**2)
            + coefs2[6]
        )

        return j0 + j2 * (2 - self.gJ) / self.gJ
