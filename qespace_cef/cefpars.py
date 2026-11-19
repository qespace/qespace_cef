import numpy as np

from . import constants as C
from .cefmatrices import StevensBase
from .data_containers import point_group_synonyms


class CEFpars:
    r"""
    Set of crystal-field parameters :math:`B_k^q` for Stevens operators.

    The names follow the operators: ``'B20'`` multiplies
    :math:`\hat{O}_2^0`, ``'B4m3'`` multiplies :math:`\hat{O}_4^{-3}`, for
    every operator that :class:`~qespace_cef.cefmatrices.StevensBase`
    implements (k = 0 ... 6). The values are stored in meV.

    The point group is stored, but does not restrict the parameters; to
    look up which parameters a point group allows, use
    :meth:`allowed_Bpars`.

    Parameters
    ----------
    Bdict : dict of str to float
        Parameter names and values, in ``units``.
    units : {'meV', 'K', 'invcm'}, optional
        Unit of the values in ``Bdict``. Default ``'meV'``.
    pointGroup : str or int, optional
        Any point-group name that
        :func:`~qespace_cef.data_containers.point_group_synonyms` accepts;
        for a name with several matches the first is kept. Default ``'C1'``.

    Attributes
    ----------
    Bpars : dict of str to float
        Parameter names and values in meV.
    symmetry : numpy.void
        Record ``(HM, Sch, lattice, id)`` of the point group: the
        Hermann-Mauguin and Schoenflies symbols, the lattice name and an
        integer id.

    Raises
    ------
    ValueError
        If ``units`` is unknown.
    KeyError
        If a name in ``Bdict`` has no implemented operator.
    NameError
        If ``pointGroup`` is unknown.

    Warnings
    --------
    P16 (AUDIT.md), open question: parameters of odd rank (k = 1, 3, 5)
    are accepted, though parity removes odd k from the CEF Hamiltonian
    within one J multiplet. Whether to accept them is open (AUDIT C3).

    Examples
    --------
    >>> import qespace_cef as qc
    >>> cfp = qc.CEFpars(dict(B20=1, B40=0.1, B6m1=-0.3))
    >>> print(cfp)  # doctest: +NORMALIZE_WHITESPACE
    <CEFpars: symmetry=('1', 'C1', 'triclinic', 0),
        B20 = 1.0 meV,
        B40 = 0.1 meV,
        B6m1 = -0.3 meV
    >
    """

    Bpars: dict[str, float]
    symmetry: np.void

    def __init__(self, Bdict: dict[str, float], units: str = "meV", pointGroup: str | int = "C1"):
        # Check units and assign conversion factors
        unitConversions = {"meV": 1.0, "K": C.K2meV, "invcm": C.invcm2meV}
        if units not in unitConversions:
            raise ValueError(
                f'The desired unit "{units}" not in the list of implemented ones: {unitConversions}'
            )

        # Assign parameters to internal dictionary
        self.Bpars = dict()
        for Bname, Bvalue in Bdict.items():
            Oname = Bname.replace("B", "O")
            # AUDIT P16: accepts odd k (1, 3, 5), which parity forbids in one J multiplet.
            if Oname not in StevensBase.implementedOps():
                raise KeyError(f"Operator corresponding to {Bname} is not implemented.")

            self.Bpars[Bname] = float(Bvalue) * unitConversions[units]

        # Symmetry
        self.symmetry = point_group_synonyms(pointGroup)[0]

    @classmethod
    def allowed_Bpars(cls, pointGroupName: str | int) -> list[str]:
        """Names of the parameters allowed (non-zero) in a point group.

        Taken from the table in the McPhase manual,
        https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html,
        which assumes z along the main symmetry axis. Relations between
        parameters (e.g. cubic ``B44 = 5 B40``) are not included.

        Parameters
        ----------
        pointGroupName : str or int
            Any point-group name that
            :func:`~qespace_cef.data_containers.point_group_synonyms`
            accepts; for a name with several matches the first is used.

        Returns
        -------
        list of str
            Parameter names such as ``'B20'``. For C1, all 49 implemented
            operators, including k = 0, 1, 3, 5.

        Raises
        ------
        NameError
            If the point group is unknown.

        Examples
        --------
        >>> import qespace_cef as qc
        >>> qc.CEFpars.allowed_Bpars("4mm")
        ['B20', 'B40', 'B44', 'B60', 'B64']
        """
        symmetry_Sch = point_group_synonyms(pointGroupName)[0]["Sch"]

        # List of allowed parameters is compiled based on McPhase manual
        # https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html

        r1 = ["B" + Oname[1:] for Oname in StevensBase.implementedOps()]
        r2 = [
            "B20",
            "B22",
            "B2m2",
            "B40",
            "B42",
            "B4m2",
            "B43",
            "B4m3",
            "B44",
            "B4m4",
            "B60",
            "B62",
            "B6m2",
            "B63",
            "B6m3",
            "B64",
            "B6m4",
            "B66",
            "B6m6",
        ]
        r3 = ["B20", "B22", "B40", "B42", "B44", "B60", "B62", "B64", "B66"]
        r4 = ["B20", "B40", "B44", "B4m4", "B60", "B64", "B6m4"]
        r5 = ["B20", "B40", "B44", "B60", "B64"]
        r6 = ["B20", "B40", "B43", "B4m3", "B60", "B63", "B6m3", "B66", "B6m6"]
        r7 = ["B20", "B40", "B43", "B60", "B63", "B66"]
        r8 = ["B20", "B40", "B60", "B66", "B6m6"]
        r9 = ["B20", "B40", "B60", "B66"]
        r10 = ["B40", "B44", "B60", "B62", "B64", "B66"]
        r11 = ["B40", "B44", "B60", "B64"]

        allowed_Bpars = {
            "C1": r1,
            "C2": r2,
            "Cs": r2,
            "C2h": r2,
            "C2v": r3,
            "D2": r3,
            "D2h": r3,
            "C4": r4,
            "S4": r4,
            "C4h": r4,
            "D4": r5,
            "C4v": r5,
            "D2d": r5,
            "D4h": r5,
            "C3": r6,
            "S6": r6,
            "D3": r7,
            "C3v": r7,
            "D3d": r7,
            "C6": r8,
            "C3h": r8,
            "C6h": r8,
            "D6": r9,
            "C6v": r9,
            "D3h": r9,
            "D6h": r9,
            "T": r10,
            "Th": r10,
            "Td": r11,
            "O": r11,
            "Oh": r11,
        }

        return allowed_Bpars[symmetry_Sch]

    def __str__(self):
        """Point group and parameters in meV, one per line."""
        ret = f"<CEFpars: symmetry={self.symmetry}"
        for Bname, Bvalue in self.Bpars.items():
            ret += f",\n\t{Bname} = {Bvalue:.4} meV"

        return ret + "\n>"
