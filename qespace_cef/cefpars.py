from . import constants as C
from .cefmatrices import StevensBase
from .data_containers import point_group_synonyms


class CEFpars:
    r"""
    ## Class representing set of crystal field parameters

    Values of crystal field parameters are stored in meV, but input in other units is possible.

    Symmetry handling was removed, as the value of that is dubious.

    It simplifies the creation of the CF parameter sets considering the point group of the ion.
    The symmetry restrictions follow from https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
    It also allows to look up the symmetry restrictions and nice printing of parameters.
    For the lookup table between notations on point groups see: https://en.wikipedia.org/wiki/Point_group

    Special care was taken for the cubic space group, as its symmetry imposes some algebraic
    relations between the non-zero elements.

    Attributes:
    -----------
        Bpars : dict[str, float]
            Bname:Bvalue pairs.
        symmetry: NamedTuple(Hermann-Mauguin, Schoenflies, lattice_names, id)
            Synonyms of point group symmetry.

    Examples:
    ---------
        >>> Bdict = dict(B20=1, B52=0.1, B6m1=-0.3)
        >>> cfp = qespace_cef.CEFpars(Bdict)
        >>> print(cfp)
            <CEFpars: symmetry=('2', 'C2', 'monoclinic', 1),
                B20 = 1.0 meV,
                B52 = 0.1 meV,
                B6m1 = -0.3 meV
            >
    """

    Bpars: dict[str, float]
    symmetry: str

    def __init__(self, Bdict: dict[str, float], units: str = "meV", pointGroup: str = "C1"):
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
            if Oname not in StevensBase.implementedOps():
                raise KeyError(f"Operator corresponding to {Bname} is not implemented.")

            self.Bpars[Bname] = Bvalue * unitConversions[units]

        # Symmetry
        self.symmetry = point_group_synonyms(pointGroup)[0]

    @classmethod
    def allowed_Bpars(cls, pointGroupName: str) -> list[str]:
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
        ret = f"<CEFpars: symmetry={self.symmetry}"
        for Bname, Bvalue in self.Bpars.items():
            ret += f",\n\t{Bname} = {Bvalue:.4} meV"

        return ret + "\n>"
