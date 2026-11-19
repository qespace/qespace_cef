"""Stevens operators as separate functions: the first implementation.

Kept as a reference for ``StevensBase``. Cross-checked against the McPhase
manual, https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node132.html
"""

from numpy import dot, eye
from numpy.linalg import matrix_power as mp

from qespace_cef.cefmatrices import J_minus, J_plus, J_z


def O20(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return 3 * Jz * Jz - JJ * E


def O22(J):
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return 0.5 * (mp(Jplus, 2) + mp(Jminus, 2))


def O2m2(J):
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return -0.5j * (mp(Jplus, 2) - mp(Jminus, 2))


def O40(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return 35 * mp(Jz, 4) + (25 - 30 * JJ) * mp(Jz, 2) + E * JJ * (3 * JJ - 6)


def O42(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    E = eye(J2p1)

    # helper matrices:
    M_1 = 7 * mp(Jz, 2) - E * (JJ + 5)
    M_2 = mp(Jplus, 2) + mp(Jminus, 2)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O4m2(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)
    E = eye(J2p1)

    # helper matrices:
    M_1 = 7 * mp(Jz, 2) - E * (JJ + 5)
    M_2 = mp(Jplus, 2) - mp(Jminus, 2)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O43(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = Jz
    M_2 = mp(Jplus, 3) + mp(Jminus, 3)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O4m3(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = Jz
    M_2 = mp(Jplus, 3) - mp(Jminus, 3)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O44(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return 0.5 * (mp(Jplus, 4) + mp(Jminus, 4))


def O4m4(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return -0.5j * (mp(Jplus, 4) - mp(Jminus, 4))


def O60(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    E = eye(J2p1)

    return (
        231 * mp(Jz, 6)
        + mp(Jz, 4) * (735 - 315 * JJ)
        + mp(Jz, 2) * (105 * JJ**2 - 525 * JJ + 294)
        + E * (-5 * JJ**3 + 40 * JJ**2 - 60 * JJ)
    )


def O62(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 33 * mp(Jz, 4) - mp(Jz, 2) * (18 * JJ + 123) + E * (JJ**2 + 10 * JJ + 102)
    M_2 = mp(Jplus, 2) + mp(Jminus, 2)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O6m2(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 33 * mp(Jz, 4) - mp(Jz, 2) * (18 * JJ + 123) + E * (JJ**2 + 10 * JJ + 102)
    M_2 = mp(Jplus, 2) - mp(Jminus, 2)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O63(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 11 * mp(Jz, 3) - Jz * (59 + 3 * JJ)
    M_2 = mp(Jplus, 3) + mp(Jminus, 3)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O6m3(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 11 * mp(Jz, 3) - Jz * (59 + 3 * JJ)
    M_2 = mp(Jplus, 3) - mp(Jminus, 3)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O64(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 11 * mp(Jz, 2) - E * (JJ + 38)
    M_2 = mp(Jplus, 4) + mp(Jminus, 4)

    return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O6m4(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    E = eye(J2p1)
    Jz = J_z(J)
    Jplus = J_plus(J)
    Jminus = J_minus(J)

    # helper matrices:
    M_1 = 11 * mp(Jz, 2) - E * (JJ + 38)
    M_2 = mp(Jplus, 4) - mp(Jminus, 4)

    return -0.25j * (dot(M_1, M_2) + dot(M_2, M_1))


def O66(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return 0.5 * (mp(Jplus, 6) + mp(Jminus, 6))


def O6m6(J):
    JJ = J * (J + 1)
    J2p1 = int(2 * J + 1)

    Jplus = J_plus(J)
    Jminus = J_minus(J)

    return -0.5j * (mp(Jplus, 6) - mp(Jminus, 6))
