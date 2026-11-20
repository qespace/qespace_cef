"""Stevens operators as defined in CONV-0001, used as a reference.

Built from the formulas and the table of polynomials ``P_lm`` copied into
``decisions/conventions/CONV-0001-CEF-parameters.md`` from the McPhase manual
(M. Rotter, node132,
https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node132.html),
with ``X = J(J+1)`` and ``J± = J_x ± i J_y``:

    O_l0   = P_l0
    O_lm   =  1/4 [ (J+^m + J-^m) P_lm + P_lm (J+^m + J-^m) ]   (m > 0)
    O_l,-m = -i/4 [ (J+^m - J-^m) P_lm + P_lm (J+^m - J-^m) ]   (m > 0)

Only the even ranks l = 2, 4, 6 are tabulated there. The angular-momentum
matrices are built here, not taken from the package, in the package's basis
order ``|J>, |J-1>, ..., |-J>`` (the first component is m = +J).
"""

import numpy as np

# P_lm(z, X) as in CONV-0001, z the diagonal of J_z.
POLYNOMIALS = {
    (2, 0): lambda z, X: 3 * z**2 - X,
    (2, 1): lambda z, X: z,
    (2, 2): lambda z, X: np.ones_like(z),
    (4, 0): lambda z, X: 35 * z**4 - (30 * X - 25) * z**2 + 3 * X**2 - 6 * X,
    (4, 1): lambda z, X: 7 * z**3 - (3 * X + 1) * z,
    (4, 2): lambda z, X: 7 * z**2 - X - 5,
    (4, 3): lambda z, X: z,
    (4, 4): lambda z, X: np.ones_like(z),
    (6, 0): lambda z, X: (
        231 * z**6
        - (315 * X - 735) * z**4
        + (105 * X**2 - 525 * X + 294) * z**2
        - 5 * X**3
        + 40 * X**2
        - 60 * X
    ),
    (6, 1): lambda z, X: 33 * z**5 - (30 * X - 15) * z**3 + (5 * X**2 - 10 * X + 12) * z,
    (6, 2): lambda z, X: 33 * z**4 - (18 * X + 123) * z**2 + X**2 + 10 * X + 102,
    (6, 3): lambda z, X: 11 * z**3 - (3 * X + 59) * z,
    (6, 4): lambda z, X: 11 * z**2 - X - 38,
    (6, 5): lambda z, X: z,
    (6, 6): lambda z, X: np.ones_like(z),
}

# Every (l, m) of CONV-0001: l = 2, 4, 6 and m = -l ... l.
LM = [(l, m) for l in (2, 4, 6) for m in range(-l, l + 1)]


def _angular_momentum(J):
    """J_z diagonal (2J+1,) and J+ (2J+1, 2J+1) in the basis |J>, ..., |-J>."""
    z = J - np.arange(int(round(2 * J)) + 1)  # m = J, J-1, ..., -J
    # <m+1|J+|m> = sqrt(J(J+1) - m(m+1)); |m+1> sits one row above |m>.
    plus = np.diag(np.sqrt(J * (J + 1) - z[1:] * (z[1:] + 1)), k=1)
    return z, plus


def stevens(J, l, m):
    """The Stevens operator O_lm(J) of CONV-0001, a (2J+1, 2J+1) complex matrix.

    Parameters
    ----------
    J : float
        Total angular momentum quantum number, 1/2, 1, 3/2, ...
    l : int
        Rank, 2, 4 or 6.
    m : int
        Component, -l ... l; m < 0 is the sine-type operator.
    """
    z, plus = _angular_momentum(J)
    P = np.diag(POLYNOMIALS[l, abs(m)](z, J * (J + 1))).astype(complex)
    if m == 0:
        return P
    up = np.linalg.matrix_power(plus, abs(m))
    down = up.T  # J- = (J+)^T, real
    if m > 0:
        S = up + down
        return 0.25 * (S @ P + P @ S)
    A = up - down
    return -0.25j * (A @ P + P @ A)
