"""Tests of ``qespace_cef.Ion`` against Hund's rules and the form-factor limit."""

import numpy as np
import pytest

import qespace_cef as qc

# Hund's-rule ground terms (S, L, J) of the trivalent rare-earth ions,
# e.g. J. Jensen and A. R. Mackintosh, Rare Earth Magnetism (Clarendon, 1991),
# table 1.1.
HUND = {
    "Ce": (1 / 2, 3, 5 / 2),
    "Pr": (1, 5, 4),
    "Nd": (3 / 2, 6, 9 / 2),
    "Sm": (5 / 2, 5, 5 / 2),
    "Tb": (3, 3, 6),
    "Dy": (5 / 2, 5, 15 / 2),
    "Ho": (2, 6, 8),
    "Er": (3 / 2, 6, 15 / 2),
    "Tm": (1, 5, 6),
    "Yb": (1 / 2, 3, 7 / 2),
}


@pytest.mark.parametrize("name", HUND)
def test_J_and_lande_factor(name):
    """J from Hund's rules; gJ = 3/2 + [S(S+1) - L(L+1)] / [2 J(J+1)]."""
    S, L, J = HUND[name]
    ion = qc.Ion(name)
    assert ion.J == J
    assert ion.gJ == pytest.approx(1.5 + (S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1)))


@pytest.mark.parametrize("name", HUND)
def test_form_factor_at_zero_q(name):
    """f(Q = 0) = 1 within the 1 % accuracy of the tabulated <j0> fit."""
    assert qc.Ion(name).mff(0.0) == pytest.approx(1.0, abs=0.01)


# Sm3+ is left out: with gJ = 2/7 the <j2> term dominates and its dipole form
# factor is known to rise again at larger Q.
@pytest.mark.parametrize("name", [n for n in HUND if n != "Sm"])
def test_form_factor_decreases(name):
    """The form factor falls from Q = 0 to Q = 6 per angstrom."""
    f = qc.Ion(name).mff(np.linspace(0, 6, 61))
    assert np.all(np.diff(f) < 0)


def test_name_is_case_insensitive():
    assert qc.Ion("ce").J == qc.Ion("CE").J == qc.Ion("Ce").J


def test_unknown_ion_raises():
    with pytest.raises(KeyError):
        qc.Ion("Eu")
