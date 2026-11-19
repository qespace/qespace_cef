"""Tests of ``qespace_cef.observables``.

Only behaviour that AUDIT.md does not question is tested here; the known
problems (P1-P4, P9-P11) get their own tests when they are fixed.
"""

import qespace_cef as qc


def test_magnetization_returns_python_floats():
    """``magnetization`` returns a tuple of three Python floats (CLAUDE.md 5.7.3)."""
    cef = qc.CEFion(qc.Ion("Ce"), qc.CEFpars(dict(B20=1.0)))
    M = qc.magnetization(cef, 2.0, [0.0, 0.0, 1.0])
    assert len(M) == 3 and all(type(m) is float for m in M)
