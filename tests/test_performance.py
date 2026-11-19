"""Performance stubs: time the main computations on the largest rare-earth J.

Skipped by default; run them with ``uv run pytest -m perf -s``. They print the
time and assert only that the result is sensible: no time limits yet. Real
benchmarks will come from real applications.
"""

import time

import numpy as np
import pytest

import qespace_cef as qc
from qespace_cef.cefmatrices import StevensBase

pytestmark = pytest.mark.perf

J_MAX = 8.0  # Ho3+
HO_PARS = dict(B20=0.1, B22=0.05, B40=1e-3, B43=2e-3, B44=1e-3, B60=1e-5, B63=2e-5, B66=1e-5)


def _timed(label, func, repeat):
    """Call ``func`` ``repeat`` times, print the mean time and return the last result."""
    start = time.perf_counter()
    for _ in range(repeat):
        result = func()
    mean = (time.perf_counter() - start) / repeat
    print(f"\n{label}: {mean * 1e3:.3f} ms per call ({repeat} calls)")
    return result


def _ho():
    return qc.CEFion(qc.Ion("Ho"), qc.CEFpars(HO_PARS))


def test_stevens_operators():
    """Build every Stevens operator for J = 8."""
    names = [n[1:] for n in StevensBase.implementedOps()]

    def build():
        ops = StevensBase(J_MAX)
        return [ops[n] for n in names]

    matrices = _timed("StevensBase(8), all operators", build, repeat=20)
    assert len(matrices) == 49 and matrices[-1].shape == (17, 17)


def test_diagonalisation():
    """Build and diagonalise the Ho3+ Hamiltonian."""
    cef = _timed("CEFion(Ho) incl. diagonalisation", _ho, repeat=50)
    assert len(cef.energies) == 17


def test_susceptibility():
    """Susceptibility on 1000 temperatures."""
    cef = _ho()
    T = np.linspace(1, 300, 1000)
    chi = _timed(
        "susceptibility(magnetization), 1000 T",
        lambda: qc.susceptibility(cef, T, [0, 0, 1], method="magnetization"),
        repeat=1,
    )
    assert chi.shape == T.shape and np.all(np.isfinite(chi))


def test_spectral_weight():
    """Powder spectral weight at 1000 values of |Q|."""
    cef = _ho()
    Qs = np.linspace(0.1, 6, 1000)

    def run():
        return [qc.spectral_weight(cef, 10.0, [Q, 0, 0], "powder") for Q in Qs]

    out = _timed("spectral_weight(powder), 1000 Q", run, repeat=1)
    assert len(out) == 1000 and out[0][1].shape == (17, 17)
