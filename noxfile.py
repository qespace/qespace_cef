"""Test matrix: numpy versions against the supported Python versions.

Run every session with ``uv run nox``, one with ``uv run nox -s "tests(3.12-numpy-latest)"``,
and list them with ``uv run nox -l``. The environments are made with uv; their
folder is ``.nox`` unless the environment variable ``QESPACE_NOX_ENVDIR``
names a shorter one (long paths can break numpy's DLL loading on Windows).
"""

import os

import nox

nox.options.default_venv_backend = "uv"
nox.options.envdir = os.environ.get("QESPACE_NOX_ENVDIR", ".nox")

# (Python, numpy) pairs that have wheels. "latest" lets the resolver pick the
# newest numpy for that Python.
MATRIX = [
    ("3.12", "latest"),
    ("3.13", "latest"),
    ("3.14", "latest"),
]


@nox.session
@nox.parametrize(("python", "numpy"), MATRIX, ids=[f"{p}-numpy-{n}" for p, n in MATRIX])
def tests(session, numpy):
    """The test suite and the docstring examples with one numpy version."""
    session.install("-e", ".", "pytest", "numpy" if numpy == "latest" else f"numpy=={numpy}")
    session.run(
        "python",
        "-c",
        "import numpy, sys; print('numpy', numpy.__version__, 'Python', sys.version.split()[0])",
    )
    session.run("pytest", "-q", "-p", "no:cacheprovider")
    session.run(
        "pytest",
        "-q",
        "-p",
        "no:cacheprovider",
        "--doctest-modules",
        "qespace_cef",
        "--ignore=qespace_cef/symbolic",
    )
