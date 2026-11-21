# Developing qespace_cef

Commands for working on the package. The rules for agents are in
[`CLAUDE.md`](CLAUDE.md) and the decisions in [`decisions/`](decisions/)
(ADRs and CONVs); the code is described in
[`ARCHITECTURE.md`](ARCHITECTURE.md); the open problems are in
[`AUDIT.md`](AUDIT.md), the next step in [`PLAN.md`](PLAN.md) and feature
ideas in [`ROADMAP.md`](ROADMAP.md). Everything goes through `uv`; `python`
is not expected on the PATH.

## Environment

```bash
uv sync                      # the package and the dev group into .venv
uv sync --group symbolic     # also sympy, for qespace_cef/symbolic
```

Dependencies: `numpy>=1.26` at run time; the groups `dev` (pytest, ruff, nox,
sphinx, furo) and `symbolic` (sympy) are for development only. After
changing them, run `uv lock`; `uv.lock` is tracked. Ask the maintainer before
adding a dependency (CLAUDE.md 4.2).

To try something without touching `.venv`:

```bash
uv run --no-project --isolated --with-editable . --with pytest pytest
```

## Checks

```bash
uv run pytest                     # tests; perf stubs and known problems excluded
uv run pytest -m perf -s          # the performance stubs, with timings
uv run pytest --doctest-modules qespace_cef --ignore=qespace_cef/symbolic
uv run ruff check .
uv run ruff format --check .      # `uv run ruff format .` to apply
uv run nox                        # before a release: numpy x Python matrix
```

`uv run pytest` reports the known problems from `AUDIT.md` as `xfailed`.
They are strict: a fixed problem shows as a failure (XPASS) until its
marker is removed.

### numpy and Python matrix

`noxfile.py` runs the tests and the doctests in fresh uv environments:

| Python | numpy |
|---|---|
| 3.12 | 1.26.0 (the floor in `pyproject.toml`), 2.0.2 (the first numpy 2), newest |
| 3.13 | newest |
| 3.14 | newest |

```bash
uv run nox -l                                  # list the sessions
uv run nox -s "tests(3.12-numpy-1.26.0)"      # one session
```

The environments go to `.nox/`. If numpy fails to load its DLLs because the
path is too long, point them to a short folder:
`QESPACE_NOX_ENVDIR=C:/tmp/nox uv run nox`. The interpreters come from
`uv python list`. When the floor changes, change the matrix with it.

## Documentation

```bash
uv run python docs/build.py       # Sphinx, into docs/_build/html
```

## Releases

The version is set in `pyproject.toml` only, and `__version__` reads it from
the package metadata. As done for 1.3.0:

1. On a branch: `uv version --bump minor` (or `patch`, `major`), and the
   `Unreleased` section of `CHANGELOG.md` renamed to the new version.
2. After the maintainer's review, merged into `main`.
3. The maintainer tags `main` as `vX.Y.Z` and pushes `main` and the tag.
