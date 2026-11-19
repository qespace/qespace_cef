"""Build the HTML documentation into ``docs/_build/html``.

Usage, from the repository root::

    uv run python docs/build.py
"""

import sys
from pathlib import Path

from sphinx.cmd.build import main

DOCS = Path(__file__).resolve().parent

if __name__ == '__main__':
    sys.exit(main(['-M', 'html', str(DOCS), str(DOCS / '_build')]))
