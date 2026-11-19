"""
Library for crystal field calculations.
"""

from importlib.metadata import version as _version

# Import main classes
from .cefion import CEFion
from .cefpars import CEFpars
from .ion import Ion

# Import the functions
from .observables import *

# The version is set in pyproject.toml only.
__version__ = _version("qespace_cef")
del _version

# Original copyright:

# Copyright 2014-2018 Petr Čermák, Jan Zubáč and Karel Pajskr
# This file is part of qespace_cef.
# qespace_cef is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# qespace_cef is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# <http://www.gnu.org/licenses/>.

# Further reworks by Michal Stekiel
