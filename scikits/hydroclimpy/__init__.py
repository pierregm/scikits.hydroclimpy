"""
Hydroclimpy

A collection of tools to manipulate hydro-climatologic data

:author: Pierre GF Gerard-Marchant 
:contact: pierregmcode_AT_gmail_DOT_com
"""

import scikits.timeseries as ts
from scikits.timeseries import *


from scikits.hydroclimpy.version import __version__

import core
from core import *
import enso
import io
import lib



__all__ = ['core', 'enso', 'io', 'lib', 'stats']
__all__.extend(core.__all__)
__all__.extend(ts.__all__)

from numpy.testing import Tester
test = Tester().test