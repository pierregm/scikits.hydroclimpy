"""
Generic collection of objects and functions used by several submodules.
"""


from .. import _config

# Import modules from scikits.hydroclimpy.core
import base
from base import *
import tools
from tools import *
import ts_addons
from ts_addons import *
import stats_addons

__all__ = []
__all__.extend(tools.__all__)
__all__.extend(ts_addons.__all__)
__all__.extend(base.__all__)

from numpy.testing import Tester
test = Tester().test