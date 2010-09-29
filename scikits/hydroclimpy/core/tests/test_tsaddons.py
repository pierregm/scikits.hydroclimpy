"""
Tests suite for climpy.core.base

:author: Pierre GERARD-MARCHANT
:contact: pierregm_AT_uga_DOT_edu

"""


import warnings
warnings.simplefilter("ignore", DeprecationWarning)

import numpy as np
import numpy.ma as ma
import scikits.timeseries as ts

from numpy.testing import *
from numpy.ma.testutils import *


from scikits.hydroclimpy import force_reference, apply_on_fields




class TestGenericFunctions(TestCase):
    #
    def test_force_reference(self):
        mseries = ts.time_series(np.arange(24),
                                 start_date=ts.Date('M','2001-01'))
        aseries = ts.time_series([1,2,3], start_date=ts.Date('A', '2001-01'))
        #
        mtest = force_reference(aseries, mseries)
        assert_equal(mtest.freq, ts.check_freq('M'))
        assert_equal(mtest.dates[[0,-1]], mseries.dates[[0,-1]])
        assert_equal(mtest, [1]*12+[2]*12)
        mtest = force_reference(aseries, mseries, ma.sum)
        assert_equal(mtest, [1]*12+[2]*12)
        #
        atest = force_reference(mseries, aseries)
        assert_equal(atest.freq, ts.check_freq('A'))
        assert_equal(atest.dates[[0,-1]], aseries.dates[[0,-1]])
        assert_equal(atest, ma.array([5.5, 17.5, 0], mask=[0,0,1]))
        atest = force_reference(mseries, aseries, ma.sum)
        assert_equal(atest, ma.array([66, 210, 0], mask=[0,0,1]))

    def test_apply_on_fields_scalar(self):
        "Test apply_on_fields : scalar output"
        adtype = [('fi', int), ('ff', float)]
        a = np.array([(1, 1.), (2, 2.), (3, 3.)], dtype=adtype)
        func = ma.sum
        test = apply_on_fields(a, func)
        control = np.array([(6, 6.)], dtype=adtype)
        assert_equal(test, control)

    def test_apply_on_fields_scalar_wmask(self):
        "Test apply_on_fields scalar output w/ mask"
        adtype = [('fi', int), ('ff', float)]
        a = ma.array([(1, 1.), (2, 2.), (3, 3.)],
                     mask=[(0, 1), (0, 1), (0, 1)],
                     dtype=adtype)
        func = ma.sum
        test = apply_on_fields(a, func)
        control = ma.array((6, 6.), mask=(0, 1), dtype=adtype)
        assert_equal(test, control)

    def test_apply_on_fields_array(self):
        "Test apply_on_fields: array output"
        adtype = [('fi', int), ('ff', float)]
        a = ma.array([(1, 1.), (2, 2.), (3, 3.)],
                     mask=[(0, 0), (0, 1), (0, 0)],
                     dtype=adtype)
        func = ma.cumsum
        test = apply_on_fields(a, func)
        control = ma.array([(1, 1.), (3, -1), (6, 4.)],
                           mask=[(0, 0), (0, 1), (0, 0)],
                           dtype=adtype)
        assert_equal(test, control)

    def test_apply_on_fields_series(self):
        "Test apply_on_fields w/ time_series"
        adtype = [('fi', int), ('ff', float)]
        a = ts.time_series([(1, 1.), (2, 2.), (3, 3.)],
                           mask=[(0, 0), (0, 1), (0, 0)],
                           dtype=adtype,
                           start_date=ts.now('M'))
        func = ma.cumsum
        test = apply_on_fields(a, func)
        control = ts.time_series([(1, 1.), (3, -1), (6, 4.)],
                                 mask=[(0, 0), (0, 1), (0, 0)],
                                 dtype=adtype,
                                 start_date=ts.now('M'))
        assert_equal(test, control)
        self.failUnless(isinstance(test, ts.TimeSeries))
        assert_equal(test.dates, control.dates)



#------------------------------------------------------------------------------
if __name__ == "__main__":
    run_module_suite()