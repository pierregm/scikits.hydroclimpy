"""
Test suite for :mod:`scikits.hydroclimpy.io`
"""

from numpy.ma.testutils import *

import scikits.timeseries as ts
from scikits.hydroclimpy.io.usgs import load_usgs_flows


class TestUSGS(TestCase):
    #
    def test_load_usgs_flows(self):
        "Test loadind USGS water data."
        #
        F02217900 = load_usgs_flows('02217900')
        data_192810 = [180, 180, 180, 180, 180, 180, 180, 180, 180, 180,
                       180, 180, 180, 180, 180, 180, 230, 300, 350, 300,
                       230, 200, 250, 270, 240, 200, 180, 180, 180, 180, 180]
        data_194912 = [276, 267, 267, 267, 262, 264, 269, 302, 291, 279,
                       310, 307, 304, 358, 425, 425, 371, 336, 336, 371,
                       336, 333, 336, 333, 317, 336, 410, 425, 384, 336, 330]
        assert(isinstance(F02217900, ts.TimeSeries))
        assert_equal(F02217900.start_date, ts.Date('D', '1928-10-01'))
        assert_equal(F02217900.end_date, ts.Date('D', '1949-12-31'))
        dates = F02217900.dates
        assert_equal(F02217900[(dates.year == 1928) & (dates.month == 10)],
                     data_192810)
        assert_equal(F02217900[(dates.year == 1949) & (dates.month == 12)],
                     data_194912)
        #
        F02217500 = load_usgs_flows('02217500')
        data_190110 = [ 570, 1000, 1200, 900, 650, 580, 540, 520, 510, 530,
                        499,  535,  535, 572, 535, 535, 499, 499, 499, 499,
                        499,  499,  499, 499, 499, 499, 464, 499, 499, 499, 499]
        assert(isinstance(F02217500, ts.TimeSeries))
        assert_equal(F02217500.start_date, ts.Date('D', '1901-10-01'))
        dates = F02217500.dates
        assert_equal(F02217500[(dates.year == 1901) & (dates.month == 10)],
                     data_190110)


################################################################################
if __name__ == '__main__':
    run_module_suite()

