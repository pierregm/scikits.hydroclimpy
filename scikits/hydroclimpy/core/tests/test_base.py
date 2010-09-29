"""
Tests suite for climpy.core.base

:author: Pierre GERARD-MARCHANT
:contact: pierregm_AT_uga_DOT_edu

"""


import warnings
warnings.simplefilter("ignore", DeprecationWarning)

import numpy as np
import numpy.ma as ma
from numpy.ma import masked

from numpy.ma.testutils import *

import scikits.timeseries as ts

from scikits.hydroclimpy import ReferencedSeries


class TestReferencedSeries(TestCase):
    #
    def setUp(self):
        series = ReferencedSeries(np.random.rand(360),
                                  start_date=ts.Date('M','1970-01-01'))
        series[1:-1:10] = masked
        refperiod = ts.DateArray([ts.Date('M','1980-01-01'),
                                  ts.Date('M','1990-01-01')], freq='M')
        self.series = series
        self.refperiod = refperiod
    #
    def test_set_refperiod(self):
        "Test setting the reference period"
        series = self.series.copy()
        assert_equal(series.refperiod, series._dates[[0,-1]])
        series.refperiod = self.refperiod
        assert_equal(series.refperiod, self.refperiod.tovalue())
    #
    def test_keep_refperiod(self):
        "Tests keeping the reference period."
        (_series, _refperiod) = (self.series, self.refperiod)
        series = _series.copy()
        series.refperiod = _refperiod
        assert_equal(series.refperiod, _refperiod.tovalue())
        
        piece = series[120:]
        assert_equal(piece.refperiod, _refperiod.tovalue())
        piece.refperiod = _refperiod + 60
        assert_equal(series.refperiod, _refperiod.tovalue())
        assert_equal(piece._mask, _series._mask[120:])
    #
    def test_check_refperiod(self):
        "Test that the reference period is always a DateArray"
        series = ReferencedSeries(np.arange(10),
                                  start_date=ts.Date('M', '2001-01'))
        refperiod = series._optinfo['reference_period']
        self.failUnless(refperiod is None)
        series.refperiod = None
        refperiod = series._optinfo['reference_period']
        self.failUnless(isinstance(refperiod, ts.DateArray))
        assert_equal(refperiod.tovalues(), series.dates[[0, -1]].tovalues())
        series.refperiod = None
        assert_equal(refperiod, series.dates[[0, -1]])
        self.failUnless(isinstance(refperiod, ts.DateArray))
    #
    def test_period_statistics(self):
        "Tests the period statistics"
        (_series, _refperiod) = (self.series, self.refperiod)
        series = _series.copy()
        series.refperiod = _refperiod
        control = series[_refperiod[0]:_refperiod[-1]+1]
        #
        assert_equal(series.period_average(), control.mean())
        assert_equal(series.period_std(), control.std(ddof=1))
        #
        anomalies = series.period_anomalies()
        self.failUnless(isinstance(anomalies, ReferencedSeries))
        assert_equal(anomalies, series - control.mean())
        assert_equal(anomalies.mask, series.mask)
        assert_equal(anomalies.refperiod, series.refperiod)
    #
    def test_pickling(self):
        "Tests pickling/unpickling"
        import cPickle
        (_series, _refperiod) = (self.series, self.refperiod)
        series = _series.copy()
        series.refperiod = _refperiod
        pickled = cPickle.loads(cPickle.dumps(series))
        assert_equal(pickled, series)
        assert_equal(pickled.refperiod, _refperiod)

    def test_convert(self):
        "Test convert on ReferencedSeries"
        mseries = self.series.copy()
        mseries.refperiod = self.refperiod
        aseries = mseries.convert('A')
        assert_equal(aseries.refperiod, self.refperiod.asfreq('A'))
        dseries = mseries.convert('D')
        assert_equal(dseries.refperiod, self.refperiod.asfreq('D'))

    def test_get_refperiod(self):
        "Check that the reference period has always the proper frequency"
        a = ReferencedSeries(np.arange(2000, 2010),
                             start_date=ts.Date('Y','2000'),
                             refperiod=(ts.Date('Y','2003'),
                                        ts.Date('A','2006')))
        assert(isinstance(a.refperiod, ts.DateArray))
        assert_equal(a.refperiod.freq, a.dates.freq)
        a = a.convert('M')
        assert_equal(a.refperiod.freq, a.dates.freq)



#------------------------------------------------------------------------------
if __name__ == "__main__":
    run_module_suite()