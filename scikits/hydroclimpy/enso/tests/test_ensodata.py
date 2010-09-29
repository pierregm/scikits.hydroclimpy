"""
Suite test for :mod:`scikits.hydroclimpy.enso.ensodata`
"""

import numpy as np
import numpy.ma as ma

from numpy.ma.testutils import *

import scikits.timeseries as ts
from scikits.timeseries import time_series

from scikits.hydroclimpy.enso import climate_series, load_oni, load_jma


class TestJMA(TestCase):
    #
    def setUp(self):
        "Setting common information"
        self.JMA_std = load_jma('standard')
        self.JMA_coaps = load_jma('coaps')
        self.JMA_slide = load_jma('sliding')
    #
    def test_jma_refperiod(self):
        "Check that we have the proper frequency for the reference period"
        refperiod = self.JMA_std.reference_period
        assert_equal(refperiod.freq, self.JMA_std._dates.freq)
        refperiod = self.JMA_coaps.reference_period
        assert_equal(refperiod.freq, self.JMA_coaps._dates.freq)
        refperiod = self.JMA_slide.reference_period
        assert_equal(refperiod.freq, self.JMA_slide._dates.freq)


class TestONI(TestCase):
    #
    def setUp(self):
        "Setting common information"
        try:
            from BeautifulSoup import BeautifulSoup, SoupStrainer
        except ImportError:
            self.indices = None
            return
        # Load the file as a tree, but only take the SST table (border=1)
        from urllib import urlopen
        url = "http://www.cpc.noaa.gov/products/analysis_monitoring/"\
              "ensostuff/ensoyears.shtml"
        url = urlopen(url)
        table = BeautifulSoup(url.read(),
                              parseOnlyThese=SoupStrainer("table", border=1))
        # Separate it by rows, but skip the first one (the header)
        years = []
        indices = []
        color = dict(red=+1, white=0, blue=-1)
        deft = [(None,'color:white')]
        for row in table.findAll("tr")[1:]:
            cols = row.findAll('td')
            years.append(int(cols.pop(0).strong.string))
            indices.append([color[getattr(_.span, 'attrs', deft)[0][-1].split(':')[-1]]
                            for _ in cols])
        start_date = ts.Date('M', year=years[0], month=1)
        self.indices = time_series(np.array(indices).ravel(),
                                   start_date=start_date)
    #
    def test_ONI_standard(self):
        "Test the 'standard' option of load_oni"
        if self.indices is None:
            return
        ONI = load_oni('standard')
        ONI.set_indices(full_year=False, minimum_size=5, reference_season=None)
        assert_equal(*ts.align_series(ONI.indices, self.indices))
    #
    def test_set_ensoindices(self):
        "Test setting indices"
        ONI = load_oni('standard')
        mseries = climate_series(np.random.rand(120),
                                 start_date=ts.Date('M','1991-01'),
                                 ensoindicator=ONI)
        control = ONI['1991-01':'2001-01']
        assert_equal(mseries.ensoindices, control.indices)
        test = mseries.set_ensoindices(full_year=True).copy()
        self.failUnless(mseries.ensoindicator.optinfo['full_year'])
        assert_equal(test,
                     control.set_indices(full_year=True))
        self.failUnless(control.optinfo['full_year'])
        assert_equal(test, mseries.ensoindices)
        assert_equal(mseries.set_ensoindices(lag=6), test.tshift(-6))
        #
        mseries.set_ensoindices(full_year=True, lag=0)
        series = ts.lib.backward_fill(mseries.convert('D'))
        series_m = series.convert('M', ma.mean)
        assert_equal(test, series_m.ensoindices)
    #
    def test_oni_refperiod(self):
        "Check that we have the proper frequency for the reference period"
        ONI = load_oni('standard')
        refperiod = ONI.refperiod
        assert(isinstance(refperiod, ts.DateArray))
        assert_equal(refperiod.freqstr, 'M')


################################################################################
if __name__ == '__main__':
    run_module_suite()


