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

from scipy.stats.mstats import mquantiles, mode

import scikits.hydroclimpy as hc
from scikits.hydroclimpy.core import Cluster, ReferencedSeries,\
                                     force_reference
from scikits.hydroclimpy.enso import ENSOIndicator, ClimateSeries,\
                                     ClimateRecords,\
                                     apply_on_phase, climate_series




class TestENSOIndicator(TestCase):
    #
    def setUp(self):
        "Initializes"
        data = [
        -0.49,-0.38,-0.26, 0.54, 0.95, 0.60, 0.93, 0.80, 0.52, 0.15,-0.50,-0.58,
        -0.58,-0.70,-0.54,-0.80,-0.75, 0.33,-0.18,-0.55,-0.77,-0.73,-0.58,-0.68,
        -0.62, 0.75, 0.72, 0.50, 0.84, 0.64, 0.77,-0.03, 0.71, 0.17, 0.64,-0.12,
        -0.17, 0.34, 0.57,-0.57,-0.77,-0.56,-0.51,-0.80, 0.19,-0.83, 0.53, 0.72,
         0.77, 0.64, 0.57, 0.51, 0.51, 0.67, 0.70, 0.40, 0.29, 0.23, 0.13, 0.00]
        refperiod = ts.DateArray([ts.Date('M','1980-01-01'),
                                  ts.Date('M','1990-01-01')], freq='M')
        series = ENSOIndicator(data,
                               start_date=ts.Date('M','1980-01-01'),
                               refperiod=refperiod)
        self.series = series
        self.refperiod = refperiod
    #
    def test_thresholds(self):
        "Tests thresholds"
        series = self.series
        # Setting a threshold w/ a single value should fail
        try:
            series.thresholds = 10
        except TypeError:
            pass
        # So should setting with more than 2 values
        try:
            series.thresholds = (-1,0,1)
        except TypeError:
            pass
        # Or when the values can't be converted to float
        try:
            series.thresholds = ('A','B')
        except ValueError:
            pass
    #
    def test_pickling(self):
        "Tests survival of pickled info"
        import cPickle
        series = self.series.copy()
        altinfo = dict(full_year=False, minimum_size = 5,
                       starting_month = 10, ensotag = '(trial)',
                       thresholds=(-0.5,0.5), reference_period=self.refperiod)
        series.optinfo.update(altinfo)
        pickled = cPickle.loads(cPickle.dumps(series))
        assert_equal(pickled._optinfo, altinfo)
        assert_equal(pickled.thresholds, altinfo['thresholds'])
    #
    def test_caching(self):
        "Tests caching"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        series.minimum_size = 6
        series.set_indices()
        cached = series._cachedmonthly
        # w/ slice
        test = series[:12]
        assert_equal(test.indices, cached['indices_monthly'][:12])
        # w/ reshape
        test = series.reshape(-1,12)
        assert_equal(test.indices, cached['indices_monthly'].reshape(-1,12))
    
    #
    def test_clustering(self):
        "Tests clustering values"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        base = np.zeros(len(series), dtype=int)
        base[(series._data >= 0.5)] = +1
        base[(series._data <= -0.5)] = -1
        klust = Cluster(base, 0)
        assert_equal(np.concatenate(series._clustered_values().clustered),
                     base)
    #
    def test_set_monthly_indices(self):
        "Tests setting monthly indices"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is None)
        midx = series.set_monthly_indices(minimum_size=6)
        assert_equal(midx,
                     [ 0, 0, 0,+1,+1,+1,+1,+1,+1, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is not None)
        midx = series.set_monthly_indices(minimum_size=5)
        assert_equal(midx,
                     [ 0, 0, 0,+1,+1,+1,+1,+1,+1, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,
                       0, 0, 0,-1,-1,-1,-1,-1, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
        midx = series.set_monthly_indices(minimum_size=5, reference_season='NDJ')
        assert_equal(midx,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
    #
    def test_set_annual_indices(self):
        "Tests setting annual indices"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        midx = series.set_annual_indices(minimum_size=5, reference_season='NDJ')
        assert_equal(midx,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1,+1,+1,+1, 0, 0,])
    #
    def test_ensoindices(self):
        "Test ENSO indices"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        series._cachedmonthly = {}
        # Accessing indices when nothing has been initialized should fail
        try:
            midx = series.ensoindices
        except ValueError:
            pass
        #
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is None)
        midx = series.set_indices(full_year=False, minimum_size=5,
                                  reference_season='NDJ')
        assert_equal(midx,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is not None)
        midx = series.ensoindices
        assert_equal(midx,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
        series.minimum_size = 6
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is None)
        assert_equal(series.indices,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,]
                     )
    #
    def test_enso_indices_w_lag(self):
        "Tests setting monthly indices"
        series = self.series
        series.thresholds = (-0.5, +0.5)
        # lag = 0
        midx = series.set_monthly_indices(minimum_size=6)
        assert_equal(midx,
                     [ 0, 0, 0,+1,+1,+1,+1,+1,+1, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
        # lag = 3
        midx = series.set_monthly_indices(minimum_size=6, lag=3)
        assert_equal(midx,
                     [-9,-9,-9, 0, 0, 0,+1,+1,+1,+1,+1,+1,
                       0,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,
                      -1,-1,-1,-1,+1,+1,+1,+1,+1,+1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0,+1,+1,+1,+1,+1,+1,+1,+1,+1, 0, 0,])
        # lag = 12
        midx = series.set_monthly_indices(minimum_size=6, lag=12)
        assert_equal(midx,
                     [-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,
                       0, 0, 0,+1,+1,+1,+1,+1,+1, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,])
    #
    def test_convert(self):
        series = self.series
        series.thresholds = (-0.5, +0.5)
        series.minimum_size = 5
        _cached = series._cachedmonthly.get('indices_monthly', None)
        self.failUnless(_cached is None)
        control = [ 0, 0, 0,+1,+1,+1,+1,+1,+1, 0,-1,-1,
                   -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                   -1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,
                    0, 0, 0,-1,-1,-1,-1,-1, 0, 0,+1,+1,
                   +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,]
        control = ts.time_series(control, dates=series._dates)
        assert_equal(series.indices, control)
        # Convert to daily
        dseries = series.convert('D')
        dcontrol = ts.lib.backward_fill(control.convert('D'))
        assert_equal(dseries.indices, dcontrol)
        #
        control = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                   +1,+1,+1,+1,+1,+1,+1,+1,+1,+1, 0, 0,]
        control = ts.time_series(control, dates=series._dates)
        assert_equal(dseries.set_indices(full_year=True, reference_season='NDJ'),
                     ts.lib.backward_fill(control.convert('D')))
    #
    def test_reshape(self):
        series = self.series
        series.thresholds = (-0.5, +0.5)
        series.minimum_size = 5
        indices = series.indices
        _series = series.reshape(-1,1)
        assert_equal(_series.shape, (len(series),1))
        assert_equal(_series.indices.shape, (len(series),1))
    #
    def test_sort(self):
        series =self.series
        series.thresholds = (-0.5, +0.5)
        series.minimum_size = 5
        indices = series.indices
        idx = series.argsort()
        _series = ma.sort(series)
        assert_equal(_series, series[idx])
        assert_equal(_series.indices, indices[idx])




class TestClimateSeries(TestCase):
    #
    def setUp(self):
        "Initializes"
        data = [
        -0.49,-0.38,-0.26, 0.54, 0.95, 0.60, 0.93, 0.80, 0.52, 0.15,-0.50,-0.58,
        -0.58,-0.70,-0.54,-0.80,-0.75, 0.33,-0.18,-0.55,-0.77,-0.73,-0.58,-0.68,
        -0.62, 0.75, 0.72, 0.50, 0.84, 0.64, 0.77,-0.03, 0.71, 0.17, 0.64,-0.12,
        -0.17, 0.34, 0.57,-0.57,-0.77,-0.56,-0.51,-0.80, 0.19,-0.83, 0.53, 0.72,
         0.77, 0.64, 0.57, 0.51, 0.51, 0.67, 0.70, 0.40, 0.29, 0.23, 0.13, 0.00]
        refperiod = ts.DateArray([ts.Date('M','1980-01-01'),
                                  ts.Date('M','1990-01-01')], freq='M')
        ensoi = ENSOIndicator(data, start_date=ts.Date('M','1980-01-01'))
        ensoi.thresholds = (-0.5,0.5)
        self.ensoi = ensoi

    def test_creation(self):
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        local_ensoi = series.ensoindicator
        assert_equal(local_ensoi, ensoi)
        assert_equal(local_ensoi.optinfo, ensoi.optinfo)

    def test_creation_shifted(self):
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates-24,
                               ensoindicator = ensoi)
        local_enso = series.ensoindicator
        assert_equal(local_enso._mask,
                     [1]*24+[0]*(len(ensoi)-24))
        assert_equal(local_enso.dates, series.dates)
        assert_equal(local_enso.compressed(), ensoi[:(len(ensoi)-24)])
        assert_equal(local_enso.optinfo, ensoi.optinfo)

    def test_indices(self):
        "Test indices"
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        # Try without initializiation: it should fail
        try:
            indc = series.ensoindices
        except ValueError:
            pass
        # So, initialize it:
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        indc = series.ensoindices
        assert_equal(indc,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                      -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                      +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,])
    #
    def test_indices_convert(self):
        "Test the conversion of ensoindices from one frequency to another."
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        control = ts.time_series([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,
                                  -1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,
                                  -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,
                                  +1,+1,+1,+1,+1,+1,+1, 0, 0, 0, 0, 0,],
                                  dates=ensoi.dates)
        assert_equal(series.ensoindices, control)
        # Conversion 'M' to 'D'
        dseries = series.convert('D')
        assert_equal(dseries.ensoindices,
                     ts.lib.backward_fill(control.convert('D')))
        # Conversion 'M' to 'A'
        aseries = series.convert('A', func=ma.mean)
        assert_equal(aseries.ensoindices,
                     mode(control.convert('A'), axis=1)[0].squeeze())
    #
    def test_phases(self):
        "Test indices"
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        # Test cold/neutral
        cold = series.cold
        assert_equal(cold._mask,
                     [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                       0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
                       0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,])
        neutral = series.neutral
        assert_equal(neutral._mask,
                     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                       1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
                       1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,])
        warm = series.warm
        assert_equal(warm._mask,
                     [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,])
    #
    def test_pickling(self):
        "Tests pickling"
        import cPickle
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        pickled = cPickle.loads(cPickle.dumps(series))
        assert_equal(series, pickled)
        assert_equal(pickled.ensoindicator, pickled.ensoindicator)
        assert_equal(series.ensoindices, pickled.ensoindices)
    #
    def test_apply_on_phase(self):
        "Test apply on phase"
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        func = ma.mean
        _cold = func(series[(series.ensoindices==-1).filled()])
        _neut = func(series[(series.ensoindices==0).filled()])
        _warm = func(series[(series.ensoindices==+1).filled()])
        _glob = func(series)
        control = ma.array([(_cold, _neut, _warm, _glob)],
                           dtype=[('cold',float),('neutral',float),
                                  ('warm',float),('global',float)])
        control['warm'] = masked
        result = apply_on_phase(series, func)
    #
    def test_reshape(self):
        "Test reshape"
        ensoi = self.ensoi
        series = ClimateSeries(np.random.rand(len(ensoi)),
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        series.set_ensoindices(minimum_size=5, reference_season='NDJ')
        #
        series = series.reshape(1, -1)
        assert_equal(series.ensoindicator.shape, series.shape)
        self.failUnless(series.ensoindicator.shape != ensoi.shape)
        assert_equal(series.ensoindicator, ensoi.reshape(1, -1))
        #
        series.shape = (-1, 1)
        assert_equal(series.ensoindicator.shape, series.shape)
        self.failUnless(series.ensoindicator.shape != ensoi.shape)
        assert_equal(series.ensoindicator, ensoi.reshape(-1, 1))
    #
    def test_getitem_with_string(self):
        "Test getitem w/ a string"
        ensoi = self.ensoi
        series = ClimateSeries(zip(np.random.rand(len(ensoi)),
                                   np.random.rand(len(ensoi)),),
                               dtype=[('A', float), ('B', float)],
                               dates=ensoi._dates,
                               ensoindicator = ensoi)
        # With field
        series_A = series['A']
        assert_equal(series_A.ensoindicator, ensoi)
        assert_equal(series['B'].ensoindicator, ensoi)
        # With date
        series_first = series['1980-01':'1981-01']
        assert_equal(series_first, series[:12])
    #
    def test_ensoindices_from_daily_indicator(self):
        "Try setting ensoindices from a daily indicator"
        ensoi = self.ensoi
        ensoi.set_monthly_indices(minimum_size=5, reference_season='NDJ')
        #
        bseries = climate_series(range(24),
                                 start_date=ts.Date('M', '1982-01'),
                                 ensoindicator=ensoi)
        control_indices = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
        control_indices = ts.time_series(control_indices,
                                         start_date=bseries.dates[0])
        # Basic check
        assert_equal(bseries.ensoindices, control_indices)
        # Define daily series
        dseries = climate_series(range(730),
                                 start_date=ts.Date('D', '1982-01-01'),
                                 ensoindicator=ensoi)
        mseries = climate_series(range(24),
                                 start_date=ts.Date('M', '1982-01'),
                                 ensoindicator=dseries.ensoindicator)
        assert_equal(mseries.ensoindices, bseries.ensoindices)
    #
    def test_ensoindices_on_irregular_series(self):
        "Try setting ensoindices on a non-regular series"
        ensoi = self.ensoi
        ensoi.set_monthly_indices(minimum_size=5, reference_season='NDJ')
        #
        bseries = climate_series(range(24),
                                 start_date=ts.Date('M', '1982-01'),
                                 ensoindicator=ensoi)
        control_indices = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
        control_indices = ts.time_series(control_indices,
                                         start_date=bseries.dates[0])
        tseries = bseries[::2]
        assert_equal(tseries.ensoindices, control_indices[::2])
        mseries = climate_series(tseries, ensoindicator=bseries.ensoindicator)
        assert_equal(mseries.ensoindices, bseries.ensoindices[::2])
        mseries = climate_series(tseries, ensoindicator=bseries.ensoindicator)
    #
    def test_ensoindices_on_empty_series(self):
        "Try on an empty series"
        ensoi = self.ensoi
        ensoi.set_monthly_indices(minimum_size=5, reference_season='NDJ')
        bseries = climate_series(range(24),
                                 start_date=hydro.Date('M', '1982-01'),
                                 ensoindicator=ensoi)
        bseries = bseries[[]]
        assert(bseries.ensoindicator is None)



class TestClimateRecords(TestCase):
    #
    def setUp(self):
        "Initializes"
        ndtype = [('lin',float),('rand',float)]
        dlin = np.linspace(0,10,120)
        drnd = np.random.rand(120)
        data = np.array(zip(dlin, drnd), dtype=ndtype)
        dates = ts.date_array(start_date=ts.now('M')-120, length=120, freq='M')
        enso = ENSOIndicator(np.random.rand(120) + np.linspace(-1,1,120), 
                             dates=dates,
                             thresholds=(-0.5,0.5),
                             full_year='False', refseason='NDH', minsize=5)
        cdat = data.view(ClimateRecords)
        cdat._dates = dates
        cdat.ensoindicator = enso
        self.dlin = dlin
        self.cdat=cdat
        self.enso=enso
    #
    def test_attribute_access_as_key(self):
        "Tests accessing an attribute as item"
        (cdat, dlin, enso) = (self.cdat, self.dlin, self.enso)
        cdlin = cdat['lin']
        assert_equal(cdlin, dlin)
        self.failUnless(isinstance(cdlin, ClimateSeries))
        assert_equal(cdlin._dates, cdat._dates)
        assert_equal(cdlin.optinfo, cdat.optinfo)
        assert_equal(cdlin.ensoindicator, enso)
    #
    def test_attribute_access_as_attribute(self):
        "Tests accessing an attribute as attribute"
        (cdat, dlin, enso) = (self.cdat, self.dlin, self.enso)
        cdlin = cdat.lin
        assert_equal(cdlin, dlin)
        self.failUnless(isinstance(cdlin, ClimateSeries))
        assert_equal(cdlin._dates, cdat._dates)
        assert_equal(cdlin.optinfo, cdat.optinfo)
        assert_equal(cdlin.ensoindicator, enso)
    #
    def test_ensoindices(self):
        "Tests the ensoindices"
        (cdat, enso) = (self.cdat, self.enso)
        assert_equal(cdat.ensoindices, enso.ensoindices)



#------------------------------------------------------------------------------
if __name__ == "__main__":
    run_module_suite()