"""
Test suite for :mod:`scikits.hydroclimpy.io`
"""

import numpy as np
import numpy.ma as ma
from numpy.ma.testutils import *

import scikits.hydroclimpy as hydro
from scikits.hydroclimpy.lib import evaporation as evapo


class TestEvaporation(TestCase):
    #
    def setUp(self):
        self.defaults = ['Hamon', 'Hansen', 'Hargreaves', 'Kharrufa',
                         'Makkink', 'PenmanMonteith', 'PenmanMonteithASCE',
                         'PriestleyTaylor', 'Thornthwaite', 'Turc']
        self.latitude = 40


    def test_daily_random(self):
        "Test some PET models on some random daily data"
        (freq, base) = ('D', 365)
        wave = -10 * np.cos(np.arange(120) * (2.*np.pi)/base)
        tmin = hydro.time_series(5 + wave + np.random.rand(120) * 3,
                                 start_date=hydro.Date(freq, '1990-01-01'))
        tmax = hydro.time_series(10 + wave + np.random.rand(120) * 3,
                                 start_date=hydro.Date(freq, '1990-01-01'))
        ma.putmask(tmax, (tmax < tmin), tmin+1)
        freeze = ((tmin + tmax)/2. <= 0)
        PET = evapo.PotentialEvapoTranspiration(tmin, tmax, self.latitude)
        for d in self.defaults:
            err_msg = "Error w/ model '%s' : the PET should be masked" % d
            P = ma.masked_values(getattr(PET, d).__call__(), 0)
            if d in ['Hamon', 'Kharrufa', 'Turc']:
                assert_equal(ma.getmaskarray(P), freeze, err_msg=err_msg)
            else:
                assert(not ma.getmaskarray(P).any())


    def test_monthly_random(self):
        "Test some PET models on some random monthly data"
        (freq, base) = ('M', 12)
        wave = -10 * np.cos(np.arange(120) * (2.*np.pi)/base)
        tmin = hydro.time_series(5 + wave + np.random.rand(120) * 3,
                                 start_date=hydro.Date(freq, '1990-01-01'))
        tmax = hydro.time_series(10 + wave + np.random.rand(120) * 3,
                                 start_date=hydro.Date(freq, '1990-01-01'))
        ma.putmask(tmax, (tmax < tmin), tmin+1)
        freeze = ((tmin + tmax)/2. <= 0)
        PET = evapo.PotentialEvapoTranspiration(tmin, tmax, self.latitude)
        for d in self.defaults:
            err_msg = "Error w/ model '%s' : the PET should be masked" % d
            P = ma.masked_values(getattr(PET, d).__call__(), 0)
            if d in ['Hamon', 'Kharrufa', 'Thornthwaite', 'Turc']:
                assert_equal(ma.getmaskarray(P), freeze, err_msg=err_msg)
            else:
                assert(not ma.getmaskarray(P).any())
            assert((P.filled(0) >= 0).all())
