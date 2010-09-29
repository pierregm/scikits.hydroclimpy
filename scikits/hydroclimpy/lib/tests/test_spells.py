"""
Test suite for the spells modules

@author: pierregm
"""


from numpy.ma.testutils import *
import scikits.hydroclimpy as hydro
from scikits.hydroclimpy.lib.spells import dry_spells, wet_spells

class TestSpells(TestCase):
    #
    def setUp(self):
        data = [1, 4, 3, 2, 0, 1, 0, 0, 0, 0,   3, 3,
                3, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0.1, 5,]
        self.rainfall = hydro.time_series(data,
                                          start_date=hydro.Date('M', '2001-01'))
        self.data = data
    #
    def test_dry_spells(self):
        "Test dry_spells"
        rainfall = self.rainfall
        dates = rainfall.dates
        control_dates_onstart = dates[[4, 6, 13, 22]]
        control_dates_onend = dates[[4, 9, 20, 22]]
        spells = [1, 4, 8, 1]
        #
        _test = dry_spells(rainfall, threshold=0.254)
        assert_equal(_test, spells)
        assert_equal(_test.dates, control_dates_onstart)
        _test = dry_spells(rainfall, threshold=0.254, start=False)
        assert_equal(_test, spells)
        assert_equal(_test.dates, control_dates_onend)
    #
    def test_dry_spells_on_list(self):
        "Test dry spells on list"
        data = self.data
        spells = [1, 4, 8, 1]
        _test = dry_spells(data, threshold=0.254)
        assert_equal(_test, spells)
        self.failUnless(getattr(_test, 'dates', None) is None)
        assert(isinstance(_test, np.ndarray))
    #
    def test_wet_spells(self):
        "Test wet_spells"
        rainfall = self.rainfall
        dates = rainfall.dates
        control_dates_onstart = dates[[0, 5, 10, 21, 23]]
        control_dates_onend = dates[[3, 5, 12, 21, 23]]
        durations = [4, 1, 3, 1, 1]
        intensities = [10, 1, 9, 7, 5]
        #
        _test = wet_spells(rainfall, threshold=0.254)
        assert_equal(_test['durations'], durations)
        assert_equal(_test['intensities'], intensities)
        assert_equal(_test.dates, control_dates_onstart)
        _test = wet_spells(rainfall, threshold=0.254, start=False)
        assert_equal(_test['durations'], durations)
        assert_equal(_test['intensities'], intensities)
        assert_equal(_test.dates, control_dates_onend)
        #
        _test = wet_spells(rainfall, threshold=0.254, mode='durations')
        assert_equal(_test, durations)
        _test = wet_spells(rainfall, threshold=0.254, mode='d')
        assert_equal(_test, durations)
        _test = wet_spells(rainfall, threshold=0.254, mode='intensities')
        assert_equal(_test, intensities)
        _test = wet_spells(rainfall, threshold=0.254, mode='i')
        assert_equal(_test, intensities)
    #
    def test_wet_spells_on_list(self):
        "Test wet_spells"
        data = self.data
        durations = [4, 1, 3, 1, 1]
        intensities = [10, 1, 9, 7, 5]
        _test = wet_spells(data, threshold=0.254)
        assert_equal(_test['durations'], durations)
        assert_equal(_test['intensities'], intensities)
    #
    def test_wet_and_dry_spells(self):
        "Test wet and dry spells"
        data = self.data + [-3]
        rainfall = hydro.time_series(data, mask=24*[0]+[1],
                                     start_date=hydro.Date('M', '2001-01'))
        threshold = 0.254
        dspells = dry_spells(rainfall, threshold=threshold)
        wspells = wet_spells(rainfall, threshold=threshold)
        assert_equal(dspells.sum()+wspells.durations.sum(), rainfall.count())

