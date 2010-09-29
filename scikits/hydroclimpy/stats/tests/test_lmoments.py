"""
Test suite for the lmoments module
"""


from numpy.ma.testutils import *

import scipy.stats.distributions as dist
import scikits.hydroclimpy.stats.lmoments as lmoments
from scikits.hydroclimpy.stats.extradistributions import \
    gennorm, glogistic, kappa, pearson3, wakeby



class TestLMoments(TestCase):
    #
    def test_lmoments(self):
        "Try getting L moments"
        #
        self.failUnlessRaises(ValueError, lmoments, np.arange(10), -1)
        self.failUnlessRaises(ValueError, lmoments, np.arange(10), 4., -1., -1)
    #
    def test_lstats_expon(self):
        "Try lstats on the exponential distribution"
        _exp = [ 1.        ,  0.5       ,  0.33333333,  0.16666667,  0.1       ,
                 0.06666667,  0.04761905,  0.03571429,  0.02777778,  0.02222222,
                 0.01818182,  0.01515152,  0.01282051,  0.01098901,  0.00952381,
                 0.00833333,  0.00735294,  0.00653595,  0.00584795,  0.00526316]
        assert_almost_equal(dist.expon.lstats(20), np.array(_exp))
        _exp = ( 4.        ,  0.5       ,  0.33333333,  0.16666667,  0.1)
        assert_almost_equal(dist.expon.lstats(5, loc=3., scale=1.),
                            np.array(_exp))
        assert_almost_equal(dist.expon(3., 1.).lstats(5), np.array(_exp))
    #
    def test_lstats_gamma(self):
        "Try lstats on the gamma distribution"
        _gam = [ 1.,  0.5, 0.33333334, 0.16666667]
        assert_almost_equal(dist.gamma.lstats(4, 1), np.array(_gam))
        assert_almost_equal(dist.gamma(1.,).lstats(4), np.array(_gam))
        _gam =  [10.0, 1.76197052, 0.10350336, 0.12585956]
        assert_almost_equal(dist.gamma.lstats(4, 10.), np.array(_gam))
        assert_almost_equal(dist.gamma(10.,).lstats(4), np.array(_gam))
    #
    def test_lstats_gaussian(self):
        "Try lstats on Gaussian"
        _nor = [ 0.        ,  0.56418958,  0.        ,  0.12260172,  0.        ,
                 0.04366115,  0.        ,  0.02184314,  0.        ,  0.0129635 ,
                 0.        ,  0.00852962,  0.        ,  0.00601389,  0.        ,
                 0.00445558,  0.        ,  0.00342643,  0.        ,  0.00271268]
        assert_almost_equal(dist.norm.lstats(20), _nor)
        assert_almost_equal(dist.norm.lstats(20, 0., 1.), _nor)
        _nor = [ 1.        ,  1.69256875,  0.        ,  0.12260172,  0.        ,
                 0.04366115,  0.        ,  0.02184314,  0.        ,  0.0129635 ,
                 0.        ,  0.00852962,  0.        ,  0.00601389,  0.        ,
                 0.00445558,  0.        ,  0.00342643,  0.        ,  0.00271268]
        assert_almost_equal(dist.norm(1., 3.).lstats(20), _nor)
        assert_almost_equal(dist.norm.lstats(20, 1., 3.), _nor)
    #
    def test_lstats_gev(self):
        "Try lstats on the GEV distribution"
        _gev = [-0.5       ,  0.75      , -0.62962963,  0.39814815, -0.27444444,
                 0.20148148, -0.1547997 ,  0.12304422, -0.10040417,  0.08365667,
                -0.0708951 ,  0.0609311 , -0.05299155,  0.04655511, -0.04125945,
                 0.03684606, -0.03312638,  0.02996009, -0.027241  ,  0.02488765]
        assert_almost_equal(dist.genextreme.lstats(20, 2, loc=0, scale=1),
                            _gev)
        assert_almost_equal(dist.genextreme(2, 0, 1).lstats(20), _gev)
        _gev = [ 0.        ,  0.5       , -0.33333333,  0.16666667, -0.1       ,
                 0.06666667, -0.04761905,  0.03571429, -0.02777778,  0.02222222,
                -0.01818182,  0.01515152, -0.01282051,  0.01098901, -0.00952381,
                 0.00833334, -0.00735302,  0.00653646, -0.00585061,  0.00527522]
        assert_almost_equal(dist.genextreme.lstats(20, 1, loc=0, scale=1),
                            _gev)
        assert_almost_equal(dist.genextreme(1, 0, 1).lstats(20), _gev)
        _gev = [ 0.        ,  0.69314718,  0.169925  ,  0.15037499,  0.05586835,
                 0.05811002,  0.02762426,  0.03055638,  0.01646503,  0.01878466,
                 0.01093282,  0.01269731,  0.00778983,  0.00914836,  0.00583332,
                 0.00690104,  0.00453268,  0.00538917,  0.00362408,  0.00432388]
        assert_almost_equal(dist.genextreme.lstats(20, 0, loc=0, scale=1),
                            _gev)
        assert_almost_equal(dist.genextreme(0, 0, 1).lstats(20), _gev)
    #
    def test_lstats_gumbel(self):
        "Try lstats on Gumbel_r"
        _gum = [ 0.57721566,  0.69314718,  0.169925  ,  0.15037499,  0.05586835,
                 0.05811002,  0.02762426,  0.03055638,  0.01646503,  0.01878466,
                 0.01093282,  0.01269731,  0.00778983,  0.00914836,  0.00583332,
                 0.00690104,  0.00453268,  0.00538917,  0.00362408,  0.00432388]
        assert_almost_equal(dist.gumbel_r.lstats(20), _gum)
        assert_almost_equal(dist.gumbel_r.lstats(20, 0., 1.), _gum)
        _gum = [ 2.73164699,  2.07944154,  0.169925  ,  0.15037499,  0.05586835,
                 0.05811002,  0.02762426,  0.03055638,  0.01646503,  0.01878466,
                 0.01093282,  0.01269731,  0.00778983,  0.00914836,  0.00583332,
                 0.00690104,  0.00453268,  0.00538917,  0.00362408,  0.00432388]
        assert_almost_equal(dist.gumbel_r.lstats(20, 1, 3), _gum)
        assert_almost_equal(dist.gumbel_r(1., 3.).lstats(20), _gum)
    #
    def test_lstats_genpareto(self):
        "Try LMT on Generalized Pareto"
        _gpa = [ 2.        ,  1.33333333,  0.6       ,  0.42857143,  0.33333333,
                 0.27272727,  0.23076923,  0.2       ,  0.17647059,  0.15789474,
                 0.14285714,  0.13043478,  0.12      ,  0.11111111,  0.10344828,
                 0.09677419,  0.09090909,  0.08571429,  0.08108108,  0.07692308]
        assert_almost_equal(dist.genpareto.lstats(20, 0.5, loc=0, scale=1),
                            _gpa)
        _gpa = assert_almost_equal(dist.genpareto.lstats(20, 0.5), _gpa)
        _gpa = [  6.66666667e-01,   2.66666667e-01,   1.42857143e-01,
                  4.76190476e-02,   2.16450216e-02,   1.16550117e-02,
                  6.99300699e-03,   4.52488688e-03,   3.09597523e-03,
                  2.21141088e-03,   1.63452109e-03,   1.24223602e-03,
                  9.66183575e-04,   7.66283525e-04,   6.17970585e-04,
                  5.05612296e-04,   4.18935903e-04,   3.51000351e-04,
                  2.97000297e-04,   2.53536839e-04]
        assert_almost_equal(dist.genpareto.lstats(20, -0.5, loc=0, scale=1),
                            _gpa)
    #
    def test_lmparams_expon(self):
        "Test parameter estimation from L_moments for the exponential distrib."
        _expon = dist.expon
        (loc, scl) = (0., 1.)
        assert_almost_equal(_expon.lmparams(_expon(loc, scl).lstats(4)),
                                            np.array((loc, scl)))
        (loc, scl) = (0.5, 0.5)
        assert_almost_equal(_expon.lmparams(_expon(loc, scl).lstats(4)),
                                            np.array((loc, scl)))
        (loc, scl) = (-0.5, 2.)
        assert_almost_equal(_expon.lmparams(_expon(loc, scl).lstats(4)),
                                            np.array((loc, scl)))
    #
    def test_lmparams_gamma(self):
        "Test parameter estimation from L_moments for the gamma distribution"
        _gamma = dist.gamma
        (loc, scl, shp) = (0., 1., 1.)
        assert_almost_equal(_gamma.lmparams(_gamma(shp, loc, scl,).lstats(4)),
                                            np.array((loc, scl, shp)),
                            6)
        (loc, scl, shp) = (0., 0.5, 1.)
        assert_almost_equal(_gamma.lmparams(_gamma(shp, loc, scl).lstats(4)),
                                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0., 2.0, 1.)
        assert_almost_equal(_gamma.lmparams(_gamma(shp, loc, scl).lstats(4)),
                                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0, 2, 1)
        assert_almost_equal(_gamma.lmparams(_gamma(shp, loc, scl).lstats(4)),
                                            np.array((loc, scl, shp)), 6)
    #
    def test_lmparams_gaussian(self):
        "Test parameter estimation from L_moments for the normal distribution"
        _norm = dist.norm
        (loc, scl) = (0., 1.)
        assert_almost_equal(_norm.lmparams(_norm(loc, scl).lstats(4)),
                            np.array((loc, scl)))
        (loc, scl) = (0.5, 0.5)
        assert_almost_equal(_norm.lmparams(_norm(loc, scl).lstats(4)),
                            np.array((loc, scl)))
        (loc, scl) = (-0.5, 2.)
        assert_almost_equal(_norm.lmparams(_norm(loc, scl).lstats(4)),
                            np.array((loc, scl)))
        #
        moms = (0., 1.)
        assert_almost_equal(_norm.lmparams(moms), np.array((0., np.pi**0.5)))
        self.failUnlessRaises(ValueError, _norm.lmparams, (1., 0.))
    #
    def test_lmparams_gev(self):
        "Test parameter estimation from L_moments for the GEV distribution"
        _genextreme = dist.genextreme
        (loc, scl, shp) = (0., 1., 1.)
        lmoms = _genextreme(shp, loc, scl,).lstats(4)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0., 0.5, 1.)
        lmoms = _genextreme(shp, loc, scl).lstats(4)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0., 2.0, 1.)
        lmoms = _genextreme(shp, loc, scl).lstats(4)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0, 2, 1)
        lmoms = _genextreme(shp, loc, scl).lstats(4)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        #
        lmoms = (-1., 1., 3.)
        #self.failUnlessRaises(ValueError, pel_genextreme, lmoms)
        lmoms = (-1., 1., 1./3)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array([-1.96741479,  1.09472348, -0.2393871 ]))
        lmoms = (-1., 1., -0.9)
        assert_almost_equal(_genextreme.lmparams(lmoms),
                            np.array([ 0.02398541,  0.15472353,  4.10559021]))
        lmoms = (-1., 1., -0.98)
        assert_almost_equal(np.array(_genextreme.lmparams(lmoms)),
                            [1.02682873e-02, 3.18386104e-03, 6.55440863])
    #
    def test_lmparams_gumbel(self):
        "Test parameter estimation from L_moments for the Gumbel distribution"
        _gumbel_r = dist.gumbel_r
        (loc, scl) = (0., 1.)
        lmoms = _gumbel_r(loc, scl).lstats(4)
        assert_almost_equal(_gumbel_r.lmparams(lmoms), np.array((loc, scl)), 6)
        (loc, scl) = (0., 0.5)
        lmoms = _gumbel_r(loc, scl).lstats(4)
        assert_almost_equal(_gumbel_r.lmparams(lmoms), np.array((loc, scl)), 6)
        (loc, scl) = (0., 2.0)
        lmoms = _gumbel_r(loc, scl).lstats(4)
        assert_almost_equal(_gumbel_r.lmparams(lmoms), np.array((loc, scl)), 6)
        (loc, scl) = (0, 2)
        lmoms = _gumbel_r(loc, scl).lstats(4)
        assert_almost_equal(_gumbel_r.lmparams(lmoms), np.array((loc, scl)), 6)
        #
        lmoms = (-1., 1., 3.)
        #self.failUnlessRaises(ValueError, pel_genextreme, lmoms)
        lmoms = (-1., 1., 1./3)
        assert_almost_equal(_gumbel_r.lmparams(lmoms),
                            np.array([-1.83274618,  1.44269504]))
    #
    def test_lmparams_genpareto(self):
        "Test parameter estimation from L_moments for the GEV distribution"
        _genpareto = dist.genpareto
        (loc, scl, shp) = (0.0, 1., 0.999)
        lmoms = _genpareto(shp, loc, scl,).lstats(4)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.0, 0.5, -0.999)
        lmoms = _genpareto(shp, loc, scl).lstats(4)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.5, 2.0, 0.5)
        lmoms = _genpareto(shp, loc, scl).lstats(4)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0, 2, 0)
        lmoms = _genpareto(shp, loc, scl).lstats(4)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        #
        lmoms = (-1., 1., 3.)
        #self.failUnlessRaises(ValueError, pel_genpareto, lmoms)
        lmoms = (-1., 1., 1./3)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array([-3., 2., 0.]))
        lmoms = (-1., 1., -0.9)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array([-40., 1482., -37.]))
        lmoms = (-1., 1., +0.9)
        assert_almost_equal(_genpareto.lmparams(lmoms),
                            np.array([-2.10526316,  0.11634349, +0.89473684]))
        lmoms = (+1., +1., +0.7)
        assert_almost_equal(np.array(_genpareto.lmparams(lmoms)),
                            np.array([-0.35294118,  0.47750865, +0.64705882]))
    #
    def test_lmparams_glogistic(self):
        "Test parameter estimation from L_moments for the generalized logistic distribution"
        (loc, scl, shp) = (0.0, 1., 0.999)
        lmoms = glogistic(shp, loc, scl,).lstats(4)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.0, 0.5, -0.999)
        lmoms = glogistic(shp, loc, scl).lstats(4)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.5, 2.0, 0.5)
        lmoms = glogistic(shp, loc, scl).lstats(4)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0, 2, 0)
        lmoms = glogistic(shp, loc, scl).lstats(4)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        #
        lmoms = (-1., 1., 3.)
        #self.failUnlessRaises(ValueError, pel_glogistic, lmoms)
        lmoms = (-1., 1., 1./3)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array([-1.51901997,  0.82699334, -0.33333333]))
        lmoms = (-1., 1., -0.9)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array([-0.01032489, 0.1092924, 0.9]))
        lmoms = (-1., 1., +0.9)
        assert_almost_equal(glogistic.lmparams(lmoms),
                            np.array([-1.98967511, 0.1092924, -0.9]))
        lmoms = (+1., +1., +0.7)
        assert_almost_equal(np.array(glogistic.lmparams(lmoms)),
                            np.array([ 0.09697573,  0.36788301, -0.7]))
    #
    def test_lmparams_gennorm(self):
        "Test parameter estimation from L_moments for the generalized normal distribution"
        (loc, scl, shp) = (0.0, 1., 0.999)
        lmoms = gennorm(shp, loc, scl,).lstats(4)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.0, 0.5, -0.999)
        lmoms = gennorm(shp, loc, scl).lstats(4)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0.5, 2.0, 0.5)
        lmoms = gennorm(shp, loc, scl).lstats(4)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        (loc, scl, shp) = (0, 2, 0)
        lmoms = gennorm(shp, loc, scl).lstats(4)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array((loc, scl, shp)), 6)
        #
        lmoms = (-1., 1., 3.)
        #self.failUnlessRaises(ValueError, pel_gennorm, lmoms)
        lmoms = (-1., 1., 1./3)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array([-1.57345421, 1.44332159, -0.7010033]))
        lmoms = (-1., 1., -0.9)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array([ 0.03455279, 0.0987164, 2.58243705]))
        lmoms = (-1., 1., +0.9)
        assert_almost_equal(gennorm.lmparams(lmoms),
                            np.array([-2.03455279, 0.0987164, -2.58243705]))
        lmoms = (+1., +1., +0.7)
        assert_almost_equal(np.array(gennorm.lmparams(lmoms)),
                            np.array([0.014923134, 0.54923500, -1.66231215]))
    #
    def test_lmparams_wakeby(self):
        "Test parameter estimation from L_moments for the Wakeby distribution"
        (e, a, b, c, d) = (0.0, 1.0,  1.0, 1.0, 0.99)
        lmoms = wakeby(b, c, d, e, a,).lstats(5)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array((e, a, b, c, d)), 6)
        (e, a, b, c, d) = (0.0, 0.5, -0.5, 1.0, 0.99)
        lmoms = wakeby(b, c, d, e, a).lstats(5)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array((e, a, b, c, d)), 6)
        (e, a, b, c, d) = (0.5, 2.0,  0.5, 1.0, 0.5)
        lmoms = wakeby(b, c, d, e, a).lstats(5)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array((e, a, b, c, d)), 6)
        (e, a, b, c, d) = (0.0, 2.0,  0.5, 0.5, 0.99)
        lmoms = wakeby(b, c, d, e, a).lstats(5)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array((e, a, b, c, d)), 6)
        #
        lmoms = (-1., 1., 3., 1., 1.)
        #self.failUnlessRaises(ValueError, pel_wakeby, lmoms)
        lmoms = (-1., 1., 1./3, 0.999, 0.999)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array([-3., 2., 0., 0., 0.]))
        lmoms = (-1., 1., 0., 0., 0.)
        assert_almost_equal(wakeby.lmparams(lmoms),
                            np.array([-4., 6., 1., 0., 0.]))

    def test_lmparams_kappa(self):
        "Test parameter estimation from L_moments for the Wakeby distribution"
        (loc, scale, h, k) = (0.0, 1.0,  1.0, 1.0)
        lmoms = kappa(h, k, loc, scale).lstats(4)
        assert_almost_equal(kappa.lmparams(lmoms),
                            np.array((loc, scale, h, k)), 6)
        (loc, scale, h, k) = (0.0, 0.5, -0.5, 0.99)
        lmoms = kappa(h, k, loc, scale).lstats(4)
        assert_almost_equal(kappa.lmparams(lmoms),
                            np.array((loc, scale, h, k)), 6)
        (loc, scale, h, k) = (0.5, 2.0, 1.0, 0.5)
        lmoms = kappa(h, k, loc, scale).lstats(4)
        assert_almost_equal(kappa.lmparams(lmoms),
                            np.array((loc, scale, h, k)), 6)
        (loc, scale, h, k) = (0.0, 2.0,  0.5, 0.99)
        lmoms = kappa(h, k, loc, scale).lstats(4)
        assert_almost_equal(kappa.lmparams(lmoms),
                            np.array((loc, scale, h, k)), 6)



if __name__ == "__main__":
    run_module_suite()

