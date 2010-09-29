"""
Test suite for the lmoments module
"""


from numpy.ma.testutils import *

import scipy.stats.distributions as dist
from scikits.hydroclimpy.stats.extradistributions import \
    gennorm, glogistic, kappa, logseries, pearson3, wakeby, ztnbinom


class TestKappa(TestCase):
    "Test the Kappa distribution"
    #
    def test_cdf(self):
        "Test the CDF of the Kappa distribution"
        (k, h) = (1., 1.)
        control = (0., 0.25, 0.5, 1.)
        assert_almost_equal(kappa.cdf([0., 0.25, 0.5, 1.], k, h), control)
        assert_almost_equal(kappa(k, h).cdf([0., 0.25, 0.5, 1.]), control)
        (k, h) = (0.1, 1.)
        control = [0.0, 0.22367038, 0.401263061, 0.65132156]
        assert_almost_equal(kappa(k, h).cdf([0., 0.25, 0.5, 1.]), control)
        (k, h) = (1., 0.1)
        control = [0.34867844, 0.45858234, 0.59873694, 1.0]
        assert_almost_equal(kappa(k, h).cdf([0., 0.25, 0.5, 1.]), control)
        (k, h) = (0.5, 1.)
        control = [0, 0.234375, 0.4375, 0.75]
        assert_almost_equal(kappa(k, h).cdf([0., 0.25, 0.5, 1.]), control)
        (k, h) = (0., 1.)
        control = [0., 0.221199217, 0.393469340, 0.63212055883]
        assert_almost_equal(kappa(k, h).cdf([0., 0.25, 0.5, 1.]), control)
    #
    def test_ppf(self):
        "Test the PPF of the Kappa distribution"
        (k, h) = (1., 1.)
        control = (0.01, 0.25, 0.5, 0.99)
        assert_almost_equal(kappa.ppf([0.01, 0.25, 0.5, 0.99], k, h), control)
        assert_almost_equal(kappa(k, h).ppf([0.01, 0.25, 0.5, 0.99]), control)
        (k, h) = (0.1, 1.)
        control = [0.01004529, 0.28358342, 0.669670085, 3.69042656]
        assert_almost_equal(kappa.ppf([0.01, 0.25, 0.5, 0.99], k, h), control)
        assert_almost_equal(kappa(k, h).ppf([0.01, 0.25, 0.5, 0.99]), control)
        (k, h) = (1., 0.1)
        control = [-2.6904266, -0.29449437, 0.33032992, 0.9899547]
        assert_almost_equal(kappa(k, h).ppf([0.01, 0.25, 0.5, 0.99]), control)



class TestGLogistic(TestCase):
    "Test the alternate generaized logistic distribution"
    #
    def test_cdf(self):
        "Try the CDF of the glogistic distribution"
        c = 1.
        _glo = [0.5, 0.57142857, 0.66666667, 1.0]
        assert_almost_equal(glogistic(c).cdf([0., 0.25, 0.5, 1.0]), _glo)
        c = 0.1
        _glo =  [0.5, 0.56295858, 0.62549377, 0.74146659]
        assert_almost_equal(glogistic(c).cdf([0., 0.25, 0.5, 1.0]), _glo)
        c = 0.
        _glo =  [0.5, 0.5621765009, 0.622459331, 0.7310585786]
        assert_almost_equal(glogistic(c).cdf([0., 0.25, 0.5, 1.0]), _glo)
    #
    def test_ppf(self):
        "Try the PPF of the glogistic distribution"
        c = 0.
        _glo = [-6.906754779, -1.098612289, 0.0, 6.906754779]
        assert_almost_equal(glogistic(c).ppf([0.001, 0.25, 0.5, 0.999]), _glo)
        c = 0.1
        _glo =  [-9.95062699, -1.16123174, 0.0, 4.98762620]
        assert_almost_equal(glogistic(c).ppf([0.001, 0.25, 0.5, 0.999]), _glo)
        c = 1.
        _glo =  [-998.0, -2.0, 0.0, 0.998998998999]
        assert_almost_equal(glogistic(c).ppf([0.001, 0.25, 0.5, 0.999]), _glo)



class TestGenNorm(TestCase):
    "Test the generalized normal distribution"
    #
    def test_cdf(self):
        "Try the CDF of the Generalized Normal distribution"
        k = 0.
        ctr =  [0.5, 0.59870633, 0.69146246, 0.84134475]
        assert_almost_equal(gennorm(k).cdf([0., .25, .5, 1.]), ctr)
        k = 0.1
        ctr =  [0.5, 0.59993470, 0.69600089, 0.853968136]
        assert_almost_equal(gennorm(k).cdf([0., .25, .5, 1.]), ctr)
        k = 1.
        ctr =   [0.5, 0.613204943, 0.755891404, 1.0]
        assert_almost_equal(gennorm(k).cdf([0., .25, .5, 1.]), ctr)
    #
    def test_ppf(self):
        "Try the PPF of the generalized normal distribution"
        k = 0.
        ctr = [-3.0902323062, -0.674489750, 0.0, 3.090232306]
        assert_almost_equal(gennorm(k).ppf([0.001, 0.25, 0.5, 0.999]), ctr)
        k = 0.1
        ctr = [-3.6209401242, -0.697756729, 0.0, 2.658362852]
        assert_almost_equal(gennorm(k).ppf([0.001, 0.25, 0.5, 0.999]), ctr)
        k = 1.
        ctr =  [-20.98218398, -0.96303108, 0.0, 0.95450861]
        assert_almost_equal(gennorm(k).ppf([0.001, 0.25, 0.5, 0.999]), ctr)



class TestWakeby(TestCase):
    "Test the Wakeby distribution"
    #
    def test_ppf(self):
        "Test the PPF of the Wakeby distribution"
        quantiles = [0.01, 0.25, 0.5, 0.99]
        (b, g, d) = (1, 1, 1)
        _wak = [0.02010101, 0.58333333, 1.5, 99.98999999999]
        assert_almost_equal(wakeby(b, g, d).ppf(quantiles), _wak)
        (b, g, d) = (.1, .1, .1)
        _wak = [0.011050826, 0.312769430, 0.741443547, 4.275319748]
        assert_almost_equal(wakeby(b, g, d).ppf(quantiles), _wak)
        (b, g, d) = (0, 0, 0)
        _wak = [0.010050336, 0.28768207, 0.693147181, 4.605170186]
        assert_almost_equal(wakeby(b, g, d).ppf(quantiles), _wak)
        (b, g, d) = (1., 0, 1.)
        assert_almost_equal(np.isnan(wakeby(b, g, d).ppf(quantiles)),
                            (True, True, True, True))
        (b, g, d) = (1., 0, 0.)
        assert_almost_equal(wakeby(b, g, d).ppf(quantiles), quantiles)
        (b, g, d) = (-1., 1., 2.)
        _wak = [0.0202530354, 0.72222222, 2.5, 5098.5]
    #
    def test_cdf(self):
        "Test the CDF of the Wakeby distribution"
        x = [-1., 0, 0.25, .5, 1.]
        (b, g, d) = (1., 1, 1.)
        _wak = [0.0, 0.0, 0.117217781, 0.219223594, 0.38196601]
        assert_almost_equal(wakeby(b, g, d).cdf(x), _wak)
        (b, g, d) = (0.1, 0.1, 0.1)
        _wak = [0.0, 0.0, 0.204994002, 0.3707061778, 0.61100326]
        assert_almost_equal(wakeby(b, g, d).cdf(x), _wak)
        (b, g, d) = (0., 0., 0.)
        _wak = [0.0, 0.0, 0.221199217, 0.3934693403, 0.632120559]
        assert_almost_equal(wakeby(b, g, d).cdf(x), _wak)
        (b, g, d) = (1., 0, 0.)
        _wak = [0.0, 0.0, 0.25, 0.50, 1.0]
        assert_almost_equal(wakeby(b, g, d).cdf(x), _wak)
        (b, g, d) = (-1., 1., 2.)
        _wak = [0.0, 0.0, 0.108194188, 0.1909830056, 0.310102051]
        assert_almost_equal(wakeby(b, g, d).cdf(x), _wak)
        (b, g, d) = (1., 0, 1.)
        assert_almost_equal(np.isnan(wakeby(b, g, d).cdf(x)),
                            (True, True, True, True, True))



class TestPEIII(TestCase):
    "Test the Pearson III distribution"
    def test_pdf(self):
        "Try the pdf for the pearson3 distribution"
        (g, m, s) = (1., 0., 1.)
        assert_almost_equal(pearson3(g).pdf([0., .25, .5, 1.]),
                            [0.3907336, 0.3374358, 0.2807478, 0.1784702])
        assert_almost_equal(pearson3(g, m, s).pdf([0., .25, .5, 1.]),
                            [0.3907336, 0.3374358, 0.2807478, 0.1784702])
        (g, m, s) = (-1., 0., 1.)
        assert_almost_equal(pearson3(g, m, s).pdf([0., .25, .5, 1.]),
                            [0.3907336, 0.4315709, 0.4480836, 0.3608941])
        (g, m, s) = (-1., -1., 0.5)
        assert_almost_equal(pearson3(g, m, s).pdf([0., .25, .5, 1.]),
                            [0.000000, -1.812188, -39.408299, -2329.521068], 6)
        (g, m, s) = (+np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).pdf([0., .25, .5, 1.]),
                            [0.10564050, 0.09362841, 0.08359751, 0.06777973], 6)
        (g, m, s) = (-np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).pdf([0., .25, .5, 1.]),
                            [0.1056405, 0.1203153, 0.1387192, 0.1953728], 6)
    #
    def test_cdf(self):
        "Try the cdf for the pearson3 distribution"
        (g, m, s) = (1., 0., 1.)
        assert_almost_equal(pearson3(g).cdf([0., .25, .5, 1.]),
                            [0.5665299, 0.6577040, 0.7349741, 0.8487961])
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            [0.5665299, 0.6577040, 0.7349741, 0.8487961])
        (g, m, s) = (-1., 0., 1.)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            [0.4334701, 0.5366327, 0.6472319, 0.8571235])
        (g, m, s) = (-1., -1., 0.5)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            [1.000000, 1.000000, 1.000000, 1.000000], 6)
        (g, m, s) = (+np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            [0.7003103, 0.7251716, 0.7472890, 0.7849331], 6)
        (g, m, s) = (-np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            [0.2996897, 0.3278695, 0.3601563, 0.4422018], 6)
    #
    def test_cdf_normal_case(self):
        (g, m, s) = (0., 0., 1.)
        assert_almost_equal(pearson3(g).cdf([0., .25, .5, 1.]),
                            dist.norm().cdf([0., .25, .5, 1.]))
        assert_almost_equal(pearson3(g).pdf([0., .25, .5, 1.]),
                            dist.norm().pdf([0., .25, .5, 1.]))
        assert_almost_equal(pearson3(g).ppf([0.001, 0.25, 0.50, 0.999]),
                            dist.norm().ppf([0.001, 0.25, 0.50, 0.999]),)
        assert_almost_equal(pearson3(g).sf([0., .25, .5, 1.]),
                            dist.norm().sf([0., .25, .5, 1.]))
        assert_almost_equal(pearson3(g).isf([0.001, 0.25, 0.50, 0.999]),
                            dist.norm().isf([0.001, 0.25, 0.50, 0.999]))
    #
    def test_ppf(self):
        "Try the ppf for the pearson3 distribution"
        quantiles = (0.001, 0.250000001, 0.50, 0.999)
        (g, m, s) = (1., 0., 1.)
        assert_almost_equal(pearson3(g).ppf(quantiles),
                            [-1.7857238, -0.7323399, -0.1639696, 4.5311204])
        assert_almost_equal(pearson3(g, m, s).ppf(quantiles),
                            [-1.7857238, -0.7323399, -0.1639696, 4.5311204])
        qtl = pearson3(g, m, s).ppf(quantiles)
        assert_almost_equal(pearson3(g, m, s).cdf(qtl), quantiles)
        (g, m, s) = (-1., 0., 1.)
        assert_almost_equal(pearson3(g, m, s).ppf(quantiles),
                            [-4.5311204, -0.5547137, 0.1639696, 1.7857238])
        qtl = pearson3(g, m, s).ppf(quantiles)
        assert_almost_equal(pearson3(g, m, s).cdf(qtl), quantiles)
        (g, m, s) = (-1., -1., 0.5)
        assert_almost_equal(pearson3(g, m, s).ppf(quantiles),
                            [-3.2655602, -1.2773569, -0.9180152, -0.1071381], 6)
        qtl = pearson3(g, m, s).ppf(quantiles)
        assert_almost_equal(pearson3(g, m, s).cdf(qtl), quantiles)
        (g, m, s) = (+np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).ppf(quantiles),
                            [-2.000000, -1.877859, -1.263601, 22.986482], 6)
        qtl = pearson3(g, m, s).ppf(quantiles)
        assert_almost_equal(pearson3(g, m, s).cdf(qtl), quantiles)
        (g, m, s) = (-np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).ppf(quantiles),
                            [-22.9864825, -0.5326628, 1.2636014, 1.9999999], 6)
        qtl = pearson3(g, m, s).ppf(quantiles)
        assert_almost_equal(pearson3(g, m, s).cdf(qtl), quantiles)
    #
    def test_sf(self):
        "Try the sf for the pearson3 distribution"
        (g, m, s) = (1., 0., 1.)
        assert_almost_equal(pearson3(g).sf([0., .25, .5, 1.]),
                            1.-pearson3(g).cdf([0., .25, .5, 1.]))
        (g, m, s) = (-1., 0., 1.)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            1.-pearson3(g, m, s).sf([0., .25, .5, 1.]))
        (g, m, s) = (-1., -1., 0.5)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            1.-pearson3(g, m, s).sf([0., .25, .5, 1.]), 6)
        (g, m, s) = (+np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            1.-pearson3(g, m, s).sf([0., .25, .5, 1.]), 6)
        (g, m, s) = (-np.pi, 0., +np.pi)
        assert_almost_equal(pearson3(g, m, s).cdf([0., .25, .5, 1.]),
                            1.-pearson3(g, m, s).sf([0., .25, .5, 1.]), 6)
    #
    def test_isf(self):
        "Try the isf for the pearson3 distribution"
        quantiles = (0.001, 0.250000001, 0.50, 0.999)
        (g, m, s) = (1., 0., 1.)
        qtl = pearson3(g, m, s).isf(quantiles)
        assert_almost_equal(pearson3(g, m, s).sf(qtl), quantiles)
        (g, m, s) = (-1., 0., 1.)
        qtl = pearson3(g, m, s).isf(quantiles)
        assert_almost_equal(pearson3(g, m, s).sf(qtl), quantiles)
        (g, m, s) = (-1., -1., 0.5)
        qtl = pearson3(g, m, s).isf(quantiles)
        assert_almost_equal(pearson3(g, m, s).sf(qtl), quantiles)
        (g, m, s) = (+np.pi, 0., +np.pi)
        qtl = pearson3(g, m, s).isf(quantiles)
        assert_almost_equal(pearson3(g, m, s).sf(qtl), quantiles)
        (g, m, s) = (-np.pi, 0., +np.pi)
        qtl = pearson3(g, m, s).isf(quantiles)
        assert_almost_equal(pearson3(g, m, s).sf(qtl), quantiles)



class TestExtraDiscrete(TestCase):
    #
    def test_basic(self):
        (n, p) = (.3, .15)
        v = np.arange(0, 21)
        vals = dist.nbinom(n, p).cdf(v)
        control = [dist.nbinom(n, p).pmf(np.arange(_+1)).sum() for _ in v]
        assert_almost_equal(vals, control)
    #
    def test_ztnbinom_pmf(self):
        "Test ztnbinom.pmf"
        (n, p) = (.3, .15)
        v = np.arange(1, 21)
        vals = ztnbinom(n, p).pmf(v)
        ctrl = dist.nbinom(n, p).pmf(v)/(1.-dist.nbinom(n, p).pmf(0))
        assert_almost_equal(vals, ctrl)
        #
        vals = ztnbinom(n, p).cdf(v)
        control = [ztnbinom(n, p).pmf(np.arange(1, _+1)).sum() for _ in v]
        assert_almost_equal(vals, control)
        #
        assert_equal(ztnbinom(n, p).ppf(vals), v)
        #
        (n, p) = (1.5, .15)
        v = np.arange(1, 21)
        vals = ztnbinom(n, p).pmf(v)
        ctrl = dist.nbinom(n, p).pmf(v)/(1.-dist.nbinom(n, p).pmf(0))
        assert_almost_equal(vals, ctrl)
        #
        vals = ztnbinom(n, p).cdf(v)
        control = [ztnbinom(n, p).pmf(np.arange(1, _+1)).sum() for _ in v]
        assert_almost_equal(vals, control)
        #
        assert_equal(ztnbinom(n, p).ppf(vals), v)
    #
    def test_logseries_pmf(self):
        "Test logseries.pmf"
        p = 0.15
        v = np.arange(1, 21)
        #
        vals = logseries(p).cdf(v)
        control = [logseries(p).pmf(np.arange(1, _+1)).sum() for _ in v]
        assert_almost_equal(vals, control)
        #
        control = [0.3171178, 0.4677487, 0.5631483, 0.6311206, 0.6827794,
                   0.7236761, 0.7569776, 0.7846595, 0.8080353, 0.8280216,
                   0.8452826, 0.8603140, 0.8734953, 0.8851232, 0.8954332,
                   0.9046156, 0.9128257, 0.9201920, 0.9268217, 0.9328049]
        assert_almost_equal(logseries(0.95).cdf(v), control)
        #
        control = [0.7213475, 0.9016844, 0.9617967, 0.9843388, 0.9933556,
                   0.9971127, 0.9987228, 0.9994273, 0.9997403, 0.9998812,
                   0.9999453, 0.9999746, 0.9999882, 0.9999945, 0.9999974,
                   0.9999988, 0.9999994, 0.9999997, 0.9999999, 0.9999999]
        assert_almost_equal(logseries(0.5).cdf(v), control)





if __name__ == "__main__":
    run_module_suite()

