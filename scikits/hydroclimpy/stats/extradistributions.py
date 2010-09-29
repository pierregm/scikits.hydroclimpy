"""
.. currentmodule:: scikits.hydroclimpy.stats.extradistributions


Conventions
===========

* :math:`F` is the cumulative distribution function (CDF).
* :math:`f` is the probability distribution function (PDF).
* :math:`x = \cfrac{x - \\xi}{\\alpha}` is a reduced variable,
  where :math:`\\xi` is a location parameter and :math:`\\alpha` is a scale
  parameter.


Generalized normal
==================

The generalized normal distribution is defined by its cumulative distribution:

.. math::
   F(x) = \Phi{\left[ -k^{-1} \log\{1 - kx\} \\right]}

where :math:`\Phi` is the normal CDF.

* When :math:`k=0`, this expression reduces to the standard normal distribution
  with location parameter :math:`\\xi` and scale parameter :math:`\\alpha`.



Generalized logistic
====================

The generalized logistic is defined by its CDF:

.. math::
   F(x;k) = \\frac{1}{1 + \left[1 - kx\\right]^{1/k}}

where :math:`k` is a shape parameter.

.. The PDF is given by :math:`f(x) = `
The PPF is given by :math:`x = k^{-1} [1 - \{(1 - F)/F\}^{k}]`.



Kappa
=====

The Kappa distribution is defined by its CDF:

.. math::
   F(x; a, b) = \left[1-h\{1-kx\}^{1/k}\\right]^{1/h}

where :math:`k` and :math:`h` are two scale parameters.

The PPF is given by :math:`x = k^{-1} [1 - \{(1 - F)^{h}/h\}^{k}]`



Wakeby
======

The Wakeby distribution is defined by the transformation:

.. math::
     \cfrac{x - \\xi}{\\alpha} = \\frac{1}{\\beta} \left[1 - (1-U)^{\\beta}\\right] -
                                 \\frac{\gamma}{\delta} \left[1 - (1-U)^{-\delta} \\right]

where :math:`U` is a standard uniform random variable, and where
:math:`\\beta`, :math:`\gamma` and :math:`\delta` are three shape parameters.
That is, the above equation defines the percent point function for the Wakeby 
distribution.

The cumulative distribution function is computed by numerically inverting
the percent point function.
The probability density function is then found by using the following relation
(given on page 46 of Johnson, Kotz, and Balakrishnan):

.. math::
   f(x) = \cfrac{[1 - F(x)]^{\delta+1}}{%
                 \\alpha [1 - F(x)]^{\\beta+\delta} + \gamma}

The parameters :math:`\\beta`, :math:`\gamma` and :math:`\delta` are shape
parameters; :math:`\\xi` is a location parameter and :math:`\\alpha` a scale
parameter.
With three shape parameters, the Wakeby distrbution can model a wide variety
of shapes.

The following restrictions apply to the parameters of this distribution:

* :math:`\\delta < 1.`
* :math:`\\beta + \delta \geq 0`
* either :math:`\\beta + \delta > 0` or :math:`\\beta = \gamma = \delta = 0`
* :math:`\gamma \geq 0`. If :math:`\gamma > 0`, then :math:`\delta > 0`
* :math:`\\alpha + \gamma \geq 0`.

The domain of the Wakeby distribution is 
:math:`[\\xi ; +\inf)` for :math:`\delta \geq 0` and :math:`\gamma > 0` and
:math:`[\\xi ; \\xi + \\alpha (1/\\beta - \gamma/\delta))`
for :math:`\delta < 0` or :math:`\gamma = 0`.



Notes about adding distributions to :mod:`scipy.stats.distributions`.
=====================================================================

In scipy, a distribution (`distrib`) is implemented as a generic class 
(`distrib_gen`) which is a subclass of either :class:`rv_continuous` or
:class:`rv_distrib`.
This `distrib_gen` class is a catch-all to define the basic methods, such
as :meth:`_pdf`, :meth:`_cdf` or :math:`_ppf`.


Attributes
----------

.. attribute:: a

   The lower bound of the validity domain.
   If the lower bound is independent of the input parameters, it can be set
   when creating a new instance of the distribution.
   If the lower bound depends on the parameters, it can be set
   with the :meth:`_argcheck` method.


.. attribute:: b

   The upper bound of the validity domain.
   If the lower bound is independent of the input parameters, it can be set
   when creating a new instance of the distribution.
   If the lower bound depends on the parameters, it can be set
   with the :meth:`_argcheck` method.


Methods
-------

.. method:: _argcheck

   Check the validity of the extra scale parameters.
   By default, the parameters must be strictly positive.
   If this is not the case, a generic condition ``(k==k)`` can be used instead.


"""


import numpy as np
from numpy import exp, floor, log, power, sqrt, where
import numpy.random as mtrand

import scipy.special as special
import scipy.optimize as optimize
import scipy.stats.distributions as dist
from scipy.stats.distributions import argsreduce, valarray

import _lmoments

arr = np.asarray

__all__ = ['gennorm_gen', 'gennorm', 'glogistic_gen', 'glogistic',
           'kappa_gen', 'kappa',
           'logseries_gen', 'logseries',
           'pearson3_gen', 'pearson3',
           'wakeby_gen', 'wakeby',
           'ztnbinom_gen', 'ztnbinom']


##--- Continuous distributions -------------------------------------------------

class kappa_gen(dist.rv_continuous):
    """
    The CDF is given by 

.. math::
   F(x; a, b) = \left[1-h\{1-kx\}^{1/k}\\right]^{1/h}
    """
    def _argcheck(self, k, h):
        k = np.asarray(k)
        h = np.asarray(h)
        # Upper bound
        self.b = where(k <= 0, np.inf, 1. / k)
        # Lower bound
        self.a = where(h > 0,
                       where(k == 0, 0., (1 - h ** (-k)) / k),
                       where(k < 0, 1. / k, -np.inf))
        return (k == k) | (h == h)

    def _cdf(self, x, k, h):
        y = where(k == 0, exp(-x), (1 - k * x) ** (1. / k))
        return where(h == 0, exp(-y), (1. - h * y) ** (1. / h))
    def _pdf(self, x, k, h):
        y = (1 - k * x) ** (1. / k - 1.)
        y *= self._cdf(x, k, h) ** (1. - h)
        return y
    def _ppf(self, q, k, h):
        y = where(h == 0, -log(q), (1. - q ** h) / h)
        y = where(k == 0, -log(y), (1. - y ** k) / k)
        return y

kappa = kappa_gen(name='kappa', shapes='k,h', extradoc="""

""")


class glogistic_gen(dist.rv_continuous):
    """
The CDF is given by

.. math::
   F(x;k) = \\frac{1}{1 + \left[1 - kx\\right]^{1/k}}
    """
    #
    numargs = 1
    #
    def _argcheck(self, k):
        return (k == k)
    def _cdf(self, x, k):
        u = where(k == 0, exp(-x), (1. - k * x) ** (1. / k))
        return 1. / (1. + u)
    def _pdf(self, x, k):
        u = where(k == 0, exp(-x), (1. - k * x) ** (1. / k))
        return u ** (1. - k) / (1. + u) ** 2
    def _ppf(self, q, k):
        F = q / (1. - q)
        return where(k == 0, log(F), (1 - F ** (-k)) / k)

glogistic = glogistic_gen(name='glogistic', shapes='k', extradoc="""

""")


class gennorm_gen(dist.rv_continuous):
    """
    The CDF is given by 

.. math::
   F(x) = \Phi{\left[ -k^{-1} \log\{1 - kx\} \\right]}
    """
    #
    numargs = 1
    #
    def _argcheck(self, k):
        return (k == k)
    #
    def _cdf(self, x, k):
        y = where(k == 0, x, -np.log(1. - k * x) / k)
        return 0.5 * (1 + special.erf(y * np.sqrt(0.5)))
    #
    def _pdf(self, x, k):
        u = where(k == 0, x, -log(1. - k * x) / k)
        return exp(k * u - u * u / 2.) / np.sqrt(2 * np.pi)
    #
    def _ppf(self, q, k):
        u = dist._norm_ppf(q)
        return where(k == 0, u, (1. - exp(-k * u)) / k)

gennorm = gennorm_gen(name='gennorm', shapes='k', extradoc="""

""")


class wakeby_gen(dist.rv_continuous):
    """
The Wakeby distribution is defined by the transformation:
(x-xi)/a = (1/b).[1 - (1-U)^b] - (c/d).[1 - (1-U)^(-d)]

    """
    #
    def _argcheck(self, b, c, d):
        b = np.asarray(b)
        c = np.asarray(c)
        d = np.asarray(d)
        check = where(b + d > 0,
                      where(c == 0, d == 0, True),
                      (b == c) & (c == d) & (d == 0))
        np.putmask(check, c > 0, d > 0)
        np.putmask(check, c < 0, False)
        return check
    #
    def _ppf(self, q, b, c, d):
        z = -np.log(1. - q)
        u = where(b == 0, z, (1. - exp(-b * z)) / b)
        v = where(d == 0, z, (1. - exp(d * z)) / d)
        return u - c * v
    #
    def _cdf(self, x, b, c, d):
        if hasattr(x, '__iter__'):
            cdf = np.array([_lmoments.cdfwak(_, parameters)
                            for (_, parameters) in zip(x, zip(b, c, d))])
        else:
            cdf = _lmoments.cdfwak(x, (b, c, d))
        return cdf
    #
    def _pdf(self, x, b, c, d):
        t = (1. - self._cdf(x, b, c, d))
        f = t ** (d + 1) / (t ** (b + d) + c)
        return f


wakeby = wakeby_gen(name='wakeby', shapes='beta, gamma, delta', extradoc="""

""")



class pearson3_gen(dist.gamma_gen):
    """
The Pearson III is a particular case of the Gamma distribution.

    """
    #
    def _fix_loc_scale(self, args, mu, scale=1):
        N = len(args)
        # Get mu and scale (sigma) from args
        if N == 2 and mu is None:
            mu = args[-1]
        elif N == 3 and scale is None:
            (mu, scale) = args[-2:]
        gamma = args[0]
        # Set the default
        if scale is None:
            scale = 1.
        if mu is None:
            mu = 0.0
        # Transforms (gamma, mu, sigma) into (alpha, xi, beta)
        sigma = scale
        if gamma == 0:
            return ((0.,), mu, sigma)
        xi = mu - 2 * sigma / gamma
        beta = sigma * abs(gamma) / 2.
        beta = sigma * gamma / 2.
        alpha = 4. / gamma ** 2.
        # Returns the resul: the first argument must be a tuple
        return ((alpha,), xi, beta)


    def pdf(self, x, *args, **kwds):
        """
        Probability density function at x of the given RV.

        Parameters
        ----------
        x : array-like
            quantiles
        arg1, arg2, arg3,... : array-like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array-like, optional
            location parameter (default=0)
        scale : array-like, optional
            scale parameter (default=1)

        Returns
        -------
        pdf : array-like
            Probability density function evaluated at x

        """
        (loc, scale) = map(kwds.get, ['loc', 'scale'])
        (args, loc, scale) = self._fix_loc_scale(args, loc, scale)
        (x, loc, scale, shape) = map(arr, (x, loc, scale, args[0]))
        x = (x - loc) * 1.0 / scale
        scale = np.abs(scale)
        # 
        isgamma = (shape > 0) & (scale != 0)
        isnorm = (shape == 0) & (scale > 0)
        ispe3 = (isgamma | isnorm)
        indom = (x > self.a) & (x < self.b)
        valid = ispe3 & indom
        #
        output = np.zeros(np.shape(valid), 'd')
        np.putmask(output, (1 - ispe3) * np.array(indom, bool), self.badvalue)
        (x, shape, scale) = argsreduce(valid, *((x, shape, scale,)))
        np.place(output, (valid & isgamma), self._pdf(x, shape) / scale)
        np.place(output, (valid & isnorm), dist._norm_pdf(x) / scale)
        if output.ndim == 0:
            return output[()]
        return output

    def cdf(self, x, *args, **kwds):
        """
        Cumulative distribution function at x of the given RV.

        Parameters
        ----------
        x : array-like
            quantiles
        arg1, arg2, arg3,... : array-like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array-like, optional
            location parameter (default=0)
        scale : array-like, optional
            scale parameter (default=1)

        Returns
        -------
        cdf : array-like
            Cumulative distribution function evaluated at x

        """
        (loc, scale) = map(kwds.get, ['loc', 'scale'])
        (args, loc, scale) = self._fix_loc_scale(args, loc, scale)
        (x, loc, scale, shape) = map(arr, (x, loc, scale, args[0]))
        x = (x - loc) * 1.0 / scale
        # 
        isgamma = (shape > 0) & (scale != 0)
        isnorm = (shape == 0) & (scale > 0)
        ispe3 = (isnorm | isgamma)
        indomain = (x > self.a) & (x < self.b)
        toolarge = (x >= self.b)
        valid = ispe3 & indomain
        output = np.zeros(np.shape(valid), 'd')
        np.place(output, (1 - ispe3) * (indomain == indomain), self.badvalue)
        np.place(output, toolarge, 1.0)
        if any(valid):  #call only if at least 1 entry
            (x, shape) = argsreduce(valid, *((x,) + (shape,)))
            vals = self._cdf(x, shape)
            np.place(output, (valid & isgamma),
                     np.where(scale > 0, vals, 1. - vals))
            np.place(output, (valid & isnorm), dist._norm_cdf(x))
        if output.ndim == 0:
            return output[()]
        return output


    def sf(self, x, *args, **kwds):
        """
        Survival function (1-cdf) at x of the given RV.

        Parameters
        ----------
        x : array-like
            quantiles
        arg1, arg2, arg3,... : array-like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array-like, optional
            location parameter (default=0)
        scale : array-like, optional
            scale parameter (default=1)

        Returns
        -------
        sf : array-like
            Survival function evaluated at x

        """
        (loc, scale) = map(kwds.get, ['loc', 'scale'])
        (args, loc, scale) = self._fix_loc_scale(args, loc, scale)
        (x, loc, scale, shape) = map(arr, (x, loc, scale, args[0]))
        x = (x - loc) * 1.0 / scale
        #
        isgamma = (shape > 0) & (scale != 0)
        isnorm = (shape == 0) & (scale > 0)
        ispe3 = (isgamma | isnorm)
        indom = (x > self.a) & (x < self.b)
        toosmall = ispe3 & (x <= self.a)
        valid = ispe3 & indom
        output = np.zeros(np.shape(valid), 'd')
        np.place(output, (1 - ispe3) * (indom == indom), self.badvalue)
        np.place(output, toosmall, 1.0)
        (x, shape) = argsreduce(valid, *((x, shape,)))
        vals = self._sf(x, shape)
        np.place(output, (valid & isgamma), np.where(scale > 0, vals, 1. - vals))
        np.place(output, (valid & isnorm), 1. - dist._norm_cdf(x))
        if output.ndim == 0:
            return output[()]
        return output

    def ppf(self, q, *args, **kwds):
        """
        Percent point function (inverse of cdf) at q of the given RV.

        Parameters
        ----------
        q : array-like
            lower tail probability
        arg1, arg2, arg3,... : array-like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array-like, optional
            location parameter (default=0)
        scale : array-like, optional
            scale parameter (default=1)

        Returns
        -------
        x : array-like
            quantile corresponding to the lower tail probability q.

        """
        (loc, scale) = map(kwds.get, ['loc', 'scale'])
        (args, loc, scale) = self._fix_loc_scale(args, loc, scale)
        (q, loc, scale, shape) = map(arr, (q, loc, scale, args[0]))
        q = np.where(scale > 0, q, 1. - q)
        #
        isgamma = (shape > 0) & (scale != 0)
        isnorm = (shape == 0) & (scale > 0)
        ispe3 = (isgamma | isnorm)
        indom = (q > 0) & (q < 1)
        islarge = (q == 1) & ispe3
        valid = ispe3 & indom
        output = valarray(np.shape(valid), value=self.a * np.abs(scale) + loc)
        np.place(output, (1 - ispe3) + (1 - indom) * (q != 0.0), self.badvalue)
        np.place(output, islarge, self.b * np.abs(scale) + loc)
        if any(valid):  #call only if at least 1 entry
            (q, shape, scale, loc) = argsreduce(valid, *((q, shape, scale, loc)))
            np.place(output, (valid & isgamma), self._ppf(q, shape) * scale + loc)
            np.place(output, (valid & isnorm), dist._norm_ppf(q) * scale + loc)
        if output.ndim == 0:
            return output[()]
        return output

    def isf(self, q, *args, **kwds):
        """
        Inverse survival function at q of the given RV.

        Parameters
        ----------
        q : array-like
            upper tail probability
        arg1, arg2, arg3,... : array-like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array-like, optional
            location parameter (default=0)
        scale : array-like, optional
            scale parameter (default=1)

        Returns
        -------
        x : array-like
            quantile corresponding to the upper tail probability q.

        """
        (loc, scale) = map(kwds.get, ['loc', 'scale'])
        (args, loc, scale) = self._fix_loc_scale(args, loc, scale)
        (q, loc, scale, shape) = map(arr, (q, loc, scale, args[0]))
        q = np.where(scale > 0, q, 1. - q)
        #
        isgamma = (shape > 0) & (scale != 0)
        isnorm = (shape == 0) & (scale > 0)
        ispe3 = (isgamma | isnorm)
        indom = (q > 0) & (q < 1)
        islarge = (q == 1) & ispe3
        valid = ispe3 & indom
        output = valarray(np.shape(valid), value=self.b)
        #place(output,(1-cond0)*(cond1==cond1), self.badvalue)
        np.place(output, (1 - ispe3) * (indom == indom) + (1 - indom) * (q != 0.0),
                 self.badvalue)
        np.place(output, islarge, self.a)
        if np.any(valid):  #call only if at least 1 entry
            goodargs = argsreduce(valid, *((q,) + args + (scale, loc)))  #PB replace 1-q by q
            (q, shape, scale, loc) = argsreduce(valid, *((q, shape, scale, loc)))
            np.place(output, (valid & isgamma), self._isf(q, shape) * scale + loc)
            np.place(output, (valid & isnorm), -dist._norm_ppf(q))
        if output.ndim == 0:
            return output[()]
        return output

pearson3 = pearson3_gen(name='Pearson III', shapes='gamma', extradoc="""

""")



##--- Discrete distributions ---------------------------------------------------


# Negative binomial
class ztnbinom_gen(dist.rv_discrete):
    def _rvs(self, n, pr):
        return mtrand.negative_binomial(n, pr, self._size)
    def _argcheck(self, n, pr):
        return (n >= 0) & (pr >= 0) & (pr <= 1)
    def _pmf(self, x, n, pr):
        coeff = exp(special.gammaln(n + x) - special.gammaln(x + 1) - special.gammaln(n))
        return coeff * power(pr, n) * power(1 - pr, x) / (1 - power(pr, n))
    def _cdf(self, x, n, pr):
        k = floor(x)
        vals = (special.betainc(n, k + 1, pr) - special.betainc(n, 1, pr))
        return vals / (1. - pr ** n)
    def _sf_skip(self, x, n, pr):
        #skip because special.nbdtrc doesn't work for 0<n<1
        k = floor(x)
        return special.nbdtrc(k, n, pr)
    def _stats(self, n, pr):
        qr = 1. - pr
        qrn = 1. - power(pr, n)
        nqr = n * qr
        mu = nqr / (pr * qrn)
        mc2 = mu * (1. + nqr) / pr
        mc3 = (nqr + n * (3 * n + 1) * qr ** 2 + nqr ** 3) / (pr ** 3 * qrn)
        mc4 = (nqr * (6 - 6 * pr + pr ** 2) + (11 - 4 * pr) * nqr ** 2 + 6 * nqr ** 3 + nqr ** 4)
        mc4 /= (pr ** 4 * qrn)
        #
        var = m2 = mc2 - mu * mu
        m3 = mc3 - 3 * mc2 * mu + 2 * mu ** 3
        m4 = mc4 - 4 * mc3 * mu + 6 * mc2 * mu ** 2 - 3 * mu ** 4
        #
        g1 = m3 / m2 ** (3 / 2.)
        g2 = m4 / m2 ** 2
        return mu, var, g1, g2
ztnbinom = ztnbinom_gen(a=1, shapes="n,pr", name='ztnbinom',
                        longname="A zero-truncated negative binomial",
                        extradoc="""

Zero-truncated negative binomial distribution

nbinom.pmf(k,n,p) = choose(k+n-1,n-1) * p**n * (1-p)**k / (1-p**n)
for k >= 1.
"""
                    )



class logseries_gen(dist.rv_discrete):
    def _argcheck(self, pr):
        return (pr > 0) & (pr < 1)
    def _pmf(self, x, pr):
        a = -1. / log(1. - pr)
        return a * pr ** x / x
    def _stats(self, p):
        a = -1. / log(1. - p)
        mu = a * p / (1. - p)
        var = mu * (1. / (1. - p) - mu)
        ap = a * p
        g1 = (1 + p - 3 * ap + 2 * ap * ap) / sqrt(ap * (1 - ap) ** 3)
        g2 = (1 + 4 * p + p ** 2 - 4 * ap * (1 + p) + (6 - 3 * ap) * ap * ap) / (ap * (1 - ap) ** 2)
        return (mu, var, g1, g2)
    def fit(self, distrib):
        (mu, var) = (distrib.mean(), distrib.var())
        pr_0 = 1. - mu / (mu * mu + var)
        def mle(pr, x):
            F = x * np.log(pr) - np.log(x) - np.log(-np.log(1. - pr))
            return - np.sum(F, axis=0)
        return optimize.fmin(mle, pr_0, args=(distrib.ravel(),), disp=0)
logseries = logseries_gen(a=1, name='logseries', shapes="pr", extradoc="""

Log-series distribution

   logseries.pmf(k,p) = -1/log(1-p) * p**k / k
   for k in {1,...,n}
""")





