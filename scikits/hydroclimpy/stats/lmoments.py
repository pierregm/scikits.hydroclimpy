"""
L-moments
=========

Let :math:`F` be the Cumulative Distribution Function (CDF) of a continuous 
random variable :math:`x`.
Hosking (1990) defined **L-moments** to be the quantities

.. math::
   \lambda_r = E[x P_{r-1}^{*}\{F(x)\}]

where :math:`P_{r}^{*}(\cdot)` is the :math:`r^\mathrm{th}` shifted Legendre
polynomial.

The first L-moment :math:`\lambda_1` is called L-location and is the mean of the
distribution.
The second L-moment :math:`\lambda_2` is called the L-scale, a value always
positive.


**L-moments ratios** (hereafter abridged as *L-ratios*) are defined as
the quantities :math:`\\tau_r = \lambda_r/\lambda_2` for any integer
:math:`r \geq 3`.
In particular, the quantities :math:`\\tau_3` and :math:`\\tau_4` are 
the L-skewness and L-kurtosis respectively.
Another L-moment ratio, the so-called `coefficient of L-variation`,
is the variable :math:`\\tau = \lambda_2/\lambda_1`.
L-ratios satisfy :math:`|\\tau_r| < 1` for all :math:`r \geq 3`.
Individual :math:`\\tau_r` values can have tighter bounds: for example,
the L-kurtosis :math:`\\tau_4` satisfies
:math:`1/4 (5\\tau_3^2 - 1) \leq \\tau_4 < 1`.

In many applications, a combination of the first two L-moments and L-ratios 
:math:`\\tau_r` for :math:`r \geq 3` is more useful than only L-moments.
For convenience, such a combination is refered to as **L-statistics**.



Probability Weighted Moments
============================

Given a random variable :math:`x` with a cumulative distribution function
:math:`F`, the Probability Weighted Moments (PWMs) are defined to be

.. math::
   M(p, r, s) = E[x^p \{F(x)\}^r \{1-F(x)\}^{s}]

Particularly useful cases are the probability weighted moments:

.. math::
   \\alpha_r = M(1,0,r) = E[x^p \{1-F(x)\}^{s}]

.. math::
   \\beta_r = M(1,r,0) = E[x^p \{F(x)\}^r]


L-moments and probability weighted moments are related by

.. math::
   \lambda_{r+1} = \sum_{k=0}^{n}{p_{r,k}^{*} \\beta_k}

where

.. math::
   p_{r,k}^{*} = (-1)^{r-k} \\binom{r}{k} \\binom{r+k}{k}

The foregoing quantities are for a probability distribution, but in practice can
often be estimated from a finite sample.
Let :math:`x_{(1:n)} \leq x_{(2:n)} \leq \cdots \leq x_{(n:n)}` bet the ordered
sample.
Let

.. math::
   a_r = n^{-1} \sum_{j=1}^{n}{\cfrac{(n-j)(n-j-1) \ldots (n-j-r+1)}{
                                      (n-1)(n-2) \ldots (n-r)} x_{(j:n)}}

.. math::
   b_r = n^{-1} \sum_{j=1}^{n}{\cfrac{(n-j)(n-j-1) \ldots (n-j-r+1)}{
                                      (n-1)(n-2) \ldots (n-r)} x_{(j:n)}}

.. math::
   \\ell_{r+1} = \sum_{k=0}^{n}{p_{r,k}^{*} b_{k}}

Then, :math:`a_r`, :math:`b_r` and :math:`\\ell_r` are unbiased estimators of
:math:`\\alpha_r`, :math:`\\beta_r` and :math:`\lambda_r` respectively.
The estimator :math:`t_r=\\ell_r/\\ell_2` of :math:`\\tau_r` is consistent but
not unbiased.
Landwehr et al. (1979a, 1979b) also considered 'plotting-positions' estimators
of probability weighted moments.
Let :math:`p_{(j:n)}` be a plotting position, i.e. a distribution free parameter
of :math:`F(x_{(j:n)})`.
Then, :math:`\\alpha_r`, :math:`\\beta_r`, :math:`\\lambda_r` and :math:`\\tau_r`
are estimated respectively by

.. math::
   \\tilde{\\alpha}_{r} = n^{-1} \sum_{j=1}^{n}{(1-p_{(j:n)})^{r} x_{(j:n)}}

.. math::
   \\tilde{\\beta}_{r} = n^{-1} \sum_{j=1}^{n}{p_{(j:n)}^{r} x_{(j:n)}}

.. math::
   \\tilde{\\lambda}_{r} = n^{-1} \sum_{j=1}^{n}{P_{r-1}^{*}(p_{(j:n)}) x_{(j:n)}}

.. math::
   \\tilde{\\tau}_r = \\tilde{\lambda}_{r} / \\tilde{\lambda}_{2}

Plotting position estimators are consistent but can have large bias and cannot
be generally recommended (Hosking & Wallis 1995; Hosking & Wallis 1997).
Hosking (1990) and Hosking & Wallis (1997) give expositions of the theory of
L-moments and L-moment ratios.
Hosking & Wallis give, for many distributions in common use, expressions for the
L-moments of the distributions and algorithms for estimating the parameters of
the distributions by equating sample and population L-moments (the 'method of
L-moments').



Errr ???
========


From (2) and (3), :math:`\\ell_r` is a linear combination of the ordered sample
values :math:`x_{(1:n)}, \ldots, x_{(n:n)}`, and we can write:

.. math::
   \ell_r = n^{-1} \sum_{j=1}^{n}{w_{(j:n)}^{(r)} x_{(j:n)}}

The weights :math:`w_{(j:n)}^{(r)}` are related to the 'discrete Legendre
polynomials' defined by Neuman and Schonbach (1974).
In Neuman and Schonbach's notation, the weight :math:`w_{(j:n)}^{(r)}` 
is the discrete Legendre polynomial :math:`(-1)^r P_{r-1}(j-1, n-1)`.
The weights can be computed as :math:`w_{(j:n)}^{(1)}=1`, 
:math:`w_{(j:n)}^{(2)}=2(j-1)/(n-1)-1`, and, for :math:`r \geq 3`,

.. math::
   (r - 1)(n - r + 1) w_{(j:n)}^{(r)} = (2r - 3)(2j - n - 1) w_{(j:n)}^{(r-1)} -
                                        (r - 2)(n + r - 2) w_{(j:n)}^{(r-2)} 

This recurrence relation is equivalent to equation (9)
of Neuman and Schonbach (1974).

Routine SAMLMU evaluates sample L-moments using this recurrence,
and achieves additional
efficiency by making use of the symmetry in the weights, w(r)
:math:`w_{(j:n)}^{(r)} = (-1)^{r-1} w_{(n+1-j:n)}^{(r)}`.
This routine is generally equivalent to SAMLMR, but is particularly suitable
when highorder L-moment ratios are to be calculated. Routine SAMLMU is accurate
to 8 decimal places for about the first 50 L-moment ratios, whereas SAMLMR
achieves similar accuracy for only the first 12 L-moment ratios.
The routines require approximately the same amount of computing time, but SAMLMU
can be a little slower than SAMLMR when NMOM is large or odd.



Implementation
==============

Three new methods are introduced for the
:class:`~scipy.stats.distributions.rv_continuous` class:


.. automethod:: scipy.stats.distributions.rv_continuous.lmoments
.. automethod:: scipy.stats.distributions.rv_continuous.lratios
.. automethod:: scipy.stats.distributions.rv_continuous.lstats


.. autofunction:: lmoments
.. autofunction:: PWM



References
==========

.. [Hosking (1990)] Hosking, J.R.M. (1990).
   L-moments analysis and estimation of distributions using linear combinations
   of order statistics:
   Journal of the Royal Statistical Society, Series B, vol. 52, p. 105--124.

.. [Hosking (1996)] Hosking, J.R.M. (1996).
   FORTRAN routines for use with the method of L-moments: Version 3,
   IBM Research Report RC20525,
   T.J. Watson Research Center, Yorktown Heights, New York.

.. [Hosking and Wallis (1997)] Hosking, J.R.M. and Wallis, J.R. (1997).
   Regional frequency analysis: An approach based on Lmoments:
   Cambridge University Press.


"""



import warnings
import numpy as np

import _lmoments
import scipy.special as special
import scipy.stats.distributions as dist

import extradistributions as extradist
from extradistributions import gennorm_gen, glogistic_gen, kappa_gen,\
                               pearson3_gen, wakeby_gen,\
                               gennorm, glogistic, kappa, pearson3, wakeby



def powseries(a, coeffs):
    """
    Performs a linear combination of the coefficients ``coeffs`` 
    with powers of ``a``

    Parameters
    ----------
    a : array
        Input array
    coeffs : sequence
        Sequence of coefficients, with entry `i` corresponding to exponent `i`.
    """
    if len(coeffs) == 1.:
        n = np.size(a)
        if n > 1:
            return [coeffs[0]] * n
        return coeffs[0]
    output = 0
    a = np.asarray(a)
    for c in coeffs[-1:0:-1]:
        output += c
        output *= a
    output += coeffs[0]
    return output



def lmoments(x, nmom=4, a=0, b=0, mode='stats'):
    """
    Sample L-moments of a data array.

    Parameters
    ----------
    x : array
        Input data array
    nmom : integer, optional
        Number of L-moments to return
    a : float, optional
        First parameter for computing plotting position
    b : float, optional
        Second parameter for computing plotting position
    mode : {'stats', 'moments', 'ratios'}, optional
        Select the kind of information to output.
        
        * ``moments`` , returns the sample L-moments :math:`\\ell_r`.
        * ``ratios`` : returns the sample L-ratios :math:`t, t_2, t_3`
        * ``stats`` (default) : returns a combination of L-moments and L-ratios
          (:math:`\\ell_1, \\ell_2, t_3...)

    Returns
    -------
    xmom : ndarray
        Information.

    Notes
    -----
    For unbiased estimates (of the lambda's) set :math:`a=b=0`.
    Otherwise, plotting-position estimators are used, based on the plotting
    position :math:`(j+a)/(n+b)` for the `j`'th smallest of `n` observations.
    For example, :math:`a=-0.35` and :math:`b=0.0` yields the estimators
    recommended by Hosking et al. (1985, Technometrics) for the GEV distribution.
    """
    x = np.sort(np.asarray(x))
    if (a == 0) and (b == 0):
        (xmom, ierr) = _lmoments.samlmu(x, nmom)
    else:
        (xmom, ierr) = _lmoments.samlmr(x, nmom, a, b)
    if ierr == -1:
        raise ValueError("Invalid parameter nmom")
    elif ierr == -2:
        raise ValueError("Invalid plotting-position parameters")
    elif ierr == -3:
        warnings.warn("All data input are equal")
    xmom = np.asarray(xmom)
    if mode[0].lower() == 'm':
        if nmom >= 2:
            xmom[2:] /= xmom[1]
    elif mode[1].lower() == 'r':
        if nmom >= 1:
            xmom[1] /= xmom[0]
    return xmom



def PWM(x, nmom=4, a=0, b=0, kind=2):
    """
    Probability weighted moments of a data array.

    Parameters
    ----------
    x : sequence
        Input data array.
    nmom : integer, optional
        Number of L-moments to return
    a : float, optional
        First parameter of plotting position.
    b : float, optional
        Second parameter of plotting positions.
    kind : integer, optional
        Specifies which kind of PWMs are to be found.
        
        * 1: :math:`\\alpha_r = E[X.\{1-F(X)\}^r]`
        * 2: :math:`\\beta_r = E[X.\{F(X)\}^r]`

    Returns
    -------
    xmom : array of length nmom. 
        Contains the sample probability weighted moments.

    Notes
    -----
    * For unbiased estimates set :math:`a=b=0`. Otherwise, plotting-position
      estimators are used, based on the plotting position :math:`(j+a)/(n+b)`  
      for the :math:`j`'th smallest of N observations.
      For example, :math:`a=-0.35` and :math:`b=0.0` yields the estimators
      recommended by Hosking et al. (1985, Technometrics)
      for the GEV distribution.
      """
    x = np.sort(np.asarray(x))
    (pwm, ierr) = _lmoments.sampwm(x, nmom, a, b, kind)
    if ierr == -1:
        errmsg = "Invalid parameter nmom: should be < min(20, %i), got %i"
        raise ValueError(errmsg % (len(x), nmom))
    elif ierr == -2:
        raise ValueError("Invalid plotting-position parameters")
    elif ierr == -4:
        raise ValueError("Invalid parameter kind: should be 1 or 2.")
    return pwm


#-------------------------------------------------------------------------------
# Additional distributions
#------------------------------------------------------------------------------
# Parameters to L-moments

# moment from definition
def _lmomg(self, m, *args):
    "Compute the mth L-moment with a Legendre polynomial."
    P_r = special.sh_legendre(m-1)
    func = lambda x : self._ppf(x, *args) * P_r(x)
    return integrate.quad(func, 0, 1)[0]
dist.rv_continuous._lmomg = _lmomg

def _lmomg_frozen(self, nmom):
    return self.dist._lmomg(nmom, *self.args, **self.kwds)
dist.rv_frozen._lmomg = _lmomg_frozen



def _lstats(self, nmom, *args, **kwargs):
    """
    Calculate L-statistics from the parameters of the distribution.

    Parameters
    ----------
    nmom : int
        Number of L-moments to return.
    arg1, arg2, arg3,... : array-like
        The shape parameter(s) for the distribution (see docstring of the
        instance object for more information)
    loc : array-like, optional
        Location parameter (default=0).
    scale : array-like, optional
        Scale parameter (default=1).

    Returns
    -------
    xmom : sequence
        A sequence of ``nmom`` L-statistics.
        The first two items are :math:`\lambda_1` and :math:`\lambda_2` .
        The following items are :math:`\\tau_r` for :math:`r \geq 3` .
    """
    lmoms = np.array([self._lmomg(_, *args, **kwargs)
                      for _ in range(1, nmom+1)])
    if nmom > 2:
        lmoms[2:] /= lmoms[1]
    return lmoms
dist.rv_continuous.lstats = _lstats

def _lstats_frozen(self, nmom):
    return self.dist.lstats(nmom, *self.args, **self.kwds)
dist.rv_frozen.lstats = _lstats_frozen



def _lmoments_(self, nmom, *args, **kwargs):
    """
    Calculate L-moments from the parameters of the distribution.

    Parameters
    ----------
    nmom : int
        Number of L-moments to return.
    arg1, arg2, arg3,... : array-like
        The shape parameter(s) for the distribution (see docstring of the
        instance object for more information)
    loc : array-like, optional
        Location parameter (default=0).
    scale : array-like, optional
        Scale parameter (default=1).

    Returns
    -------
    xmom : sequence
        A sequence of `nmom` L-moments :math:`\lambda_r`
        for :math:`1 \geq r \geq nmom`.
    """
    lstats = self.lstats(nmom, *args, **kwargs)
    if nmom > 2:
        lstats[2:] *= lstats[1]
    return lstats
dist.rv_continuous.lmoments = _lmoments_

def _lmoments_frozen(self, nmom):
    return self.dist.lmoments(nmom, *self.args, **self.kwds)
dist.rv_frozen.lmoments = _lmoments_frozen



def _lratios(self, nmom, *args, **kwargs):
    """
    Calculate L-moment ratios from the parameters of the distribution.

    Parameters
    ----------
    nmom : int
        Number of L-moments to return.
    arg1, arg2, arg3,... : array-like
        The shape parameter(s) for the distribution (see docstring of the
        instance object for more information)
    loc : array-like, optional
        Location parameter (default=0).
    scale : array-like, optional
        Scale parameter (default=1).

    Returns
    -------
    xmom : sequence
        A sequence of :math`nmom-1` L-ratios.
        The first item is the coefficient of L-variation :math:`\\tau`.
        The following items are :math:`\\tau_r` for :math:`r \geq 3`.
    """
    lstats = self.lstats(nmom, *args, **kwargs)
    if nmom > 1:
        lstats[1] /= lstats[0]
    return lstats[1:]
dist.rv_continuous.lratios = _lratios

def _lratios_frozen(self, nmom):
    return self.dist.lratios(nmom, *self.args, **self.kwds)
dist.rv_frozen.lratios = _lratios_frozen




def _list_parameters(self, *args, **kwargs):
    "reorganize the parameters of a distribution as (loc, scale, shapes)"
    (loc, scale) = (kwargs.get('loc'), kwargs.get('scale'))
    (args_, loc, scale) = self._fix_loc_scale(args, loc, scale)
    parameters = [loc, scale] + list(args_)
    return tuple(parameters)
dist.rv_generic._list_parameters = _list_parameters



class _lstats_direct(object):
    """
    Directly calculates L-statistics from the parameters of the distribution.

    Parameters
    ----------
    nmom : int
        Number of L-statistics to return.
    loc : float, optional
        Location parameter.
    scale : float, optional
        Scale parameter.

    Returns
    -------
    xmom : sequence
        A sequence of `nmom` L-statistics.
    """
    #
    def __init__(self, func, dist=None):
        self.func = func
        self.dist = dist
    #
    def __get__(self, obj, objtype=None):
        self.dist = obj
        return self
    #
    def __call__(self, nmom, *args, **kwargs):
        parameters = self.dist._list_parameters(*args, **kwargs)
        (xmom, ierr) = self.func(tuple(parameters), nmom)
        if ierr == -1:
            errmsg = "Invalid number of moments."
            raise ValueError(errmsg % nmom)
        elif ierr == -2:
            errmsg = "Invalid input parameters"
            raise ValueError(errmsg)
        return xmom

dist.expon_gen.lstats = _lstats_direct(_lmoments.lmrexp)
dist.gamma_gen.lstats = _lstats_direct(_lmoments.lmrgam)
dist.genextreme_gen.lstats = _lstats_direct(_lmoments.lmrgev)
extradist.glogistic_gen.lstats = _lstats_direct(_lmoments.lmrglo)
extradist.gennorm_gen.lstats = _lstats_direct(_lmoments.lmrgno)
dist.genpareto_gen.lstats = _lstats_direct(_lmoments.lmrgpa)
dist.gumbel_r_gen.lstats = _lstats_direct(_lmoments.lmrgum)
extradist.kappa_gen.lstats = _lstats_direct(_lmoments.lmrkap)
dist.norm_gen.lstats = _lstats_direct(_lmoments.lmrnor)
extradist.wakeby_gen.lstats = _lstats_direct(_lmoments.lmrwak)


def _lstats_pe3(self, nmom, *args, **kwargs):
    """
    Calculate L-moments from the parameters of the distribution.

    Parameters
    ----------
    nmom : int
        Number of L-moments to return.
    loc : float, optional
        Location parameter.
    scale : float, optional
        Scale parameter.

    Returns
    -------
    xmom : sequence
        A sequence of `nmom` L-moments.
    """
    # Check the arguments
    (mu, sigma) = (kwargs.get('loc', 0), kwargs.get('scale', 1))
    N = len(args)
    if N == 1:
        gamma = args[0]
    elif N == 2:
        (gamma, mu) = args
    elif N >= 3:
        (gamma, mu, sigma) = args[:3]
    (xmom, ierr) = _lmoments.lmrpe3((mu, sigma, gamma), nmom)
    if ierr == -1:
        errmsg = "Invalid number of moments."
        raise ValueError(errmsg % nmom)
    elif ierr == -2:
        errmsg = "Invalid input parameters"
        raise ValueError(errmsg)
    return xmom
extradist.pearson3_gen.lstats = _lstats_pe3


#------------------------------------------------------------------------------
# L-moments to parameters


def _lmparams_generic(self, lmoms):
    errmsg = "Parameter estimation from L-statistics is not available "\
             "for the current distribution."
    raise NotImplementedError(errmsg)
dist.rv_continuous.lmparams = _lmparams_generic

def _lmparams_frozen(self, lmoms):
    return self.dist.lmparams(lmoms)
dist.rv_frozen.lmparams = _lmparams_frozen


def _validate_lstats(nmom, lstats):
    " Private function: check the first nmom L-statistics."
    errmsg = "There should be at least %s L-moments to estimate "\
             "the parameters of the distribution."
    try:
        n = len(lstats)
    except TypeError:
        raise ValueError(errmsg % nmom)
    # Make sure we have enough L-moments
    if n < nmom:
        raise ValueError(errmsg % n)
    lstats = lstats[:nmom]
    if lstats[1] <= 0:
        errmsg = "Invalid second L-moment: must be strictly positive."
        raise ValueError(errmsg)
    if nmom >= 3:
        for i in range(2, nmom):
            if abs(lstats[i]) >= 1.:
                errmsg = "Invalid %s L-moment : should be in ]-1;1["
                raise ValueError(errmsg % (i+1))
        if nmom >= 4:
            lbound = (5*lstats[2]**2 - 1)/4.
            if lstats[3] < lbound:
                errmsg = "Invalid 4th L-moment: should be larger than %s"
                raise ValueError(errmsg % lbound)
    return lstats


def _lmparams_expon(lstats):
    "(needs to be overwritten)"
    (loc, scale) = _validate_lstats(2, lstats)
    scale *= 2
    loc -= scale
    return (loc, scale)
dist.expon_gen.lmparams = staticmethod(_lmparams_expon)


def _lmparams_gamma(lstats):
    "(needs to be overwritten)"
    (lmom_1, lmom_2) = _validate_lstats(2, lstats)
    if lmom_1 <= lmom_2:
        errmsg = "Invalid first L-moment: "\
                 "must be larger than the second L-moment."
        raise ValueError(errmsg)
    #
    (a1, a2, a3) = (-0.3080, -0.05812, 0.01765)
    (b1, b2, b3, b4) = (0.7213, -0.5947, -2.1817, 1.2113)
    cv = lmom_2*1./lmom_1
    if cv < 0.5:
        t = np.pi * cv * cv
        alpha = (1 + a1*t)/(t*(1. + (a2 + a3*t)*t))
    else:
        t = 1. - cv
        alpha = t*(b1 + b2*t)/(1. + (b3 + b4*t)*t)
    scale = lmom_1/alpha
    return (0, scale, alpha)
dist.gamma_gen.lmparams = staticmethod(_lmparams_gamma)



def _lmparams_genextreme(lstats):
    "(needs to be overwritten)"
    (small, eps, maxit) = (1e-5, 1e-6, 20)
    (lmom_1, lmom_2, tau_3) = _validate_lstats(3, lstats)
    #
    if tau_3 > 0:
        (c0, c1, c2, c3) = (-1., 1.59921491, -0.48832213,  0.01573152)
        (d0, d1, d2) = (1., -0.64363929,  0.08985247)
        z = (1. - tau_3)
        g = (c0 + z*(c1 + z*(c2 + z*c3)))/(d0 + z*(d1 + z*d2))
        if abs(g) < small:
            (scale, shape) = (lmom_2/np.log(2.), 0)
            loc = lmom_1 - dist._EULER * scale
            return (loc, scale, shape)
    else:
        (a0, a1, a2, a3, a4) = ( 0.28377530, -1.21096399, -2.50728214,
                                -1.13455566, -0.07138022)
        (b0, b1, b2, b3) = (1. , 2.06189696, 1.31912239, 0.25077104)
        g =  (a0 + tau_3*(a1 + tau_3*(a2 + tau_3*(a3 + tau_3*a4))))
        g /= (b0 + tau_3*(b1 + tau_3*(b2 + tau_3*b3)))
        if tau_3 < -0.8:
            if (tau_3 <= -0.97):
                g = 1. - np.log2(1. + tau_3)
            to = (tau_3 + 3)/2.
            for i in range(maxit):
                x2 = 2**(-g)
                x3 = 3**(-g)
                xx2 = (1. - x2)
                xx3 = (1. - x3)
                t = xx3/xx2
                d = (xx2 * x3 * np.log(3.) - xx3 * x2 * np.log(2)) / (xx2 * xx2)
                gold = g
                g -= (t - to)/d
                if (abs(g - gold) <= eps*g):
                    break
    shape = g
    gam = np.exp(special.gammaln(1.+g))
    scale = lmom_2 * g/(gam * (1. - 2**(-g)))
    loc = lmom_1 - scale*(1. - gam)/g
    return (loc, scale, shape)
dist.genextreme_gen.lmparams = staticmethod(_lmparams_genextreme)



def _lmparams_glogistic(lstats):
    "(needs to be overwritten)"
    (lmom_1, lmom_2, tau_3) = _validate_lstats(3, lstats)
    g = -tau_3
    #
    if (abs(g) <= 1e-6):
        return (lmom_1, lmom_2, 0.)
    #
    gg = g * np.pi/np.sin(g*np.pi)
    a = lmom_2/gg
    return (lmom_1 - a*(1.-gg)/g, a, g)
extradist.glogistic_gen.lmparams = staticmethod(_lmparams_glogistic)



def _lmparams_gennorm(lstats):
    "(needs to be overwritten)"
    small = 1e-8
    (lmom_1, lmom_2, tau_3) = _validate_lstats(3, lstats)
    #
    if abs(tau_3) >= 0.95:
        errmsg = "Invalid third L-moment : must be smaller than 0.95"
        raise ValueError(errmsg)
    if abs(tau_3) <= small:
        (loc, scale, shape) = (lmom_1, lmom_2*np.sqrt(np.pi), 0.)
    else:
        tt = tau_3 * tau_3
        (a0, a1, a2, a3) = (np.sqrt(np.pi/3)*2., -3.6544371, 1.8396733,
                            -0.20360244)
        (b0, b1, b2, b3) = (1., -2.0182173, 1.2420401, -0.21741801)
        g = -tau_3 * (a0 + tt*(a1 + tt*(a2+ tt*a3)))
        g /= (b0 + tt*(b1 + tt*(b2 + tt*b3)))
        #g = -tau_3*powseries(tt, a)/powseries(tt, b)
        e = np.exp(0.5*g*g)
        scale = lmom_2 * g / (e * special.erf(0.5*g))
        loc = lmom_1 + scale * (e - 1.)/g
        shape = g
    return (loc, scale, shape)
extradist.gennorm_gen.lmparams = staticmethod(_lmparams_gennorm)



def _lmparams_genpareto(lstats):
    "(needs to be overwritten)"
    (lmom_1, lmom_2, tau_3) = _validate_lstats(3, lstats)
    shape = (1. - 3.*tau_3)/(1. + tau_3)
    scale = (1. + shape)*(2. + shape)*lmom_2
    loc = lmom_1 - scale/(1.+shape)
    return (loc, scale, -shape)
dist.genpareto_gen.lmparams = staticmethod(_lmparams_genpareto)



def _lmparams_gumbel_r(lstats):
    "(needs to be overwritten)"
    (lmom_1, lmom_2) = _validate_lstats(2, lstats)
    scale = lmom_2 / np.log(2.)
    loc = lmom_1 - dist._EULER * scale
    return (loc, scale)
dist.gumbel_r_gen.lmparams = staticmethod(_lmparams_gumbel_r)



def _lmparams_kappa(lstats):
    "(needs to be overwritten)"
    (l1, l2, t3, t4) = _validate_lstats(4, lstats)
    if abs(t4) >= 1.:
        errmsg = "Invalid fourth L-moment : should be in ]-1;1["
        raise ValueError(errmsg)
    elif (t4 <= (5.*t3*t3 - 1.)/4.):
        errmsg = "Invalid fourth L-moment"
        raise ValueError(errmsg)
    elif (t4 >= (5.*t3*t3 + 1.)/6.):
        errmsg = "The third and fourth moments (tau-3, tau-4) lies above"\
                 " the generalized logistic line, which suggests"\
                 " that L-moments are not consistent with any Kappa"\
                 " distribution with h > -1."
        raise ValueError(errmsg)
    ((loc, scale, k, h), ierr) = _lmoments.pelkap((l1, l2, t3, t4))
    if (ierr == 3):
        errmsg = "Iteration failed to converge."
        raise ValueError(errmsg)
    elif (ierr == 4):
        errmsg = "Unable to make progress from current point in iteration."
        raise ValueError(errmsg)
    elif (ierr == 5):
        errmsg = "Iteration encountered numerical difficulties - "\
                 "overflow would have been likely to occur."
        raise ValueError(errmsg)
    elif (ierr == 6):
        errmsg = "Iteration for h and k converged, but overflow would have "\
                 "occurred when calculating xi and alpha."
        raise ValueError(errmsg)
    return (loc, scale, k, h)
extradist.kappa_gen.lmparams = staticmethod(_lmparams_kappa)



def _lmparams_norm(lstats):
    "(needs to be overwritten)"
    (loc, scale) = _validate_lstats(2, lstats)
    scale *= np.sqrt(np.pi)
    return (loc, scale)
dist.norm_gen.lmparams = staticmethod(_lmparams_norm)


def _lmparams_pearson3(lstats):
    "needs to be overwritten"
    # Method : Rational approximation is used to express the shape parameter
    # ...of the Gamma ditribution as a function of tau_3.
    # The relative accuracy of the approximation is better than 3e-5.
    small = 1e-6
    #
    (c0, c1, c2, c3) = (1., 0.2906 ,  0.18820,  0.0442)
    (d0, d1, d2, d3, d4, d5, d6) = ( 1.00000, 0.36067, -0.59567,  0.25361,
                                    -2.78861, 2.56096,-0.77045)
    #
    (lmom_1, lmom_2, tau_3) = _validate_lstats(3, lstats)
    t3 = abs(tau_3)
    # Zero skewness
    if abs(t3) <= small:
        return (lmom_1, lmom_2*np.sqrt(np.pi), 0.)
    #
    if t3 < 1./3:
        t = 3*np.pi*t3*t3
        a = (1. + c1*t)/((1. + t*(c2 + t*c3))*t)
    else:
        t = 1. - t3
        a = t*(d1 + t*(d2 + t*d3))/(1. + t*(d4+t*(d5 + t*d6)))
    sqrt_a = np.sqrt(a)
    beta = np.sqrt(np.pi) * lmom_2 
    beta *= np.exp(special.gammaln(a) - special.gammaln(a + 0.5))
    mu = lmom_1
    sigma = beta*sqrt_a
    alpha = 2./sqrt_a
    if tau_3 < 0:
        alpha = -alpha
    return (mu, sigma, alpha)
extradist.pearson3_gen.lmparams = staticmethod(_lmparams_pearson3)



def _lmparams_wakeby(lstats):
    "(needs to be overwritten)"
    (lmom_1, lmom_2, tau_3, tau_4, tau_5) = _validate_lstats(5, lstats)
    for (tau, lab) in zip((tau_4, tau_5), ('fourth', 'fifth')):
        if abs(tau) >= 1.:
            errmsg = "Invalid %s moment:must be in ]-1;+1[" % lab
            raise ValueError(errmsg)
    (lambda_1, lambda_2) = (lmom_1, lmom_2)
    (lambda_3, lambda_4, lambda_5) = lmom_2 * np.array([tau_3, tau_4, tau_5])
    #
    # Estimate n1, n2, n3, c1, c2, c3 when (xi != 0.)
    n1 = +3.*lambda_2 - 25.*lambda_3 +  32.*lambda_4
    n2 = -3.*lambda_2 +  5.*lambda_3 +   8.*lambda_4
    n3 = +3.*lambda_2 +  5.*lambda_3 +   2.*lambda_4
    c1 = +7.*lambda_2 - 85.*lambda_3 + 203.*lambda_4 - 125.*lambda_5
    c2 = -7.*lambda_2 + 25.*lambda_3 +   7.*lambda_4 -  25.*lambda_5
    c3 = +7.*lambda_2 +  5.*lambda_3 -   7.*lambda_4 -   5.*lambda_5

    # Estimates b and d
    _a = n2*c3 - c2*n3
    _b = n1*c3 - c1*n3
    _c = n1*c2 - c1*n2
    
    discrim = _b*_b - 4.*_a*_c
    if discrim >= 0.:
        discrim **= 0.5
        root_1 = 0.5*(-_b + discrim)/_a
        root_2 = 0.5*(-_b - discrim)/_a
        b =  max(root_1, root_2)
        d = -min(root_1, root_2)
        if d < 1.:
            # Estimates a, c, xi
            a =  (1.+b)*(2.+b)*(3.+b)*((1.+d)*lambda_2 - (3.-d)*lambda_3)
            a /= (4.*(b+d))
            c = -(1.-d)*(2.-d)*(3.-d)*((1.-b)*lambda_2 - (3.+b)*lambda_3)
            c /= (4.*(b+d))
            xi = lambda_1 - a/(1.+b) - c/(1.-d)
            # Are the parameters valid ?
            if (c >= 0.) & (a + c >= 0.):
                return (xi, a, b, c, d)
    # Estimates can only be found by setting xi=0
    xi = 0
    n1 = 4.*lambda_1 - 11.*lambda_2 + 9.*lambda_3
    n2 = -lambda_2 + 3.*lambda_3
    n3 = +lambda_2 + lambda_3
    c1 = 10.*lambda_1 - 29.*lambda_2 + 35.*lambda_3 - 16.*lambda_4
    c2 = -lambda_2 + 5.*lambda_3 -4.*lambda_4
    c3 = +lambda_2 - lambda_4
    _a = n2*c3 - c2*n3
    _b = n1*c3 - c1*n3
    _c = n1*c2 - c1*n2
    discrim = _b*_b - 4.*_a*_c
    if discrim >= 0.:
        discrim **= 0.5
        root_1 = 0.5*(-_b + discrim)/_a
        root_2 = 0.5*(-_b - discrim)/_a
        b =  max(root_1, root_2)
        d = -min(root_1, root_2)
        if d < 1.:
            a = +(1.+b)*(2.+b)/(b+d)*(lambda_1 - (2.-d)*lambda_2)
            c = -(1.-d)*(2.-d)/(b+d)*(lambda_1 - (2.+b)*lambda_2)
            if (c >= 0.) & (a + c >= 0.):
                return (xi, a, b, c, d)
    # Can't find valid estimates even w/ xi=0: fit a generalized Pareto instead.
    d = -(1.-3.*tau_3)/(1.+tau_3)
    c = (1.-d)*(2.-d)*lmom_2
    b = 0
    a = 0
    xi = lmom_1 - c/(1.-d)
    if d > 0:
        return (xi, a, b, c, d)
    return (xi, c, -d, 0, 0)
extradist.wakeby_gen.lmparams = staticmethod(_lmparams_wakeby)

###############################################################################
## Log-Pearson
## Note: OK, it's tricky: one has to use the moments
## of the DECIMAL LOG of the data as inputs
#def plp3(x,mu=0,sigma=1,gamma=1,parameters=None):
#    """CDF for Log-Pearson type III distribution (decimal log)
#    Input:
#        x               : sample distribution
#        mu,sigma,gamma : mean, std & skew parameters OF THE DECIMAL LOG!"""
#    if parameters: (mu,sigma,gamma) = parameters
#    if not hasattr(x,'__getitem__'):
#        result = cdfpe3(log10(x),(mu,sigma,gamma))
#    else:
#        n = len(F)
#        result = array(map(cdfpe3,log10(x),[(mu,sigma,gamma)]*n))
#    return result
##---
#def dlp3(x,mu=0,sigma=1,gamma=1,parameters=None):
#    """PDF for Log-Pearson type III distribution (decimal log)
#    Input:
#        x               : sample distribution
#        mu,sigma,gamma : mean, std & skew parameters OF THE DECIMAL LOG!"""
#    if parameters: (mu,sigma,gamma) = parameters
#    if abs(gamma) > 1e-12:
#        (alpha,beta,xi) = (4./gamma**2, 1./2.*sigma*abs(gamma), mu-2*sigma/float(gamma))
#        u = (log10(x)-xi)/float(beta)
#        denom = (beta * exp(dlgama(alpha)) * (x*log(10)))
#        f = 1./denom * abs(u)**(alpha-1.) * exp(-abs(u))
#    else:
#        u = (log10(x)-mu)/float(sigma)
#        f = 1./(sigma*sqrt(2*pi)*x*log(10.)) * exp(-0.5*u**2)
#    return f
##---
#def qlp3(F,mu=0,sigma=1,gamma=1,parameters=None):
#    """Quantile function of Log-Pearson type III Distribution (decimal log)
#    Input:
#        F               : CDF value
#        mu,sigma,gamma : mean, std & skew parameters OF THE DECIMAL LOG!"""
#    if parameters: (mu,sigma,gamma) = parameters
#    if not hasattr(F,'__getitem__'):
#        result = quape3(log10(F),(mu,sigma,gamma))
#    else:
#        n = len(F)
#        result = pow( repeat([10,],n), array(map(quape3,F,[(mu,sigma,gamma)]*n)) )
#    return result
##---
#def mlp3(lmom):
#    """L-moments estimators of the Log-Generalized Pearson III distribution parameters.
#    Input:
#        lmom : Array of L-moments (min. dim: 3) computed on the DECIMLA LOG data
#    Output:
#        para : (mu,sigma,gamma) parameters """
#    return pelpe3(lmom[0:3])
##----
#def xtlp3(mu=0,sigma=1,gamma=1,parameters=None,tol=1e-6):
#    """Domain of the Log-PearsonIII distribution : !!!DECIMAL LOG!!!
#    Input  : parameters
#    Output : (xmin,xmax)"""
#    if parameters: (mu,sigma,gamma) = parameters
#    (beta,xi) = (0.5*sigma*abs(gamma), mu-2*sigma/float(gamma))
#    huge = 740. #745.13225
#    if gamma>0: xthrg = (10**(xi-beta*huge+tol),1e20)
#    else: xthrg = (-1e20,10**(xi-beta*huge+tol))
#    return xthrg
################################################################################
#
#
################################################################################
## Cluster Analysis
##----------------------------
#def clusteragglo(x,n,method=1):
#    """Cluster analysis by any of several agglomerative hierarchical methods
#    Input:
#        method : clustering method. 
#                 1 : single-link clustering
#                 2 : complete-link clustering
#                 3 : ward's procedure
#        x      : data array : x(i,j) contains the j'th attribute of the i'th data point
#        n      : nb of data points
#    Outputs:
#      merge    : array of dimension (2,n). merge(1,i) and merge(2,i)
#                 are the labels of the clusters merged at the i'th
#                 stage.  merge(1,n) and merge(2,n) are not used.
#        disp   : array of length n.  disp(i) is a measure of the
#                 within-cluster dispersion after the i'th merge.
#                 dispersion is defined differently for each method:
#                 see below.  disp(n) is not used.
#
#Agglomerative hierarchical clustering: general description.
#Initially there are N clusters, each containing one data point,
#labeled 1 through N in the same order as the data points.  At each
#stage of clustering, two clusters are merged.  Their labels are saved
#in the MERGE array.  The smaller of the two labels is used as the
#label of the merged cluster.  After the Mth stage of clustering
#there are N-M clusters.  To find which data points belong to which
#clusters, use routine CLUINF.
#
#Single-link clustering: the distance between two clusters A and B is
#defined to be the minimum of the Euclidean distances between pairs of
#points with one point in A and one in B.  At each stage, the two
#clusters separated by the smallest distance are merged.  The square
#of this distance is saved in the corresponding element of array DISP.
#
#Complete-link clustering: the distance between two clusters A and B
#is defined to be the maximum of the Euclidean distances between pairs
#of points with one point in A and one in B.  At each stage, the two
#clusters separated by the smallest distance are merged.  The square
#of this distance is saved in the corresponding element of array DISP.
#DISP(I) is therefore the largest squared Euclidean distance between
#two points that are in the same cluster after the Ith merge.
#
#Ward's procedure: at each stage, the clusters that are merged are
#chosen to minimize the within-cluster sum of squared deviations of
#each attribute about the cluster mean.  This sum of squares is saved
#in the corresponding element of array DISP."""
#    iwork = zeros(n,'d')
#    nw = n*(n-1)/2+1
#    work = zeros(nw,'d')
#    (merge,disp) = cluagg(method,x,n,iwork,work)
#    return (merge,disp)
##----------------------------
#def clusterinfo(ncluster,merge):
#    """Obtains information about clusters arising from agglomerative
#hierarchical clustering
#
#Agglomerative hierarchical clustering procedures typically produce a
#list of the clusters merged at each stage of the clustering.  This
#routine uses this list to construct arrays that explicitly show
#which cluster a given data point belongs to, and which data points
#belong to a given cluster.
#
#Inputs
#    nclust : number of clusters
#    n      : number of data points
#    merge  : array of dimension (2,n). merge(1,i) and merge(2,i)
#             identify the clusters merged at the i'th step.
#             This is the array merge returned by routine cluagg,
#             and should be left unchanged after exit from that
#             routine.
#Outputs
#    iassgn : array of length n. its i'th element is the number
#             of the cluster to which the i'th data point belongs.
#    list   : array of length n. Contains the data points in
#             cluster 1, followed by the data points in cluster 2,
#             etc.  Data points in each cluster are listed in
#             increasing order.  The last data point in each
#             cluster is indicated by a negative number.
#             See the example below.
#    num    : array of length nclust.  number of data points in
#             each cluster.
#
#Cluster numbers used in arrays iassgn, list and num range from 1 to
#nclust.  they are arbitrary, but are uniquely defined: cluster 1
#contains data point 1, cluster m (m.ge.2) contains data point j,
#where j=merge(2,n-m).
#
#Example of the list array.
#Suppose that there are 8 data points  and 3 clusters,
#and that the elements of the list array are
# 1, -4, 3, 6, -8, 2, 5, -7.
# Then the clusters are as follows: cluster 1 contains points 1 and 4;
# cluster 2 contains points  3, 6 and 8; cluster 3 contains points 2, 5 and 7."""
#    (iassgn,datalist,num) = cluinf(nclust,merge)
#    return (iassign,datalist,num)
##----------------------------
#def clusterkm(x,nclust,iassign,maxit=1):
#    """Cluster analysis by the k-means algorithm
#
#Inputs:
#    x     : array of dimension (nx,natt).  x(i,j) should
#            contain the j'th attribute for the i'th data point.
#   nx     : the first dimension of array x, as declared in the calling program.
#   natt   : number of attributes for each data point
#   nclust : number of clusters
#   n      : number of data points
#   iassgn : Array of length n.
#            On entry, should contain the initial assignment of sites to clusters.
#            On exit, contains the final assignment.  The i'th element of the array
#            contains the label of the cluster to which the i'th data point belongs.
#            Labels must be between 1 and nclust, and each of the values 1 through nclust
#            must occur at least once.
#    maxit  : maximum number of iterations for the k-means clustering algorithm
#Outputs:
#  list    : array of length n. Contains the data points in cluster 1,
#            followed by the data points in cluster 2, etc.
#            Data points in each cluster are listed in increasing order.
#            The last data point in each cluster is indicated by a negative number.
#  num     : array of length nclust.  number of data points in each cluster.
#  ss      : within-group sum of squares of the final clusters.
#    """
#    n = len(iassign)
#    iwork = zeros(nclust*3,'l')
#    nw = (n+nclust)*(natt+1)+2*nclust+1
#    rw = zeros(nw,'f') 
#    datalist,num = clukm(x,nclust,iassgn,ss,maxit,iwork,rw)
#    (datalist,num,ss) = clukm(x,nclust,iassgn,maxit,iwork,rw)
#    return (iassign,datalist,nsum,ss)
##------------------------------------------------------------------------------
## Regional Analysis
#def reglmr(xmom,weight,nmom=4):
#    """Regional weighted average of L-moments
#Inputs
#    nmom   : number of l-moments to be found (default=4).
#    xmom   : array of dimension (nxmom,nsite). x(i,j) contains
#             the i'th l-moment ratio for site j.
#    weight : array of length nsite. contains the weights to be applied to each site.
#Outputs    :    
#    rmom   : array of length nmom. on exit, contains the regional weighted average l-moment ratios.
#    """
#    return reglmr(xmom,weight)
#
#def regstat(names,length,xmom,prob,nsim=1,a=0,b=0,seed=123,kprint=1,kout=1):
#    """calculates three statistics useful in regional frequency analysis
#
# Discordancy measure, d(i), for individual sites in a region.
#    Large values might be used as a flag to indicate potential errors
#    in the data at the site.  "Large" might be 3 for regions with 15
#    or more sites, but less (exact values in array dc1) for smaller
#    regions.
#
# Heterogeneity measures, h(j), for a region based upon either:-
#    j=1: the weighted s.d. of the l-cvs or
#    j=2: the average distance from the site to the regional average
#         on a graph of l-cv vs. l-skewness
#    j=3: the average distance from the site to the regional average
#         on a graph of l-skewness vs. l-kurtosis
#    In practice h(1) is probably sufficient.  A value greater than
#    (say) 1.0 suggests that further subdivision of the region should
#    be considered as it might improve quantile estimates.
#
# Goodness-of-fit measures, z(k), for 5 candidate distributions:
#    k=1: generalized logistic
#    k=2: generalized extreme value
#    k=3: generalized normal (lognormal)
#    k=4: pearson type iii (3-parameter gamma)
#    k=5: generalized pareto
#    Provided that the region is acceptably close to homogeneous,
#    the fit may be judged acceptable at 10% significance level
#    if z(k) is less than 1.645 in absolute value.
#
# For further details see J.R.M. Hosking and J.R. Wallis (1997),
# "Regional frequency analysis: an approach based on l-moments",
# cambridge university press, chapters 3-5.
#
# parameters of routine:
#    names  : Character*12 array of length nsites. Site names.
#    length : array of length nsites. record lengths at each site.
#    xmom   : Array of dimension (5,nsites), containing the first 5 sample
#             L-moments for each site, in the order mean, L-cv, L-skewness,
#             L-kurtosis, T-5, i.e xmom(i,j) contains the i'th l-moment for site j.
#             NB xmom(2,.) contains l-cv, not the usual l-2!
#    a, b   : Parameters of plotting position.
#             Note: a and b should be the same as the values used
#             to calculate the moments in the xmom array.
#    seed   : Seed for random number generator. Should be a whole number
#             in the range 2d0 to 2147483647d0.
#    nsim   : Number of simulated worlds for heterogeneity and
#             goodness-of-fit tests.
#             Note: nsim=0 will force return at completion of
#                 outlier test.  nsim=1 will suppress calculation of
#                 h and z statistics, but parameter and quantile
#                 estimates will be found.
#    Prob   : Array of length nprob.  probabilities for which
#             quantiles are to be calculated.
#    kprint : Output flag. should be set to
#                 0  to suppress output
#                 1  to print output
#    kout   : Channel to which output is directed
#Outputs:    
#    rmom   : Array of length 5. On exit, contains the regional
#             weighted average L-moment ratios.
#    d      : Array of length nsites. on exit, contains the
#                 discordancy measure (d statistic) for each site.
#    vobs   : Array of length 3. On exit, contains the regional
#             observed values of 3 heterogeneity statistics:
#                 (1) weighted s.d. of L-cvs;
#                 (2) average of L-cv/L-skew distances;
#                 (3) average of L-skew/L-kurtosis distances.
#    vbar   : Array of length 3. On exit, contains the mean of the
#             simulated values of the 3 heterogeneity statistics.
#    vsd    : Array of length 3. On exit, contains the s.d. of the
#             simulated values of the 3 heterogeneity statistics.
#    h      : Array of length 3. On exit, contains heterogeneity
#             measures (h statistics), i.e. h=(vobs-vbar)/vsd.
#    z      : Array of length 5. On exit, contains goodness-of-fit
#             measures (z statistics) for 5 distributions:
#                 (1) gen. logistic, (2) gen. extreme value,
#                 (3) gen. normal, (4) pearson type iii,
#                 (5) gen. pareto.
#     para   : Array of dimension (5,6). On exit, if nsim.ge.1,
#              contains parameters of growth curves fitted by the
#              above 5 distributions, plus wakeby.
#    """
#    nsites = xmom.shape[1]
#    nprob = len(prob)
#    rmom = zeros(5,'d')
#    d = zeros(nsites,'d')
#    vobs = zeros(3,'d'); vbar = zeros(3,'d'); vsd= zeros(3,'d')
#    h = zeros(5,'d'); z = zeros(5,'d')
#    para = zeros((5,6),'d')
#    regtst(nsites,names,length,xmom,a,b,seed,nsim,nprob,prob,
#           kprint,kout,rmom,d,vobs,vbar,vsd,h,z,para)
#    return (rmom, d,
#            reshape(concatenate((vobs,vbar,vsd)),(3,-1)),
#            reshape(concatenate((h,z,ravel(para))),(5,8)))


if __doc__ is not None:
    #
    _pel_doc_template = """
    Estimates the parameters of a %(name)s distribution from its L-moments.

    Parameters
    ----------
    lmoments : sequence
        A sequence of L-moments. Must be at least of size %(size)s.

    Returns
    -------
    parameters : tuple
        A tuple %(params)s
"""
    pelparams = {2: "(location, scale)",
                 3: "(location, scale, shape)"}
    _lmparams_expon.__doc__ = _pel_doc_template % dict(name="exponential",
                                                       size=2,
                                                       params=pelparams[2])
    _lmparams_gamma.__doc__ = _pel_doc_template % dict(name="Gamma", size=2,
                                                       params=pelparams[2])
    _lmparams_norm.__doc__ = _pel_doc_template % dict(name="normal", size=2,
                                                      params=pelparams[2])
    _lmparams_genextreme.__doc__ = _pel_doc_template % dict(name="GEV", size=3,
                                                            params=pelparams[3])
    _lmparams_gennorm.__doc__ = _pel_doc_template % dict(name="generalized normal",
                                                         size=3,
                                                         params=pelparams[3])
    _lmparams_genpareto.__doc__ = _pel_doc_template % dict(name="Pareto",
                                                           size=3,
                                                           params=pelparams[3])
    _lmparams_glogistic.__doc__ = _pel_doc_template % dict(name="logistic",
                                                           size=3,
                                                           params=pelparams[3])
    _lmparams_gumbel_r.__doc__ = _pel_doc_template % dict(name="Gumbel", size=2,
                                                          params=pelparams[2])
    _lmparams_pearson3.__doc__ = _pel_doc_template % dict(name="Pearson type III",
                                                          size=3,
                                                          params=pelparams[3])
    _lmparams_wakeby.__doc__ = _pel_doc_template % dict(name="Wakeby", size=5,
                                                        params="(location, "\
                                                               " scale,"\
                                                               " b, c, d)")


from numpy.ma.testutils import *

import scipy.special as special
import scipy.integrate as integrate




if __name__ == "__main__":
    if 1:
        print "normal"
        print "Generic:", [dist.norm()._lmomg(_) for _ in (1,2,3,4)]
        print "Standard:", dist.norm().lmoments(4)
    if 1:
        print "expon"
        print "Generic:", [dist.expon()._lmomg(_) for _ in (1,2,3,4)]
        print "Standard:", dist.expon().lmoments(4)
    if 1:
        reorg = lambda (m, s, g): (g, m, s)
        lmoms = (lmom_1, lmom_2, tau_3) = (0., 1., 0.)
        params = extradist.pearson3.lmparams(lmoms)
        assert_almost_equal(np.array(params), (0.000000, 1.772454, 0.000000))
        assert_almost_equal(pearson3(*reorg(params)).lmoments(3),
                            np.array(lmoms),)
        #
        lmoms = (lmom_1, lmom_2, tau_3) = (0., 1., 0.5)
        params = extradist.pearson3.lmparams(lmoms)
        assert_almost_equal(np.array(params), (0.000000, 2.299931, 3.079345), 6)
        assert_almost_equal(pearson3(*reorg(params)).lmoments(3),
                            np.array(lmoms), 5)
        #
        lmoms = (lmom_1, lmom_2, tau_3) = (-0.5, 0.5, 0.5)
        params = extradist.pearson3.lmparams(lmoms)
        assert_almost_equal(np.array(params), (-0.5000, 1.149966, +3.079345), 6)
        assert_almost_equal(pearson3(*reorg(params)).lstats(3),
                            np.array(lmoms), 5)

