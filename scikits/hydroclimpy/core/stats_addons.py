"""
Generic statistics functions, with support to masked arrays.

.. autofunction:: qqcalc

.. autofunction:: biweight



"""

__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"



__all__ = ['biweight',
           'qqcalc'
          ]

import itertools as itools

import numpy as np
import numpy.ma as ma

import scipy.stats.distributions as ssd
import scipy.stats.mstats as mstats
from scipy.stats.mstats import mquantiles, plotting_positions,\
                               hdquantiles, hdmedian
__all__ += ['mquantiles', 'plotting_positions']
__all__ += ['hdquantiles', 'hdmedian']



def qqcalc(data, distrib=ssd.norm, alpha=.4, beta=.4):
    """
    Returns the theoretical quantiles from an empirical distribution.
    
    Parameters
    ----------
    data : array
        Input data
    distrib : {norm, function}, optional
        Theoretical distribution used to compute the expected quantiles.
        If None, use a normal distribution.
        Otherwise, ``distrib`` must have a :meth:`.ppf` method.
    alpha : {float}, optional
        Coefficient for the computation of plotting positions
    beta : {float}, optional
        Coefficient for the computation of plotting positions.

    """
    pp = mstats.plotting_positions(data, alpha=alpha, beta=beta)
    qq = ma.fix_invalid(distrib.ppf(pp))
    qq._mask = pp._mask
    return qq

#------------------------------------------------------------------------------
def biweight(x, cst):
    """
    Computes the biweight average and midvariance for a given 1D array.
    Returns a tuple (biweight mean, biweight variance).

    Parameters
    ----------
    x: {ndarray}
        Input Array
    cst : {float} 
        Parameter controlling how outliers are censored.

    Notes
    -----
    The function is restricted to 1D data only.
    """
    assert (x.ndim == 1, "1D array only !")
    xmed = ma.median(x, 0)
    manom = x - xmed
    mad = ma.median(ma.absolute(manom))
    u_i = (manom/float(cst*mad))
    u_i *= ma.less_equal(ma.absolute(u_i), 1.).astype(float)
    w_i = (1-u_i**2)
    if ma.count(w_i) > 0:
        biw_m = xmed + ma.sum(manom * w_i**2)/ma.sum(w_i**2)
    else:
        biw_m = xmed
    biw_sd = ma.sqrt(ma.count(x)*ma.sum(manom**2 * w_i**4))
    biw_sd *= 1./ma.absolute(ma.sum(w_i * (1-5*u_i**2)))
    return (biw_m, biw_sd.item())
