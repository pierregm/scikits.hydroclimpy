"""

The :mod:`~scikits.hydroclimpy.lib.spi` module defines a single function, 
:func:`~scikits.hydroclimpy.lib.spi.spi` to compute the standardized precipitation
index.

.. autofunction:: spi

"""


import numpy as np
import numpy.ma as ma

import scikits.hydroclimpy
from scikits.hydroclimpy import TimeSeries, Date, DateArray,\
                                date_array, time_series, const as _c

from scikits.timeseries.lib.moving_funcs import mov_sum

import scipy.stats.distributions as ssd


def _get_gamma_cdf(aseries, condition):
    """
    Returns the CDF values for aseries.

    Parameters
    ----------
    aseries : TimeSeries
        Annual series of data (one column per period)
    condition : TimeSeries
        Period mask.
    """
    # Mask the months for which no precipitations were recorded
    aseries_ = ma.masked_values(aseries, 0)
    # Get the proportion of 0 precipitation for each period (MM/WW)
    pzero = 1. - aseries_.count(axis=0) / aseries.count(axis=0).astype(float)
    # Mask outside the reference period
    aseries_._mask |= condition._data
    meanrain = aseries_.mean(axis=0)
    aleph = ma.log(meanrain) - ma.log(aseries_).mean(axis=0)
    alpha = (1. + ma.sqrt(1.+4./3*aleph)) / (4.*aleph)
    beta = meanrain/alpha
    # Get the Gamma CDF (per month)
    gcdf = pzero + (1.-pzero) * ssd.gamma.cdf(aseries,alpha,scale=beta)
    return gcdf

def _get_gaussian_cdf(aseries, condition):
    """
    Returns the CDF values for aseries.

    Parameters
    ----------
    aseries : TimeSeries
        Annual series of data (one column per period).
    condition : TimeSeries
        Period mask.
    """
    # Mask the months for which no precipitations were recorded
    aseries_ = ma.masked_values(aseries, 0)
    # Mask outside the reference period
    aseries_._mask |= condition._data
    meanrain = aseries_.mean(axis=0)
    stdrain = aseries_.std(axis=0, ddof=1)
    # Get the Normal CDF (per month)
    gcdf = ssd.norm.cdf(aseries,loc=meanrain,scale=stdrain)
    return gcdf

_cdfdict = dict(normal=_get_gaussian_cdf,
                gaussian=_get_gaussian_cdf,
                gamma=_get_gamma_cdf,
                default=_get_gamma_cdf)


def spi(rainfall, span=3, freq='monthly', 
        start_reference=None, end_reference=None, 
        distribution='gamma'):
    """
    Computes the standardized precipitation index for the precipitation series
    at the given span.

    Parameters
    ----------
    rainfall : TimeSeries
        Series of monthly cumulative precipitation.
    span : {3, int} optional
        Span, in months.
    freq : {'monthly','weekly'} optional
        Frequency
    start_reference : {None, Date} optional
        Starting date of the reference period. 
        If None, use the first date of the series.
    end_reference : {None, Date} optional
        Ending date of the reference period. 
        If None, use the last date of the series.
    distribution : {'gamma','normal'} optional
        Type of distribution for the data of each period.
    """
    if not isinstance(rainfall, TimeSeries):
        raise TypeError,\
              "The precipitation input should be a valid TimeSeries object!"
    _freq = freq.upper()[0]
    if _freq not in ['M','W']:
        raise ValueError,"Invalid frequency! Should be in ('Monthly','Weekly')"\
                         " (got '%s' instead)" % freq
    
    rain = rainfall.convert(_freq, func=ma.sum)
    dates = rain.dates
    # Convert to cumulative over span (mm) and as a Nx12 or Nx53 array... 
    arain = mov_sum(rain, span).convert('A')
    # Define the reference period ...................
    if start_reference is None:
        startd = dates[0]
    elif not isinstance(start_reference, Date):
        startd = Date(_freq, start_reference)
    else:
        startd = start_reference.asfreq(_freq)
    if end_reference is None:
        endd = dates[-1]
    elif not isinstance(end_reference, Date):
        endd = Date(_freq, end_reference)
    else:
        endd = end_reference.asfreq(_freq)
    condition = time_series((dates < startd) | (dates > endd), 
                            dates=dates).convert('A')
    # Get the cumulative distribution
    try:
        cdffunc = _cdfdict[distribution]
    except KeyError:
        raise ValueError("Unrecognized distribution '%s':"\
                         "Should be in %s)" % (distribution,
                                               ','.join(_cdfdict.keys())))
    cdf = cdffunc(arain, condition)
    # Get the corresponding normal scores
    zscores = ma.fix_invalid(ssd.norm.ppf(cdf))
    zscores._mask |= ma.getmaskarray(arain)
    # Retransform to a regular time series    
    if _freq == 'M':
        kwargs = {'month':1}
    else:
        kwargs = {'week':1}
    kwargs.update({'year':dates[0].year})
    spi = time_series(zscores.ravel(),
                      start_date=Date(_freq,**kwargs),
                      ).view(type(rainfall))
    return spi
################################################################################
