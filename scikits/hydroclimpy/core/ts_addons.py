"""
.. currentmodule:: scikits.hydroclimpy.core.ts_addons


:class:`~scikits.timeseries.TimeSeries` manipulation
====================================================

The following functions accept basic :class:`~scikits.timeseries.TimeSeries` 
objects as inputs.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   deseasonalize
   force_reference
   mask_inside_period
   mask_outside_period
   periodize
   convert_and_fill
   apply_on_fields



Convenience functions
=====================

.. autosummary::
   :toctree: generated/

   get_maskoptions

"""
import numpy as np
from numpy import ndarray

import numpy.ma as ma
from numpy.ma import masked

import scikits.timeseries as ts
from scikits.timeseries import Date, DateError, TimeSeries, \
                               check_freq
from scikits.timeseries import const as _c
from scikits.timeseries.lib import backward_fill


FR_ANNSTART = [
_c.FR_ANNDEC, _c.FR_ANNJAN, _c.FR_ANNFEB, _c.FR_ANNMAR, _c.FR_ANNAPR, _c.FR_ANNMAY,
_c.FR_ANNJUN, _c.FR_ANNJUL, _c.FR_ANNAUG, _c.FR_ANNSEP, _c.FR_ANNOCT, _c.FR_ANNNOV,
              ]

__all__ = ['FR_ANNSTART',
           'apply_on_fields',
           'convert_and_fill',
           'deseasonalize',
           'force_reference',
           'get_maskoptions',
           'mask_inside_period', 'mask_outside_period',
           'periodize']


def get_maskoptions(**options):
    """
    Processes a dictionary of optional parameters and returns a tuple
    of (options for masked arrays, remaining options).

    """
    remaining = options
    maskoptions = dict(fill_value=remaining.pop('fill_value', None),
                       hard_mask=remaining.pop('hard_mask', False),
                       )
    return (maskoptions, remaining)


def _mask_period(data, period=None, start_date=None, end_date=None,
                 inside=True, include_edges=True, inplace=False):
    """
    Returns a series masked where dates fall outside the selection period,
    as well as where data are initially missing (masked).

    Parameters
    ----------
    data : Timeseries
        Data to process
    period : {None, sequence}, optional
        A sequence of (starting date, ending date).
    start_date : {None, string/Date }, optional
        Starting date. If None, uses the first date of the series.
    end_date : {None, string/Date }, optional
        Ending date. If None, uses the last date of the series.
    inside : {True, False}, optional
        Whether the dates inside the range should be masked.
        If not, masks outside.
    include_edges : {True, False}, optional
        Whether the starting and ending dates should be masked.
    inplace : {True, False}, optional
        Whether the data mask should be modified in place.
        If not, returns a new :class:`~scikits.timeseries.TimeSeries` object.

    """
    data = ma.array(data, subok=True, copy=not inplace)
    if not isinstance(data, TimeSeries):
        raise ValueError("Data should be a valid TimeSeries!")
    dates = data._dates
    if dates.ndim == 1:
        dates_lims = dates[[0, -1]]
    else:
        dates_lims = dates.ravel()[[0, -1]]
    # Check the period .....................
    if period is not None:
        if isinstance(period, (tuple, list, ndarray)):
            (start_date, end_date) = (period[0], period[-1])
        else:
            (start_date, end_date) = (period, start_date)
    # Check the starting date ..............
    if start_date is None:
        start_date = dates_lims[0]
    elif isinstance(start_date, str):
        start_date = Date(data.freq, string=start_date)
    elif not isinstance(start_date, Date):
        raise DateError("The starting date should be a valid Date object!")
    # Check the ending date ................
    if end_date is None:
        end_date = dates_lims[-1]
    elif isinstance(end_date, str):
        end_date = Date(data.freq, string=end_date)
    elif not isinstance(end_date, Date):
        raise DateError("The ending date should be a valid Date object!")
    # Constructs the selection mask .........
    dates = data.dates
    if inside:
        if include_edges:
            selection = (dates >= start_date) & (dates <= end_date)
        else:
            selection = (dates > start_date) & (dates < end_date)
    else:
        if include_edges:
            selection = (dates <= start_date) | (dates >= end_date)
        else:
            selection = (dates < start_date) | (dates > end_date)
    data[selection] = masked
    return data


def mask_inside_period(data, start_date=None, end_date=None,
                       include_edges=True, inplace=False):
    """
    Masks a series where the dates fall inside a given period.

    Parameters
    ----------
    data : Timeseries
        Data to process
    start_date : {None, string/Date }, optional
        Starting date of the masking period.
        If None, uses the first date of the series.
    end_date : {None, string/Date }, optional
        Ending date of the masking period.
        If None, uses the last date of the series.
    include_edges : {True, False}, optional
        Whether the starting and ending dates should be masked.
    inplace : {True, False}, optional
        Whether the data mask should be modified in place.
        If not, returns a new :class:`~scikits.timeseries.TimeSeries` object.

    """
    kwargs = dict(start_date=start_date, end_date=end_date,
                  inside=True, include_edges=include_edges, inplace=inplace)
    return _mask_period(data, **kwargs)


def mask_outside_period(data, start_date=None, end_date=None,
                       include_edges=True, inplace=False):
    """
    Masks a series where the dates fall outside a given period.

    Parameters
    ----------
    data : Timeseries
        Data to process
    start_date : {None, string/Date }, optional
        Starting date. If None, uses the first date of the series.
    end_date : {None, string/Date }, optional
        Ending date. If None, uses the last date of the series.
    include_edges : {True, False}, optional
        Whether the starting and ending dates should be masked.
    inplace : {True, False}, optional
        Whether the data mask should be modified in place.
        If not, returns a new :class:`~scikits.timeseries.TimeSeries` object.

    """
    kwargs = dict(start_date=start_date, end_date=end_date,
                  inside=False, include_edges=include_edges, inplace=inplace)
    return _mask_period(data, **kwargs)



def deseasonalize(series, season, normalize=True, unbiased=True):
    """
    Deseasonalize a time series by substracting the seasonal average and
    dividing by the seasonal standard deviation (if needed).

    Returns a new time series, with the same frequency as the original.

    Parameters
    ----------
    series : TimeSeries
        Input time series.
    season : string/integer
        Frequency of the season.
        For example, use 'M' or 3000 for 12 seasons per year (one per month).
    normalize : {True, False}, optional
        Whether to divide by the seasonal standard deviation.
    unbiased : {True, False}, optional
        Whether the biased (False) or unbiased (True) estimate of the standard
        deviation should be computed.
    """
    if not isinstance(series, TimeSeries):
        raise TypeError("The input should be a valid TimeSeries!")
    season = ts.check_freq(season)
    output = series.copy()
    # 1.Monthly .............
    if ((season // _c.FR_MTH) == 1):
        input = series.asfreq(season).months
        nseason = 12
    # 2. Quarterly ..........
    elif ((season // _c.FR_QTR) == 1):
        input = series.asfreq(season).quarters
        nseason = 4
    # 3. Weekly .............
    elif ((season // _c.FR_WK) == 1):
        input = series.asfreq(season).weeks
        nseason = 53
    #........................
    else:
        raise NotImplementedError("Frequency currently supported: MTH, QTR")
    #..................................
    ddof = int(unbiased)
    for i in range(1, nseason + 1):
        selection = (input == i)
        selected = series[selection]
        output[selection] = current = selected.anom()
        if normalize:
            current /= selected.std(ddof=ddof)
    return output



def force_reference(series, ref, func=ma.mean):
    """
    Force a timeseries to match a reference, in terms of frequency and first
    and last dates.

    Parameters
    ----------
    series : TimeSeries
        Input time series
    ref : TimeSeries
        Reference time series.
    func : function
        Function used for the conversion of series to ref.freq, if necessary.

    """
    assert(isinstance(ref, TimeSeries) and isinstance(series, TimeSeries),
           "Inputs should be valid TimeSeries object ! (got '%s' and '%s')" % \
           (type(ref), type(series)))
    (sfreq, rfreq) = (series.freq, ref.freq)
    if sfreq > rfreq:
        series = series.convert(rfreq, func=func)
    elif sfreq < rfreq:
        series = backward_fill(series.convert(rfreq))
    if (series.dates[0] != ref.dates[0]) or (series.dates[-1] != ref.dates[-1]):
        return ts.adjust_endpoints(series, ref.dates[0], ref.dates[-1])
    return series



def convert_and_fill(series, freq, func=None, position='END', *args, **kwargs):
    """
    Converts a series from one frequency to another, by manipulating both the
    `data` and `dates` attributes.

    This function replicates the standard :func:`scikits.timeseries.convert`
    function.
    When converting to a higher frequency, the result is automatically forward-
    filled (if ``position=='START'``) or backward-filled (if ``position='END'``)

    The function accepts the same input parameters and has the same signature as
    the standard :func:`scikits.timeseries.convert` function.
    Please refer to the docstring of this latter for more information.

    See Also
    --------
    convert
        Original function
    """
    if not isinstance(series, TimeSeries):
        raise TypeError("The input should be a valid TimeSeries object!")
    result = ts.convert(series, freq, func=func, position=position,
                        *args, **kwargs)
    if freq > series._dates.freq:
        _series = result._series
        if position.upper()[0] == 'S':
            _series[:] = ts.lib.forward_fill(_series)
        else:
            _series[:] = ts.lib.backward_fill(_series)
    return result



def apply_on_fields(series, func, *args, **kwargs):
    """
    Apply the function ``func`` to each field of ``series``.

    Parameters
    ----------
    series : array-like
        A (subclass) of ndarray.
    func : function
        Function to apply.
    args : var
        Additional arguments to ``func``.
    kwargs : var
        Additional optional parameters to ``func``
    """
    names = series.dtype.names
    if names is None:
        return func(series, *args, **kwargs)
    results = [func(series[f], *args, **kwargs) for f in names]
    rdtype = [(f[0], r.dtype) for (f, r) in zip(series.dtype.descr, results)]
    rdata = [ma.getdata(r) for r in results]
    rmask = [ma.getmaskarray(r) for r in results]
    isarray = isinstance(results[0], np.ndarray)
    if isarray:
        output = ma.array(zip(*rdata), mask=zip(*rmask),
                          dtype=rdtype).view(type=type(series))
        output._update_from(series)
    else:
        output = ma.array(tuple(rdata), mask=tuple(rmask), dtype=rdtype)
    return output




#TODO : The following code hasn't been thoroughly checked yet...

def periodize(series, freq=None, func=None, *args, **kwargs):
    """
    Convert a series to an annual frequency.

    The initial series must have a monthly or quarterly frequency.
    If this is not the case, a temporary series is created from the conversion
    of the initial series to the ``freq`` parameter, using ``func`` as
    conversion function.

    Parameters
    ----------
    series : TimeSeries
        Series to convert.
    freq : var, optional
        A valid frequency specifier (as a string or integer), or a 3-letter 
        string corresponding to the first quarter of the year.
        Use this parameter if the series does not have an adequate frequency.
    func : {None,function}, optional
        Function controlling how data sharing the same new dates should be
        manipulated.
        This function should handle masked values appropriately.
    args : {extra arguments for func parameter}, optional
        Additional mandatory parameters of the ``func`` function.
    kwargs : {extra keyword arguments for func parameter}, optional
        Additional optional keyword parameters of the ``func`` function.

    Returns
    -------
    periodized_series : TimeSeries
        Annual T:class:`~scikits.timeseries.TimeSeries`, with an additional key
        ``period_names`` key in the :attr:`optinfo` dictionary. 
        The name of the periods depend
        on the frequency of the initial series or on the ``freq`` parameter.

    Examples
    --------
    To group a daily series by months, with accumulation :
       
       >>> periodize(series,"M", func=ma.sum)
       
    To group a daily series by quarters with accumulation, the first quarter
    starting in October:
       
       >>> periodize(series, 'OND', func=ma.sum)
       >>> periodize(series, 'Q-E-SEP', func=ma.sum)
    
    To group a daily series by quarters with accumulation, the first quarter
    starting in December:
       
       >>> periodize(series, 'NDJ', func=ma.sum)
       >>> periodize(series, 'Q-E-NOV', func=ma.sum)

    .. todo::
       This function is vestigial and has never been properly cleaned up.
       Caution is required...

    """
    invalidmsg = "The periodize function requires the input series "\
                 "to have a monthly or quarterly frequency. (got '%s' instead)"
    period_convert = {_c.FR_MTH: _c.FR_ANN, _c.FR_QTR: _c.FR_ANN,
                      _c.FR_QTREDEC: _c.FR_ANNDEC, _c.FR_QTRSJAN: _c.FR_ANNJAN,
                      _c.FR_QTREJAN: _c.FR_ANNJAN, _c.FR_QTRSFEB: _c.FR_ANNFEB,
                      _c.FR_QTREFEB: _c.FR_ANNFEB, _c.FR_QTRSMAR: _c.FR_ANNMAR,
                      _c.FR_QTREMAR: _c.FR_ANNMAR, _c.FR_QTRSAPR: _c.FR_ANNAPR,
                      _c.FR_QTREAPR: _c.FR_ANNAPR, _c.FR_QTRSMAY: _c.FR_ANNMAY,
                      _c.FR_QTREMAY: _c.FR_ANNMAY, _c.FR_QTRSJUN: _c.FR_ANNJUN,
                      _c.FR_QTREJUN: _c.FR_ANNJUN, _c.FR_QTRSJUL: _c.FR_ANNJUL,
                      _c.FR_QTREJUL: _c.FR_ANNJUL, _c.FR_QTRSAUG: _c.FR_ANNAUG,
                      _c.FR_QTREAUG: _c.FR_ANNAUG, _c.FR_QTRSSEP: _c.FR_ANNSEP,
                      _c.FR_QTRESEP: _c.FR_ANNSEP, _c.FR_QTRSOCT: _c.FR_ANNOCT,
                      _c.FR_QTREOCT: _c.FR_ANNOCT, _c.FR_QTRSNOV: _c.FR_ANNNOV,
                      _c.FR_QTRENOV: _c.FR_ANNNOV, _c.FR_QTRSDEC: _c.FR_ANNDEC,
                      }
    periodnames = {_c.FR_MTH:['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
                   _c.FR_QTR:['JFM', 'AMJ', 'JAS', 'OND'],
_c.FR_QTREJAN:['FMA', 'MJJ', 'ASO', 'NDJ'], _c.FR_QTRSJAN:['FMA', 'MJJ', 'ASO', 'NDJ'],
_c.FR_QTREFEB:['MAM', 'JJA', 'SON', 'DJF'], _c.FR_QTRSFEB:['MAM', 'JJA', 'SON', 'DJF'],
_c.FR_QTREMAR:['AMJ', 'JAS', 'OND', 'JFM'], _c.FR_QTRSMAR:['AMJ', 'JAS', 'OND', 'JFM'],
_c.FR_QTREAPR:['MJJ', 'ASO', 'NDJ', 'FMA'], _c.FR_QTRSAPR:['MJJ', 'ASO', 'NDJ', 'FMA'],
_c.FR_QTREMAY:['JJA', 'SON', 'DJF', 'MAM'], _c.FR_QTRSMAY:['JJA', 'SON', 'DJF', 'MAM'],
_c.FR_QTREJUN:['JAS', 'OND', 'JFM', 'AMJ'], _c.FR_QTRSJUN:['JAS', 'OND', 'JFM', 'AMJ'],
_c.FR_QTREJUL:['ASO', 'NDJ', 'FMA', 'MJJ'], _c.FR_QTRSJUL:['ASO', 'NDJ', 'FMA', 'MJJ'],
_c.FR_QTREAUG:['SON', 'DJF', 'MAM', 'JJA'], _c.FR_QTRSAUG:['SON', 'DJF', 'MAM', 'JJA'],
_c.FR_QTRESEP:['OND', 'JFM', 'AMJ', 'JAS'], _c.FR_QTRSSEP:['OND', 'JFM', 'AMJ', 'JAS'],
_c.FR_QTREOCT:['NDJ', 'FMA', 'MJJ', 'ASO'], _c.FR_QTRSOCT:['NDJ', 'FMA', 'MJJ', 'ASO'],
_c.FR_QTRENOV:['DJF', 'MAM', 'JJA', 'SON'], _c.FR_QTRSNOV:['DJF', 'MAM', 'JJA', 'SON'],
_c.FR_QTREDEC:['JFM', 'AMJ', 'JAS', 'OND'], _c.FR_QTRSDEC:['JFM', 'AMJ', 'JAS', 'OND'],
                  }
    #
    if not isinstance(series, TimeSeries):
        raise TypeError("The input should be a valid TimeSeries object!")
    #
    if freq is None:
        freq = series._dates.freq
        if not (freq in period_convert):
            raise ValueError(invalidmsg % series._dates.freqstr)
    else:
        if isinstance(freq, basestring):
            freq = freq.upper()
            try:
                freq = check_freq(freq)
            except ValueError:
                months = 'JFMAMJJASONDJF'
                try:
                    idx = months.index(freq)
                except ValueError:
                    raise ValueError(invalidmsg % series._dates.freqstr)
                freq = freqlist[idx]
        else:
            freq = check_freq(freq)
        if not (freq in period_convert):
            raise ValueError(invalidmsg % series._dates.freqstr)
        if func is None:
            raise ValueError("When converting to an intermediary frequency, "
                             "the conversion function cannot be None!")
        series = series.convert(freq, func=func, *args, **kwargs)
    converted = series.convert(period_convert[freq])
    #
    names = periodnames[freq]
    output = converted
    output.optinfo['period_names'] = names
#    ndtype=[(_,converted.dtype.char) for _ in names]
#    mdtype=[(_,bool) for _ in names]
#    output = converted.__class__(np.fromiter(converted._data, dtype=ndtype),
#                                 dates=converted._dates,
#                                 mask=np.fromiter(ma.getmaskarray(converted),
#                                                  dtype=mdtype))
#    output._update_from(converted)
    return output


