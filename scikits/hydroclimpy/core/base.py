"""
.. currentmodule:: scikits.hydroclimpy.core.base


:class:`ReferencedSeries` objects
=================================

Overview
--------

Defines a :class:`~scikits.timeseries.TimeSeries` object with a given
reference period. 

The reference period consists simply of a starting date and an ending date.
These dates do not have to be part of the series.

The reference period plays a role in statistical analysis.
For example, it is customary in climate science to define anomalies as the 
deviations from the average calculated on a given reference period chosen for
its representativity.



Construction
~~~~~~~~~~~~

To construct a :class:`ReferencedSeries` object, the simplest method is 
to directly call the class constructor with the proper parameters.

However, the recommended way is to use the :func:`referenced_series` 
factory function.



.. _refseries_attributes:

Specific attibutes
------------------

.. attribute:: ReferencedSeries.refperiod

   Returns the reference period, as a :class:`~scikits.timeseries.DateArray`.
   This attribute is writable.
   Setting the reference period to :const:`None` will force the period to
   match the first and last dates (in chronological order).


.. attribute:: ReferencedSeries.optinfo

   Dictionary of optional information.
   This attribute is read-only. 
   However, the different entries of the dictionary can be modified.



.. _refseries_methods:

Specific methods
----------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ~ReferencedSeries.period_average
   ~ReferencedSeries.period_mean
   ~ReferencedSeries.period_var
   ~ReferencedSeries.period_std
   ~ReferencedSeries.period_anomalies
   ~ReferencedSeries.period_standardize



API
---

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ReferencedSeries
   referenced_series


"""


import numpy as np
from numpy import ndarray

import numpy.ma as ma
from numpy.ma import nomask, getmaskarray

import scikits.timeseries as ts
from scikits.timeseries import Date, DateArray, DateError, TimeSeries, \
                               time_series

from ts_addons import get_maskoptions, mask_outside_period


__all__ = ['ReferencedSeries',
           'referenced_series']

_docstr = dict(
out="""out : ndarray
        Alternative output array in which to place the result. It must have
        the same shape as the expected output but the type will be cast if
        necessary.
    """,
dtype="""dtype : type
        Type to use in computing the standard deviation. For arrays of
        integer type the default is float32, for arrays of float types it
        is the same as the array type.
    """,
axis="""axis : integer, optional
        Axis along which to perform the operation. The default is to perform
        the operation on a flattened view of the input series.
    """,
ddof="""ddof : integer, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N-ddof``.
    """,
period="""period : {None, tuple}, optional
        Reference period of the computation, as a tuple (starting date, ending date)
        If None, use the reference period of the input series.
    """,
)



class ReferencedSeries(TimeSeries, object):
    """
    Defines a time series with a given reference period. 
    
    The reference period consists simply of a starting date and an ending date.
    These dates do not have to be part of the series.

    The reference period plays a role in statistical analysis.
    For example, it is customary in climate science to define anomalies as the 
    deviations from the average calculated on a given reference period chosen for
    its representativity.
    
    
    Parameters
    ----------
    data : array-like
        Data
    dates : {None, sequence}, optional
        Corresponding dates.
    mask : {nomask, sequence}, optional
        Mask.
    refperiod : {None, tuple}, optional
        Reference period, as a tuple (starting date, ending date).
        If None, the reference period is set to the whole range of dates.
    freq : {None, string or integer}, optional
        A valid frequency specification.
    start_date : {None, Date}, optional
        Starting date of the series.
        This parameter is only useful if ``dates`` is None.
    autosort : {True, False}, optional
        Whether to sort the series in chronological order.
    dtype : {None, dtype}, optional.
        Datatype of the series.
        If None, the dtype of ``data`` is used.
    copy : {False, True}, optional
        Whether copy the data (``True``) or just link to it (``False``).
    \*\*optional_parameters** :
        All the parameters recognized by :class:`~numpy.ma.MaskedArray`
        are also recognized.

    See Also
    --------
    :class:`scikits.timeseries.TimeSeries`
        Base class of :class:`ReferencedSeries`

    """
    # Dirty trick: we need the following lines to stop pylint's whining
    _dates = np.array([], dtype=int)
    _optinfo = {}

    def __new__(cls, data, dates=None, mask=nomask, refperiod=None,
                freq=None, start_date=None, autosort=True,
                dtype=None, copy=False, **options):
        (maoptions, options) = get_maskoptions(**options)
        maoptions.update(dict(copy=copy, dtype=dtype))
        _data = ts.time_series(data, dates=dates, mask=mask, freq=freq,
                               start_date=start_date, autosort=autosort,
                               **maoptions).view(cls)
        # Set the reference period
        if refperiod is None or tuple(refperiod) == (None, None):
            # Don't call  refperiod yet, in case we come from a pickle
            _data._optinfo['reference_period'] = refperiod
        else:
            # OK, here we can call refperiod
            _data.refperiod = refperiod
        return _data


    optinfo = property(fget=lambda self:self._optinfo,
                       doc="Dictionary of optional information")


    def _set_refperiod(self, period=None):
        "Sets the period of reference."
        _optinfo = self._optinfo
        # Set the refperiod to the range of current dates
        _refperiod = _optinfo.get('reference_period', None)
        if period is None:
            if (_refperiod is None):
                _optinfo['reference_period'] = self._dates[[0, -1]]
            return
        if not isinstance(period, (tuple, list, ndarray)):
            msg = "The reference period should be a tuple "\
                  "(start_date,end_date)."
            raise ValueError(msg)
        # Check the starting and ending dates ...
        dates = self._dates
#        if dates.ndim == 1:
#            dates_lims = dates[[0,-1]]
#        else:
#            dates_lims = dates.ravel()[[0,-1]]
        (start_date, end_date) = (period[0], period[-1])
        if isinstance(start_date, str):
            start_date = Date(self.freq, string=start_date)
        elif not isinstance(start_date, (Date, DateArray)):
            raise DateError, "The starting date should be a valid Date object!"\
                             " ( got %s instead)" % (start_date.__class__)
        #
        if isinstance(end_date, str):
            end_date = Date(self.freq, string=end_date)
        elif not isinstance(end_date, Date):
            raise DateError("The ending date should be a valid Date object!"\
                            " ( got %s instead)" % (end_date.__class__))
        _optinfo['reference_period'] = DateArray((start_date, end_date),
                                                 freq=self.freq)
        return

    def _get_refperiod(self):
        "Returns the period of reference."
        _optinfo = self._optinfo
        _dates = self._dates
        refperiod = _optinfo.get('reference_period', None)
        if refperiod is None:
            _optinfo['reference_period'] = DateArray((_dates.start_date,
                                                      _dates.end_date),
                                                     freq=self.freq)
        else:
            _optinfo['reference_period'] = DateArray(refperiod,
                                                     freq=self._dates.freq)
        return _optinfo['reference_period']

    refperiod = property(
                fget=_get_refperiod, fset=_set_refperiod,
                doc="Reference period  (alias to :attr:`reference_period`)."
                )
    reference_period = property(fget=_get_refperiod, fset=_set_refperiod,
                                doc="Reference period.")

    def convert(self, freq, func=None, position='END', *args, **kwargs):
        "(this docstring should be overwritten)"
        result = super(ReferencedSeries, self).convert(freq, func=func,
                                                       position=position,
                                                       *args, **kwargs)
        _optinfo = result._optinfo
        refperiod = _optinfo.get('reference_period', None)
        if refperiod is not None:
            _optinfo['reference_period'] = refperiod.asfreq(freq)
        return result
    convert.__doc__ = TimeSeries.convert.__doc__


    def period_average(self, axis=None, dtype=None, out=None, period=None):
        """
    Returns the series averaged over the reference period, along the given axis.

    Parameters
    ----------
    %(axis)s
    %(dtype)s
    %(out)s
    %(period)s
        """
        if period is None:
            period = self.refperiod
        elif not isinstance(period, (tuple, list, ndarray)):
            msg = "Period should be a tuple (starting date, ending date)!"
            raise ValueError, msg
        refdata = mask_outside_period(self, period[0], period[1],
                                      include_edges=False)
        return refdata.mean(axis=axis, dtype=dtype, out=out)
    period_mean = period_average
    period_average.__doc__ = ((period_average.__doc__ or "") % _docstr) or None


    def period_var(self, axis=None, dtype=None, out=None, ddof=1,
                   period=None):
        """
    Returns the standard deviation over a reference period.
        
    Parameters
    ----------
    %(axis)s
    %(dtype)s
    %(out)s
    %(ddof)s
    %(period)s

    Warnings
    --------
    The default ``ddof`` is 1: 
    by default, the method returns the unbiased estimate of the variance.

        """
        if period is None:
            period = self.refperiod
        elif not isinstance(period, (tuple, list, ndarray)):
            msg = "Period should be a tuple (starting date, ending date)!"
            raise ValueError, msg
        refdata = mask_outside_period(self, period[0], period[1],
                                      include_edges=False)
        return refdata.var(axis=axis, dtype=dtype, out=out, ddof=ddof)
    period_var.__doc__ = ((period_var.__doc__ or "") % _docstr) or None


    def period_std(self, axis=None, dtype=None, out=None, ddof=1,
                   period=None):
        """
    Returns the standard deviation over a reference period.
        
    Parameters
    ----------
    %(axis)s
    %(dtype)s
    %(out)s
    %(ddof)s
    %(period)s

    Warnings
    --------
    The default ``ddof`` is 1:
    by default, the method returns the unbiased estimate of the standard deviation.

        """
        if period is None:
            period = self.refperiod
        elif not isinstance(period, (tuple, list, ndarray)):
            msg = "Period should be a tuple (starting date, ending date)!"
            raise ValueError, msg
        refdata = mask_outside_period(self, period[0], period[1],
                                      include_edges=False)
        return refdata.std(axis=axis, dtype=dtype, out=out, ddof=ddof)
    period_std.__doc__ = ((period_std.__doc__ or "") % _docstr) or None


    def period_anomalies(self, axis=None, dtype=None, period=None):
        """
    Returns data anomalies (deviations from average), where the average
    is defined on a reference period, along a given axis.

    Parameters
    ----------
    %(axis)s
    %(dtype)s
    %(period)s
        
        """
        if period is None:
            period = self.refperiod
        elif not isinstance(period, (tuple, list, ndarray)):
            msg = "Period should be a tuple (starting date, ending date)!"
            raise ValueError(msg)
        period_mean = self.period_mean(axis=axis, dtype=dtype, period=period)
        if not axis:
            return self - period_mean
        else:
            return self - ma.expand_dims(period_mean, axis)
    period_anomalies.__doc__ = ((period_anomalies.__doc__ or "") % _docstr) or\
                               None


    def period_standardize(self, axis=None, dtype=None, ddof=1, period=None):
        """
    Standardizes data by substracting the average over a reference period,
    and by dividing by the standard deviation over the same period.
    

    Parameters
    ----------
    %(axis)s
    %(dtype)s
    %(ddof)s
    %(period)s

    Warnings
    --------
    The default ``ddof`` is 1: 
    by default, the method returns the unbiased estimate of the standard deviation.
    
        """
        if period is None:
            period = self.refperiod
        elif not isinstance(period, (tuple, list, ndarray)):
            msg = "Period should be a tuple (starting date, ending date)!"
            raise ValueError, msg
        refdata = mask_outside_period(self, period[0], period[1],
                                      include_edges=False)
        refavg = refdata.mean(axis=axis, dtype=dtype)
        refstd = refdata.std(axis=axis, dtype=dtype, ddof=ddof)
        if not axis:
            result = (self - refavg) * 1. / refstd
        else:
            result = (self - ma.expand_dims(refavg)).astype(float)
            result /= ma.expand_dims(refstd)
        return result
    period_standardize.__doc__ = ((period_standardize.__doc__ or "") % _docstr)\
                                 or None

    #... Pickling

    def __getstate__(self):
        """
    
    Returns the internal state of the TimeSeries, for pickling purposes.
    
    """
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self._data.tostring(),
                 getmaskarray(self).tostring(),
                 self._fill_value,
                 self._dates.shape,
                 np.asarray(self._dates).tostring(),
                 self.freq,
                 self.optinfo
                 )
        return state
    #    
    def __setstate__(self, state):
        """
    Restores the internal state of the TimeSeries, for pickling purposes.
    The argument ``state`` is typically the output of the ``__getstate__`` output,
    and is a 5-element tuple:
    
        - class name
        - a tuple giving the shape of the data
        - a typecode for the data
        - a binary string for the data
        - a binary string for the mask.
        """
        (ver, shp, typ, isf, raw, msk, flv, dsh, dtm, frq, infodict) = state
        super(ReferencedSeries, self).__setstate__(state)
        self._dates.__setstate__((ver, dsh, np.dtype(int), isf, dtm, frq))
        self._dates.freq = frq
        # Update the _optinfo dictionary
        self._optinfo.update(infodict)
        # Fix the reference period
        refperiod = infodict['reference_period']
        if refperiod is not None:
            self._optinfo['reference_period'] = DateArray(refperiod, freq=frq)




def referenced_series(data, dates=None, start_date=None, freq=None, mask=nomask,
                      refperiod=None, dtype=None, copy=False, fill_value=None,
                      keep_mask=True, hard_mask=False):
    """
    Creates a :class:`ReferencedSeries` object.

    Parameters
    ----------
    data : array_like
        Data portion of the array. Any data that is valid for constructing a
        MaskedArray can be used here. May also be a TimeSeries object.
    dates : {None, DateArray}, optional
        A sequence of dates corresponding to each entry.
        If None, the dates will be constructed as a DateArray with the same
        length as ``data``, starting at ``start_date`` with frequency ``freq``.
    start_date : {Date}, optional
        Date corresponding to the first entry of the data (index 0).
        This parameter must be a valid Date object, and is mandatory if ``dates``
        is None and if ``data`` has a length greater or equal to 1.
    freq : {freq_spec}, optional
        A valid frequency specification, as a string or an integer.
        This parameter is mandatory if ``dates`` is None.
    refperiod : {None, tuple}, optional
        Reference period of the object, as a tuple (starting date, ending date).
        If None, the reference period is set to the whole range of dates.

    Notes
    -----
    * All other parameters that are accepted by the :func:`numpy.ma.array`
      function in the :mod:`numpy.ma` module are also accepted by this function.

    """
    series = time_series(data, dates=dates, start_date=start_date, freq=freq,
                         mask=mask, dtype=dtype, copy=copy,
                         fill_value=fill_value, keep_mask=keep_mask,
                         hard_mask=hard_mask).view(ReferencedSeries)
    series.refperiod = refperiod
    return series
