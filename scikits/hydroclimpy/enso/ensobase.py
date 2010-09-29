# -*- coding: utf-8 -*-
"""
Introduces new objects to associate ENSO information with time series.

A more complete description of the modul is available in the documentation.

.. currentmodule:: scikits.hydroclimpy.enso



Functions
=========




Manipulating ENSO information
-----------------------------

.. autosummary::
   :toctree: generated/

   set_ensoindicator
   group_indices
   grouped_indices_limits
   
   apply_on_phase


:author: Pierre GF Gerard-Marchant
:contact: pierregmcode_AT_gmail_DOT_com
"""
#pylint: disable-msg=E1101
import warnings
warnings.simplefilter("ignore", DeprecationWarning)


import numpy as np
from numpy import ndarray

import numpy.ma as ma
from numpy.ma import nomask, masked
from numpy.ma.mrecords import MaskedRecords

from scipy.stats.mstats import mquantiles, mode

import scikits.timeseries as ts
from scikits.timeseries import TimeSeries, TimeSeriesRecords, DateArray, \
                               adjust_endpoints, fill_missing_dates
from scikits.timeseries.lib.interpolate import forward_fill, backward_fill

import scikits.hydroclimpy.core
from ..core import Cluster, ReferencedSeries, \
                   convert_and_fill, flatten_sequence, get_maskoptions


_c = ts.const


__all__ = ['ENSOIndicator', 'ClimateSeries', 'InvalidENSOError',
           'ClimateRecords',
           'apply_on_phase', 'climate_series', 'climate_records',
           'group_indices', 'grouped_indices_limits',
           'set_ensoindicator']



class InvalidENSOError(Exception):
    """
    Defines an exception for issues with ENSO indicators.
    """
    pass



def code2months(code):
    """
    Transforms a character string to a list of months indices.

    The string should consist of uppercase letters, one letter per month, in
    chronological order.

    Examples
    --------
    >>> code2months('OND')
    array([10,11,12])
    >>> code2months('NDJ')
    array([11,12,1])
    >>> code2months('SOD')
    ValueError

    """
    strg = 'JFMAMJJASONDJFMAMJJASOND'
    if not isinstance(code, basestring):
        raise TypeError("Unrecognized string ('%s')" % code)
    code = code.upper()
    idx = strg.find(code)
    if idx != -1:
        months = np.arange(len(code)) + idx + 1
    else:
        raise ValueError("Unrecognized string ('%s')" % code)
    months[months > 12] -= 12
    return tuple(months)


def months2code(*monthlist):
    """

    Transforms a list of months indices into a three character code

    Examples
    --------
    >>> months2code([10,11,12])
    'OND'

    """
    strg = 'JFMAMJJASONDJFMAMJJASOND'
    return ''.join([str(strg[mm - 1]) for mm in flatten_sequence(monthlist)])



class ENSOIndicator(ReferencedSeries, object):
    """
    Define ENSO indicators objects, as a subclass of 
    :class:`~scikits.hydroclimpy.ReferencedSeries`.

    ENSO indices are calculated with the :meth:`set_indices` method.
    Values of the indicator meeting some thresholds condition are clustered, and
    indices are allocated depending on the size (duration) of these clusters.
    The thresholds and minimum duration are respectively set through the attributes
    :attr:`thresholds` and :attr:`minimum_size` respectively.

    The indices are typically calculated on a monthly basis and cached for convenience
    in the :attr:`_cachedmonthly` private attribute.
    This attribute is a dictionary with keys:

       * ``'indices_monthly'`` for indices allocated with each month.
       * ``'indices_annual'`` for indices allocated with period of 12 months.

    Note that the indices cached in :attr:`_cachedmonthly` have always a monthly
    frequency.

    Monthly indices converted to the frequency of the indicator, with the same
    date range and shape as the indicator, are stored in the private attribute
    :attr:`_cachedcurrent`.


    Parameters
    ----------
    data : array
        Raw indicators of ENSO phases.
    dates : {None, array-like}, optional
        Corresponding dates.
    mask : {nomask, array-like}, optional
        mask
    refperiod : {None, tuple/DateArray}, optional
        Reference period, as a tuple or :class:`DateArray`
        (`starting date`, `ending date`).
    start_date: {None, Date}, optional
        Staring date of the series (if ``dates`` is None).
    ensotag : {'', string}, optional
        Short string acting as 'Memo' field.
    thresholds : {None, tuple}, optional
        Tuple (lower threshold, upper threshold) for the definition of ENSO phases.
        If None, the thresholds are set to the lower and upper quartiles.
    minsize : {None, int}, optional
        Minimum size for the clusters of consecutive data.
    fill_value : {float}, optional
        Code for missing values in the raw data file.

    """
    _optinfo = {}

    def __new__(cls, data, dates=None, mask=nomask,
                freq=None, refperiod=None, start_date=None,
                ensotag='(undefined)', thresholds=None, minsize=None,
                dtype=None, copy=False, fill_value= -9999, **optinfo):
        """

    Parameters
    ----------
    data : array
        Raw indicators of ENSO phases.
    dates : {None, array-like}, optional
        Corresponding dates.
    mask : {nomask, array-like}, optional
        mask
    refperiod : {None, tuple/DateArray}, optional
        Reference period, as a tuple or DateArray (starting date, ending date).
    start_date: {None, Date}, optional
        Staring date of the series (if `dates` is None).
    ensotag : {'', string}, optional
        Short string acting as 'Memo' field.
    thresholds : {None, tuple}, optional
        Tuple (lower threshold, upper threshold) for the definition of ENSO phases.
        If None, the thresholds are set to the lower and upper quartiles.
    minsize : {None, int}, optional
        Minimum size for the clusters of consecutive data.
    fill_value : {float}, optional
        Code for missing values in the raw data file.

    """
        _data = ReferencedSeries(data, dates=dates, mask=mask, freq=freq,
                                 start_date=start_date, refperiod=refperiod,
                                 dtype=dtype, copy=copy, fill_value=fill_value)
        _data = _data.view(cls)
        #
        if thresholds is not None:
            _data.thresholds = thresholds
        else:
            _data._optinfo['thresholds'] = (None, None)
        #
        if minsize is not None:
            _data._optinfo['minimum_size'] = int(minsize)
        #
        refseason = optinfo.get('reference_season', None)
        if refseason is not None:
            refseason = code2months(str(refseason))
            _data._optinfo['reference_season'] = refseason
        _data._optinfo['ensotag'] = ensotag
        return _data

    def _update_from(self, obj):
        self._cachedmonthly = dict(indices_annual=None,
                                   indices_monthly=None,)
        self._cachedmonthly.update(getattr(obj, '_cachedmonthly', {}))
        self._cachedclustered = getattr(obj, '_cachedclustered', None)
        self._cachedcurrent = getattr(obj, '_cachedcurrent', None)
        TimeSeries._update_from(self, obj)


    def __array_finalize__(self, obj):
        ReferencedSeries.__array_finalize__(self, obj)
#        if self._cachedcurrent is not None:
#            self._cachedcurrent.shape = self.shape
        return


    def __getitem__(self, indx):
        result = super(ENSOIndicator, self).__getitem__(indx)
#        if np.isscalar(result):
#            _dates = self._dates
#            date = TimeSeries(_dates, dates=_dates).__getitem__(indx)
#            result = ENSOIndicator(result, dates=date, refperiod=self.refperiod)
        if not np.isscalar(result):
            if self._cachedcurrent is not None:
                result._cachedcurrent = self._cachedcurrent.__getitem__(indx)
        return result

    def __getslice__(self, i, j):
        result = super(ENSOIndicator, self).__getslice__(i, j)
        if self._cachedcurrent is not None:
            result._cachedcurrent = self._cachedcurrent.__getslice__(i, j)
        return result

    def __shortdesc__(self):
        msg = "<ENSO Indicator - %s [%s-%s]>"
        _bstr = DateArray(self.refperiod, freq=self.freq).tostring()
        return msg % (self.ensotag, _bstr[0], _bstr[1])

    def __repr__(self):
        shortdesc = self.__shortdesc__()
        str_array = super(ENSOIndicator, self).__repr__()
        return "%s\n%s" % (shortdesc, str_array)


    # Thresholding ..............................
    def _get_thresholds(self):
        "Retrieves the indicator thresholds for the definition of ENSO phases."
        _optinfo = self.optinfo
        if tuple(_optinfo['thresholds']) == (None, None):
            thresholds = mquantiles(self._series, (.25, .75), axis=None)
            _optinfo['thresholds'] = thresholds
        return _optinfo['thresholds']
    #....
    def _set_thresholds(self, newthresholds=None):#low,high):
        "Defines the indicator thresholds for the definition of ENSO phases."
        _optinfo = self.optinfo
        if (newthresholds is not None):
            try:
                (low, high) = newthresholds
            except:
                raise TypeError("The input thresholds must be given as a "\
                                "sequence (low, high)")
            if low > high:
                (low, high) = (high, low)
            thresholds = (float(low), float(high))
        else:
            thresholds = mquantiles(self._series, (.25, .75), axis=None)
        if thresholds != _optinfo.get('thresholds', None):
            self._cachedmonthly = {}
            self._cachedcurrent = None
        _optinfo['thresholds'] = thresholds
    #....
    thresholds = property(fget=_get_thresholds,
                          fset=_set_thresholds,
    doc="""
    Characteristic thresholds for the definition of ENSO phases.
    By default, the thresholds correspond to the lower and upper quartiles.

    .. note:: Modifying this value resets any cached indices previously computed.
    """)


    def _get_full_year(self):
        "Returns the full_year status flag."
        return self._optinfo.get('full_year', None)
    #....
    def _set_full_year(self, full_year):
        "Sets the full_year status flag to a new boolean value."
        full_year = bool(full_year)
        self.optinfo['full_year'] = bool(full_year)
        return
    #....
    full_year = property(fget=_get_full_year, fset=_set_full_year)


    def _get_ensotag(self):
        "Returns the tag of the ENSO Indicator."
        return self._optinfo.get('ensotag', '(undefined)')
    #....
    def _set_ensotag(self, tag):
        "Sets the ENSO indicator tag to a new string value."
        self._optinfo['ensotag'] = str(tag)
        return
    #....
    ensotag = property(fget=_get_ensotag, fset=_set_ensotag,
                       doc="Tag describing the ENSO indicator.")


    def _get_minimum_size(self):
        return self.optinfo.get('minimum_size', None)
    #....
    def _set_minimum_size(self, minsize):
        minsize = int(minsize)
        optinfo = self.optinfo
        if minsize != optinfo.get('minimum_size', None):
            optinfo['minimum_size'] = minsize
            self._cachedmonthly = {}
            self._cachedcurrent = None
        return
    #...
    minimum_size = property(fget=_get_minimum_size, fset=_set_minimum_size,
    doc="""
    Minimum numbers of consecutive months for the definition of ENSO indices.

    .. note:: Modifying this value resets any cached indices previously computed.
    """)


    def _get_refseason(self):
        return self.optinfo.get('reference_season', None)
    #....
    def _set_refseason(self, refseason):
        optinfo = self.optinfo
        if isinstance(refseason, basestring):
            refseason = code2months(refseason)
        elif isinstance(refseason, (tuple, list)):
            refseason = [int(_) for _ in refseason]
        elif not isinstance(refseason, ndarray):
            refseason = None
        if tuple(refseason) != tuple(optinfo.get('reference_season', ())):
            optinfo['reference_season'] = refseason
            self._cachedmonthly = {}
            self._cachedcurrent = None
        return
    #....
    reference_season = property(fget=_get_refseason, fset=_set_refseason,
    doc="""
    Reference season for the definition of ENSO indices.

    .. note:: Modifying this value resets any cached indices previously computed.
    """)
    refseason = property(fget=_get_refseason, fset=_set_refseason,
                         doc="Alias to :attr:`reference_season`.")

    # Monthly indicators ..................................
    def _clustered_values(self):
        """
    Returns the indices of the ENSO indicator grouped by values.

        """
        if self._cachedclustered is None:
            # Mask data if less than 6 consecutive months with the same index
            # Remmbr: `longenough` is a MASK, so True means "NOT in a cluster"
            _monthly = ma.zeros(self.size, dtype=int)
            _monthly._mask = self._mask
            _monthly[(self >= self.thresholds[-1]).filled(False)] = +1
            _monthly[(self <= self.thresholds[0]).filled(False)] = -1
            self._cachedclustered = Cluster(_monthly.filled(), increment=0)
        return self._cachedclustered


    def _set_monthly_indices(self, minimum_size=None, reference_season=None):
        """
    Sets the ENSO indices per groups of months.
    
    In a first step, each month is first qualified as +1 (-1) if the
    corresponding value of the indicator is greater (less) than the upper (lower)
    threshold, or as 0 otherwise.
    Consecutive values are then grouped according to their value.
    If the size of a group is less than the ``minimum_size`` parameters, 
    the value of each element of this group is set to 0.
    
    In a third step, groups of size larger than ``minimum_size`` are checked
    whether they cover the reference season.
    If they don't, the value of each of their elements is then also set to 0.

    Parameters
    ----------
    minimum_size : {None, int}, optional
        Minimum size for the groups of consecutive values.
        If None, defaults to ``self.minimum_size``.
    reference_season : {None, string or sequence}, optional
        Reference season.
        If None, defaults to ``self.reference_season``.
        
    Notes
    -----
    * If the ``reference_season`` parameter is a string, it should be correspond 
      to a series of consecutive months. For example, 'JFM', 'OND', 'NDJF' are
      all valid reference seasons (Jan. to Mar., Oct. to Dec., Nov. to Feb. 
      respectively), but 'JAM', 'SND' or 'XYZ' are not.
    * If the ``reference_season`` parameter is a sequence, it should be a sequence
      of integers between 1 (Jan.) and 12 (dec.)

        """
        # Get the cluster size ........
        _minimum_size = self.optinfo.get('minimum_size', None)
        minsize = minimum_size or _minimum_size
        if minsize != _minimum_size:
            self.minimum_size = minsize
        # Check the season ............
        _refseason = self.optinfo.get('reference_season', None)
        refseason = reference_season or _refseason
        if tuple(_refseason or ()) != tuple(refseason or ()):
            self.refseason = refseason
        refseason = self.optinfo.get('reference_season', None)
        # Flatten the area (it should be 1D anyway), but save the initial shape
        if self._varshape != ():
            raise NotImplementedError("The ENSO indicator should be 1D!")
        inishape = self.shape
        self.shape = -1
        # Check that the frequency is actually monthly
        if self.freq == 'M':
            klust = self._clustered_values()
            (_mask, _dates, months) = (self._mask, self._dates, self.months)
        else:
            # OK, so force to monthly for now...
            _self = self.convert('M', func=ma.mean)
            _self._optinfo = self._optinfo
            klust = _self._clustered_values()
            (_mask, _dates, months) = (_self._mask, _self._dates, _self.months)
        _indices = np.concatenate(klust.clustered)
        # Find the clusters larger than minimum_size
        _slices = klust.slices[klust.sizes >= minsize]
        # Find the large clusters that fall after min_month
        valid = np.zeros(months.shape, bool)
        if refseason is not None:
            for s in _slices:
                valid[s] = np.logical_and.reduce([(months[s] == m).any()
                                                  for m in refseason])
        else:
            for s in _slices:
                valid[s] = True
        _indices[np.logical_not(valid)] = 0
        _indices = _indices.view(TimeSeries)
        _indices._mask = _mask
        _indices._dates = _dates
        # Cache the results............
        self._cachedmonthly['indices_monthly'] = _indices
        # And restore the shape
        self.shape = inishape
        return _indices


    def _set_annual_indices(self, minimum_size=None, reference_season=None):
        """
    Sets the ENSO indices per periods of 12 months, starting at the first 
    month of the reference season if any, otherwise at October.

    The same steps are followed as for :meth:`set_monthly_indices`.

    Parameters
    ----------
    minimum_size : {None, int}, optional
        Minimum size for the groups of consecutive values.
        If None, defaults to :attr:`minimum_size`.
    reference_season : {None, string or sequence}, optional
        Reference season.
        If None, defaults to :attr:`reference_season`.

    See Also
    --------
    :meth:`set_monthly_indices`
        Sets the ENSO indices for each month.

    """
        # Get the monthly indices .....
        _monthly = self.set_monthly_indices(minimum_size=minimum_size,
                                            reference_season=reference_season)
        # Make sure we reset the full_year flag to True (we lost it w/ set_monthly
        self.full_year = True
        # Get the annual indices
        refseason = self.refseason
        if refseason:
            _annual = _monthly[self.months == refseason[0]]
            refseason = months2code(refseason)
        else:
            _annual = _monthly[self.months == 10]
        _annual = adjust_endpoints(forward_fill(fill_missing_dates(_annual)),
                                   self._dates[0], self._dates[-1])
        # Cache the results ...........
        self._cachedmonthly['indices_annual'] = _annual
        return _annual


    def _fix_cachedcurrent(self, ensoindices):
        """
    Private function to convert cached ENSO indices.
    This function is used to ensure that:
    * the frequency of the indices matches the frequency of the indicator.
    * the date range of the indices matches the date range of the indicator.
    * the shape of the indices matches the shape of the indicator.
        """
        _currentfreq = ensoindices._dates.freq
        _freq = self._dates.freq
        # Check the frequency, and convert if needed
        if _currentfreq > _freq:
            # To a lower frequency ('Q'/'A')
            if self.ndim == 2:
                conversionfunc = None
            else:
                conversionfunc = lambda x: mode(x)[0].squeeze()
            ensoindices = ensoindices.convert(_freq, func=conversionfunc)
        elif _currentfreq < _freq:
            # To a higher frequency (eg, 'D')
            ensoindices = backward_fill(ensoindices.convert(_freq))
        # Reset to the original frequency if needed...
        (start, end) = self._dates.flat[[0, -1]]
        if tuple(ensoindices._dates.flat[[0, -1]]) != (start, end):
            ensoindices = ts.adjust_endpoints(ensoindices,
                                              start_date=start, end_date=end)
        # Reset the shape of the indices
        if ensoindices.shape != self.shape:
            ensoindices.shape = self.shape
        return ensoindices


    def set_indices(self, full_year=None, reference_season=None,
                    minimum_size=None, lag=0):
        """
    Sets the ENSO indices for the ENSO indicator.

    Results are cached for convenience: if the method is called without arguments,
    the monthly indices are returned, provided they had been already computed.

    Parameters
    ----------
    full_year : {None, False, True}, optional
        Whether the ENSO indices are set for a period of 12 months (True) or by
        valid episodes (False).
        If None, defaults to :attr:`full_year`.
    reference_season : {None, string or sequence}, optional
        Reference season.
        If None, defaults to :attr:`reference_season`.
    minimum_size : {None, int}, optional
        Minimum size for the groups of consecutive values.
        If None, defaults to :attr:`minimum_size`.
    lag : {integer}, optional
        Number of months of lag for ENSO indices. For example,
        if lag=6, the ENSO phase starting in Oct. 2001 is applied starting on 
        Apr. 2002.
    
    Notes
    -----
    * If ``full_year`` is False, this method calls :meth:`set_monthly_indices`.
      Otherwise, it calls :meth:`set_annual_indices`.
    * Results are cached for convenience.
      The cached results are always with a monthly frequency, and converted to
      the current frequency ``freq`` as needed, with the command
      
         >>> scikits.timeseries.lib.backward_fill(_cached.convert(freq, func=ma.mean))
      
    * The cached indices have a monthly frequency by default.
      When the ENSOIndicator is converted to a lower frequency with a given
      conversion function, the ENSO indices are converted using
      :func:`scipy.stats.mstats.mode` as conversion function.
    
    """
        optinfo = self._optinfo
        # Check the full_year flag: reset _cachedcurrent if it changed....
        if full_year is None:
            full_year = optinfo.get('full_year', False)
        else:
            full_year = bool(full_year)
            if optinfo.get('full_year', False) != full_year:
                optinfo['full_year'] = full_year
                self._cachedcurrent = None
        # Check the reference season .....................................
        if reference_season is None:
            reference_season = optinfo.get('reference_season', None)
        else:
            optreference_season = optinfo.get('reference_season', ())
            if tuple(reference_season) != tuple(optreference_season):
                self.refseason = reference_season
        # Check the minimum_size..........................................
        if minimum_size is None:
            minimum_size = optinfo.get("minimum_size", None)
        else:
            self.minimum_size = minimum_size
        # Check the current lag ..........................................
        if lag != optinfo.get('current_lag', 0):
            optinfo['current_lag'] = lag
            self._cachedcurrent = None
        # Define what we need to get......................................
        if full_year:
            _cached = self._cachedmonthly.get('indices_annual', None)
            func = self._set_annual_indices
        else:
            _cached = self._cachedmonthly.get('indices_monthly', None)
            func = self._set_monthly_indices
        #...............
        if _cached is None:
            if (minimum_size is None) and (reference_season is None):
                err_msg = "The ENSO indices have not been initialized!\n"\
                          "Please define a reference season with `self.refseason`"\
                          " and a minimum cluster size with `self.minsize`"
                raise ValueError(err_msg)
            _cached = func(reference_season=reference_season,
                           minimum_size=minimum_size)
            # Make sure we didn't zap full_year
            optinfo['full_year'] = full_year
        # Check whether we need to take a lag into account
        if lag:
            _cached = ts.tshift(func(reference_season=reference_season,
                                     minimum_size=minimum_size),
                                - lag, copy=False)
        # Fix the frequency/date range/shape of _cached
        _cached = self._fix_cachedcurrent(_cached)
        self._cachedcurrent = _cached
        return _cached

    def get_indices(self):
        """
    Returns the indices currently cached in :attr:'_cachedcurrent'.

    If the ENSO indices are not set (:attr:'_cachedcurrent' is None), they are
    set with the :meth:'set_indices' method.
        """
        if self._cachedcurrent is not None:
            self._cachedcurrent = self._fix_cachedcurrent(self._cachedcurrent)
            return self._cachedcurrent
        return self.set_indices()

    ensoindices = property(fget=get_indices, fset=set_indices,
    doc="""
    ENSO indices :
    
    * +1: El Ni\~no
    *  0: Neutral
    * -1: La Ni\~na
    """)
    indices = property(fget=get_indices, fset=set_indices,
                       doc="Alias to :attr:`ensoindices`")


    def set_monthly_indices(self, reference_season=None,
                            minimum_size=None, lag=0):
        """
    Sets the ENSO indices per groups of months.
    
    In a first step, each month is first qualified as +1 (-1) if the
    corresponding value of the indicator is greater (less) than the upper (lower)
    threshold, or as 0 otherwise.
    Consecutive values are then grouped according to their value.
    If the size of a group is less than the ``minimum_size`` parameters, 
    the value of each element of this group is set to 0.
    
    In a third step, groups of size larger than ``minimum_size`` are checked
    whether they cover the reference season.
    If they don't, the value of each of their elements is then also set to 0.

    Parameters
    ----------
    minimum_size : {None, int}, optional
        Minimum size for the groups of consecutive values.
        If None, defaults to ``self.minimum_size``.
    reference_season : {None, string or sequence}, optional
        Reference season.
        If None, defaults to ``self.reference_season``.
    lag : {integer}, optional
        Number of months of lag for ENSO indices. For example,
        if lag=6, the ENSO phase starting in Oct. 2001 is applied starting on 
        Apr. 2002.
        
    Notes
    -----
    * If the ``reference_season`` parameter is a string, it should be correspond 
      to a series of consecutive months. For example, 'JFM', 'OND', 'NDJF' are
      all valid reference seasons (Jan. to Mar., Oct. to Dec., Nov. to Feb. 
      respectively), but 'JAM', 'SND' or 'XYZ' are not.
    * If the ``reference_season`` parameter is a sequence, it should be a sequence
      of integers between 1 (Jan.) and 12 (dec.)

        """
        return self.set_indices(full_year=False,
                                reference_season=reference_season,
                                minimum_size=minimum_size, lag=lag)


    def set_annual_indices(self, minimum_size=None, reference_season=None, lag=0):
        """
    Sets the ENSO indices per periods of 12 months, starting at the first 
    month of the reference season if any, otherwise at October.

    The same steps are followed as for :meth:`set_monthly_indices`.

    Parameters
    ----------
    minimum_size : {None, int}, optional
        Minimum size for the groups of consecutive values.
        If None, defaults to :attr:`minimum_size`.
    reference_season : {None, string or sequence}, optional
        Reference season.
        If None, defaults to :attr:`self.reference_season`.
    lag : {integer}, optional
        Number of months of lag for ENSO indices. For example,
        if lag=6, the ENSO phase starting in Oct. 2001 is applied starting on 
        Apr. 2002.

    See Also
    --------
    :meth:`set_monthly_indices`
        Sets the ENSO indices for each month.

    """
        return self.set_indices(full_year=True,
                                reference_season=reference_season,
                                minimum_size=minimum_size, lag=lag)


    def convert(self, freq, func=None, position='END', *args, **kwargs):
        """
    %s

    Notes
    -----
    * If the ENSOIndicator has a :attr:'_cachedcurrent' attribute which is not
      None, this latter is fixed with :meth:'_fix_cachedcurrent'.
        """ % TimeSeries.convert.__doc__
        result = super(ENSOIndicator, self).convert(freq,
                                                    func=func,
                                                    position=position,
                                                    *args, **kwargs)
        _cached = self._cachedcurrent
        if _cached is not None:
            result._cachedcurrent = result._fix_cachedcurrent(_cached)
        return result

#    # Strengths .................................
#    def _set_annual_strength(self):
#        "Returns the strength of the old-style ENSO events."
#        # Get the month limits (three months before, 6 months after)
#        min_m = self.optinfo['starting_month'] - 3
#        max_m = (self.optinfo['starting_month'] + 5) % 12
#        # Get the corresponding year limits
#        (start_y, end_y) = self.years[[0,-1]]
#        (first_m, last_m) = self.months[[0,-1]]
#        if first_m <= 7:
#            start_y = start_y - 1
#        if last_m >= 3:
#            end_y = end_y + 1
#        # Extend the series from the min_month of first year to the max_month of last year
#        start_date = Date(freq='M', year=start_y, month=min_m)
#        end_date = Date(freq='M', year=end_y, month=max_m)
#        J2M = ts.adjust_endpoints(self, start_date=start_date, end_date=end_date,)
#        # Reshape to a (-1,9) array
#        months = J2M.months
#        J2M = J2M[(months >= min_m) | (months <= max_m)].reshape(-1,9)
#        # Define some shortcuts ......
#        j2m = J2M._series
#        (tlow, thigh) = self.thresholds
#        # Get the annual indices ......
#        _annual = np.zeros(len(j2m), dtype=int)
#        repvalue = j2m[:,4].filled(0)
#        _annual[repvalue >= thigh] = 1
#        _annual[repvalue <= tlow] = -1       
#        (high, med, low) = ((_annual > 0), (_annual == 0), (_annual < 0)) 
#        # Compute the averages, following Hanley et al. (2003)
#        averaged = ma.array(np.zeros(len(j2m), dtype=float), nomask)
#        averaged[high] = ma.masked_less(j2m[[high]],0).mean(-1)
#        averaged[low] = ma.masked_greater(j2m[[low]],0).mean(-1)
#        averaged[med] = j2m[[med]].mean(-1)
#        # Compute the strength
#        _strength = TimeSeries(averaged, copy=True, mask=nomask,
#                               dates=J2M._dates[:,4], freq='M')
#        _strength[(averaged < 0).filled(False)] /= abs(float(tlow))
#        _strength[(averaged >= 0).filled(False)] /= abs(float(thigh))  
#        self._cachedinfo['strength_annual'] = _strength
#        return _strength
#    #................................................................
#    def strength(self):        
#        """
#    Returns the strength of the new-style ENSO events.
#
#    """
#        # Get the monthly indices .....
#        indices = self.indices()
#        values = self._series
#        # Reclusterize
#        clustered_indices = Cluster(indices.filled(0), 0)
#        c_slices = clustered_indices.slices
#        c_sizes = clustered_indices.sizes
#        # Get the averages per cluster .
#        c_averages = np.fromiter((values[k].mean() for k in c_slices), 
#                                 dtype=float)
#        c_values = np.fromiter((k[0] for k in clustered_indices.clustered), 
#                               dtype=int)
#        c_strength = np.concatenate([v.repeat(s) 
#                                     for (v,s) in zip(c_averages,c_sizes)])
#        # Fix the dates
#        _dates = self._dates
#        c_strength = c_strength.view(ClimateSeries)
#        _newdatesvalues = [DateArray(_dates[k]).tovalue() for k in c_slices]
#        c_strength._dates = DateArray(np.concatenate(_newdatesvalues), freq='M')
#        c_strength._cachedinfo = self._cachedinfo
#        c_strength.ensoindices = self.ensoindices
#        self._cachedinfo['strength_monthly'] = c_strength
#        return c_strength



class _csmethod(object):
    """
    Wrappers to TimeSeries methods.
    """
    __extradoc__ = """
    Notes
    -----
    When called, returns a new :class:`ClimateSeries` object, with the new series 
    the result of the method applied on the original series.
    If ``onindices`` is True, the same operation is performed on the ``indices``.

    """
    #
    def __init__ (self, methodname, onindices=True):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self.__name__ = methodname
        self._onindices = onindices
        if __doc__:
            self.__doc__ = getattr(TimeSeries, methodname).__doc__ + \
                           _csmethod.__extradoc__
        self.obj = None
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, *args, **kwargs):
        "Execute the call behavior."
        _name = self.__name__
        instance = self.obj
        func_series = getattr(super(ClimateSeries, instance), _name)
        result = func_series(*args, **kwargs)
        if instance.ensoindicator is not None and self._onindices:
            func_enso = getattr(instance.ensoindicator, _name)
            result.ensoindicator = func_enso(*args, **kwargs)
#        result._update_from(instance)
        return result


class ClimateSeries(ReferencedSeries, object):
    """
    TimeSeries objects with ENSO information.
    
    Parameters
    ----------
    data : array-like
        Data
    dates : {None, sequence}, optional
        Corresponding dates
    mask : {nomask, sequence}, optional
        Mask
    ensoindicator : {None, ENSOIndicator}, optional
        ENSO Indicator.
    refperiod : {None, tuple}, optional
        Reference period, as a tuple (starting date, ending date).
        If None, the reference period is set to the whole range of dates.
    freq : {None, string or integer}, optional
        A valid frequency specification
    start_date : {None, Date}, optional
        Starting date of the series.
        This parameter is only useful if ``dates`` is None.
    dtype : {None, dtype}, optional
        Datatype of the series.
        If None, the dtype of ``data`` is used.
    copy : {False, True}, optional
        Whether copy the data (True) or just link to it (False).

    See Also
    --------
    :class:`ENSOIndicator`
        Class to manipulate ENSO information.

    """

    def __new__(cls, data, dates=None, mask=nomask, ensoindicator=None,
                refperiod=None, freq=None, start_date=None, length=None,
                dtype=None, copy=False, **options):
        (tsoptions, options) = get_maskoptions(**options)
        tsoptions.update(dict(copy=copy, dtype=dtype))
        _data = ReferencedSeries(data, dates=dates, mask=mask, freq=freq,
                                 refperiod=refperiod, start_date=start_date,
                                 length=length, **tsoptions).view(cls)
        _data.ensoindicator = ensoindicator
        return _data


    def __array_finalize__(self, obj):
        ReferencedSeries.__array_finalize__(self, obj)
        ensoindicator = self.ensoindicator
        if ensoindicator is not None:
            ensoindicator.shape = self.shape

    def __array_wrap__(self, obj, context=None):
        result = super(ClimateSeries, self).__array_wrap__(obj, context)
        result.ensoindicator = self.ensoindicator
        return result

    def __getitem__(self, indx):
        obj = super(ClimateSeries, self).__getitem__(indx)
        if (self.ensoindicator is not None) and isinstance(obj, ClimateSeries):
            if isinstance(indx, basestring) and (indx in self.dtype.names):
                obj.ensoindicator = self.ensoindicator
            else:
                obj.ensoindicator = self.ensoindicator[indx]
        return obj


    def __repr__(self):
        """
        Calculates the repr representation, using masked for fill if it is
enabled. Otherwise fill with fill value.
"""
        desc = """\
climate_series(
 %(data)s,
               dates =
 %(time)s,
               freq  = %(freq)s)
"""
        desc_short = """\
climate_series(%(data)s,
              dates = %(time)s,
              freq  = %(freq)s)
"""
        if np.size(self._dates) > 2 and self.is_valid():
            _dates = self._dates[[0, -1]]
            timestr = "[%s ... %s]" % (str(_dates[0]), str(_dates[-1]))
        else:
            timestr = str(self.dates)

        if self.ndim <= 1:
            return desc_short % {'data': str(self._series),
                                 'time': timestr,
                                 'freq': self.freqstr, }
        return desc % {'data': str(self._series),
                       'time': timestr,
                       'freq': self.freqstr, }


    #--- Array methods
    ravel = _csmethod('ravel')
#    convert = _csmethod('convert')

    def reshape(self, *args, **kwargs):
        result = super(ClimateSeries, self).reshape(*args, **kwargs)
        ensoi = result.ensoindicator
        if ensoi is not None:
            result.ensoindicator = ensoi.reshape(*args, **kwargs)
        return result

    def convert(self, freq, func=None, position='END', *args, **kwargs):
        """
    (this docstring should be overwritten)
        """
        result = ReferencedSeries.convert(self, freq, func=func,
                                          position=position, *args, **kwargs)
        if self.ensoindicator is not None:
            ensoindicator = self.ensoindicator.convert(freq, func=func,
                                                       position=position,
                                                       *args, **kwargs)
            result.ensoindicator = ensoindicator
        return result

    #--- Specific methods

    def _get_ensoindicator(self):
        return self._optinfo.get('ensoindicator', None)
    #
    def set_ensoindicator(self, indicator):
        if indicator is not None:
            if not isinstance(indicator, ENSOIndicator):
                raise TypeError("The input is not a valid ENSOIndicator ! "\
                                "(got '%s' instead)" % type(indicator))
            # Reset the frequency
            (sfreq, ifreq) = (self._dates.freq, indicator.freq)
            if ifreq != sfreq:
                if sfreq < ifreq:
                    indicator = indicator.convert(sfreq, func=ma.mean)
                elif sfreq > ifreq:
                    indicator = backward_fill(indicator.convert(sfreq))
            # Set the ENSO indices to the default (if we don't have any...)
            if indicator._cachedcurrent is None:
                indicator.set_indices()
            # Reset the dates
            dates = self._dates
            if dates.isfull():
                if len(dates) > 0:
                    (start, end) = dates.flat[[0, -1]]
                    if tuple(indicator._dates.flat[[0, -1]]) != (start, end):
                        indicator = ts.adjust_endpoints(indicator, start, end)
                else:
                    indicator = None
            else:
                indicator = indicator[self._dates]
        self._optinfo['ensoindicator'] = indicator
    #
    ensoindicator = property(fget=_get_ensoindicator, fset=set_ensoindicator)


    def _get_ensoindices(self):
        ensoindicator = self.optinfo.get('ensoindicator', None)
        if ensoindicator is None:
            errmsg = "The ENSO indicator hasn't been initialized yet!"
            raise InvalidENSOError(errmsg)
        return ensoindicator.ensoindices
    #
    def set_ensoindices(self, **parameters):
        """
    Sets the ENSO indices from the current ENSO indicator.

    Returns
    -------
    ensoindices : TimeSeries
        Series of ENSO indices : +1 for El Ni\~no, 0 for Neutral, -1 for La Ni\~na.

        """
        ensoindicator = self.optinfo.get('ensoindicator', None)
        if ensoindicator is None:
            errmsg = "The ENSO indicator hasn't been initialized yet!"
            raise InvalidENSOError(errmsg)
        return ensoindicator.set_indices(**parameters)
    #
    ensoindices = property(fget=_get_ensoindices, fset=set_ensoindices,
                           doc="ENSO indices : "\
                               "(+1: El Ni\~no; 0:Neutral; -1:La Ni\~na)")
    #
    def _phase(self, index):
        """
    Gets the subset of the series falling during the given episode.

    Parameters
    ----------
    index : {-1,0,+1}
        Index corresponding to the desired ENSO phase:

        * -1 : La Ni\~na
        *  0 : Neutral
        * +1 : El Ni\~no

        """
        indices = self.ensoindices
        result = self.copy().view(TimeSeries)
        result[(indices != index).filled(True)] = masked
        return result
    #
    cold = property(fget=lambda self:self._phase(-1),
                    doc="Data falling during La Ni\~na episodes, as a TimeSeries.")
    neutral = property(fget=lambda self:self._phase(0),
                       doc="Data falling during Neutral episodes, as a TimeSeries.")
    warm = property(fget=lambda self:self._phase(+1),
                    doc="Data falling during El Ni\~no episodes, as a TimeSeries.")

climate_series = ClimateSeries


def set_ensoindicator(series, indicator):
    """

    Transforms a :class:`~scikits.timeseries.TimeSeries` into a 
    :class:`~scikits.hydroclimpy.enso.ClimateSeries`, by setting an ENSO indicator.

    Parameters
    ----------
    series : TimeSeries
        The :class:`~scikits.timeseries.TimeSeries` object to transform.
    indicator : ENSOIndicator
        The :class:`~scikits.hydroclimpy.enso.ENSOIndicator` object we want to set.

    Raises
    ------
    TypeError
        If ``series`` is not a :class:`~scikits.timeseries.TimeSeries`,
        or ``indicator`` not an :class:`~scikits.hydroclimpy.enso.ENSOIndicator`.

    Returns
    -------
    climateseries
        A :class:`~scikits.hydroclimpy.enso.ClimateSeries` object.

    """
    if not isinstance(indicator, ENSOIndicator):
        return TypeError("Argument is not a valid ENSO indicator !")
    if not isinstance(series, TimeSeries):
        raise TypeError("Time series should be a valid TimeSeries!")
    _series = series.view(ClimateSeries)
    _series.ensoindicator = indicator
    _series.refperiod = getattr(series, 'refperiod', None)
    return _series


def group_indices(ensoindices):
    grouped = Cluster(ensoindices._data, 0)
    d_indices = {}
    for k in (-1, 0, +1):
        d_indices[k] = grouped.slices[[i
                                       for (i, v) in enumerate(grouped.clustered)
                                       if v[0] == k]]
    return d_indices


def grouped_indices_limits(dates, ensoindices, **optinfo):
    dict_indices = group_indices(ensoindices)
    dict_limits = {}
    _dates = np.empty(len(dates) + 1, int)
    _dates[:-1] = dates.tovalue()
    _dates[-1] = _dates[-2]
    for (k, idx) in dict_indices.iteritems():
        dict_limits[k] = [(s, e) for (s, e) in zip(_dates[[i.start for i in idx]],
                                                 _dates[[i.stop for i in idx]])]
    return dict_limits



class ClimateRecords(TimeSeriesRecords, ClimateSeries):
    """
    TimeRecords object with climate information.
    As for any :class:`~numpy.recarray` object, individual fields can be accessed
    as items or directly as attributes.

    Parameters
    ----------
    data : array-like
        Data
    dates : {None, sequence}, optional
        Corresponding dates
    mask : {nomask, sequence}, optional
        Mask
    ensoindicator : {None, ENSOIndicator}, optional
        ENSO Indicator.
    refperiod : {None, tuple}, optional
        Reference period, as a tuple (starting date, ending date).
        If None, the reference period is set to the whole range of dates.
    freq : {None, string or integer}, optional
        A valid frequency specification
    start_date : {None, Date}, optional
        Starting date of the series.
        This parameter is only useful if ``dates`` is None.
    dtype : {None, dtype}, optional
        Datatype of the series.
        If None, the dtype of ``data`` is used.
    copy : {False, True}, optional
        Whether copy the data (True) or just link to it (False).

    See Also
    --------
    :class:`ClimateSeries`
        :class:`~scikits.timeseries.TimeSeries` with climate information.

    :class:`ENSOIndicator`
        Class to manipulate ENSO information.

    :class:`TimeSeriesRecords`
        :class:`~numpy.ma.mrecords.MaskedRecords` with support for time-indexing.

    """
    #
    def __array_finalize__(self, obj):
        TimeSeriesRecords.__array_finalize__(self, obj)
        self._ensoindicator = getattr(obj, '_ensoindicator', None)
    #
    def __getattribute__(self, attr):
        getattribute = MaskedRecords.__getattribute__
        _dict = getattribute(self, '__dict__')
        if attr in (ndarray.__getattribute__(self, 'dtype').names or []):
            obj = getattribute(self, attr).view(ClimateSeries)
            obj._dates = _dict['_dates']
            obj._ensoindicator = getattribute(self, '_ensoindicator')
            obj._optinfo = getattribute(self, '_optinfo')
            return obj
        return getattribute(self, attr)
    #
    def __getitem__(self, item):
        obj = TimeSeriesRecords.__getitem__(self, item)
        getattribute = MaskedRecords.__getattribute__
        _dict = getattribute(self, '__dict__')
        if isinstance(obj, TimeSeries):
            obj = obj.view(ClimateSeries)
            obj._ensoindicator = getattribute(self, '_ensoindicator')
            obj._optinfo = getattribute(self, '_optinfo')
        return obj


def climate_records(data, dates=None, mask=nomask, ensoindicator=None,
                    refperiod=None, freq=None, start_date=None, length=None,
                    dtype=None, copy=False, **options):
    """
    Creates a :class:`~scikits.hydroclimpy.enso.ClimateSeriesRecords` object.

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
    mask : {nomask, sequence}, optional
        Mask.  Must be convertible to an array of booleans with
        the same shape as data: True indicates a masked (eg.,
        invalid) data.
    dtype : {dtype}, optional
        Data type of the output.
        If dtype is None, the type of the data argument (`data.dtype`) is used.
        If dtype is not None and different from `data.dtype`, a copy is performed.
    copy : {False, True}, optional
        Whether to copy the input data (True), or to use a reference instead.
        Note: data are NOT copied by default.
    fill_value : {var}, optional
        Value used to fill in the masked values when necessary.
        If None, a default based on the datatype is used.
    keep_mask : {True, boolean}, optional
        Whether to combine mask with the mask of the input data,
        if any (True), or to use only mask for the output (False).
    hard_mask : {False, boolean}, optional
        Whether to use a hard mask or not.
        With a hard mask, masked values cannot be unmasked.


    Notes
    -----
    * All other parameters that are accepted by the :func:`numpy.ma.array`
      function in the :mod:`numpy.ma` module are also accepted by this function.

    * The date portion of the time series must be specified in one of the
      following ways:

       + specify a TimeSeries object for the ``data`` parameter.

       + pass a :class:`~scikits.timeseries.DateArray` 
         for the ``dates`` parameter.

       + specify a ``start_date`` (a continuous 
         :class:`~scikits.timeseries.DateArray` will be automatically
         constructed for the dates portion).

       + specify just a frequency (for :class:`~scikits.timeseries.TimeSeries` 
         of size zero).

    """
    series = climate_series(data, dates=dates, mask=mask,
                             ensoindicator=ensoindicator, refperiod=refperiod,
                             freq=freq, start_date=start_date,
                             dtype=dtype, copy=copy, **options)
    return series.view(ClimateRecords)

#-------------------------------------------------------------------------------

def apply_on_phase(series, func, *args, **kwargs):
    u"""
    
    Applies the function `func` to each phase
    of the :class:`~scikits.hydroclimpy.enso.ClimateSeries` `series`.
    
    Parameters
    ----------
    series : ClimateSeries object
        Input climate data. The ENSO indices must have been defined by setting
        a :class:`~scikits.hydroclimpy.enso.ENSOIndicator` object.
    func : function
        Function to apply.
    args : {None, sequence}, optional
        Mandatory arguments of the function ``func``.
    kwargs : {None, dictionary}, optional
        Optional parameters of the function ``func``.

    Returns
    -------
    result
        A structured :class:`~numpy.ma.MaskedArray` of results, with a four fields:

        * ``global`` for the result of the function on the whole series.
        * ``cold`` for the result of the function on La Niña episodes
        * ``neutral`` for the result of the function on Neutral episodes
        * ``warm`` for the result of the function on El Niño episodes.

    """
    if series.ensoindices is None:
        raise AttributeError("ENSO indices should be defined for the input series!")
    #
    _glob = np.asarray(func(series, *args, **kwargs))
    names = ('cold', 'neutral', 'warm', 'global')
    result = ma.empty(_glob.shape, dtype=[(_, _glob.dtype) for _ in names])
    result['global'] = _glob
    for attr in names[:-1]:
        result[attr] = func(getattr(series, attr), *args, **kwargs)
    return result
