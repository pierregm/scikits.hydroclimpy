"""
.. currentmodule:: scikits.hydroclimpy.plotlib.flowplots

This module introduces some classes and functions to help plotting hydrographs
and flow duration curves.
When available, ENSO information can be taken into account.


Plotting hydrographs
====================

:class:`Hydrograph` class
-------------------------

.. autoclass:: Hydrograph
   :show-inheritance:
   :members:


:func:`hydrograph` constructor function
---------------------------------------

.. autofunction:: hydrograph



Plotting flow duration curves (FDCs)
====================================

.. autofunction:: plot_fdc

"""
__author__ = "Pierre GF Gerard-Marchant"


import numpy as np
import numpy.ma as ma

from scipy.stats.distributions import norm

import matplotlib
from matplotlib.pyplot import figure, gca, setp, rcParamsDefault
from matplotlib.collections import LineCollection


from scikits.timeseries import TimeSeries
from scikits.timeseries.lib.plotlib import TimeSeriesFigure, TimeSeriesPlot
import scikits.timeseries.const as _c

from scikits.hydroclimpy import FR_ANNSTART
from scikits.hydroclimpy.enso import ClimateSeries, InvalidENSOError

import mpl_addons
from mpl_addons import set_normal_ticks, set_normal_limits



#####--------------------------------------------------------------------------
#---- --- Hydrograph ---
#####--------------------------------------------------------------------------
# """
# Defines a subclass of
# :class:`~scikits.timeseries.lib.plotlib.TimeSeriesFigure`
# designed to plot a hydrograph along with one hyetograph.
# 
# The figure is dived into two :class:`~scikits.timeseries.lib.plotlib.TimeSeriesPlot`.
# The top plot corresponds to the hyetograph (the time evolution of the 
# precipitation).
# The bottom plot corresponds to the hydrograph (the time evolution of the
# streamflow).
# 
# The two subplots share the same x-axis.
# Therefore, any modification on one subplot will be reflected on the other
# subplot.
# 
# 
# Parameters
# ----------
# hyetodata : {None, TimeSeries}
#     Time series of precipitation.
# hydrodata : {None, TimeSeries}
#     Time series of streamflow.
# figsize : {None, tuple}
#     Size of the figure, as a tuple (width, height) in inches.
#     If None, defaults to rc figure.figsize.
# dpi : {None, int}, optional
#     Resolution in dots per inches. 
#     If None, defaults to rc figure.dpi
# facecolor : {None, string}, optional
#     Background color.
#     If None, defaults to rc figure.facecolor.
# edgecolor : {None, string}, optional
#     Border color.
#     If None, defaults to rc figure.edgecolor.
# linewidth : {float, None}
#     Width of the patch edge line.
# frameon : {True, False}
#     Whether to draw the frame around the figure.
# 
# Attributes
# ----------
# hydro : ClimateSeriesPlot
#     Handle of the subplot corresponding to the hydrograph (bottom subplot).
# hydrodata : TimeSeries
#     Stream flows time series.
# hyeto : ClimateSeriesPlot
#     Handle of the subplot corresponding to the hyetograph (top subplot).
# hyetodata : TimeSeries
#     Precipitation time series.
# 
# 
# Raises
# ------
# TypeError
#     If the hyetodata or the hydrodata are not valid
#     :class:`~scikits.timeseries.TimeSeries` objects.
# ValueError
#     If the hyetodata and the hydrodata do not have the same frequency.
# 
# """


class Hydrograph(TimeSeriesFigure):
    """
    Defines a subclass of
    :class:`~scikits.timeseries.lib.plotlib.TimeSeriesFigure`
    designed to plot a hydrograph along with one hyetograph.
    
    The figure is dived into two :class:`~scikits.timeseries.lib.plotlib.TimeSeriesPlot`.
    The top plot corresponds to the hyetograph (the time evolution of the 
    precipitation).
    The bottom plot corresponds to the hydrograph (the time evolution of the
    streamflow).
    
    The two subplots share the same x-axis.
    Therefore, any modification on one subplot will be reflected on the other
    subplot.
    
    
    Parameters
    ----------
    hyetodata : {None, TimeSeries}
        Time series of precipitation.
    hydrodata : {None, TimeSeries}
        Time series of streamflow.
    figsize : {None, tuple}
        Size of the figure, as a tuple (width, height) in inches.
        If None, defaults to rc figure.figsize.
    dpi : {None, int}, optional
        Resolution in dots per inches. 
        If None, defaults to rc figure.dpi
    facecolor : {None, string}, optional
        Background color.
        If None, defaults to rc figure.facecolor.
    edgecolor : {None, string}, optional
        Border color.
        If None, defaults to rc figure.edgecolor.
    linewidth : {float, None}
        Width of the patch edge line.
    frameon : {True, False}
        Whether to draw the frame around the figure.


    Raises
    ------
    TypeError
        If the hyetodata or the hydrodata are not valid
        :class:`~scikits.timeseries.TimeSeries` objects.
    ValueError
        If the hyetodata and the hydrodata do not have the same frequency.

    """
    #
    def __init__(self,
                 hyetodata=None,
                 hydrodata=None,
                 figsize=None,
                 dpi=None,
                 facecolor=None,
                 edgecolor=None,
                 linewidth=1.0, # the default linewidth of the frame
                 frameon=True, # whether or not to draw the figure frame
                 subplotpars=None, # default to rc
                 ):
        # Initial checks
        self._hyetodata = hyetodata
        self._hydrodata = hydrodata
        #
        kwargs = dict(figsize=figsize, dpi=dpi, facecolor=facecolor,
                      edgecolor=edgecolor, frameon=frameon, linewidth=linewidth,
                      subplotpars=subplotpars, series=None)
        TimeSeriesFigure.__init__(self, **kwargs)
        # Add axes
        self.hyeto = self.add_tsplot(211, series=hyetodata)
        self.hyeto.set_position([0.07, 0.67, 0.90, 0.23])
        self.hydro = self.add_tsplot(212, series=hydrodata, sharex=self.hyeto)
        self.hydro.set_position([0.07, 0.12, 0.90, 0.53])


    def _get_hyetodata(self):
        return self._hyetodata
    #
    def _set_hyetodata(self, hyetodata):
        if hyetodata is not None:
            if not isinstance(hyetodata, TimeSeries):
                raise TypeError("The precipitation series should be a valid "\
                                "TimeSeries object ! "\
                                "(got <%s> instead" % type(hyetodata))
            _hydrodata = self._hydrodata
            if _hydrodata is not None:
                if hyetodata._dates.freq != _hydrodata._dates.freq:
                    raise ValueError("Incompatible frequencies ! "\
                                     "The precipitation series should have a "\
                                     "frequency of '%s', got '%s' instead " % \
                                     (_hydrodata.freqstr, hyetodata.freqstr))
            self._hyetodata = hyetodata
            self.hyeto._series = hyetodata
    #
    hyetodata = property(fget=_get_hyetodata,
                         fset=_set_hyetodata,
                         doc="Precipitation series.")


    def _get_hydrodata(self):
        return self._hydrodata
    #
    def _set_hydrodata(self, hydrodata):
        if hydrodata is not None:
            if not isinstance(hydrodata, TimeSeries):
                raise TypeError("The streamflow series should be a valid "\
                                "TimeSeries object ! "\
                                "(got <%s> instead" % type(hydrodata))
            _hyetodata = self._hyetodata
            if _hyetodata is not None:
                if hydrodata._dates.freq != _hyetodata._dates.freq:
                    raise ValueError("Incompatible frequencies ! "\
                                     "The streamflow series should have a "\
                                     "frequency of '%s', got '%s' instead " % \
                                     (_hyetodata.freqstr, hydrodata.freqstr))
            self._hydrodata = hydrodata
            self.hydro._series = hydrodata
    #
    hydrodata = property(fget=_get_hydrodata,
                         fset=_set_hydrodata,
                         doc="Streamflow series.")



    def plot_hyeto(self, *args, **kwargs):
        """
    Plots the hyetograph on the top subplot.
    
    This method accepts the same arguments and optional parameters as a standard
    :class:`~scikits.timeseries.lib.plotlib.TimeSeriesPlot`
        """
        kwargs.update({'c':'#000099'})
        hyeto = self.hyeto.tsplot(self.hyetodata, *args, **kwargs)
        setp(self.hyeto.get_xticklabels(), visible=False)
        setp(self.hyeto.get_xticklabels(), visible=False)
        setp(self.hyeto.get_xticklabels(minor=True), visible=False)
        return hyeto
    plot_hyetograph = plot_hyeto


    def plot_hydro(self, *args, **kwargs):
        """
    Plots the hydrograph on the bottom subplot.

    This method accepts the same arguments and optional parameters as a standard
    :class:`~scikits.timeseries.lib.plotlib.TimeSeriesPlot`
        """
        return self.hydro.tsplot(self.hydrodata, **kwargs)
    plot_hydrograph = plot_hydro


    def plot_default(self, **kwargs):
        """
    Plots both the hyetograph (top plot) and the hydrograph (bottom plot).

    This method accepts the same arguments and optional parameters as a standard
    :class:`~scikits.timeseries.lib.plotlib.TimeSeriesPlot`
        """
        hyeto = self.plot_hyeto(**kwargs)
        hydro = self.plot_hydro(**kwargs)
        return (hyeto, hydro)


    def set_ensoindices(self, *args, **kwargs):
        """
    Sets the ENSO indices of both the precipitation and streamflows series.
        """
        self.hyeto._series.ensoindices = self.hyetodata.set_ensoindices(*args, **kwargs)
        self.hydro._series.ensoindices = self.hydrodata.set_ensoindices(*args, **kwargs)

#
#    def set_ensoindicator(self, ensoindicator):
#        self.hyetodata.set_ensoindicator(ensoindicator)
#        self.hydrodata.set_ensoindicator(ensoindicator)
#    set_ensoindicator.__doc__ = ClimateSeriesPlot.set_ensoindicator.__doc__
#
#
#    def plot_enso_background(self, *args, **kwargs):
#        self.hyeto.plot_enso_background(*args, **kwargs)
#        self.hydro.plot_enso_background(*args, **kwargs)
#        self.hyeto.set_xticklabels([])
#    plot_enso_background.__doc__ = ClimateSeriesPlot.plot_enso_background.__doc__


    def set_datelimits(self, start_date=None, end_date=None):
        """
    Sets the date limits of the plot to ``start_date`` and ``end_date``.
    The dates can be given as :class:`~scikits.timeseries.Date` objects, strings
    or integers.

    Parameters
    ----------
    start_date : {var}
        Starting date of the plot.
        If None, the current left limit is used.
    end_date : {var}
        Ending date of the plot.
        If None, the current right limit is used.

    Notes
    -----
    In practice, we set the date range of the hyetograph.
    Because the hydrograph and the hyetograph share the same x axis,
    the modification are automatically reported to the former graph.
        """
        self.hyeto.set_datelimits(start_date=start_date, end_date=end_date)
        setp(self.hyeto.get_xticklabels(), visible=False)
        setp(self.hyeto.get_xticklabels(minor=True), visible=False)


#................................................
def hydrograph(hyetodata, hydrodata, hyetoargs=None, hydroargs=None, num=None,
               figsize=None, dpi=None, facecolor=None, edgecolor=None,
               frameon=True, subplotpars=None):
    """
    Creates a new :class:`Hydrograph` object.

    Parameters
    ----------
    hyetodata : TimeSeries object
        Rainfall input data.
    hydrograph : TimeSeries object
        Streamflow data
    hyetoargs : dictionary, optional
        Arguments for the hyetograph plot.
    hydroargs : dictionary, optional
        Arguments for the hydrograph plot.
    num : {None, int}, optional
        Number of the figure.
        If None, a new figure is created and ``num`` is incremented.
    figsize : {None, tuple}
        Size of the figure, as a tuple (width, height) in inches.
        If None, defaults to rc figure.figsize.
    dpi : {None, int}, optional
        Resolution in dots per inches. 
        If None, defaults to rc figure.dpi
    facecolor : {None, string}, optional
        Background color.
        If None, defaults to rc figure.facecolor.
    edgecolor : {None, string}, optional
        Border color.
        If None, defaults to rc figure.edgecolor.
    linewidth : {float, None}
        Width of the patch edge line.
    frameon : {True, False}
        Whether to draw the frame around the figure.

    """
    hyetoargs = dict(hyetoargs or {})
    hydroargs = dict(hydroargs or {})
    figargs = dict(num=num, figsize=figsize, dpi=dpi, facecolor=facecolor,
                   frameon=frameon, FigureClass=Hydrograph,
                   subplotpars=subplotpars)
    figargs.update(hyetodata=hyetodata, hydrodata=hydrodata)

    fig = figure(**figargs)
    #
    fig.plot_hydro(**hydroargs)
    fig.plot_hyeto(**hyetoargs)
    #
    return fig





#######--------------------------------------------------------------------------
##---- --- Flow Duration Curves ---
#######--------------------------------------------------------------------------
#
##................................................
#def _plotfdc1d(subplot, series, **kwargs):
#    """Plots the flow duration curve for the 1D series series on the subplot plot.
#    """
#    scale = kwargs.pop('scale','log')
#    idx = series.argsort(fill_value=-9999)[::-1]
#    ydata = series[idx]
#    xdata = _getxdata(series, scale)
#    subplot.plot(xdata, ydata, **kwargs) 
#    
##................................................    
#def _getxdata(series, scale, exceedence=True):
#    """Computes the plotting positions, as (i/(n+1)) where n is the number of 
#    unmasked data of the 1D-series 'series'. 
#    If scale='log', the plotting positions are transformed into normal positions.
#    """ 
#    n = series.count()
#    xdata = ma.array(numpy.empty_like(series), mask=True)
#    xdata[:n] = np.linspace(1./(n+1), 1.-1./(n+1), n)
#    if scale == 'log':
#        xdata = norm.ppf(xdata)
#    if not exceedence:
#        xdata = xdata[::-1]
#    return xdata   
#
##..............................................................................
#def plot_ENSOcomparison_fdc(subplot, series, scale='log', **kwargs):
#    """Plots the comparison of flow duration curves of a single 1d series 
#    for the different ENSO phases.
#    """
#    if not hasattr(series, 'ensoindices') or series.ensoindices is None:
#        raise AttributeError("Undefined ENSO indices!")
#    if series.ndim > 1:
#        raise ValueError("The input series should be 1D only!")
#    #........................
#    xlims = kwargs.pop('xlims',(0.0005,0.9995))
#    #........................
#    scale = scale[:3].lower()   
#    if scale not in  ['lin', 'log']:
#        raise ValueError("Unrecognized option '%' for scale: "\
#                         "should be in ['lin','log'])")
#    #........................
#    _plotfdc1d(subplot, series, 
#               ls='-', c=ENSOlines[-9999], label='Global', scale=scale)
#    for i in (-1,0,1):
#        (attr, lab) = (ENSOattr[i], ENSOlabels[i])
#        _plotfdc1d(subplot, getattr(series,attr)(), 
#                   c=ENSOlines[i], label=lab, scale=scale, **kwargs)
#    #........................
#    _setfdclims(subplot, series, scale=scale, xlims=xlims)
#    return subplot
#
#...............................................................................
def plot_fdc(series, multimode=True, plot_enso=False,
             starting_month=None, lag=6,
             scale='log', xmin=0.0005, xmax=0.9995, ax=None, **kwargs):
    """
    Plots one or several flow duration curves (FDCs) for the series.
    
    The input series should be 1D or 2D.
    By default, if the series is 1D, one curve only will be plotted, whereas if
    the series is 2D, a curve will be plotted for each line of the series.
    
    A 1D series can also be converted into an annual series
    with the :keyword:`starting_month` parameter. 
    In that case, ``starting_month`` should be an integer between 1 and 12
    precising the month at which the 12-month period should start.
    For example, to plot the FDCs for each water year (usually from April
    to the following March), use ``starting_month=4``.
    
    When ``enso=True``, ENSO phases are plotted with different colors. 
    When the series is 2D or if it has been converted to an annual frequency, 
    the ENSO indices are defined with the ``full_year=True`` option, where an ENSO
    episode lasts at least 12 consecutive months.

    Parameters
    ----------
    series : TimeSeries
        Flow data.
    ax : {None, :class:`matplotlib.axes.Axes`}, optional
        Subplot where to plot the flow duration curves.
        If None, use the current plot.
    multimode : {True, False}, optional
        Whether to interpret a 2D input series as several series or a single one.
    starting_month : {None, integer}, optional
        First month of each year.
        If None, plots the global flow duration curve.
        Otherwise, ``starting_month`` must be an integer between 1 and 12,
        corresponding to the first month of the water year
        (usually, 4 for April).
    plot_enso : {True, False}, optional
        Whether to plot each ENSO phase with a different color.
    lag : {integer}, optional
        Number of months of lag for the definition of ENSO indices. For example,
        if lag=6, the ENSO phase starting in Oct. 2001 is applied starting on 
        Apr. 2002.
        If None, use a lag computed as the time difference between ``starting_month``
        and the first month of the reference season of the ENSO indicator (or October
        if undefined).
    scale : {'log','lin'}, optional
        String indicating whether the x-axis is in log (``'log'``) or linear
        (``'lin'``) scale.
        If ``'log'``, each plotting position is expressed as a Gaussian pdf.
    other parameters :
        The parameters recognized by the :func:`matplotlib.pyplot.plot` function are 
        also recognized.

    Raises
    ------
    TypeError
        If ``plot_enso=True`` but the series is not a
        :class:`~scikits.hydroclimpy.enso.ClimateSeries`.

    ValueError
        * If ``starting_month`` is not between 1 and 12.
        * If ``starting_month`` is defined but the initial series is not 1D.
    """
    if ax is None:
        ax = gca()
    # Make sure we have at most a 2D series ...............
    if series.ndim > 2:
        raise ValueError("The input series should be 2D at most!")
    # Get the ENSO indicator associated w/ the series (if any)
    ensoindicator = getattr(series, 'ensoindicator', None)
    # Check the starting month ............................
    if starting_month is not None:
        # Make sure we have an integer between 1 and 12
        starting_month = int(starting_month)
        if (starting_month < 1) or (starting_month > 12):
            errmsg = "The starting month should be between 1 (Jan.) and "\
                     "12 (Dec.)! (got %s instead)" % starting_month
            raise ValueError(errmsg)

    # Check whether we need to plot the ENSO information ..
    if plot_enso is True:
        # Make sure we have some ENSO information .........
        if ensoindicator is None:
            errmsg = "No ENSO information is associated with the input series."
            raise InvalidENSOError(errmsg)
        # Reset the indices if we have a starting_month ...
        if starting_month is not None:
            if lag is None:
                refmonth = (ensoindicator.reference_season or [10, ])[0]
                lag = (starting_month + 12 - refmonth) % 12
            series.set_ensoindices(full_year=True, lag=lag)
        else:
            # Make sure that the indices are already set
            series.set_ensoindices()
        # Load the default marker colors ..................
        from scikits.hydroclimpy.plotlib.ensotools import ENSOlines, \
                                                          ENSOmarkers, \
                                                          ENSOlabels
    # No ENSO information to plot : get basic lines & markers
    else:
        ENSOlines = {'G':'#cccccc'}
        ENSOmarkers = {'G':'#cccccc'}
    # Check whether we are in multimode or not ............
    ## 1D input
    if series.ndim == 1:
        # Convert to annual if needed
        if starting_month:
            multimode = True
            series = series.convert(FR_ANNSTART[starting_month - 1], func=None)
        else:
            multimode = False
        _series = series.view(ma.MaskedArray)
    ## 2D input
    else:
        #  w/ starting month
        if starting_month is not None:
            errmsg = "The input series should be 2D! (got %s instead)"
            raise ValueError(errmsg % str(series.shape))
        # w/o multimode
        if not multimode:
            _series = series.view(ma.MaskedArray).ravel()
    # Get the number of valid data per year (ie, per row)
    n = _series.count(axis= -1)
    # Get the xdata .........
    scale = scale[:3].lower()
    if scale == 'lin':
        if multimode:
            xdata = [np.linspace(1. / (nx + 1), 1 - 1. / (nx + 1), nx)
                     for nx in n]
        else:
            xdata = np.linspace(1. / (n + 1), 1 - 1. / (n + 1), n)
#            xdata = ma.empty(len(series), dtype=float)
#            xdata[:n] = np.linspace(1. / (n + 1), 1 - 1. / (n + 1), n)
    elif scale == 'log':
        if multimode:
            xdata = [norm.ppf(np.linspace(1. / (nx + 1), 1 - 1. / (nx + 1), nx))
                     for nx in n]
        else:
            xdata = norm.ppf(np.linspace(1. / (n + 1), 1 - 1. / (n + 1), n))
#            xdata = ma.empty(len(series), dtype=float)
#            xdata[:n] = norm.ppf(np.linspace(1. / (n + 1), 1 - 1. / (n + 1), n))
    else:
        raise ValueError("Unrecognized option '%' for scale: "\
                         "should be in ['lin','log'])")
    # Get some defaults .....
    if multimode:
        lwdefault = 0.8
        zorderdefault = 3
        colordefault = ENSOlines['G']
    else:
        lwdefault = 2
        zorderdefault = 10
        colordefault = 'k'
    marker = kwargs.pop('marker', 'o')
    markersize = kwargs.get('markersize', kwargs.get('ms', 3))
    lw = kwargs.pop('linewidth', kwargs.pop('lw', lwdefault))
    zorder = kwargs.pop('zorder', zorderdefault)
    color = kwargs.pop('color', kwargs.pop('c', colordefault))

    # Multi-mode : one line per year ......................
    if multimode:
        if plot_enso:
            ensoindices = series.ensoindices
            if ensoindices.ndim > 1:
                ensoindices = ensoindices[:, 0]
            # ENSO mode : different colors for different phases
#            eidx = series.ensoindices._data
#            # Take the first column if it's 2D
#            if eidx.ndim > 1:
#                eidx=eidx[:,0]
            for(i, attr) in zip((-1, 0, 1), ('cold', 'neutral', 'warm')):
                key = attr[0].upper()
                label = ENSOlabels[key]
                ydata = series[ensoindices == i]
                ydata = [np.sort(_).compressed()[::-1] for _ in ydata]
#                ydata = np.sort(getattr(series, attr).compressed())[::-1]
                points = [zip(x, y) for (x, y) in zip(xdata, ydata)]
                collec = LineCollection(points,
                                        label=ENSOlabels[key],
                                        color=ENSOlines[key],
                                        zorder=zorder, linewidth=lw)
                ax.add_collection(collec, autolim=True)
        else:
            ydata = [np.sort(y.compressed())[::-1] for y in _series]
            points = [zip(x, y) for (x, y) in zip(xdata, ydata)]
            label = kwargs.pop('label', None)
            collec = LineCollection(points, label=label, linewidth=lw,
                                    colors=ENSOlines['G'])
            ax.add_collection(collec, autolim=True)
    # One line for the while dataset ......................
    else:
        ydata = ma.sort(series.compressed(), endwith=False)[::-1]
        points = [zip(xdata, ydata._series)]
        label = kwargs.pop('label', 'none')
        collec = LineCollection(points, label=label, linewidth=lw,
                                colors=color, zorder=zorder)
        ax.add_collection(collec, autolim=True)
        # If we need to add some colors
        if plot_enso and marker:
            for attr in ('cold', 'neutral', 'warm'):
                key = attr[0].upper()
                label = ENSOlabels[key]
                color = ENSOmarkers[key]
                #ydata = ma.sort(getattr(series, attr), endwith=False)[::-1]
                current = getattr(ydata, attr)._series
                _fdc = ax.plot(xdata, current, ls='', lw=0,
                               marker=marker, ms=markersize,
                               mfc=color, mec=color,
                               label=label, zorder=zorder)
    #........................
    set_normal_limits(ax, xmin=xmin, xmax=xmax, scale=scale)
    ax.set_ylim(_series.min(), _series.max())
    return ax


if __doc__ is not None:
    hydrograph.__doc__ = hydrograph.__doc__ % mpl_addons._doc_parameters
#    Hydrograph.__doc__ = Hydrograph.__doc__ % doc_parameters
