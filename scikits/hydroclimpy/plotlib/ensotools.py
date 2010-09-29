# -*- coding: utf-8 -*-
"""
.. currentmodule:: scikits.hydroclimpy.plotlib.ensotools

This module defines new classes and functions to plot series with ENSO information.

Classes and Functions
=====================

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ClimateFigure
   ClimateSeriesPlot
   ENSOPhaseComparisonPlot
   csfigure
   csplot
   add_csplot

_____


Module-wide Data
================

.. data:: ENSOcolors

   Defines a dictionary storing the color codes corresponding to each ENSO phase.
   
   ``markers`` : dictionary
      Color codes for matplotlib markers.
   ``lines`` : dictionary
      Color codes for matplotlib lines.
   ``polygons`` : dictionary
      Color codes for matplotlib polygons.
   ``fill`` : dictionary
      Color codes for matplotlib filled areas.


.. data:: ENSOmarkers

   Alias to :data:`ENSOcolors['markers']`.
   
   + **'G'** : Color code corresponding to global data : ``'#cccccc'``
   + **'C'** : ... to cold (La Niña) episodes:``'#6666cc'``
   + **'N'** : ... Neutral episodes :``'#669933'``
   + **'W'** : ... to warm (El Niño) episodes :``'#ffcc33'``


.. data:: ENSOlines

   Alias to :data:`ENSOcolors['lines']`.
   
   + **'G'** : Color code corresponding to global data : ``'#cccccc'``
   + **'C'** : ... to cold (La Niña) episodes : ``'#6666cc'``
   + **'N'** : ... to Neutral episodes : ``'#009900'``
   + **'W'** : ... to warm (El Niño) episodes : ``'#990000'``


.. data:: ENSOpolygons

   Alias to :data:`ENSOcolors['polygons']`.
   
   + **'G'** : Color code corresponding to global data : ``'#cccccc'``.
   + **'C'** : ... to cold (La Niña) episodes : ``'#6666cc'``
   + **'N'** : ... to Neutral episodes : ``'#ffffff'``
   + **'W'** : ... to warm (El Niño) episodes : ``'#ffcc33'``


.. data:: ENSOfill

   Alias to :data:`ENSOcolors['fill']`.
   
   + **'G'** : Color code corresponding to global data : ``'#cccccc'``
   + **'C'** : ... to cold (La Niña) episodes : ``'#6666cc'``
   + **'N'** : ... to Neutral episodes : ``'#669933'``
   + **'W'** : ... to warm (El Niño) episodes : ``'#ffcc33'``


.. data:: ENSOlabel

   Dictionary storing the labels for each ENSO phase
   
   * ``'G'`` : 'Global'
   * ``'C'`` : 'Cold'
   * ``'N'`` : 'Neutral'
   * ``'W'`` : 'Warm'

"""
__author__ = "Pierre GF Gerard-Marchant"



import numpy as np
import numpy.ma as ma

import matplotlib
from matplotlib import pyplot, _pylab_helpers
from matplotlib.axes import Subplot
from matplotlib.collections import LineCollection, BrokenBarHCollection
from matplotlib.colors import BoundaryNorm
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, Locator, FixedLocator
from matplotlib.transforms import blended_transform_factory

import scikits.timeseries
from scikits.timeseries.cseries import check_freq
import scikits.timeseries.lib.plotlib as mpl
from scikits.timeseries.lib.plotlib import TimeSeriesFigure, TimeSeriesPlot, tsplot

from scikits.hydroclimpy import Cluster, periodize
from scikits.hydroclimpy.plotlib import whiskerbox
import scikits.hydroclimpy.enso as enso
#hpl_doc_parameters = scikits.hydroclimpy.plotlib._doc_parameters

import logging
testlogging = logging.getLogger("CPL")


ENSOcolors = {
'markers':  {'G':'#cccccc', 'C':'#6666cc', 'N':'#669933', 'W':'#ffcc33'},
'lines' :   {'G':'#cccccc', 'C':'#6666cc', 'N':'#009900', 'W':'#990000'},
'polygons': {'G':'#cccccc', 'C':'#6666cc', 'N':'#ffffff', 'W':'#ffcc33'},
'fill':     {'G':'#cccccc', 'C':'#6666cc', 'N':'#669933', 'W':'#ffcc33'},
}
ENSOlabels = {'G': 'Global', 'C':'Cold', 'N':'Neutral', 'W':'Warm'}
ENSOattr = {'G': None, 'C':'cold', 'N':'neutral', 'W':'warm'}
ENSOlines = ENSOcolors['lines']
ENSOmarkers = ENSOcolors['markers']
ENSOpolygons = ENSOcolors['polygons']

ENSOmap = BoundaryNorm((-9999, -1, 0, +1),
                       ('#cccccc', '#6666cc', '#669933', '#ffcc33'))

#####--------------------------------------------------------------------------
#---- --- ClimateSeriesPlot ---
#####--------------------------------------------------------------------------
class ClimateSeriesPlot(TimeSeriesPlot, object):
    """
    Defines a subclass of :class:`~matplotlib.axes.Subplot` suitable for plotting
    :class:`~scikits.hydroclimpy.enso.ClimateSeries` objects.

    The subplot can be instantiated using the same keywords as for a standard
    subplot.
    A specific keyword :keyword:`series` is also recognized, which defines
    the underlying :class:`~scikits.hydroclimpy.enso.ClimateSeries`.
    
    
    See Also
    --------
    :class:`scikits.timeseries.lib.plotlib.TimeSeriesPlot`
    
    """

    def __init__(self, fig=None, *args, **kwargs):
        # Retrieve the series ...................
        TimeSeriesPlot.__init__(self, fig, *args, **kwargs)
        if self._series is not None:
            if isinstance(self._series, enso.ENSOIndicator):
                self.ensoindicator = self._series
            else:
                self.ensoindicator = self._series.ensoindicator
    #
    def gca(self, **kwargs):
        """
    Return the current axes, creating a new :class:`ClimateSeriesPlot` if necessary.
        """
        ax = self._axstack()
        if ax is not None:
            return ax
        return self.add_csplot(111, **kwargs)
    #        
    def set_ensoindicator(self, ensoindicator):
        """
    Sets the :class:`~scikits.hydroclimpy.enso.ENSOIndicator` of the
    :class:`ClimateSeriesPlot`.
    
    Parameters
    ----------
    ensoindicator : :class:`~scikits.hydroclimpy.enso.ENSOIndicator`
        The ENSO indicator used for the definition of ENSO phases.

    Raises
    ------
    TypeError
        If the ``ensoindicator`` parameter is not a valid
        :class:`~scikits.hydroclimpy.enso.ENSOIndicator` object.
        """
        if not isinstance(ensoindicator, enso.ENSOIndicator):
            raise TypeError("Invalid ENSO indicator !")
        self.ensoindicator = ensoindicator
        self.ensoindices_limits = enso.grouped_indices_limits(self._series)
    #
    csplot = TimeSeriesPlot.tsplot
    #        
    def plot_enso_background(self, ensoindices=None, lag=0, **optinfo):
        """
    Plots colored stripes in the background of plot to represent ENSO phases.

    Parameters
    ----------
    ensophases: {array-like}, optional
        Array of ENSO indices (``+1`` for El Niño, ``0`` for Neutral and
        ``-1`` for La Niña episodes).
        If None, the ENSO indices of the underlying series are used instead.
    
        """
        if ensoindices is None:
            series = self._series
            if series is None or not hasattr(series, 'ensoindices') or \
                series.ensoindices is None:
                raise ValueError("Undefined ENSO indices!")
            ensoindices = series.ensoindices
        if self.xdata is None:
            errmsg = "Unable to retrieve the dates of the current plot!"
            raise ValueError(errmsg)
        #
        dates = self.xdata
        clust_indices = Cluster(ensoindices.filled(0), 0)
        _dates = np.empty(len(dates) + 1, int)
        _dates[:-1] = dates
        _dates[-1] = dates[-2]
        episodes = dict([(k, zip(_dates[v[:, 0]], _dates[v[:, 1]]))
                         for (k, v) in clust_indices.grouped_limits().items()])
        #
        colors = ENSOcolors['polygons']
        for (key, idx) in {'C':-1, 'N':0, 'W':+1}.iteritems():
            colors[idx] = colors[key]
        #
        trans = blended_transform_factory(self.transData, self.transAxes)
        for (k, lim) in episodes.iteritems():
            _bbc = BrokenBarHCollection([(x + lag, y - x) for (x, y) in lim],
                                         (0, 1),
                                         facecolors=colors[k],
                                         edgecolors=colors[k],)
            _bbc.set_alpha(0.2)
            _bbc.set_transform(trans)
            self.add_collection(_bbc)

#####--------------------------------------------------------------------------
#---- --- Phases Comparisons ---
#####--------------------------------------------------------------------------

class ENSOPhaseComparisonPlot(Subplot):
    """Class for plotting the comparison of data between ENSO phases.

    The object is instantiated with the same optional parameters 
    as for a standard :class:`matplotlib.axes.Subplot`.
    
    The parameter ``series`` must be given to define the series to plot. 
    ``series`` must be a :class:`~scikits.hydroclimpy.enso.ClimateSeries`, and an
    :class:`~scikits.hydroclimpy.enso.ENSOIndicator` must have been allocated.
    The series will be transformed into an annual series with the 
    :func:`~scikits.hydroclimpy.core.periodize` function.

    """
    def __init__(self, fig=None, *args, **kwargs):
        """Initializes instance.
        
    Parameters
    ----------
    fig : {None, figure instance}
        Figure from which the plot will depend.
        If None, the current figure is selected.
    series : ClimateSeries
        Series to convert.
    freq : {var}, optional
        A valid frequency specifier (as a string or integer), or a 3-letter 
        string corresponding to the first quarter of the year.
        Use this parameter if the series does not have an adequate frequency.
    func : {None,function}, optional
        Function controlling how data sharing the same new dates should be
        manipulated.
        This function should handle masked values appropriately.
    kwargs : dictionary
        Dictionary of optional parameters.
        The same parameters as for a standard subplot instantiation are recognized.

    See Also
    --------
    scikits.hydroclimpy.periodize
        Function to convert a series into a 2D series with years as rows and
        periods as columns.
        """
        # Initialize the plot ..............
        fig = fig or kwargs.pop('fig', None)
        if fig is None:
            fig = pyplot.gcf()
        elif not isinstance(fig, Figure):
            fig = pyplot.figure(fig)
        if len(args) == 0:
            args = (111,)
        # Get the information about the series
        if not hasattr(series, 'ensoindices'):
            errmsg = "ENSO indices should be defined for the input series!"
            raise AttributeError(errmsg)
        series = kwargs.pop('series', None)
        freq = kwargs.pop('freq', None)
        func = kwargs.pop('func', None)
        self.series = periodize(series, freq=freq, func=func)
        self.cold = self.series.cold
        self.neutral = self.series.neutral
        self.warm = self.series.warm
        self.freq = check_freq(self.series._dates.freq)
        self.positions = np.arange(self.series.shape[-1])
        self.periodlabels = self.series.optinfo['period_names']
        Subplot.__init__(self, fig, *args, **kwargs)


    def barplot(self, func=None, *args, **kwargs):
        """
    Plots a bar chart comparing the phases for each period, as transformed
    by the ``func`` function.
        
    Parameters
    ----------
    func : function, optional
        Function to apply.
        By default, use the :func:`numpy.ma.mean` function.
    args : var
        Mandatory arguments of function ``func``.
    kwargs : var
        Optional arguments of function ``func``.
    
        """
        func = func or ma.mean
        width = 0.2
        colordict = ENSOcolors['fill']
        barlist = []
        pos = self.positions - 3 * width
        for (i, s) in zip(['W', -1, 0, 1],
                         [self.series, self.cold, self.neutral, self.warm]):
            pos += width
            series = ma.apply_along_axis(func, 0, s, *funcopt)
            b = self.bar(pos, series, width=width, bottom=0,
                         color=colordict[i], ecolor='k', capsize=3)
            barlist.append(b[0])
        self.barlist = barlist
        self.figure.axes.append(self)
        self.format_xaxis()
        return barlist


    def whiskerplot(self, plot_mean=False):
        """
    Plots a whiskerbox plot comparison of the different ENSO phases.

    Parameters
    ----------
    plot_mean : {False, True}, optional
        Whether to plot an horizontal bar across each box to represent the mean.

    Returns
    -------
    barlist
        The list of the boxes.
        """
        width = 0.2
        colordict = ENSOcolors['fill']
        barlist = []
        pos = self.positions - 3 * width
        #.........................
        for (i, s) in zip(['G', 'C', 'N', 'W'],
                         [self.series, self.cold, self.neutral, self.warm]):
            pos += width
            b = whiskerbox(s, self, positions=pos, width=width,
                           color=colordict[i], plot_mean=plot_mean)
            barlist.append(b[0])
        self.barlist = barlist
        self.figure.axes.append(self)
        self.format_xaxis()
        return barlist

    def legend(self, **kwargs):
        """
    Adds a legend to the plot.
        """
        if not hasattr(self, 'barlist'):
            args = (('Global', 'Cold', 'Neutral', 'Warm'),)
        else:
            args = (self.barlist, ('Global', 'Cold', 'Neutral', 'Warm'),)
        Subplot.legend(self, *args, **kwargs)

    def format_xaxis(self, **kwargs):
        """
    Formats the x axis.
        """
        self.xaxis.set_major_locator(FixedLocator(self.positions))
        self.set_xticklabels(self.periodlabels, **kwargs)
        self.set_xlim((-1, len(self.positions) + 1))



#####--------------------------------------------------------------------------
#---- --- ClimateSeries Figure ---
#####--------------------------------------------------------------------------

class ClimateFigure(TimeSeriesFigure):
    """
    A subclass of :class:`~scikits.timeseries.lib.plotlib.TimeSeriesFigure`
    designe to plot :class:`~scikits.hydroclimpy.enso.ClimateSeries` objects.
    """
    def __init__(self, **kwargs):
        TimeSeriesFigure.__init__(self, **kwargs)

def csfigure(series=None,
             num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None,
             frameon=True, subplotpars=None):
    """
    Creates a new :class:`ClimateFigure` object.

    Parameters
    ----------
    series : {None, TimeSeries object}
        Input data.
    num : {None, int}, optional
        Number of the figure.
        If None, a new figure is created and ``num`` is incremented.
    %(figsize)s
    %(dpi)s
    %(facecolor)s
    %(edgecolor)s
    %(frameon)s
    %(subplotpars)s
    """
    figargs = dict(num=num, figsize=figsize, dpi=dpi,
                   facecolor=facecolor, edgecolor=edgecolor,
                   frameon=frameon, subplotpars=subplotpars)
    figargs.update(FigureClass=ClimateFigure, series=series)
    fig = pyplot.figure(**figargs)
    return fig

def add_csplot(axes, *args, **kwargs):
    """
    Adds a :class:`ClimateSeriesPlot` to the current plot.
    
    """
    kwargs.update(SubplotClass=ClimateSeriesPlot)
    if 'series' not in kwargs.keys():
        if hasattr(axes, 'series'):
            kwargs['series'] = axes.series
        elif hasattr(axes, '_series'):
            kwargs['series'] = axes._series
        else:
            kwargs['series'] = None
    return mpl.add_generic_subplot(axes, *args, **kwargs)

Figure.add_csplot = add_csplot
Figure.add_subplot = add_csplot
#................................................

def csplot(series, *args, **kwargs):
    """
    Plots the series to the current :class:`ClimateSeriesPlot` subplot.
    If the current plot is not a :class:`ClimateSeriesPlot`, 
    a new :class:`ClimateFigure` is created.

    """
    # allow callers to override the hold state by passing hold=True|False
    b = pyplot.ishold()
    h = kwargs.pop('hold', None)
    if h is not None:
        pyplot.hold(h)
    # Get the current figure, or create one
    figManager = _pylab_helpers.Gcf.get_active()
    if figManager is not None :
        fig = figManager.canvas.figure
        if not isinstance(fig, ClimateFigure):
            fig = csfigure(series=series)
    else:
        fig = csfigure(series=series)
    # Get the current axe, or create one
    sub = fig._axstack()
    if sub is None:
        sub = fig.add_csplot(111, series=series, **kwargs)
    try:
        ret = sub.csplot(series, *args, **kwargs)
        pyplot.draw_if_interactive()
    except:
        pyplot.hold(b)
        raise
    pyplot.hold(b)
    return ret


def add_ensocomparisonplot(axes, *args, **kwargs):
    """
    Adds a :class:`ClimateSeriesPlot` to the current plot.
    
    """
    kwargs.update(SubplotClass=ENSOPhaseComparisonPlot)
    if 'series' not in kwargs.keys():
        series = getattr(axes, '_series', getattr(axes, 'series', None))
        kwargs['series'] = series
    return mpl.add_generic_subplot(axes, *args, **kwargs)


Figure.add_ensocomparisonplot = add_ensocomparisonplot


###############################################################################
#if __doc__ is not None:
#    csfigure.__doc__ = csfigure.__doc__ % hpl_doc_parameters
