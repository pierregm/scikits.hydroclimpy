"""
===============================
Generic additions to matplotlib
===============================

The :mod:`~scikits.climply.plotlib` module provides several additions to matplotlib.

Formatters
==========

Overview
--------

.. autosummary::
   :nosignatures:

   EngrFormatter
   MinimalFormatter

API
---

.. autoclass:: EngrFormatter

.. autoclass:: MinimalFormatter


Functions
=========

.. autosummary::
   :nosignatures:

   whiskerbox
   qqplot
   set_normal_ticks
   set_normal_limits


API
---

.. autofunction:: whiskerbox
.. autofunction:: qqplot


.. autofunction:: set_normal_ticks
.. autofunction:: set_normal_limits


"""
__author__ = "Pierre GF Gerard-Marchant"



import numpy as np
from numpy import ndarray
import numpy.ma as ma

import scipy.stats.distributions as ssd
import scipy.stats.mstats as mstats

from matplotlib import pyplot
from matplotlib.pyplot import FixedLocator, Formatter, FixedFormatter, \
                              ScalarFormatter
from ..core.stats_addons import qqcalc


###############################################################################

_doc_parameters = dict(
figsize="""    figsize : {None, tuple}
        Size of the figure, as a tuple (width, height) in inches.
        If None, defaults to rc figure.figsize.
    """,
dpi="""    dpi : {None, int}, optional
        Resolution in dots per inches. 
        If None, defaults to rc figure.dpi
    """,
facecolor="""    facecolor : {None, string}, optional
        Background color.
        If None, defaults to rc figure.facecolor.
    """,
edgecolor="""    edgecolor : {None, string}, optional
        Border color.
        If None, defaults to rc figure.edgecolor.
    """,
linewidth="""    linewidth : {float, None}
        Width of the patch edge line.
    """,
frameon="""    frameon : {True, False}
        Whether to draw the frame around the figure.
    """,
subplotpars="""    subplotpars : {None, var}
    A :class:`SubplotParams` instance, defaults to rc
    """,
)



day_of_year_locator = FixedLocator([  1, 16, 32, 47, 60, 75,
                                     91, 106, 121, 136, 152, 167,
                                    182, 197, 213, 228, 244, 259,
                                    274, 289, 305, 320, 335, 350, ])
day_of_year_formatter = FixedFormatter(
                        ['Jan-01', '', 'Feb-01', '', 'Mar-01', '', 'Apr-01', '',
                         'May-01', '', 'Jun-01', '', 'Jul-01', '', 'Aug-01', '',
                         'Sep-01', '', 'Oct-01', '', 'Nov-01', '', 'Dec-01', '']
                        )
monthlist = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'
             ]


##### -------------------------------------------------------------------------
#---- --- Formatters ---


class EngrFormatter(ScalarFormatter):
    """
    A variation of the standard :class:`ScalarFormatter`, using only multiples 
    of three in the mantissa.
    A fixed number of decimals can be displayed with the optional parameter 
    ``ndec``.

    Parameters
    ----------
    ndec : {None, int}, optional
        Number of decimals to be displayed.
        If ``None``, the number of decimals is defined from the current ticks.
    useOffest : {True, False}, optional
    useMathTtest : {True, False}
    strict : {False, True}
    """
    def __init__(self,
                 ndec=None, useOffset=True, useMathText=False, strict=False):
        ScalarFormatter.__init__(self, useOffset, useMathText)
        self._powerlimits = (-3, 3)
        if ndec is None or ndec < 0:
            self.format = None
        elif ndec == 0:
            self.format = "%d"
        else:
            self.format = "%%1.%if" % ndec
        self.strict = strict
    #
    def _set_orderOfMagnitude(self, mrange):
        """Sets the order of magnitude."""
        ScalarFormatter._set_orderOfMagnitude(self, mrange)
        self.orderOfMagnitude = 3 * (self.orderOfMagnitude // 3)
    #
    def _set_format(self):
        """Sets the format string to format all ticklabels."""
        # set the format string to format all the ticklabels
        locs = (np.array(self.locs) - self.offset) / 10 ** self.orderOfMagnitude
        locs += 1e-15
        sigfigs = [len(str('%1.3f' % loc).split('.')[1].rstrip('0')) \
                   for loc in locs]
        sigfigs.sort()
        if self.format is None:
            self.format = '%1.' + str(sigfigs[-1]) + 'f'
        if self._usetex or self._useMathText:
            self.format = '$%s$' % self.format


#..............................................................................
class MinimalFormatter(Formatter):
    """
    A minimal formatter: just the plain data !
    """
    def __init__(self, sigfigs=None):
        if sigfigs is None:
            self.fmt = "%f"
        else:
            self.fmt = "%.%if" % sigfigs
    #
    def __call__(self, x, pos=None):
        return str(self.fmt % x).rstrip('0')


#####--------------------------------------------------------------------------
#---- --- Functions ---


def whiskerbox(series, fsp=None, positions=None, mode='mquantiles', width=0.8,
               wisk=None, plot_mean=False, logscale=None, color=None,
               outliers=None):
    """
    Draws a whisker plot.
    The bottom and top of the boxes correspond to the lower and upper quartiles
    respectively (25th and 75th percentiles).
    

    Parameters
    ----------
    series : Sequence
        Input data. 
        If the sequence is 2D, each column is assumed to represent a different variable.
    fsp : :class:`Subplot`
        Subplot where to draw the data.
        If None, uses the current axe.
    positions : {None, sequence}, optional
        Positions along the x-axis.
        If None, use a scale from 1 to the number of columns.
    mode : {'mquantiles', 'hdquantiles'}, optional
        Type of algorithm used to compute the quantiles. 
        If 'mquantiles', use the classical form :func:`~scipy.stats.mstats.mquantiles`
        If 'hdquantiles', use the Harrell-Davies estimators of the function
        :func:`~scipy.stats.mmorestats.hdquantiles`.
    wisk : {None, float}, optional
        Whiskers size, as a multiplier of the inter-quartile range. 
        If None, the whiskers are drawn between the 5th and 95th percentiles.
    plot_mean : {False, True}, optional
        Whether to overlay the mean on the box.
    color : {None, string}, optional
        Color of the main box.
    outliers : {dictionary}, optional
        Options for plotting outliers.
        By default, the dictionary uses 
        ``dict(marker='x', ms=4, mfc='#999999', ls='')``

    """
    outliers = outliers or dict(marker='x', ms=4, mfc='#999999', mec='#999999',
                                ls='',)
    if fsp is None:
        fsp = pyplot.gca()
    if not fsp._hold:
        fsp.cla()
    # Make sure the series is a masked array 
    series = ma.array(series, copy=False, subok=False)
    # Reshape the series ...................
    if series.ndim == 1:
        series = series.reshape(-1, 1)
    elif series.ndim > 2:
        series = np.swapaxes(series, 1, -1).reshape(-1, series.shape[1])
    if positions is None:
        positions = np.arange(1, series.shape[1] + 1)
    # Get the quantiles ....................
    plist = [0.05, 0.25, 0.5, 0.75, 0.95]
    # Harrell-Davies ........
    if mode == 'hdquantiles':
        # 1D data ...........
        if series.ndim == 0:
            (qb, ql, qm, qh, qt) = mstats.hdquantiles(series.ravel(), plist)
        # 2D data ...........
        else:
            (qb, ql, qm, qh, qt) = ma.apply_along_axis(mstats.hdquantiles, 0,
                                                       series, plist)
    # Basic quantiles .......
    else:
        (qb, ql, qm, qh, qt) = mstats.mquantiles(series, plist, axis=0)
    # Get the heights, bottoms, and whiskers positions
    heights = qh - ql
    bottoms = ql
    if wisk is not None:
        hival = qh + wisk * heights
        loval = ql - wisk * heights
    else:
        (hival, loval) = (qt, qb)
    # Plot the whiskers and outliers .......
    for i, pos, xh, xl in np.broadcast(np.arange(len(positions)),
                                    positions, hival, loval):
        x = series[:, i]
        # Get high extreme ..
        wisk_h = x[(x <= xh).filled(False)]
        if len(wisk_h) == 0:
            wisk_h = qh[i]
        else:
            wisk_h = max(wisk_h)
        # Low extremes ......
        wisk_l = x[(x >= xl).filled(False)]
        if len(wisk_l) == 0:
            wisk_l = ql[i]
        else:
            wisk_l = min(wisk_l)
        fsp.plot((pos, pos), (wisk_l, wisk_h), dashes=(1, 1), c='k', zorder=1)
        fsp.plot((pos - 0.25 * width, pos + 0.25 * width), (wisk_l, wisk_l),
                 '-', c='k')
        fsp.plot((pos - 0.25 * width, pos + 0.25 * width), (wisk_h, wisk_h),
                 '-', c='k')
        # Outliers, if any...
        if outliers is not None and len(outliers) > 0:
            flh = x[(x > xh).filled(False)].view(ndarray)
            fll = x[(x < xl).filled(False)].view(ndarray)
            if len(flh) > 0 and len(fll) > 0:
                fsp.plot([pos] * (len(flh) + len(fll)), np.r_[flh, fll],
                         **outliers)
        # Plot the median....
        fsp.plot((pos - 0.5 * width, pos + 0.5 * width), (qm[i], qm[i]),
                 ls='-', c='k', lw=1.2, zorder=99)
        # Plot the mean......
        if plot_mean:
            fsp.plot((pos - 0.5 * width, pos + 0.5 * width),
                     (x.mean(), x.mean()),
                     ls=':', dashes=(1, 1), c='#000000', lw=1.1, zorder=99)
#            fsp.plot((pos,), (x.mean(),), marker='o', color=color, zorder=99)
    # Plot the boxes .......................
    bars = fsp.bar(positions - 0.5 * width, heights,
                   width=width, bottom=bottoms,
                   color=color, yerr=None, xerr=None,
                   ecolor='k', capsize=3, zorder=50)
    if logscale:
        fsp.set_yscale('log')
    return bars


def set_normal_ticks(fsp=None):
    """
    Change the ticks of the x-axis to a gaussian quantiles.
    
    Parameters
    ----------
    fsp : {None,Subplot}
        Subplot where to draw the data.
        If ``None``, uses the current axis.
    """
    xtmaj = np.array([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
#    xtmaj_str = [r"$10^{-1}$", r"$1$", r"$10$", 
#                 r"$50$", r"$90$", r"$99$", r"$99.9$"]
    xtmaj_str = [r"$10^{-1}$", 1, 10, 50, 90, 99, 99.9]
    xtmin = np.concatenate([np.linspace(0.0001, 0.001, 10),
                            np.linspace(0.001, 0.01, 10),
                            np.linspace(0.01, 0.1, 10),
                            np.linspace(0.1, 0.9, 9),
                            np.linspace(0.9, 0.99, 10),
                            np.linspace(0.99, 0.999, 10),
                            np.linspace(0.999, 0.9999, 10)])
    if fsp is None:
        fsp = pyplot.gca()
    fsp.xaxis.set_major_locator(FixedLocator(ssd.norm.ppf(xtmaj)))
    fsp.xaxis.set_minor_locator(FixedLocator(ssd.norm.ppf(xtmin)))
    fsp.set_xticklabels(xtmaj_str)


def set_normal_limits(fsp, xmin=None, xmax=None, emit=True, scale='log'):
    """
    Sets the limits of the xaxis for a normal plot.
    
    Parameters
    ----------
    xmin : {None, float}, optional
        Lower limit of the axis, as a float between 0 and 1.
        If None, use the current limit.
    xmax : {None, float}, optional
        Upper limit of the axis, as a float between 0 and 1.
        If None, use the current limit.
    emit : {True, False}, optional
        Whether notify the observers of a limit change.
    scale : {'log', 'lin'}, optional
        Whether the
    """
    if xmin is not None:
        if (xmin < 0) or (xmin > 1):
            raise ValueError("Invalid lower limit: should be between 0 and 1")
    if xmax is not None:
        if (xmax < 0) or (xmax > 1):
            raise ValueError("Invalid lower limit: should be between 0 and 1")
    #
    if scale == 'log':
        output = fsp.set_xlim(xmin=ssd.norm.ppf(xmin),
                              xmax=ssd.norm.ppf(xmax), emit=emit)
        set_normal_ticks(fsp)
    else:
        output = fsp.set_xlim(xmin=xmin, xmax=xmax, emit=emit)
    return output



def qqplot(data, distrib=ssd.norm, alpha=.4, beta=.4,
           fsp=None, plot_line=True, **kwargs):
    """
    Returns a quantile-quantile plot with theoretical quantiles in abscissae, and
    experimental quantiles in ordinates.
    The experimental quantiles are estimated through the equation:
    
    .. math::
       q_i = \\frac{i-\\alpha}{n-\\alpha-\\beta}

    where :math:`i` is the rank order statistic, :math:`n` the number of 
    unmasked data and :math:`\\alpha` and :math:`\\beta` two parameters between 0
    and 1. The default :math:`\\alpha=\\beta=0.4` gives approximately unbiased
    estimates of the quantiles.
    
    
    Parameters
    ----------
    data : array
        Input data
    distrib : {norm, function}, optional
        Theoretical distribution used to compute the expected quantiles.
        If None, use a normal distribution.
    alpha : {float} optional
        Coefficient for the computation of plotting positions.
    beta : {float} optional
        Coefficient for the computation of plotting positions.
    fsp : {None, Subplot}, optional
        Subplot where to plot. If None, use the current axe.
    plot_line : {True, False}
        Whether to compute a regression line.

    Returns
    -------
    plotted : :class:`matplotlib.lines.Line2D`
        Plotted data
    lines :  :class:`matplotlib.lines.Line2D`
        Plotted regression line
    (a,b) : tuple
        Slope and intercept of the regression line.

    Notes
    -----
    * The ``distrib`` parameter must be a function with a :meth:`.ppf` method.
    * The input data is ravelled beforehand.

    See Also
    --------
    scipy.stats.mstats.mquantiles
        Computes quantiles from a population.
    """
    data = np.ravel(data)
    qq = qqcalc(data, distrib=distrib, alpha=alpha, beta=beta)
    #
    if fsp is None:
        fsp = pyplot.gca()
    if not len(kwargs):
        kwargs.update(marker='o', c='k', ls='')
    plotted = fsp.plot(qq, data, **kwargs)
    #
    if plot_line:
        (a, b) = ma.polyfit(qq, data, 1)
        xlims = fsp.get_xlim()
        regline = np.polyval((a, b), xlims)
        lines = fsp.plot(xlims, regline, 'k:')
        return (plotted, lines, (a, b))
    return plotted

