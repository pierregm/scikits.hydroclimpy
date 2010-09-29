
.. currentmodule:: scikits.hydroclimpy.plotlib

====================
Flow duration curves
====================

A flow duration curve represents how often any given flow discharge is likely to
be equaled or exceeded.
The x axis corresponds to probabilities of exceedence, while the y axis corresponds
to stream flow discharges.

To construct a flow duration curve (FDC), stream flow discharges are ranked in 
decreasing order and associated with a given plotting position: if ``n`` is the
number of flow data, the highest flow is associated with ``1/(n+1)`` and the 
lowest flow with ``n/(n+1)``.

The function :func:`~flowplots.plot_fdc` provides a convenient
way to plot FDCs.

This function requires a time series as input, either as a 1D or a 2D series.
In that second case, each row should correspond to a given period of flow records
(usually a year).

By default, if the series is 1D, only one FDC will be plotted.
If the series is 2D, one FDC will be plotted for each row of the series.

A 1D series can also be converted to an annual frequency
using the :keyword:`starting_month` parameter.
This parameter precises the month at which the year should start.
For example, water years are usually defined as starting in April and lasting to
the following march.
To plot FDCs for each water year, use ``starting_month=4``.

>>> import numpy as np
>>> import scikits.hydroclimpy as hydro
>>> import scikits.hydroclimpy.enso as enso
>>> import scikits.hydroclimpy.io.usgs as usgs
>>> import scikits.hydroclimpy.plotlib as cpl


Let's import the streamflow data for the North Oconee River in Athens, GA.
The corresponding USGS code of the gage is '02217770'.


>>> series = usgs.load_usgs_flows('02217770')


Let's create a figure and add a subplot.

>>> fig = cpl.figure()
>>> fsp = fig.add_subplot(111)


Let's call the :func:`~flowplots.plot_fdc` function to plot the
corresponding flow duration curve.
By default, the curve is plotted on the current subplot. 
We can also select a particular subplot with the :keyword:`subplot` parameter.
We will use the additional parameters :keyword:`lw` to specify the linewidth 
of the curve, :keyword:`ls` for the linestyle, :keyword:`c` for the color, 
and :keyword:`zorder` for the priority order.


>>> cpl.plot_fdc(series, subplot=fsp, lw=2, ls='-', c='k', zorder=10)


Let's add the FDCs for each water year.
By default, the multiple curves are drawn using a linewidth of 0.8 and a grey color.


>>> cpl.plot_fdc(series, subplot=fsp, starting_month=4)


We can now finalize the plot by adding some labels and legends.

>>> fsp.set_ylabel("Flows (cfs)", fontweight='bold')
>>> fig.suptitle("Flow duration curve for the North Oconee River at Athens, GA",
                 fontweight="bold", fontsize=12)
>>> fig.savefig('athens_fdc.png')


.. plot:: source/plots/plot_fdc_noenso.py



Adding ENSO information
=======================

The :func:`~flowplots.plot_fdc` function accepts also :class:`~scikits.hydroclimpy.enso.ClimateSeries`
objects to add ENSO related information to the plot.
This option can be activated by setting the ``plot_enso`` optional boolean
flag to ``True``.
Naturally, an :class:`~scikits.hydroclimpy.enso.ENSOIndicator` must have been associated with the input series.

When the input series is 1D, individual points will plotted with the color
associated with their ENSO condition, 
provided a ``marker`` is explicitly given as an input parameter.

When the input series is 2D, each individual year will be associated 
with an ENSO condition, depending on the ``lag`` optional parameter.
By default, the ``lag`` parameter is calculated as the time difference between
the first month of the ENSO indicator reference 
(or October, if the ENSO indicator does not have any reference season) 
and the ``starting_month`` parameter.
For example, if ``starting_month`` is 4 (April) and the reference season runs
from October to December included, the lag will be 6 months, 
and each water year will be associated with the ENSO conditions prevailing 
the previous October.

In both cases, the colors corresponding to each ENSO conditions are available
through the :data:`~scikits.hydroclimpy.plotlib.ENSOcolors` dictionary.

Let us replot the previous FDCs with ENSO information.
First, we need to associate an ENSO indicator to the series.

>>> ONI = enso.load_oni()
>>> series = enso.set_ensoindicator(series, ONI)

We can now create a new figure and a new plot.

>>> fig = cpl.figure()
>>> cpl.plot_fdc(series, enso=True, marker='o')
>>> cpl.plot_fdc(series, enso=True, starting_month=4)
>>> cpl.gca().legend()
>>> fig.savefig('athens_fdc_enso.png')


.. plot:: source/plots/plot_fdc_enso.py



