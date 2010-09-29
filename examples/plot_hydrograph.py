
"""
Introduction
============


The :mod:`scikits.hydroclimpy` module was designed to automatize some basic
operations in hydrology, such as drawing hydrographs.
In this example, we describe how to draw a hydrograph for the North Oconee
River in Athens, GA.
An hydrograph is the combination of two plots:
* the hydrograph itself, which plots the evolution of discharge in time;
* the hyetograph, which plots the time evolution of rainfall for the corresponding site.

First, we need to import some basic modules: :mod:`numpy`, :mod:`numpy.ma` and
of course, :mod:`scikits.hydroclimpy` itself.
"""
import numpy as np
import numpy.ma as ma
import scikits.hydroclimpy as hydro

"""
We will need to download the rainfall data from the Athens, GA weather station.
As this station is part of the COAPS network, we need to import the
:mod:`scikits.hydroclimpy.io.coaps` module.
We will also need discharge data from the USGS site, so we must also import the
:mod:`scikits.hydroclimpy.io.usgs` module.
"""
import scikits.hydroclimpy.io.coaps as coaps
import scikits.hydroclimpy.io.usgs as usgs

import matplotlib.pyplot as pyplot
import scikits.hydroclimpy.plotlib as cpl

#
"""
Importing the rainfall information
==================================

Let's import the rainfall data for Athens, GA.
As we already know from a previous :ref:`section <examples_wweatherdata>`,
the COAPS identification code for this is 90435.
The :func:`coaps.load_coaps_data` function returns a series with a structured
dtype, but we only need to take the ``'rain'`` field.
"""

weatherdata = coaps.load_coaps_data(90435)
rainfall = weatherdata['rain']

"""
We can check the frequency of the series and its dates range:

   >>> print rainfall.freqstr
   'D'
   >>> print rainfall.dates[[0,-1]]
   [13-Jan-1944 31-Dec-2007]
"""
print rainfall.freqstr
print rainfall.dates[[0,-1]]
"""

Importing the streamflow information
====================================

Let's import the streamflows recorded on the North Oconee River in Athens, GA.
We use the :func:`~scikits.hydroclimpy.io.usgs.load_usgs_flows` function, that requires
the identification code(s) of one or several USGS streamflow gages.
The corresponding code for the station of interest is '02217770'.
"""

flowdata = usgs.load_usgs_flows('02217770')

"""
Here also we can check the frequency and range of dates of the series:

   >>> print flowdata.freqstr
   D
   >>> print flowdata.dates[[0,-1]]
   [10-Aug-2002 14-Sep-2008]


Adjusting the series
====================


The rainfall and streamflows series do not have the same length.
Let's select the overlapping region of the two series.
We use the :func:`~scikits.timeseries.adjust_endpoints` function of the
:mod:`scikits.timeseries` package, that allows us to specify common starting and
ending points to the series.
This function takes two optional parameters, ``start_date`` and ``end_date``.
We will force the ``rainfall`` series to start on the same date as the ``flowdata`` series,
and the ``flowdata`` series to end at the same date as the ``rainfall`` series.

"""
rainfall = hydro.adjust_endpoints(rainfall, start_date=flowdata.dates[0])
flowdata = hydro.adjust_endpoints(flowdata, end_date=rainfall.dates[-1])

"""
As an alternative, we could have used the :func:`~scikits.timeseries.align_series`
function, that sets a sequence of series to the same range of dates.
By default, the series are extended to cover the widest range possible.


Plotting the hydrograph
=======================

Now that we have our two series, we can plot an hydrograph.
The :mod:`scikits.climpy` module provides a function for this specific purpose,
:func:`~scikits.climpy.plotlib.hydrograph`.
This function creates a new matplotlib Figure object, with two superposed plots:

* the hyetograph, that is, the plot of rainfall along time.
* the hydrograph per se, that is, the plot of stream flows along time.

The function requires two mandatory arguments, ``hyetodata`` and ``hydrodata``,
for the precipitation and streamflows series.
It also accepts all the standard parameters of a Figure object.

"""

fig = cpl.hydrograph(rainfall, flowdata, figsize=(12,6))

"""
The handle of the hyetograph subplot can be accessed through the :attr:`hyeto` 
attribute of the figure, and the handle of hydrograph through the :attr:`hydro`
attribute.

We can set the labels of the y-axis:
"""
fig.hyeto.set_ylabel("Rainfall (mm)", fontweight='bold')
fig.hydro.set_ylabel("Flows (cfs)", fontweight='bold')
fig.suptitle("Hydrograph for the North Oconee River at Athens, GA",
             fontweight="bold", fontsize=12)
fig.savefig("athens_hydrograph.png")

"""
As the hydrograph is a subclass of :class:`~scikits.timeseries.lib.plotlib.TimeSeriesFigure`,
the ticks on the x axis of each subplot are automatically adjusted with the level
of zoom.
Moreover, as the two subplots share the same x-axis, any modification on one
subplot is reflected on the other.
As an example, let's focus on the year 2005.

"""

fig.hyeto.set_datelimits('2005-01-01', '2005-12-31')
fig.savefig("athens_hydrograph_zoomed.png")
cpl.show()

