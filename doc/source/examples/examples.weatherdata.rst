
.. _examples_weatherdata:

======================
Importing weather data
======================

The Center for Ocean-Atmospheric Prediction Studies (COAPS) at Florida State 
University (FSU) stores weather-related information recorded by a network of 
stations over three Southern states: Alabama (AL), Florida (FL) and Georgia 
(GA).
The :mod:`scikits.hydroclimpy` provides some functions to access this information.
These functions are available in the :mod:`scikit.hydroclimpy.io.coaps` module.

   >>> import scikits.hydroclimpy.io.coaps as coaps


Finding the data to import
==========================
.. currentmodule:: scikits.hydroclimpy.io

Each station of the COAPS network is identified by a 5-digit integer code.
The identification code for a station can be found with the 
:func:`coaps.ids_bystate` function.
This function takes the 2-letter code of a state (``'AL'``, ``'FL'``, ``'GA'``) 
and returns a dictionary with station names as keys and stations ids as values.
Let's get the information for Georgia.

   >>> stationdict = coaps.ids_bystate('GA')

Let's find the id code corresponding to the station located in Athens, GA.
We must find the key matching the name ``'Athens'`` and take the corresponding value.
A fast solution is to use list completion:
we can loop on the items of the ``stationdict`` dictionary, selecting only the
values for which the key matches the string ``'Athens'``.

   >>> stationid = [v for (k, v) in stationdict.items() if 'Athens' in k.capitalize()]
   >>> print stationid
   [90435]

Thus, the COAPS id code for the station at Athens WSO Airport, GA is ``90435``.



Importing the data
==================

We can now use this identification code with the :func:`coaps.load_coaps_data`
function.
This function takes a COAPS id code as input parameter.

   >>> data = coaps.load_coaps_data(90435)

The function returns a new :class:`scikits.timeseries.TimeSeries` object 
with a four-field structured datatype:

   * ``rain`` reports the recorded precipitation (in mm);
   * ``tmin`` reports the minimum recorded temperature (in Celsius);
   * ``tmax`` reports the maximum recorded temperature (in Celsius);
   * ``tobs`` reports the number of temperature records for that particular day.

To access only the rainfall information, we must take the ``'rain'`` field.

   >>> rainfall = data['rain']

The object ``rainfall`` is itself a :class:`~scikits.timeseries.TimeSeries` object.
We can check its frequency by printing its :attr:`~scikits.timeseries.TimeSeries.freqstr`:

   >>> print rainfall.freqstr
   'D'

which means that the data are available at a ``D``\aily timestep.
We can find the range of available data by simply querying its
:attr:`~scikits.timeseries.TimeSeries.dates` attribute:

   >>> print rainfall.dates[[0,-1]]
   [13-Jan-1944 31-Dec-2007]



Converting to another frequency
===============================

The :meth:`~scikits.timeseries.TimeSeries.convert` method allows us to convert
the series to another frequency.
For example, if we need to calculate the cumulative precipitation over each month,
we can simply create a new series with a monthly frequency by using the method
and setting its optional parameter ``func`` to ``ma.sum`` (so that we can
handle potential missing data).

   >>> mrainfall = rainfall.convert('M',func=ma.sum).round(1)

We can also use this method to calculate the mean monthly precipitations
over the whole range of years.
First, we must convert the series to an annual frequency, but this time we will
not specify any value for the ``func`` parameter.
In that case, the :meth:`~scikits.timeseries.TimeSeries.convert` method groups 
the data by year and returns a 2D array.

   >>> arainfall = mrainfall.convert('A')
   >>> print arainfall.shape
   (64,12)

With the previous command, we verified that the series ``arainfall`` is indeed
2D, spanning 64 years.
Each column of ``arainfall`` corresponds to a given month: the first column to
January, the second to February, and so forth.
To compute the monthly means, we just have to use the standard :meth:`mean` method,
using ``axis=0`` as parameter to compute the mean per column:

   >>> monthly_means = arainfall.mean(axis=0).round(1)
   >>> print monthly_means
   [ 115.9  112.4  130.8   94.1   99.4  103.5  122.3   90.2   94.5   77.8
      92.7   98.6]

Note that we used the :meth:`round` method to round the monthly means to the
first decimal.


