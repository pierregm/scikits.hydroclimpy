
.. currentmodule:: scikits.hydroclimpy

================
Climate analysis
================

The :mod:`scikits.hydroclimpy` module defines several objects to keep track
of climate information by associating ENSO information to a time series:

* The :class:`enso.ENSOIndicator` stores the information required for the definition
  of ENSO phases.
* The :class:`enso.ClimateSeries` class is a subclass of :class:`~scikits.timeseries.TimeSeries`
  with a special :attr:`~ClimateSeries.ensoindices` attribute.

These objects are available through the :mod:`scikits.hydroclimpy.enso` package.

   >>> import scikits.hydroclimpy.enso as enso


Setting ENSO information
========================

Let us now add some climate information to ``mrainfall``, the time series 
of precipitation for Athens, GA defined
in the :ref:`previous <examples_weatherdata>` section.

First, we need to transform ``mrainfall`` from a regular :class:`~scikits.timeseries.TimeSeries`
to a :class:`enso.ClimateSeries` object, by taking a view.

   >>> mrainfall = mrainfall.view(enso.ClimateSeries)

We must now allocate an ENSO indicator to the new object.
Two ENSO indicators are available in :mod:`~scikits.hydroclimpy`:

- The Oceanic Niño Index (ONI) is based on the 3-month running mean of the 
  sea surface temperatures (SST) anomalies recorded over the Ni\~no-3.4 region of
  the tropical Pacific Ocean (5N-5S, 120W-170W). This indicator is used by the
  US Climate Prediction Center 
  (`CPC <http://www.cpc.noaa.gov/products/monitoring_and_data/ENSO_connections.shtml>`_)
  for its climate outputs.
  To load the ONI, use the :func:`enso.load_oni` function.
  A detailed description of the optional input parameters of this function is
  given in the :ref:`page <enso_data>` describing 
  :class:`~scikits.hydroclimpy.enso.ENSOIndicator` objects.

- The Japanese Meteorological Association (JMA) provides another ENSO indicator,
  based on the SST anomalies record over the Niño-3 region of the Pacific Ocean
  (5N-5S, 150W-90W).
  To load the JMA index, use the :func:`enso.load_jma` function.

As an example, we will use the CPC definition of ENSO phases.
The corresponding :class:`~scikits.hydroclimpy.enso.ENSOIndicator` is defined as

   >>> ONI = enso.load_oni()

Note that by default, no ENSO information is available before January 1950.

   >>> print ONI.dates[[0,-1]]
   [Jan-1950 Dec-2008]



We can now set the :attr:`~scikits.hydroclimpy.enso.ClimateSeries.ensoindicator`
attribute of our ``mrainfall`` series of monthly cumulative precipitation.

   >>> mrainfall.ensoindicator = ONI
   >>> print mrainfall.ensoindicator.dates[[0,-1]]
   [Jan-1944 Dec-2007]

We verified that the ``ONI`` series is automatically adjusted to the same range
of dates as the ``mrainfall`` series.
The missing information (from January 1944 to December 1949 included) has been
automatically replaced by missing data.


Analyzing the climate information
=================================

Now that we have defined an :class:`~scikits.hydroclimpy.enso.ENSOIndicator`, 
we can access the data falling in any ENSO episode with the corresponding 
attributes:

* :attr:`~scikits.hydroclimpy.enso.ClimateSeries.cold` for La Niña episodes.
* :attr:`~scikits.hydroclimpy.enso.ClimateSeries.neutral` for Neutral episodes.
* :attr:`~scikits.hydroclimpy.enso.ClimateSeries.warm` for El Niño episodes.

These three attributes return each the initial series masked where the series
does not fall during the corresponding episodes.
For example, according to the CPC, 2000 was a year where all months were categorized
as falling in La Niña, except for the three months July, August and September that
were categorized as Neutral.
We can verify that our computations are valid by selecting the values of 
``mrainfall`` falling in 2000:

   >>> mrainfall_2K = mrainfall[(mrainfall.year == 2000)]
   >>> print mrainfall_2K.cold
   [113.8 50.5 86.6 43.2 54.9 50.3 -- -- -- 5.8 106.7 87.9]
   >>> print mrainfall_2K.neutral
   [-- -- -- -- -- -- 85.3 93.2 122.2 -- -- --]
   >>> print mrainfall_2K.warm
   [-- -- -- -- -- -- -- -- -- -- -- --]

We can see that indeed the period July-September 2000 fell during a Neutral episode,
and that the cumulative precipitation were 85.3, 93.2 and 122.2 mm respectively.
The other months of 2000 were all characterized as falling in El Niño episodes.


We can perform some analysis on the rainfall data, such as calculating the mean
monthly precipitation for each ENSO phase.
The most efficient way to compare the information is to create an empty masked
array with a flexible dtype, and to fill it as needed.
First, let us define a flexible datatype, with each phase as a field and an
additional field for the 'global' monthly means (`ie` without ENSO information).

   >>> mdtype = [('cold',float), ('neutral',float), ('warm',float), ('global',float)]

Now, let us create an empty masked array.
We already know that this array must have a length of 12 (one entry for each month)

   >>> monthly_means_enso = ma.empty(12, dtype=mdtype)

We can already fill the ``'global'`` field.

   >>> monthly_means_enso['global'] = monthly_means

Now, we can fill the other phases using a ``for`` loop

   >>> for phase in ('cold', 'neutral', 'warm'):
   >>>     mcurrent = getattr(mrainfall, phase)
   >>>     acurrent = mcurrent.convert('A')
   >>>     monthly_means_enso[phase] = acurrent.mean(axis=0).round(1)


