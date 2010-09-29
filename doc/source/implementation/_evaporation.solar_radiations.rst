.. currentmodule:: scikits.hydroclimpy.lib.evaporation


The :class:`SolarInformation` class
===================================

.. autoclass:: SolarInformation



Attributes
----------

.. autosummary::
   :nosignatures:

   ~SolarInformation.freq
   ~SolarInformation.latitude
   ~SolarInformation.earth_sun_distance
   ~SolarInformation.solar_declination
   ~SolarInformation.sunset_angle
   ~SolarInformation.daylength
   ~SolarInformation.daylighthours
   ~SolarInformation.extraterrestrial_solar_radiations


Scalar attributes
~~~~~~~~~~~~~~~~~

.. attribute:: SolarInformation.freq

   Frequency specifier of the input dates.

   .. seealso::
      :attr:`scikits.timeseries.DateArray.freq`

.. autoattribute:: SolarInformation.latitude


Series attributes
~~~~~~~~~~~~~~~~~

The following attributes are :class:`~scikits.timeseries.TimeSeries` objects,
with the same :attr:`dates` and :attr:`freq` attributes as the instance.
If the initial frequency is less than daily (for example, monthly), the values
of the series correspond to daily values average over the period.

.. autoattribute:: SolarInformation.earth_sun_distance
.. autoattribute:: SolarInformation.solar_declination
.. autoattribute:: SolarInformation.sunset_angle
.. autoattribute:: SolarInformation.daylength
.. autoattribute:: SolarInformation.daylighthours
.. autoattribute:: SolarInformation.extraterrestrial_solar_radiations




Methods
-------

.. autosummary::
   :nosignatures:

   ~SolarInformation.frac_day_of_year
   ~SolarInformation.clearsky_solar_radiations
   ~SolarInformation.solar_radiations_from_sunshinehours
   ~SolarInformation.solar_radiations_from_temperatures


.. automethod:: SolarInformation.frac_day_of_year
.. automethod:: SolarInformation.clearsky_solar_radiations
.. automethod:: SolarInformation.solar_radiations_from_sunshinehours
.. automethod:: SolarInformation.solar_radiations_from_temperatures

