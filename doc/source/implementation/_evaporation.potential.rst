
.. currentmodule:: scikits.hydroclimpy.lib.evaporation


The :class:`PotentialEvapoTranspiration` class
==============================================

.. autoclass:: PotentialEvapoTranspiration
   :show-inheritance:


Attributes
----------

Because the :class:`PotentialEvapoTranspiration` class is a subclass of 
:class:`SolarInformation`, it inherits all its attributes 
(:attr:`~SolarInformation.freq`, 
:attr:`~SolarInformation.latitude`,
:attr:`~SolarInformation.extraterrestrial_solar_radiations`...).
In addition, the class has its own attributes:


.. attribute:: PotentialEvapoTranspiration.tmin

   Minimum air temperature series [ºC], as a 
   :class:`~scikits.timeseries.TimeSeries` object.


.. attribute:: PotentialEvapoTranspiration.tmax

   Maximum air temperatures series [ºC], as a 
   :class:`~scikits.timeseries.TimeSeries` object.


.. attribute:: PotentialEvapoTranspiration.tavg

   Average air temperatures series step [ºC], as a 
   :class:`~scikits.timeseries.TimeSeries` object.

.. autoattribute:: PotentialEvapoTranspiration.Kt
.. autoattribute:: PotentialEvapoTranspiration.lambdah
.. autoattribute:: PotentialEvapoTranspiration.Delta
.. autoattribute:: PotentialEvapoTranspiration.gamma
.. autoattribute:: PotentialEvapoTranspiration.saturated_vapor_pressure

.. attribute:: PotentialEvapoTranspiration.atm_pressure
   
   Atmospheric pressure [kPa].

.. autoattribute:: PotentialEvapoTranspiration.apx_solar_radiations
.. autoattribute:: PotentialEvapoTranspiration.net_lw_radiations
.. autoattribute:: PotentialEvapoTranspiration.soil_heat_flux_density


Methods
-------

.. autosummary::

   ~PotentialEvapoTranspiration.vapor_pressure_from_humidity
   ~PotentialEvapoTranspiration.vapor_pressure_from_dewpoint
   ~PotentialEvapoTranspiration.set_atmospheric_pressure

.. automethod:: PotentialEvapoTranspiration.vapor_pressure_from_humidity
.. automethod:: PotentialEvapoTranspiration.vapor_pressure_from_dewpoint
.. automethod:: PotentialEvapoTranspiration.set_atmospheric_pressure



Potential EvapoTranspiration Models
===================================

The :mod:`~scikits.climpy.lib.evaporation` module provides different models
for the computation of reference evapotranspiration.
These models are available as methods of the :class:`PotentialEvapoTranspiration`
class.


.. autosummary::
   :toctree: generated/

   ~PotentialEvapoTranspiration.PenmanMonteith
   ~PotentialEvapoTranspiration.PenmanMonteithASCE
   ~PotentialEvapoTranspiration.PriestleyTaylor
   ~PotentialEvapoTranspiration.Makkink
   ~PotentialEvapoTranspiration.Hansen
   ~PotentialEvapoTranspiration.Turc
   ~PotentialEvapoTranspiration.Hargreaves
   ~PotentialEvapoTranspiration.Droogers

