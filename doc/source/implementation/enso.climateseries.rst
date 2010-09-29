.. currentmodule:: scikits.hydroclimpy.enso

===========================================
Adding climate information to a time series
===========================================

.. inheritance-diagram:: scikits.hydroclimpy.enso.ClimateSeries




The :class:`ClimateSeries` class
===========================================

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ClimateSeries



Attributes
----------

As a subclass of :class:`~scikits.hydroclimpy.core.base.ReferencedSeries`, 
a :class:`ClimateSeries` instance
inherits all the :ref:`attributes <refseries_attributes>` of its parents.
In addition, it has the following specific attributes.


.. attribute:: ClimateSeries.ensoindicator

   Indicator of the ENSO phases, as a :class:`ENSOIndicator`
   object.


.. attribute:: ClimateSeries.ensoindices

   Shortcut to the :attr:`ENSOIndicator.ensoindices` attribute of the :attr:`ensoindicator`.
   This attribute is a :class:`~scikits.timeseries.TimeSeries` object, with values:

   * +1 for El Niño (warm) phases;
   *  0 for Neutral phases;
   * -1 for La Niña (cold) phases.


.. attribute:: ClimateSeries.cold

   Returns a copy of the series where values that do not fall during a cold
   (La Niña) episode are masked.


.. attribute:: ClimateSeries.neutral

   Returns a copy of the series where values that do not fall during a Neutral
   episode are masked.


.. attribute:: ClimateSeries.warm

   Returns a copy of the series where values that do not fall during a warm
   (El Niño) episode are masked.



Methods
-------

As a subclass of :class:`~scikits.hydroclimpy.core.base.ReferencedSeries`, 
a :class:`ClimateSeries` instance inherits 
all the :ref:`methods <refseries_methods>` of its parents.



Construction
============

To create a new :class:`ClimateSeries` object, 
the class can be directly called with the adequate input parameters.
However, it is also possible to use the constructor function :func:`climate_series`,
which accepts the same input parameters.

Another possibility is to transform a standard 
:class:`~scikits.timeseries.TimeSeries`
or :class:`~scikits.hydroclimpy.core.base.ReferencedSeries` 
into a :class:`ClimateSeries` with the :func:`set_ensoindicator` function.


.. autosummary::
   :nosignatures:
   :toctree: generated/

   climate_series
   set_ensoindicator




:class:`ClimateRecords`
==================================

A standard :class:`ClimateSeries` supports naturally named fields.
A specific field can be accessed as a regular item.
To make fields accessible as attributes, a :class:`ClimateRecords`
object can be defined.

The easiest way to define a :class:`ClimateRecords` is to take
a view of a :class:`ClimateSeries`.
A constructpr function, :func:`climate_records` can also be used.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ClimateRecords
   climate_records



Manipulating ENSO information
=============================

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ~ensobase.group_indices
   ~ensobase.grouped_indices_limits

   ~ensobase.apply_on_phase



Additional information
======================

For convenience, the :mod:`~scikits.hydroclimpy.enso` module also defines
extensions.

.. autoexception:: InvalidENSOError
   

