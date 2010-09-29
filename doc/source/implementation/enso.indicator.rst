.. currentmodule:: scikits.hydroclimpy.enso


Representing the ENSO signal: the :class:`ENSOIndicator` class
====================================================================


.. autoclass:: ENSOIndicator
   :show-inheritance:


Attributes
----------

Public attributes
~~~~~~~~~~~~~~~~~

As a subclass of :class:`~scikits.hydroclimpy.core.base.ReferencedSeries`, 
:class:`ENSOIndicator` inherits all its attributes, 
as described in the :ref:`attributes <refseries_attributes>` section.
In addition, each instance of the class has the following specific attributes:

.. attribute:: ENSOIndicator.thresholds

   Characteristic thresholds for the definition of ENSO phases.
   By default, the thresholds correspond to the lower and upper quartiles.
   Setting the attribute to ``None`` reverts the thresholds to these defaults.

   .. note:: Modifying this value resets any cached indices previously computed.


.. attribute:: ENSOIndicator.full_year

   Boolean indicating whether indices are calculated for a whole year
   (twelve consecutive months) or until the conditions are met.


.. attribute:: ENSOIndicator.ensotag

   Short string describing the indicator.


.. attribute:: ENSOIndicator.minimum_size

   Returns the minimum numbers of consecutive months for the definition of
   ENSO indices.

   .. note:: Modifying this value resets any cached indices previously computed.


.. attribute:: ENSOIndicator.reference_season

   Reference season for the definition of ENSO indices.

   .. note:: Modifying this value resets any cached indices previously computed.
 

.. attribute:: ENSOIndicator.refseason

   Alias for :attr:`reference_season`.


.. attribute:: ENSOIndicator.ensoindices

   Indices corresponding to each ENSO phase:
     * +1 for El Niño (warm) phases;
     *  0 for Neutral phases;
     * -1 for La Niña (cold) phases.

.. attribute:: ENSOIndicator.indices

   Alias to :attr:`ensoindices`.



Private attributes
~~~~~~~~~~~~~~~~~~

The following attributes are *private*, and are not meant to be used directly.
They are described here for the sake of completion.


.. attribute:: ENSOIndicator._cachedmonthly

   Dictionary storing the values of monthly (``'indices_monthly'``) 
   or annual (``'indices_annual'``) indices.


.. attribute:: ENSOIndicator._cachedcurrent

   Dictionary storing the values of the indices for the current set of options.


.. attribute:: ENSOIndicator._cachedclustered

   Stores groups of consecutive indices, 
   as a :class:`~scikits.hydrcolimpy.core.tools.Cluster` object.


.. attribute:: ENSOIndicator._optinfo

   Dictionary storing various attribute values.



Methods
-------

In addition to all the methods inherited from :class:`~scikits.hydroclimpy.ReferencedSeries` and
described in the :ref:`refseries_attributes` section, any instance
of :class:`~scikits.hydroclimpy.enso.ensobase.ENSOIndicator` has the following specific methods.

.. .. autosummary::
..    :nosignatures:
.. 
..    ENSOIndicator.set_indices
..    ENSOIndicator.set_monthly_indices
..    ENSOIndicator.set_annual_indices

=========================================== ====================================================
:class:`~ENSOIndicator.set_indices`         Sets the indices.
:class:`~ENSOIndicator.set_monthly_indices` Sets the indices per month.
:class:`~ENSOIndicator.set_annual_indices`  Sets the indices per group of 12 consecutive months,
                                            starting at the first month of the reference.
=========================================== ====================================================

_____

.. automethod:: ENSOIndicator.set_indices
.. automethod:: ENSOIndicator.set_monthly_indices
.. automethod:: ENSOIndicator.set_annual_indices

