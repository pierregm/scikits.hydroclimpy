
.. currentmodule:: scikits.hydroclimpy.lib.evaporation


The :class:`PETEvaluation` class
=================================

For convenience, the :mod:`~scikits.hydroclimpy.lib.evaporation` module 
defines a class for the comparison of different PET models for a given site.


.. autoclass:: PETEvaluation
   :show-inheritance:


Attributes
----------

.. attribute:: PETEvaluation.estimates

   A structured :class:`~scikits.timeseries.TimeSeries` object, where each
   field stores the PET estimate for the given model.
   The models are listed in :attr:`~PETEvaluation.models`

.. attribute:: PETEvaluation.models

   A list of strings corresponding to the different PET models to analyze.
   By default, the list is 
   ``['Droogers', 'Hamon', 'Hansen', 'Hargreaves', 'Kharrufa', 'Makkink',
   'PenmanMonteith', 'PenmanMonteithASCE', 'PriestleyTaylor', 'Thornthwaite',
   'Turc']``.


.. attribute:: PETEvaluation.reference

   Name of the PET model to use as reference.
   By default, the model is ``'PenmanMonteith'``.


.. attribute:: PETEvaluation.ndtype

   A :class:`~numpy.dtype` object, equivalent to
   ``np.dtype([(m, float) for m in self.models])``.


Methods
-------

.. autosummary::
   :nosignatures:

   ~PETEvaluation.monthly_average
   ~PETEvaluation.RMSE
   ~PETEvaluation.RMSE_per_period
   ~PETEvaluation.kendalltau
   ~PETEvaluation.spearmanr


.. automethod:: PETEvaluation.monthly_average
.. automethod:: PETEvaluation.RMSE
.. automethod:: PETEvaluation.RMSE_per_period
.. automethod:: PETEvaluation.kendalltau
.. automethod:: PETEvaluation.spearmanr

