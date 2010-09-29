.. currentmodule:: scikits.timeseries


***********
Information
***********


Overview
========

About
-----

History
~~~~~~~

The :mod:`scikits.hydroclimpy` module was originally developed by 
Pierre Gérard-Marchant to analyze the influence of climate variations (ENSO) 
on streamflows in the Southeastern US.

Part of the initial code related to the manipulation of series with missing 
data.
This code was eventually incorporated in :mod:`numpy` as the :mod:`numpy.ma` 
module.
Another part of the code was related to the handling of series indexed in time,
and was merged with the :mod:`scikits.timeseries` code.


Thanks
~~~~~~

Code development was sponsored by the 
`Southeast Climate Consortium <http://www.seclimate.org/>`_ and the
`University of Georgia <http://www.engineering.uga.edu/>`_.



Legalese
========

License 
-------

The :mod:`scikits.hydroclimpy` scikit is free for both commercial and 
under the BSD_ license.


.. _BSD: http://www.opensource.org/licenses/bsd-license.php



Support
-------

No commercial support is offered for the :mod:`scikits.hydroclimpy` module.
Requests for help should be directed to 
`Pierre Gérard-Marchant <pierregmcode@gmail.com>`_



Conventions
===========

By convention, the following imports are assumed throughout the documentation::

   >>> import numpy as np
   >>> import numpy.ma as ma
   >>> import datetime
   >>> import scikits.timeseries as ts

Input variables and keywords are represented in :keyword:`this type`. 

