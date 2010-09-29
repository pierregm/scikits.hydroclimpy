.. _installing:

************
Installation
************

Dependencies
============

Requirements
------------

In order to use the :mod:`scikits.hydroclimpy` package, the following external
packages must be installed beforehand:

Python_ 2.5 or later.
   Please note that Python_ 3 is not supported yet.

setuptools_
   :mod:`scikits` is a namespace package, and as a result every :mod:`scikit` requires setuptools_ to be installed to function properly.

Numpy_ 1.3.0 or later.
   Numpy_ is a library to manipulate large arrays of numerical data.

`scikits.timeseries <timeseries>`_
   The timeseries_ scikits is an extension to Numpy_ designed to manipulate 
   series indexed in time.

SciPy_ 0.7.0 or later:
   SciPy_ is a set of Numpy_\-based tools for engineering and scientific applications.
   Some of the sub-modules of scikits.timeseries_ use SciPy_ 
   interpolation and signal functions.

.. _Python: http://www.python.org/download/
.. _setuptools: http://pypi.python.org/pypi/setuptools
.. _Numpy: http://www.scipy.org/Download
.. _timeseries: http://pytseries.sourceforge.net/
.. _SciPy: http://www.scipy.org/Download


Recommended
-----------

The following packages are optional, but strongly recommended:

matplotlib_ 0.98.0 or later:
   matplotlib_ is a Python 2D plotting library.
   :mod:`scikits.hydroclimpy` includes some extensions to matplotlib_ to plot time series.

basemap_ 0.99 or later
   basemap_ is an extension to matplotlib_ designed to plot 2D data 
   on maps with different projections.

PyTables_ 2.0 or later:
   PyTables_ is a package for managing hierarchical datasets, 
   using the `HDF5 <http://www.hdfgroup.org/HDF5/>`_ format.

Optional
--------

BeautifulSoup_ :
   BeautifulSoup_ is a HTML/XML parser designed to quickly load information 
   stored in a parse tree.
   It is used by the :mod:`scikits.hydroclimpy.enso` module to retrieve information
   about the Oceanic Niño index.

xlrd_ :
   xlrd_ is a package to extract information from Microsoft Excel spreadsheets.
   This package is used by the :mod:`scikits.hydroclimpy.io.coaps` module to
   retrieve information about the stations of the COAPS network.

.. _matplotlib: http://matplotlib.sourceforge.net
.. _PyTables: http://www.pytables.org
.. _basemap: http://matplotlib.sourceforge.net/basemap/doc/html
.. _BeautifulSoup: http://www.crummy.com/software/BeautifulSoup/
.. _xlrd: http://www.lexicon.net/sjmachin/xlrd.htm





Download
========


The :mod:`scikits.hydroclimpy` module itself is currently available through subversion only.
You can download the latest version of the source files by checking out the repository with the command::

   svn co http://svn.scipy.org/svn/scikits/trunk/hydroclimpy hydroclimpy

This command will create a :file:`hydroclimpy` folder in the current directory.
On Windows, you can also use a SVN client such as `Tortoise SVN <http://tortoisesvn.net/>`_.



Installation
============

To install the :mod:`scikits.hydroclimpy` package, run the command::

    python setup.py install

in the directory you checked out the source code to.

If you run Python 2.6, note that you should probably use the ``--user`` or
``--home`` flags to install the package locally, without interfering with
the system installation.


Configuration
=============

The package requires a configuration file (``'hydroclimpyrc'``) to be imported.
An example of configuration file is provided with the sources and is copied in the installation directory.
This file should be modified by the user according to his/her needs.

When the package is imported, the file is looked for successively 
in the following locations successively:

   * the current working directory
   * the directory reprensented by the environment variable ``HYDROCLIMPYRC``
   * the ``$HOME/.hydroclimpy`` directory
   * the ``$PYTHONPATH/hydroclimpy`` directory

.. toctree::

   configurationfile.rst

