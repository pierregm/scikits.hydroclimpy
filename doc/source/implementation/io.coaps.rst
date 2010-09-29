.. currentmodule:: scikits.hydroclimpy.io

========================================
Importing weather information from COAPS
========================================

The Center for Ocean-Atmospheric Prediction Studies (COAPS) at Florida State 
University (FSU) stores weather-related information recorded by a network of 
stations over three Southern states: Alabama (``AL``), Florida (``FL``) 
and Georgia (``GA``).

A description of the stations for each state is available by FTP.
At the time of writing (Feb. 2009), this information is stored into three 
Microsoft Excel files.
This information can be directly accessed in Python with the following 
functions, provided that the external module xlrd_ is installed.
If this module is not available, an ASCII version of the file must be available.

.. _xlrd: www.lexicon.net/sjmachin/xlrd.htm


Loading stations information
============================

The following functions take the 2-letter code of a state and return some information
about the stations in the state:

.. function:: coaps.names_bystate(two_letter_state_code)

   Returns a dictionary with stations identification codes as keys and stations
   names as values.
   The COAPS station identification code is usually a 5-digit integer.


.. function:: coaps.ids_bystate(two_letter_state_code)

   Returns a dictionary with stations names as keys and stations ids codes as
   values.
   This is the reverse of :func:`names_bystate`.


.. function:: coaps.coordinates_bystate(two_letter_state_code)

   Returns a dictionary with stations identification codes as keys and stations
   coordinates as values.
   The coordinates are given in degrees, as a tuple (longitude, latitude).


.. function:: coaps.stationinfo_bystate(two_letter_state_code)

   Returns all the information relative to the stations in the given state,
   as a :class:`numpy.recarray` object.


.. function:: coaps.stationinfo_HDF

   Returns all the information relative to the stations in any of the three
   states.
   This information is stored into a HDF5 file.
   The optional module pytables_ must be installed.

   .. _pytables: www.pytables.org/


.. function:: coaps.load_stationinfo(stationcode)

   Returns the information relative to the station whose id code is `stationcode`.
   The input parameter must be a valid 5-digit COAPS station id code.
   The result is a structured array, containing:

   * the station name (``'STATION_NAME'``);
   * the station id (``'COOP_ID'`` or ``'COOPID'``);
   * an indentifier for the kind of crop (``'CROP_ID'``);
   * the coordinates of the station, as strings (``'LATITUDE'`` and 
     ``'LONGITUDE'``) and floats (``'LATITUDE_DEC'`` and ``'LONGITUDE_DEC'``);
   * the name of the county the station is in (``'COUNTY_WITHIN'``);
   * the names of the counties the station represents (``'COUNTY_REPRESENT'``).



Loading weather information
===========================

The following function accepts the following input parameters:

* **station_id_code** : *int*
      Identification code of the COAPS station, as a 5-digit integer.
* **datadir** : *string*, optional
      Local directory where to store archived data
* **netfile_dir** : *string*, optional
      URL of the directory where to download the raw input data (ASCII).
* **headerlines** : *int*, optional
      Number of lines of the header of the raw input data.
* **nodata** : float
      Number used to represent missing data.

The default values of the optional parameters are read from the configuration
file :file:`hydroclimpyrc`.


.. function:: coaps.load_coaps_data(station_id_code, **options)

   Loads the temperatures and rainfall data for the current station.
   Returns a new :class:`~scikits.timeseries.TimeSeries` object of daily
   frequency, with a structured dtype:
   
   ``tobs`` : int
       Number of records per day.

   ``tmin`` : float [Celsius]
       Minimum recorded temperature on any given day.

   ``tmax`` : float [Celsius]
       Maximum recorded temperature on any given day.

    ``rainfall`` : float [mm]
       Total precipitation recorded on any given day.

