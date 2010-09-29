"""
.. currentmodule:: scikits.hydroclimpy.io.coaps

Utilities for importing weather data from the COAPS network.
The COAPS network covers the three states of Alabama (``AL``), Florida (``FL``) 
and Georgia (``GA``).


Stations information
====================

The following functions take the 2-letter code of a state and return some information
about the stations in the state:

.. function:: names_bystate(two_letter_state_code)

   Returns a dictionary with stations identification codes as keys and stations
   names as values.
   The COAPS station identification code is usually a 5-digit integer.


.. function:: ids_bystate(two_letter_state_code)

   Returns a dictionary with stations names as keys and stations ids codes as
   values.
   This is the reverse of :func:`names_bystate`.


.. function:: coordinates_bystate(two_letter_state_code)

   Returns a dictionary with stations identification codes as keys and stations
   coordinates as values.
   The coordinates are given in degrees, as a tuple (longitude, latitude).


.. function:: load_coaps_stationinfo(two_letter_state_code)

   Returns all the information relative to the stations' in the given state,
   as a :class:`numpy.recarray` object.

.. function:: load_coaps_networkinfo()

   Returns a structured array of the COAPS id and coordinates of all stations
   in the COAPS network.
   The output dtype is ``[('coop_id', '|S5'), ('lat', float), ('lon', float)]``.


Data Retrieval
==============

The following functions all accept the same input parameters:

* **station_id_code** : *int*
      Identification code of the COAPS station, as a 5-digit integer.
* **datadir** : *string*, optional
      Local directory where to store archived data
* **netfile_dir** : *string*, optional
      URL of the directory where to download the raw input data (ASCII).
* **headerlines** : *int*, optional
      Number of lines of the headr of the raw input data.
* **nodata** : float
      Number used to represent missing data.

The default values of the optional parameters are read from the configuration
file :file:`hydroclimpyrc`.


.. function:: load_coaps_stationdata(station_id_code, **options)

   Loads the temperatures and rainfall data for the current station.
   Returns a new :class:`~scikits.timeseries.TimeSeries` object of daily
   frequency, with a structured dtype:
   
   ``tobs`` : int
       Number of records per day.
   ``tmin`` : float [Celsius]
       Minimum recorded temperature on any given day.
   ``tmax`` : float [Celsius]
       Maximum recorded temperature on any given day.
    ``rain`` : float[mm]
       Total precipitation recorded on any given day.

.. autosummary::
   :toctree: generated/

   load_coaps_adjusted_networkdata
   load_coaps_networkdata
   load_coaps_period_networkdata
   load_coaps_stationdata


"""
__author__ = "Pierre GF Gerard-Marchant"

__all__ = ['stationinfo_bystate', 'stationinfo_HDF',
           'coordinates_bystate', 'names_bystate', 'ids_bystate',
           'load_coaps_stationinfo', 'load_coaps_networkinfo',
           'load_coaps_stationdata', 'load_coaps_networkdata']


import cPickle
import itertools
import logging
import os
import sys
import zipfile


import numpy as np

import scikits.timeseries as ts

from ..core import _config
#from ..core.tools import deprecated, deprecated_for

config = _config.configdict
coaps_config = dict(config.items('COAPS'))
coapslogger = logging.getLogger('io.coaps')
coapslogger = sys.stdout

try:
    import xlrd
    has_xlrd = True
except ImportError:
    coapslogger.write("The package 'xlrd' is not installed!")
    has_xlrd = False



def _check_archive(cfgdict):
    "Returns the name of the ZIP archive."
    archive = cfgdict['archive']
    if archive[-4:].lower() != '.zip':
        archive += '.zip'
    return archive


def _latlon_str2dec(latstr):
    """
    Converts a latitude in format DD:MM:SS (or DD MM SS ) to decimal.
    """
    if ':' in latstr:
        val = [float(_) for _ in latstr.split(':')]
    else:
        val = [float(_) for _ in latstr.split()]
    if len(val) == 2:
        (DD, MM) = val
        SS = 0
    elif len(val) == 3:
        (DD, MM, SS) = val
    else:
        raise ValueError("Unable to process the coordinate :%s" % latstr)
    return np.sign(DD) * (abs(DD) + MM / 60. + SS / 3600.)



#--- COAPS station information ------------------

def read_coaps_stationinfo_xls(filename):
    """
    Loads the data from filename (as XLS), and return a structured array.
    """
    if isinstance(filename, file):
        filename = filename.name
    book = xlrd.open_workbook(filename)
    sheet = book.sheet_by_index(0)
    # Find the header...
    row = 0
    while not (np.array(sheet.row_types(row)) == 1).all():
        row += 1
    header = [r.upper().replace(' ', '_').replace('.', '')
                  for r in sheet.row_values(row)]
    # Find the first non-null row from then
    row += 1
    while not np.nonzero(sheet.row_types(row))[0].size:
        row += 1
    # Load the data in an array ................
    data = []
    for i in range(sheet.ncols):
        data.append(np.array(sheet.col_values(i, row)))
    ndtype = [(str(name), a.dtype.str) for (name, a) in zip(header, data)]
    result = np.empty((sheet.nrows - row,), dtype=ndtype)
    for (name, a) in zip(header, data):
        result[name] = a
    if 'LONG_DEC' in header and (result['LONG_DEC'] > 0).all():
        result['LONG_DEC'] *= -1
    return result


def read_coaps_stationinfo_csv(filename):
    """
    Loads the data from filename (as CSV), and return a new structured array.
    """
    raise NotImplementedError
#    import csv
#    reader = csv.reader(open(filename), delimiter=',', 
#                        quoting=csv.QUOTE_NONNUMERIC)
#    # Find the header ....
#    row = reader.next()
#    while allstrings(row) and min([len(r) for r in row]) == 0:
#        row = reader.next()
#    assert(allstrings(row),"Cannot define a header from the current row!"\
#                           "(%s)" % str(row))
#    header = [r.upper().replace(' ','_').replace('.','') for r in row]
#    # Find the firsts non-empty record .........
#    row = reader.next()
#    while allstrings(row) and min([len(r) for r in row]) > 0:
#        row = reader.next()
#    # Load the data in an array ................
#    data = np.array([r for r in reader])
#    data = [convert_columns(col) for col in data.T]
#    data = np.core.records.fromarrays(data, names=header)
#    if hasattr(data, 'LONG_DEC')  and (data.LONG_DEC > 0).all():
#        data.LONG_DEC *= -1       
#    return data


def stationinfo_bystate(state, cfgdict=coaps_config):
    """
    Loads the information about the COAPS stations for the selected state.

    If the information is not already available (archived), reads the 
    corresponding raw data.
    Returns a numpy.recarray.
    
    Parameters
    ----------
    state : string
        Abbreviation of the state where the stations are based.
        - AL : Alabama
        - FL : Florida
        - GA : Georgia

    Notes
    -----
    The function tries to open a `XXcntylist.xls` file in the archive directory,
    with :var:`XX` the two-letter code for the state.
    If the file is not available, or if the :mod:`xlrd` module is not available,
    a :exception:`IOError` exception is raised.

    """
    # Force state to a 2-letter state code
    state = str(state).upper()[:2]
    if not state in ('AL', 'FL', 'GA'):
        errmsg = "The input should be a valid two-letter state code.\n"\
                 "('AL', 'FL', 'GA')"
        raise ValueError(errmsg)
    #
    coapsarchive = _check_archive(cfgdict)
    archfilename = '%s_stationinfo' % state
    datadir = cfgdict['datadir']
    try:
        zipf = zipfile.ZipFile(coapsarchive, 'r')
        stationinfo = cPickle.loads(zipf.read(archfilename))
    except IOError:
        zipf = zipfile.ZipFile(coapsarchive, 'w')
    except KeyError:
        zipf = zipfile.ZipFile(coapsarchive, 'a')
    else:
        return stationinfo

    # Get the COAPS info file
    if has_xlrd:
        filename = '%scntylist.xls' % state
        netdir = os.path.join(cfgdict['netfile.dir'], '..')
        netfile = os.path.join(netdir, filename)
        sourcedir = np.lib._datasource.DataSource(netdir)
        xlsfile = sourcedir.open(netfile)
        stationinfo = read_coaps_stationinfo_xls(xlsfile)
    else:
        filename = '%scntylist.csv' % state
        sourcedir = np.lib._datasource.DataSource(datadir)
        try:
            csvfile = sourcedir.open(filename)
        except IOError:
            logging.exception(":(")
            raise
        stationinfo = read_coaps_stationinfo_csv(csvfile)
    #
    zipf.writestr(archfilename, cPickle.dumps(stationinfo))
    zipf.close()
    return stationinfo


def coordinates_bystate(state, cfgdict=coaps_config):
    """
    Returns a dictionary of COOP_ID <> (Lon.,Lat.)
    """
    data = stationinfo_bystate(state, cfgdict=cfgdict)
    # Get the station ids
    try:
        coopid = data['COOP_ID']
    except ValueError:
        coopid = data['COOPID']
    strgsize = np.ceil(np.log10(max(coopid)))
    # Get the longitudes
    try:
        lon = data['LONG_DEC']
    except ValueError:
        lon = np.array([[float(_) for _ in s.split(':')]
                        for s in data['LONGITUDE']])
        lon[:, 1] /= 60.
        lon[:, 1] *= np.sign(lon[:, 0])
        lon = lon.sum(axis= -1)
    # Get the latitudes
    try:
        lat = data['LATITUDE_DEC']
    except ValueError:
        data_latitude = data['LATITUDE']
        if data_latitude[0].dtype.char in 'SU':
            lat = np.array([[float(_) for _ in s.split(':')]
                            for s in data_latitude])
            lat[:, 1] /= 60.
            lat = lat.sum(axis= -1)
        else:
            lat = data_latitude * 24
    return dict(zip(coopid.astype("|S%s" % strgsize), zip(lon, lat)))


def names_bystate(state, cfgdict=coaps_config):
    """
    Returns a dictionary (COOP id <> station name).

    Parameters
    ----------
    state : string
        Two-letter abbreviation of the state.

    See Also
    --------
    ids_bystate
        Returns a dictionary (station name <> COOP id)`
    """
    data = stationinfo_bystate(state, cfgdict=cfgdict)
    try:
        coopid = data['COOP_ID']
    except ValueError:
        coopid = data['COOPID']
    strgsize = np.ceil(np.log10(max(coopid)))
    return dict(zip(coopid.astype("|S%s" % strgsize), data['STATION_NAME']))


def ids_bystate(state, cfgdict=coaps_config):
    """
    Returns a dictionary of (station name  <> COOP id)

    Parameters
    ----------
    state : string
        Two-letter abbreviation of the state.

    See Also
    --------
    names_bystate
        Returns a dictionary (COOP id <> station name)
    """
    data = stationinfo_bystate(state, cfgdict=cfgdict)
    try:
        coopid = data['COOP_ID']
    except ValueError:
        coopid = data['COOPID']
    return dict(zip(data['STATION_NAME'], coopid.astype(int)))



def stationinfo_HDF(cfgdict=coaps_config):
    """
    Load information for the COAPS stations from a HDF5 file.
    If the file is not available, it is automatically created. 
    """
    try:
        import tables
    except ImportError:
        err_msg = "Function 'stationinfo_HDF' unavailable: "\
                  "the package pytables is not installed."
        raise NotImplementedError(err_msg)
    #----------------------------------------
    Stations = dict(name=tables.StringCol(30),
                    coop_id=tables.Int32Col(),
                    crop_id=tables.StringCol(5),
                    latitude=tables.Float64Col(),
                    longitude=tables.Float64Col(),
                    county_within=tables.StringCol(15),
                    county_represent=tables.StringCol(60))
    compressor = tables.Filters(complevel=1, complib='zlib')
    hdf5name = os.path.join(cfgdict['datadir'], "COAPS_info.hdf5")
    try:
        h5file = tables.openFile(hdf5name, 'r', filters=compressor)
    except IOError:
        coapslogger.info("Unable to find the file %s :" \
                          "Attempting creation..." % hdf5name)
        h5file = tables.openFile(hdf5name, 'w', filters=compressor)
        table = h5file.createTable("/", 'stations', Stations, "Station info")
        station = table.row
        #
        iterables = itertools.chain(load_coaps_stationinfo('AL'),
                                    load_coaps_stationinfo('FL'),
                                    load_coaps_stationinfo('GA'))
        lastid = 0
        for sinfo in iterables:
            try:
                coopid = sinfo['COOPID']
            except ValueError:
                coopid = sinfo['COOP_ID']
            if coopid != lastid:
                if lastid > 0:
                    station.append()
                lastid = coopid
                station['name'] = sinfo['STATION_NAME'].title()
                station['coop_id'] = coopid
                try:
                    station['crop_id'] = sinfo['CROPID']
                except ValueError:
                    station['crop_id'] = sinfo['CROP_ID']
                station['latitude'] = _latlon_str2dec(sinfo['LATITUDE'])
                station['longitude'] = _latlon_str2dec(sinfo['LONGITUDE'])
                if station['longitude'] > 0:
                    station['longitude'] *= -1
                station['county_within'] = sinfo['COUNTY_WITHIN'].title()
                station['county_represent'] = sinfo['COUNTY_REPRESENT'].title()
            else:
                station['county_represent'] += ', %s' % \
                                               sinfo['COUNTY_REPRESENT'].title()
        table.flush()
    return h5file


def load_coaps_stationinfo(stationcode, cfgdict=coaps_config):
    """
    Load the information relative to the station with ID code `stationcode`
    """
    code = int(stationcode)
    statenb = code // 10000
    if statenb == 9:
        statecode = 'GA'
    elif statenb == 8:
        statecode = 'FL'
    elif statenb == 1:
        statecode = 'AL'
    else:
        errmsg = "Unable to recognize station code '%s'" % code
        raise ValueError(errmsg)
    stateinfo = stationinfo_bystate(statecode, cfgdict=cfgdict)
    try:
        selector = (stateinfo['COOP_ID'] == code)
    except ValueError:
        selector = (stateinfo['COOPID'] == code)
    return stateinfo[selector]


def load_coaps_networkinfo(cfgdict=coaps_config):
    """
    Load stations information for all the stations of the COAPS network.

    Parameters
    ----------
    cfgdict : dictionary, optional
        Dictionary of options

    Returns
    -------
    astations_coords
        A structured ndarray of station coordinates, with fields 
        `coop_id`, `lat` and `lon`.
    """
    # Read the coordinates from the state files and store them in a dictionary
    stations_coords = {}
    for state in ('AL', 'FL', 'GA'):
        stations_coords.update(coordinates_bystate(state, cfgdict=cfgdict))
    # Construct a structured array for the coordinates
    station_dtype = [('coop_id', '|S5'), ('lat', float), ('lon', float)]
    astations_coords = np.empty(len(stations_coords), dtype=station_dtype)
    astations_coords['coop_id'] = sorted(stations_coords)
    (astations_coords['lon'], astations_coords['lat']) = \
                            zip(*(stations_coords[k] for k in stations_coords))
    return astations_coords



#---

def load_coaps_stationdata(station, cfgdict=coaps_config):
    """
    Loads the temperatures and rainfall data for the current COAPS station.

    Parameters
    ----------
    station : int
        COAPS station identification code.
    cfgdict : dictionary
        Dictionary of configuration options. 
        By default, the configuration options are read from the :file:`.climpyrc`
        file
    
    Returns
    -------
    stationdata : ClimateRecords
        A :class:`~scikits.climpy.core.ClimateRecords` object, with fields:
        tobs : int
            Number of temperature observations for the day.
        tmin : float
            Minimum recorded temperature (oC).
        tmax : float
            Maximum recorded temperature (oC).
        rain : float
            Recorded amount of precipitation (mm).

    """
    #  
    datadir = cfgdict['datadir']
    coapsarchive = _check_archive(cfgdict)
    #
    archfilename = "COAPS%s" % station
    try:
        zipf = zipfile.ZipFile(coapsarchive, 'r')
        stationdata = cPickle.loads(zipf.read(archfilename))
    except IOError:
        zipf = zipfile.ZipFile(coapsarchive, 'w')
    except (KeyError, ValueError):
        zipf = zipfile.ZipFile(coapsarchive, 'a')
    else:
        return stationdata
    #
    netfile = os.path.join(cfgdict['netfile.dir'], "%sascii.txt" % station)
    sourcedir = np.lib._datasource.DataSource(datadir)
    dfile = sourcedir.open(netfile)

    skip = int(cfgdict['headerlines'])
    missing_values = cfgdict['nodata']
    dateconverter = lambda y, m, d : ts.Date('D',
                                             year=int(y), month=int(m),
                                             day=int(d))
    converters = {4 : lambda s: np.round((float(s) - 32) / 1.8, 1),
                  5 : lambda s: np.round((float(s) - 32) / 1.8, 1),
                  6 : lambda s: 0.254 * float(s)}
    ndtype = [('nobs', int), ('tmin', float), ('tmax', float), ('rain', float)]
    series = ts.tsfromtxt(dfile, delimiter="\t", skiprows=skip,
                          missing=missing_values, converters=converters,
                          datecols=(0, 1, 2), dateconverter=dateconverter)
    series = series.view(ndtype)
    # Store metainfo
    metainfo = []
    dfile.seek(0)
    for i in range(skip - 2):
        row = dfile.readline().strip()
        if not len(row):
            continue
        info = tuple([_.strip() for _ in row.split(":")])
        if len(info) == 1:
            info = (info, '')
        metainfo.append(info)
    series._optinfo.update(metainfo=dict(metainfo))

    zipf.writestr(archfilename, cPickle.dumps(series))
    zipf.close()
    return series



def load_coaps_networkdata(start_date=None, end_date=None,
                           cfgdict=coaps_config):
    """
    Load data for all the COAPS stations.

    If both `start_date` or `end_date` are not None, the series are adjusted
    to the same date range defined by `start_date` and `end_date`.

    Parameters
    ----------
    start_date : Date, str, optional
        Starting date of the series.
    end_date : Date, str, optional
        Ending date of the series.
    coaps_config
        Dictionary of options

    Returns
    -------
    stationdict
        A dictionary of daily data from the COAPS network.
        Keys are COAPS station ids.

    See Also
    --------
    load_coaps_adjusted_networkdata
        Equivalent function where the series have all the same starting and
        ending dates
    """
    stationids = {}
    for state in ('AL', 'FL', 'GA'):
        stationids.update(names_bystate(state, cfgdict=cfgdict))
    # Create an empty dictionary
    stationdatadict = {}
    prompt = "Loading station information : [%5s]\r"
    kwargs = dict(start_date=start_date, end_date=end_date)
    # Loop on the station ids and fill the dictionary
    if (start_date is not None) or (end_date is not None):
        for sid in stationids:
            coapslogger.write(prompt % sid)
            coapslogger.flush()
            current = load_coaps_stationdata(sid, cfgdict=cfgdict)
            current = current.fill_missing_dates()
            stationdatadict[sid] = current.adjust_endpoints(**kwargs)
    else:
        for sid in stationids:
            coapslogger.write(prompt % sid)
            coapslogger.flush()
            stationdatadict[sid] = load_coaps_stationdata(sid, cfgdict=cfgdict)
    return stationdatadict



def load_coaps_adjusted_networkdata(start_date=None, end_date=None,
                                    cfgdict=coaps_config):
    """
    Load the data for all COAPS stations.

    The series are adjusted to the same `start_date` and `end_date`.
    If `start_date` is None, the earliest date is used.
    if `end_date` is None, the latest date is used.

    Parameters
    ----------
    start_date : Date, str, optional
        Starting date of the series.
    end_date : Date, str, optional
        Ending date of the series.

    Returns
    -------
    stationdict
        A dictionary of daily data from the COAPS network.
        Keys are COAPS station ids.

    See Also
    --------
    load_coaps_networkdata
        Equivalent function where the series are not necessarily adjusted.
    """
    datadict = load_coaps_networkdata(start_date=start_date, end_date=end_date,
                                      cfgdict=cfgdict)
    readjust = False
    # No starting date ? Get the earliest from the dataset
    if start_date is None:
        start_date = np.min([d.dates.min() for d in datadict.values()])
        readjust = True
    # No ending date ? Get the latest from the dataset
    if end_date is None:
        end_date = np.max([d.dates.max() for d in datadict.values()])
        readjust = True
    # Readjust the dataset to the same range
    if readjust:
        for (coopid, data) in datadict.items():
            _data = data.fill_missing_dates()
            datadict[coopid] = _data.adjust_endpoints(start_date=start_date,
                                                      end_date=end_date)
    return datadict


def load_coaps_period_networkdata(field, freq=None, func=None,
                                  start_date=None, end_date=None,
                                  cfgdict=coaps_config):
    """
    Load data converted the given period for all the stations

    Parameters
    ----------
    field : str
        Type of data to select.
        Must be one of ('tmin', 'tmax', 'rain')
    freq : var, optional
        Period to convert the dataset to.
    func : function, optional
        Function with which to convert the dataset.
        The function must output a 1D dataset.
    start_date : var, optional
        Starting date of the dataset.
    end_date : var, optional
        Ending date of the dataset.

    Returns
    -------
    period_networkdata
        Structured array of the converted data for all stations of the network.
        The output dtype is ``[(station_id, float)]`` for all the station ids.
    """
    # Make sure we have a valid data type
    valid_field = ('tmin', 'tmax', 'rain')
    if field not in valid_field:
        errmsg = "Invalid datatype: should be in %s.\n(got '%s')"
        raise ValueError(errmsg % (valid_field, field))
    # Check the frequency
    freq = ts.check_freq(freq)
    # Load the dictionary of adjusted data
    datadict = load_coaps_adjusted_networkdata(start_date=start_date,
                                               end_date=end_date,
                                               cfgdict=cfgdict)
    # Define the list of station ids
    coaps_ids = datadict.keys()
    # Define the output dtype
    ndtype = [("%s" % _, float) for _ in sorted(coaps_ids)]
    # Choose one series as reference and convert it
    reference = datadict[coaps_ids[0]]
    reference = reference[field].convert(freq, func=func)
    # Exit if we don't have a 1D series
    if reference.ndim != 1:
        errmsg = "Conversion error: the output dataset should be 1D.\n"\
                 "(got a %iD series instead)"
        raise TypeError(errmsg % reference.ndim)
    series = ts.time_series(np.empty(len(reference), dtype=ndtype),
                            dates=reference.dates)
    series_values = series.series
    for (id_, data) in datadict.items():
        series_values[id_] = data[field].convert(freq, func=func).series
    return series



