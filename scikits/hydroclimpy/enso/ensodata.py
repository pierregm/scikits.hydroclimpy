# -*- coding: utf-8 -*-
u"""
.. currentmodule:: scikits.hydroclimpy.enso

The :mod:`scikits.hydroclimpy.enso` module provides two functions to create 
some standard :class:`ENSOIndicator` objects
from online data repository.
For convenience, online data are downloaded the first time they are accessed, 
and then a local copy is used instead.

JMAI indices
------------

Indicators based on the data of the Japanese Meteorological Association (JMA) 
can be accessed with the :func:`load_jma` function.
The data is based on the Sea Surface Temperature (SST) anomalies recorded in
the Niño-3 region of the Pacific Ocean (5N-5S, 150W-90W).
Three different data sets can be used:

- the ``COAPS`` dataset uses a slightly narrower area than Niño-3
  (4N-4S instead of 5N-5S).
  Monthly mean SSTs were computed on boxes of 2-deg latitude x 2-deg longitude
  over the given area.
  Normal mean SSTs were computed using 1961-1990 as a reference period and are
  reported in Table.
  The SST anomalies correspond therefore to the deviation of the monthly average
  SSTs from the normal SSTs.
  A 5-month running mean filter was then applied on the data.
  From January 1949 to present, the anomalies are based on the actual data.
  Before 1949, the anomalies are based on reconstructed monthly mean SST fields,
  computed using orthogonal projection technique, and using the same normal SSTs.
- The ``Standard`` dataset is based on the data recorded over region Niño-3 by
  the JMA.
  Anomalies are computed using a fixed reference period of 1971-2000.
- The ``Sliding`` dataset is similar to the ``Standard`` dataset, except that
  anomalies are computed as the deviations from a sliding 30 year average.

.. tabularcolumns:: |l|c|l|c|
   ===== ====  ===== ====
   Month [oC]  Month [oC]
   ===== ====  ===== ====
   Jan   25.4  Jul   25.2
   Feb   26.2  Aug   24.6
   Mar   26.9  Sep   24.6
   Apr   27.1  Oct   24.6
   May   26.6  Nov   24.6
   Jun   26.1  Dec   24.9
   ===== ====  ===== ====

.. autosummary::
   :toctree: generated/

   load_jma


ONI
---

The Oceanic Niño Index (ONI) used by the Climate Prediction Center (CPC) 
for its previsions can be accessed with the :func:`~ensodata.load_oni` function.
The ONI is based on the 3-m running mean of  ERSST.v3 SST anomalies in the 
Niño-3.4 region (5N-5S, 120W-170W), using the period 1971-2000 as a reference 
for the calculation of monthly normal SSTs.


.. autosummary::
   :toctree: generated/

   load_oni

"""
__author__ = "Pierre GF Gerard-Marchant"


__all__ = ['load_jma', 'load_oni']

import zipfile
import cPickle
import warnings

import numpy as np
import numpy.ma as ma


import scikits.timeseries as ts
from scikits.timeseries import Date, time_series, date_array
from scikits.timeseries.lib.moving_funcs import cmov_average, cmov_mean


from ..core import _config
import ensobase as enso
config = _config.configdict

import logging
ensologger = logging.getLogger('ensodata')

# ...............................................

def _set_ensoindicator_options(ensoindicator, **options):
    """
    Sets some options of an :class:`~scikits.hydroclimpy.enso.ENSOIndicator`

    """
    # Define some defaults..................
    defaults = dict(ensotag='(undefined)',
                    full_year=False,
                    thresholds=(None, None),
                    minimum_size=None,
                    reference_season=None,
                    reference_period=None,
                    )
    optinfo = ensoindicator.optinfo
    # Set the defaults......................
    defaults.update([(k, v) for (k, v) in options.items() if k in defaults])
    optinfo.update(defaults)
    # Fix some issues.......................
    ensoindicator.minimum_size = optinfo['minimum_size']
    ensoindicator.refseason = optinfo['reference_season']
    try:
        ensoindicator.full_year = eval(optinfo['full_year'])
    except NameError:
        optinfo['full_year'] = False
    except TypeError:
        ensoindicator.full_year = optinfo['full_year']
    #
    refstring = (optinfo['reference_period'] or '').split(',')
    if refstring == ['']:
        ensoindicator.refperiod = None
    else:
        if len(refstring) != 2 and refstring != ['']:
            raise ValueError("The starting and ending dates of the reference "
                             "period should be separated bv a comma")
        start_date = refstring[0]
        if len(start_date) == 4:
            start_date = Date('M', year=int(start_date), month=1)
        elif len(start_date) >= 7:
            start_date = Date('M', year=int(start_date[:4]),
                                   month=int(start_date[5:7]))
        else:
            start_date = None
        end_date = refstring[1]
        if len(end_date) == 4:
            end_date = Date('M', year=int(end_date), month=12)
        elif len(end_date) >= 7:
            end_date = Date('M', year=int(end_date[:4]),
                                 month=int(end_date[5:7]))
        else:
            end_date = None
        ensoindicator.refperiod = (start_date, end_date)
    #
    opt_thresholds = optinfo['thresholds']
    if opt_thresholds != (None, None):
        if isinstance(opt_thresholds, tuple):
            ensoindicator.thresholds = opt_thresholds
        else:
            ensoindicator.thresholds = eval(opt_thresholds)
    return



def load_jma(mode='COAPS', **options):
    """
    Load the JMA 5-month running mean monthly SST anomalies over the Nino-3 zone
    (5N-5S, 150W-90W) and returns a :class:`~scikits.hydroclimpy.enso.ENSOIndicator`
    object.
    Three different modes are available:

    ``coaps``
        The data is downloaded from the 
        `COAPS website <ftp://www.coaps.fsu.edu/pub/JMA_SST_Index/jmasst1868-today.filter-5>`_.
    ``standard``
        The data is downloaded from the JMA site and corresponds to the SST
        anomalies computed using a 
        `fixed reference <http://ds.data.jma.go.jp/tcc/tcc/products/elnino/index/sstindex/base_period_7100/Nino_3/5rmean>`_
        period (1971-2000).
    ``sliding``
        The data is also downloaded from the JMA site, but corresponds this time
        to anomalies computed using a 
        `sliding 30-year rule <http://ds.data.jma.go.jp/tcc/tcc/products/elnino/index/sstindex/sliding_30year_period/Nino_3/5rmean>`_.

    By default, the configuration options are read from the configuration file.
    These options can be overwritten by the optional input parameters.


    Parameters
    ----------
    mode : {'coaps','standard','sliding'}, optional
        Mode describing the data to download.
    options : dictionary
        Optional parameters to parse to the ENSOIndicator for the definition of
        ENSO indices.
    thresholds : tuple
        Low and high temperature thresholds for the definition of El Niño or
        La Niña conditions.
        By default, the JMA uses -0.5oC and +0.5oC.
    minimum_size : int
        Minimum number of consecutive months in El Niño/La Niña conditions
        required for the definition of an episode.
        By default, the JMA use 6 consecutive months.
    reference_season : string, tuple 
        Months that must be in an episode for it to be valid.
        By default, the original JMA COAPS index uses ``'OND'`` (October to
        December included), while the newer COAPS index uses ``'NDJ'`` (November
        to January included).
    full_year : boolean
        Whether episodes are defined from Oct. 1st to the following Sep. 30th,
        or last only until the SST anomalies conditions are no longer met.
        The original JMA COAPS index uses ``full_year=True``, while the newer
        index uses ``full_year=False``.

    """
    #
    ensoarchive = dict(config.items('ENSO'))['ensoarchive']
    if ensoarchive[-4:].lower() != '.zip':
        ensoarchive += '.zip'
    # Check the mode and load the proper options
    modes = ('coaps', 'standard', 'sliding')
    mode = mode.lower()
    if mode not in modes:
        raise ValueError("Invalid '%s' mode for the JMA indices: "\
                         "it should be in %s" % (mode, modes))
    if mode == 'coaps':
        cfg = dict(config.items('ENSO.JMA.COAPS'))
    elif mode == 'standard':
        cfg = dict(config.items('ENSO.JMA.Standard'))
    elif mode == 'sliding':
        cfg = dict(config.items('ENSO.JMA.Sliding'))
    else:
        raise ValueError("Invalid '%s' mode for the JMA indices: "\
                         "it should be in %s" % (mode, modes))
    cfg.update(options)
    #
    datadir = cfg['datadir']
    netfile = cfg['netfile']
    archive = cfg['archive']
    scalefactor = float(cfg.get('scalefactor', 1.))
    #
    try:
        zipf = zipfile.ZipFile(ensoarchive, 'r')
        ensoi = cPickle.loads(zipf.read(archive))
        errmsg = "... Loading from existing archived file '%s'" % ensoarchive
        ensologger.info(errmsg)
    except (IOError, ValueError):
        zipf = zipfile.ZipFile(ensoarchive, 'w')
        ensologger.info("... Creating archive")
    except KeyError:
        zipf = zipfile.ZipFile(ensoarchive, 'a')
        ensologger.info("... Appending to archive")
    else:
        if isinstance(ensoi, enso.ENSOIndicator):
            ensoi.set_indices()
            return ensoi
    #
    sourcedir = np.lib._datasource.DataSource(datadir)
    dfile = sourcedir.open(netfile)
    #
    tmp = ts.tsfromtxt(dfile, delimiter=None, dtype=float, datecols=0, freq='A',
                       skiprows=int(cfg.get('skipline', 0)),
                       missing=cfg.get('nodata', '999'),)
    if tmp.varshape != (12,):
        msg = "Unrecognized shape: it should have been (x, 12)."
        raise TypeError(msg)
    ensoi = enso.ENSOIndicator(tmp._series.ravel(),
                               start_date=tmp._dates[0].asfreq('M', 'START'))
    ensoi *= scalefactor
    _set_ensoindicator_options(ensoi, **cfg)
    ensoi.set_indices()
    #
    # Store in the archive
    ensologger.info("... Saving to archive")
    zipf.writestr(archive, cPickle.dumps(ensoi))
    zipf.close()
    return ensoi



def load_oni(mode='standard', **options):
    """
    Loads the ONI 3-m averaged monthly SST anomalies over the Niño-3.4 region
    and returns a :class:`~scikits.hydroclimpy.enso.ENSOIndicator` object.

    Two modes are accepted as arguments:
    
    - in the ``standard`` mode, the SSTs are retrieved from the original CPC
      website_.
      Data are available from Jan. 1950 to present.
    - in the ``backup`` mode, the SSTs are retrieved from the CPC `ftp site <ftpsite>`_.
      Data are available from Jan. 1900 to present.

    .. _website : http://www.cpc.noaa.gov/products/analysis_monitoring/ensostuff/ensoyears.shtml
    .. _ftpsite : ftp://eclipse.ncdc.noaa.gov/pub/ersst/pdo/el_nino_v3.dat.


    Parameters
    ----------
    mode : {'standard','backup'}, optional
        Mode describing the data to download.
    options : dictionary
        Optional parameters to parse to the ENSOIndicator for the definition of
        ENSO indices.
    thresholds : tuple of floats, optional
        Low and high temperature thresholds for the definition of El Niño and
        La Niña conditions.
        By default, the CPC uses -0.5oC and +0.5oC.
    minimum_size : int, optional
        Minimum number of consecutive months in El Niño / La Niña conditions
        required for the definition of an episode.
        By default, the CPC use 5 consecutive months.
    reference_season : string or tuple, optional
        Months that must be in an episode for it to be valid.
        By default, the CPC uses None (no restriction on the months).
    full_year : boolean, optional
        The CPC uses ``full_year=False``.

    References
    ----------
    Xue, Y., T. M. Smith, and R. W. Reynolds, 2003: Interdecadal changes of 30-yr
    SST normals during 1871-2000. *J. Climate*, 16, 1601-1612.

    """
    # Initialization .......................
    ensoarchive = dict(config.items('ENSO'))['ensoarchive']
    if ensoarchive[-4:].lower() != '.zip':
        ensoarchive += '.zip'
    #
    mode = mode.lower()
    cfg = dict(config.items('ENSO.ONI'))
    cfg.update(options)
    try:
        from BeautifulSoup import BeautifulSoup, SoupStrainer
    except ImportError:
        warnings.warn("The module 'BeautifulSoup' is unavailable.\n"\
                      "Reverting to backup mode")
        mode = 'backup'
    #
    datadir = cfg['datadir']
    if mode == 'standard':
        netfile = cfg['netfile']
        archive = cfg['archive']
    else:
        netfile = cfg['netfile_backup']
        archive = cfg['archive_backup']
    # Try to open an existing ENSOIndicator

    ensoarchive = dict(config.items('ENSO'))['ensoarchive']
    if ensoarchive[-4:].lower() != '.zip':
        ensoarchive += '.zip'
    #
    try:
        zipf = zipfile.ZipFile(ensoarchive, 'r')
        ensoi = cPickle.loads(zipf.read(archive))
        ensologger.info("... Loading from existing archived file")
    except IOError:
        zipf = zipfile.ZipFile(ensoarchive, 'w')
        ensologger.info("... Creating archive")
    except KeyError:
        zipf = zipfile.ZipFile(ensoarchive, 'a')
        ensologger.info("... Appending to archive")
    else:
        if isinstance(ensoi, enso.ENSOIndicator):
            return ensoi
    #
    sourcedir = np.lib._datasource.DataSource(datadir)
    dfile = sourcedir.open(netfile)
    #
    #
    if mode == 'standard':
        # Load the file as a tree, but only take the SST table (border=1)
        table = BeautifulSoup(dfile.read(),
                              parseOnlyThese=SoupStrainer("table", border=1))
        # Separate it by rows, but skip the first one (the header)
        years = []
        data = []
        indices = []
        color = {'red':+1, 'white':0, 'blue':-1}
        deft = [(None, 'color:white')]
        for row in table.findAll("tr")[1:]:
            cols = row.findAll('td')
            years.append(int(cols.pop(0).strong.string))
            data.append([float(_.fetchText()[-1].string.replace('&nbsp;', '99.9'))
                         for _ in cols])
            indices.append([color[getattr(_.span, 'attrs', deft)[0][-1].split(':')[-1]]
                                for _ in cols])
        #
        start_date = Date('M', year=years[0], month=1)
        ensoi = enso.ENSOIndicator(ma.masked_values(data, 99.9).ravel(),
                                 start_date=start_date,)
#        oni.set_indices(full_year=False, minsize=5, refseason=None)
        indices = time_series(np.array(indices).ravel(), start_date=start_date)
    else:
        rawdata = np.loadtxt(dfile)
        dates = date_array([Date('M', year=yy, month=mm)
                            for (yy, mm) in rawdata[:, :2]], freq='M')
        ensoi = enso.ENSOIndicator(cmov_mean(rawdata[:, -1], 3).round(2),
                                   dates,)
    #
    _set_ensoindicator_options(ensoi, **cfg)
    ensoi.set_indices()
    #
    # Store in the archive
    zipf.writestr(archive, cPickle.dumps(ensoi))
    zipf.close()
    return ensoi

