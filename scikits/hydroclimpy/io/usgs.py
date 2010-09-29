"""
Utilities for importing USGS streamflows data.

.. autofunction:: load_usgs_flows

"""

import itertools
import os
import zipfile
import cPickle
import urllib

import numpy as np
import numpy.ma as ma

import scikits.timeseries as ts

import logging
usgslogger = logging.getLogger('usgs')

from .. import _config
config = _config.configdict

usgs_config = dict(config.items('USGS'))
data_dict = {'00060':'F'}

__all__ = ['load_usgs_flows']


def _check_archive(cfgdict):
    "Returns the name of the ZIP archive."
    archive = cfgdict['archive']
    if archive[-4:].lower() != '.zip':
        archive += '.zip'
    return archive


def load_usgs_flows(site_no, cfgdict=usgs_config):
    """
    Downloads hydrological data from the USGS site.

    Parameters
    ----------
    site_no : str
        2- to 15-digit identification number of the site whose information must
        be downloaded.
    netfile : {string}, optional
        URL of the USGS flow gage site.
        By default, the value is read from the configuration file.
    begin_date : {string}, optional
        First date to be retrieved.
        By default, the value is read from the configuration file.
    end_date : {None, string}, optional
        Last date to be retrieved.
        If None, uses the current date.
    datadir : {None, string}, optional
        Directory where to store the output.
        If None, uses the current directory.
    commentchar : {string}, optional
        Character corresponding to a comment.
    headerlines : {int}, optional
        Number of uncommented lines to skip from the beginning of the file.
    archive : {None, string}, optional
        Name of the ZIP archive where to read pre-stored data or when to store 
        data newly retrieved.
        The default is read from the configuration file.

    Notes
    -----
    The list of data codes valid for the ``data_type`` parameters are available
    on the `USGS site <http://waterdata.usgs.gov/nwis/pmcodes/pmcodes?pm_group=All+--+include+all+parameter+groups&pm_search=&format=html_table&show=parameter_group_nm&show=parameter_nm>`_

    """
    # Force the site_no to a 0-padded, 8-character string
    site_no = "%08i" % int(site_no)
    #data_type = cfgdict.get('data_type', '00060')
    data_type = '00060'
    data_letter = data_dict[data_type]
    usgsarchv = _check_archive(cfgdict)
    archfname = "%s%sD" % (data_letter, site_no)
    #
    try:
        zipf = zipfile.ZipFile(usgsarchv, 'r')
        series = cPickle.loads(zipf.read(archfname))
    except IOError:
        zipf = zipfile.ZipFile(usgsarchv, 'w')
    except KeyError:
        zipf = zipfile.ZipFile(usgsarchv, 'a')
    else:
        if series.size > 0:
            return series
        else:
            zipf = zipfile.ZipFile(usgsarchv, 'a')
    #
    options = dict(usgs_site=cfgdict['netfile'],
                   begin_date=cfgdict['begin_date'],
                   data_type=data_type,
                   site_no= site_no)
    end_date = cfgdict.get('end_date', None)
    if end_date is None:
        end_date = ts.now('D').strftime("%Y-%m-%d")
    options['end_date'] = end_date
    commd = "%(usgs_site)s/dv?format=rdb&cb_%(data_type)s=on&" + \
            "begin_date=%(begin_date)s&end_date=%(end_date)s&" + \
            "site_no=%(site_no)s"
    datadir = cfgdict.get('datadir', os.getcwd())
    # We can't use np.lib._datasource.DataSource, the downloaded file has always
    # the same name and would be cached. 
    sourcedir = np.lib._datasource.DataSource(datadir)
    dfile = urllib.urlopen(commd % options)
    #
    headers = int(cfgdict['headerlines'])
    #
    (datelist, datalist) = ([], [])
    for line in dfile.readlines():
        line = line.strip().split("\t")
        if (not line) or (line[0][0] == '#'):
            continue
        try:
            datalist.append(line[3])
        except IndexError:
            continue
        else:
            datelist.append(line[2])
    series = ts.time_series([float(_) for _ in datalist[headers:]],
                            dates=[ts.Date('D', _) for _ in datelist[headers:]],
                            freq='D')
    #
    zipf.writestr(archfname, cPickle.dumps(series))
    zipf.close()
    return series

