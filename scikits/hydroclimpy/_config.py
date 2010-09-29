"""
Configuration information for climpy.iotools

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date: 2008-05-19 16:56:16 -0400 (Mon, 19 May 2008) $
:version: $Id: _config.py 374 2008-05-19 20:56:16Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 374 $"
__date__     = '$Date: 2008-05-19 16:56:16 -0400 (Mon, 19 May 2008) $'

import os
import shutil
import sys
import tempfile
import warnings
from operator import itemgetter
from ConfigParser import ConfigParser, NoOptionError


###-----------------------------------------------------------------------------
# Initialize the logger

import logging

# Basic info
formatstr = "%(levelname)-8s : %(name)s/%(module)s.%(funcName)s: "\
            "%(message)s"
logging.basicConfig(level=logging.DEBUG,
                    format=formatstr,
                    datefmt='%m-%d %H:%M',)

## Define a handler to the sys.stderr
#console = logging.StreamHandler()
#console.setLevel(logging.DEBUG)
## Set a format which is simpler for console use
#formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
#console.setFormatter(formatter)
## Add the handler to the root logger
#logging.getLogger('').addHandler(console)






#--- --- Misc functions from matplotlib ---

def _is_writable_dir(p):
    """

    Returns `True` if `p` is a string pointing to a potentially writable
    directory, `False` otherwise.

    """
    try: 
        p + ''  # test is string like
    except TypeError: 
        return False
    try:
        t = tempfile.TemporaryFile(dir=p)
        t.write('1')
        t.close()
    except OSError: 
        return False
    else: 
        return True

def _get_home():
    """
    Find the home directory of the user, if possible, or raises a RuntimeError

    :see:  http://mail.python.org/pipermail/python-list/2005-February/263921.html

    """
    path = ''
    try:
        path = os.path.expanduser("~")
    except:
        pass
    if not os.path.isdir(path):
        for evar in ('HOME', 'USERPROFILE', 'TMP'):
            try:
                path = os.environ[evar]
                if os.path.isdir(path):
                    break
            except: 
                pass
    if path:
        return path
    else:
        raise RuntimeError('Please define the environment variable $HOME')

def _get_configdir():
    """
    Returns the string representing the configuration directory.

    default is HOME/.hydroclimpy.  you can override this with the
    HYDROCLIMPYONFIGDIR environment variable

    """
    configdir = os.environ.get('HYDROCLIMPYONFIGDIR')
    if configdir is not None:
        if not _is_writable_dir(configdir):
            errmsg = 'Could not write to HYDROCLIMPYONFIGDIR="%s"'
            raise RuntimeError(errmsg % configdir)
        return configdir

    h = _get_home()
    p = os.path.join(_get_home(), '.hydroclimpy')

    if os.path.exists(p):
        if not _is_writable_dir(p):
            errmsg = "'%s' is not a writable directory.\n"\
                     "Please set %s/.hydroclimpy to a writable directory.\n"\
                     "Alternatively, please define the environment variable "\
                     "HYDROCLIMPYCONFIGDIR to any writable directory where "\
                     "the data must be stored."
            raise RuntimeError(errmsg % (h, h))
    else:
        if not _is_writable_dir(h):
            errmsg = "Failed to create the %s/.hydroclimpy directory.\n"\
                     "Please consider setting HYDROCIMPYCONFIGDIR to "\
                     "a writable directory for configuration data."
            raise RuntimeError(errmsg % h)
        os.mkdir(p)

    return p



def _get_data_path():
    'get the path to matplotlib data'

    path = os.environ.get('HYDROCLIMPYDATA', None)
    if path is not None:
        if not os.path.isdir(path):
            errmsg = "The path in environment HYDROCLIMPYDATA is not a directory."
            raise RuntimeError(errmsg)
        return path

    path = os.path.dirname(__file__)
    if os.path.isdir(path):
        return path




def _get_hydroclimpy_filename():
    """
    Return the path to the configuration file

    Search order:

     * current working dir
     * environ var HYDROCLIMPYRC
     * HOME/.mhydroclimpy/hydroclimpyrc
     * $PYTHONPATH/hydroclimpy/matplotlibrc


    """
    # Check that we don't have an old config file lurking in the current dir
    oldname = os.path.join(os.getcwd(), '.hydroclimpyrc')
    if os.path.exists(oldname):
        sys.stderr.write("""
WARNING:
Old rc filename ".hydroclimpyrc" found in working directory and renamed to
new default file name "hydroclimpyrc" (no leading ".").""")
        sys.stderr.flush()
        shutil.move('.hydroclimpyrc', 'hydroclimpyrc')
    # Make sure we don't have another config file lurking around
    home = _get_home()
    oldname = os.path.join(home, '.hydroclimpyrc')
    if os.path.exists(oldname):
        configdir = _get_configdir()
        newname = os.path.join(configdir, 'hydroclimpyrc')
        sys.stderr.write("""
WARNING:
Old rc filename "%s" found and renamed to new default rc file name "%s".
        """ % (oldname, newname))
        shutil.move(oldname, newname)

    # Is the config file in the current directory ?
    fname = os.path.join(os.getcwd(), 'hydroclimpyrc')
    if os.path.exists(fname):
        return fname
    # Is the HYDROCLIMPYRC environment variable defined ?
    path = os.environ.get('HYDROCLIMPYRC', None)
    if path is not None:
        if os.path.exists(path):
            fname = os.path.join(path, 'hydroclimpyrc')
            if os.path.exists(fname):
                return fname
    # Is the config file in the config dir ?
    fname = os.path.join(_get_configdir(), 'hydroclimpyrc')
    if os.path.exists(fname):
        return fname

    path = _get_data_path() # guaranteed to exist or raise
    fname = os.path.join(path, 'hydroclimpyrc')
    if not os.path.exists(fname):
        warnings.warn('Could not find hydroclimpyrc; using defaults')
    return fname



#--- --- Configuration function ---

def update_config(confg):
    """
    Updates the information of higher-level sections to take lower section 
    into account.
    """
    # Check whether we need to stay local ..
    try:
        force_local = confg.getboolean('DEFAULT', 'force_local')
    except NoOptionError:
        force_local = False
    # Update the root directory ............
    if force_local:
        rootdir = os.path.join(os.getcwd(),'hydroclimpywork')
    else:
        rootdir = os.path.expanduser(confg.get('DEFAULT', 'rootdir'))
    # Get the sections .....................
    # Rank sections on the nb of . and create a dictionary section <-> level
    sect_levels = dict([(sect, len(sect.split('.')))
                        for sect in confg.sections()])
    # Get a list of sections, ordered by increasing level
    sections_ordered = [s[0] for s in sorted(sect_levels.items(),
                                             key=itemgetter(1))]
    # Update the items of the higher level sections
    for sect in sections_ordered:
        if sect_levels[sect] == 1:
            confg.set(sect, 'rootdir', rootdir)
        if sect_levels[sect] > 1:
            ancestor = '.'.join(sect.split('.')[:-1])
            for (item, value) in confg.items(ancestor):
                confg.set(sect, item, value)
    return confg

#..............................................................................

def _check_localdirectories(confg):
    """
    Checks if the `datadir` directories stored in the configuration file `confg`
    exist, and creates them if needed.

    """
    # Get the list of directory .................
    dirlist = []
    for section in confg.sections():
        rootdir = os.path.expanduser(confg.get(section, 'rootdir'))
        dirlist.append(rootdir)
        local_list = [confg.get(section, item, False,)
                        for item in ('rawdatadir', 'datadir') 
                            if confg.has_option(section, item)]
        for item in local_list:
            if item is not None:
                newpath = os.path.join(rootdir, item)
                dirlist.append(newpath)
    # Check the existence and create if needed ..
    for currentdir in dirlist:
        if not os.path.isdir(currentdir):
            msg = "Invalid directory '%s'.\nAttempting to create it..." 
            warnings.warn(msg % currentdir)
            try:
                os.makedirs(currentdir)
            except IOError:
                msg = "OSError: a file named '%s' already exists!"
                warnings.warn(msg % currentdir)
    return dirlist


#---- --- Initialization ---

def _get_hydr():
    """
    Gets the configuration file.
    
    """
    try:
        configfile = os.environ['HYDROCLIMPYCONFIG']
    except KeyError:
        configfile = os.path.join(os.getcwd(), 'hydroclimpyrc')
        try:
            open(configfile)
        except IOError:
            p = _get_configdir()
            configfile = os.path.join(p, 'hydroclimpyrc')
    return configfile

def configure(configfile=None):
    "Configure the environment."
    if configfile is None:
        configfile = _get_hydroclimpy_filename()
    config = ConfigParser()
    flag = config.read(configfile)
    if not flag:
        errmsg = "Cannot find the configuration file 'hydroclimpyrc'"
        errmsg = "Cannot find the configuration file '%s'" % configfile
        raise RuntimeError(errmsg)
    config = update_config(config)
    _check_localdirectories(config)
    return config

configdict = configure()

# Change the default verbosity level
deflevel = getattr(logging, configdict.defaults().get('verbose', 'DEBUG'))
logging.basicConfig(level=deflevel)

# Add an extra logger to the logfile (if it is defined)
# The 'logfile'entry is '' if empty : make it None by default
logfilename = configdict.defaults().get('logfile', None) or None
if logfilename is not None:
    # Save in the current directory by default (only a file name, not a path...)
    if os.path.dirname(logfilename) is None:
        logfilename = os.path.join(os.getcwd(), logfilename)
    # Define the new logfile
    logfile = logging.FileHandler(logfilename)
    #logfile = logging.TimedRotatingFileHandler(logfilename, 'D')
    formatstr = "%(asctime)s: " + formatstr
    formatter = logging.Formatter(formatstr, "%Y-%m-%d %H:%M")
    logfile.setFormatter(formatter)
    logfile.setLevel(deflevel)
    logging.getLogger('').addHandler(logfile)


