"""
:mod:`scikits.hydroclimpy`
--------------------------

A collection of tools to manipulate environmental and climatologic time series.

setuptools must be installed first. If you do not have setuptools installed
please download and install it from http://pypi.python.org/pypi/setuptools
"""

version = '0.67.1'

classifiers = ['Development Status :: 4 - Beta',
               'Intended Audience :: Science/Research',
               'Intended Audience :: Developers',
               'License :: OSI Approved :: BSD License',
               'Operating System :: Microsoft :: Windows',
               'Operating System :: POSIX',
               'Operating System :: Unix',
               'Operating System :: MacOS',
               'Programming Language :: C',
               'Programming Language :: Python',
               'Topic :: Education',
               'Topic :: Scientific/Engineering',
               'Topic :: Software Development',
              ]

distname = 'scikits.hydroclimpy'

long_description = """
The scikits.hydroclimpy module is a collection of tools for manipulating and 
plotting environmental time series of various frequencies. This package is 
an extension for scikits.timeseries, focusing on tools for the analysis of
hydroclimatologic datasets.
"""

import os
import sys
import setuptools
from numpy.distutils.core import setup, Extension


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path,
                           namespace_packages=['scikits'])

    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True,
    )

    config.add_subpackage('scikits')
    config.add_subpackage('doc')
    config.add_subpackage('examples')
    config.add_subpackage(distname)
    config.add_data_files('scikits/__init__.py')

    return config

package_data = {'': 'hydroclimpyrc'}

def setup_package():

    setup(
          install_requires=['numpy > 1.2.5',
                            'scipy >= 0.7',
                            'scikits.timeseries >= 0.91'],
          namespace_packages=['scikits'],
          packages=setuptools.find_packages(),
          test_suite = 'nose.collector',
          name = distname,
          version = version,
          description = "Environmental time series manipulation",
          long_description = long_description,
          license = "BSD",
          author = "Pierre GF GERARD-MARCHANT",
          author_email = "pierregmcode_AT_gmail_DOT_com",
          maintainer = "Pierre GF GERARD-MARCHANT",
          maintainer_email = "pierregmcode_AT_gmail_DOT_com",
          url = "http://hydroclimpy.sourceforge.net",
          classifiers = classifiers,
          platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
          configuration = configuration,
    )

    return

if __name__ == '__main__':
    setup_package()
