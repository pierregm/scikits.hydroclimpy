
import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    nxheader = join(get_numpy_include_dirs()[0], 'numpy',)
    confgr = Configuration('hydroclimpy', parent_package, top_path)
    #
    confgr.add_data_files('hydroclimpyrc','hydroclimpyrc.template')
    confgr.add_subpackage('stats')
    return confgr

if __name__ == "__main__":
    from numpy.distutils.core import setup
    #setup.update(nmasetup)
    config = configuration(top_path='').todict()
    setup(**config
          )
