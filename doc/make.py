#!/usr/bin/env python

"""
Python script for building documentation. This script was designed for
building the docs on windows, but may work on other platforms as well.

To build the docs you must have all optional dependencies for the timeseries
installed. See the installation instructions for a list of these.

Note: currently latex builds do not work because of table formats that are not
supported in the latex generation.

Usage
-----

python make.py clean
python make.py html
"""

import fileinput
import fnmatch
import glob
import os
import shutil
import sys

def check_build():
    build_dirs = ['build', 'build/doctrees', 'build/html', 'build/latex',
                  '_static', '_templates']
    for d in build_dirs:
        try:
            os.mkdir(d)
        except OSError:
            pass

def sf():
    'push a copy to the sf site'
#    shutil.copy('../CHANGELOG', 'build/html/_static/CHANGELOG')
    os.system('cd build/html; rsync -avz . pierregmsf,hydroclimpy@web.sourceforge.net:/home/groups/h/hy/hydroclimpy/htdocs/ -essh --cvs-exclude')

#def sfpdf():
#    'push a copy to the sf site'
#    os.system('cd build/latex; scp Matplotlib.pdf pierregmsf,hydroclimpy@web.sourceforge.net:/home/groups/h/hy/hydroclimpy/htdocs/')
# 
# def figs():
#     os.system('cd users/figures/ && python make.py')

def html():
    check_build()
    shutil.copy('../scikits/hydroclimpy/hydroclimpyrc', '_static/hydroclimpyrc')
    if small_docs:
        options = "-D plot_formats=\"[('png', 80)]\""
    else:
        options = ''
    if os.system('sphinx-build %s -P -b html -d build/doctrees source build/html' % options):
        raise SystemExit("Building HTML failed.")

#    figures_dest_path = 'build/html/plots'
#    if os.path.exists(figures_dest_path):
#        shutil.rmtree(figures_dest_path)
#    shutil.copytree('plots', figures_dest_path)

def latex():
    check_build()
    #figs()
    if sys.platform != 'win32':
        # LaTeX format.
        if os.system('sphinx-build -b latex -d build/doctrees source build/latex'):
            raise SystemExit("Building LaTeX failed.")

        # Produce pdf.
        os.chdir('build/latex')

        # Copying the makefile produced by sphinx...
        if (os.system('pdflatex HydroClimpy') or
            os.system('pdflatex HydroClimpy') or
            os.system('makeindex -s python.ist HydroClimpy') or
            os.system('makeindex -s python.ist HydroClimpy') or
            os.system('pdflatex HydroClimpy')):
            raise SystemExit("Rendering LaTeX failed.")

        os.chdir('../..')
    else:
        print 'latex build has not been tested on windows'

def clean():
    try:
        shutil.rmtree("build")
    except OSError:
        pass
    generated = []
    for (r, d, f) in os.walk("source"):
        for dirname in fnmatch.filter(d, 'generated'):
             dirname = os.path.join(r, dirname)
             if os.path.exists(dirname):
                 shutil.rmtree(dirname)
    # for pattern in ['doc/source/generated/*.',
    #                 'doc/plots/*.pdf',
    #                 'doc/_static/hydroclimpyrc', ]:
    #     for filename in glob.glob(pattern):
    #         if os.path.exists(filename):
    #             os.remove(filename)

def all():
    #figs()
    html()
    latex()


funcd = {
#    'figs'     : figs,
    'html'     : html,
    'latex'    : latex,
    'clean'    : clean,
    'sf'       : sf,
#    'sfpdf'    : sfpdf,
    'all'      : all,
    }


small_docs = False

# Change directory to the one containing this file
current_dir = os.getcwd()
os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

if len(sys.argv) > 1:
    if '--small' in sys.argv[1:]:
        small_docs = True
        sys.argv.remove('--small')
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are %s' % (
                    arg, funcd.keys()))
        func()
else:
    small_docs = False
    all()
os.chdir(current_dir)
