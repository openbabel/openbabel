#!/usr/bin/env python
from setuptools import setup, Extension
import os
import subprocess
import sys


__author__ = 'Noel O\'Boyle'
__email__ = 'openbabel-discuss@lists.sourceforge.net'
__version__ = '1.8'
__license__ = 'GPL'


if os.path.exists('README.rst'):
    long_description = open('README.rst').read()
else:
    long_description = '''
        The Open Babel package provides a Python wrapper to the Open Babel C++ chemistry library. Open Babel is a
        chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project
        allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state
        materials, biochemistry, or related areas. It provides a broad base of chemical functionality for custom
        development.
    '''

try:
    subprocess.check_output('pkg-config --exists openbabel-2.0'.split())
except:
    print ('pkg-config could not find an openbabel install. Please install openbabel and rerun setup.py.')
    sys.exit(1)


babel_version = subprocess.Popen('pkg-config --modversion openbabel-2.0'.split(),
    stdout=subprocess.PIPE).stdout.readline().decode('utf-8')

include_dirs = subprocess.Popen('pkg-config --variable=pkgincludedir openbabel-2.0'.split(),
    stdout=subprocess.PIPE).stdout.readline().decode('utf-8').split()

lib_dirs = subprocess.Popen('pkg-config --variable=libdir openbabel-2.0'.split(),
    stdout=subprocess.PIPE).stdout.readline().decode('utf-8').split()

#Code should be added here to ensure these directories exist, and that this version of babel is compatible
#with this python package. 

swig_opts = ['-c++', '-small', '-O', '-templatereduce', '-naturalvar']
swig_opts += ['-I%s' % i for i in include_dirs]

obextension = Extension('_openbabel',
                        ['openbabel-python.i'],
                        include_dirs=include_dirs,
                        library_dirs=lib_dirs,
                        swig_opts=swig_opts,
                        libraries=['openbabel'])


setup(name='openbabel',
      version=__version__,
      author=__author__,
      author_email=__email__,
      license=__license__,
      url='http://openbabel.org/',
      description='Python interface to the Open Babel chemistry library',
      long_description=long_description,
      zip_safe=True,
      py_modules=['openbabel', 'pybel'],
      ext_modules=[obextension],
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Environment :: Other Environment',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: OS Independent',
          'Operating System :: POSIX',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Unix',
          'Programming Language :: C++',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Software Development :: Libraries',
      ],
)
