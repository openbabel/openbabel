#!/usr/bin/env python
from distutils.core import *
import os,sys,shutil

about = """The Open Babel package provides a Python wrapper
to the Open Babel C++ chemistry library. Open Babel is a chemical
toolbox designed to speak the many languages of chemical data. It's an
open, collaborative project allowing anyone to search, convert,
analyze, or store data from molecular modeling, chemistry, solid-state
materials, biochemistry, or related areas. It provides a broad base of
chemical functionality for custom development.
"""

srcdir = os.path.dirname(__file__)

obExtension = Extension('_openbabel',
                 [os.path.join(srcdir, "openbabel-python.cpp")],
                 include_dirs=[os.path.join(srcdir, "..", "..", "include"),
                               os.path.join("..", "include")],
                 library_dirs=[os.path.join(srcdir, "..", "..", "lib"),
                               os.path.join(srcdir, "..", "..", "lib64"),
                               os.path.join("..", "lib")],
                 libraries=['openbabel']
                 )

if "build" in sys.argv:
    shutil.copyfile(os.path.join(srcdir, "pybel_py%dx.py" % sys.version_info[0]), "pybel.py")
    shutil.copyfile(os.path.join(srcdir, "openbabel.py"), "openbabel.py")

setup(name='openbabel',
      version='1.7',
      author='Noel O\'Boyle',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.org/',
      license='http://www.gnu.org/copyleft/gpl.html',
      py_modules=['openbabel','pybel'],
      ext_modules=[obExtension],
      description = 'openbabel: Python interface to the Open Babel chemistry library',
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
      long_description = about,
      )
