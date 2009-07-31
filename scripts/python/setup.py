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

def find_likely_directory():
    """Find (guess!) where Open Babel is installed.

    Order of precedence is:
      $OPENBABEL_INSTALL > ../../src > /usr/local > /usr
    """
    includedir, libdir = None, None
    name = os.environ.get("OPENBABEL_INSTALL")
    if name: # OPENBABEL_INSTALL is set
        sys.stderr.write("INFO: Using the value of $OPENBABEL_INSTALL (%s)\n" % name)
        if not os.path.isdir(name):
            sys.stderr.write("ERROR: $OPENBABEL_INSTALL (%s) is not a directory\n" % name)
        else:
            includedir = [name+"/include/openbabel-2.0"]
            libdir = [name+"/lib"]

    else: # OPENBABEL_INSTALL is not set
        sys.stderr.write("WARNING: Environment variable OPENBABEL_INSTALL is not set\n")
        sys.stderr.write("INFO: Looking for library and include files in ../../src and ../../include\n")
        if os.path.isfile("../../include/openbabel/atom.h"):
            includedir = ["../../include"]
        if os.path.isfile("../../src/.libs/libopenbabel.so") or os.path.isfile("../../src/.libs/libopenbabel.dylib"):
            libdir = ["../../src/.libs"]
        elif os.path.isfile("../../src/libopenbabel.so"):
            libdir = ["../../src"]

        if not (libdir and includedir):
            for dirname in ["/usr/local","/usr"]:
# Look for each of these directories in turn
# for the directory include/openbabel-2.0
# (This is version specific, so I may do as
# Andrew Dalke did for PyDaylight and use
# a regular expression to find the latest version of openbabel)
                if os.path.isdir(dirname+"/include/openbabel-2.0"):
                    sys.stderr.write("INFO: Setting OPENBABEL_INSTALL to %s\n" % dirname)
                    includedir = [dirname+"/include/openbabel-2.0"]
                    libdir = [dirname+"/lib"]
                    break
                
    if not (libdir and includedir):
        sys.stderr.write("ERROR: Cannot find Open Babel include or library directory\n")
    return (includedir, libdir)




OBinclude,OBlibrary = find_likely_directory()

obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=OBinclude,
                        library_dirs=OBlibrary,
                        libraries=['openbabel']
                        )

shutil.copyfile("pybel_py%d.x.py" % sys.version_info[0], "pybel.py")

setup(name='openbabel',
      version='1.5',
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

os.remove("pybel.py")
