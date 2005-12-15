#!/usr/bin/env python
from distutils.core import *
import os,sys

def find_likely_directory():
    """Find where Open Babel is installed.

    Order of precedence is:
      $OPENBABEL_LIBDIR > /usr/local/openbabel
    """
    name = os.environ.get("OPENBABEL_LIBDIR")
    if name: # OPENBABEL_LIBDIR is set
        if not os.path.isdir(name):
            sys.stderr.write("WARNING: $OPENBABEL_LIBDIR (%s) is not a directory\n" % name)
            return name
    else: # OPENBABEL_LIBDIR is not set
        sys.stderr.write("WARNING: Environment variable OPENBABEL_LIBDIR is not set")
        for dirname in ["/usr/local/openbabel"]: # Look for each of these directories in turn
            if os.path.isdir(dirname):
               sys.stderr.write("INFO: Setting OPENBABEL_LIBDIR to %s\n" % dirname)
               name = dirname
               return name

    sys.stderr.write("ERROR: Cannot find Open Babel library directory\n")
    return None
        
# Need to edit the next statement to use find_likely_directory
obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=['../../src'],
                        library_dirs=['../../src'],
                        libraries=['openbabel']
                        )

setup(name='openbabel',
      version='1.0.0',
      description='Python interface to Open Babel',
      author='Geoff Hutchison',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/',
      py_modules=['openbabel'],
      ext_modules=[obExtension])
