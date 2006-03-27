#!/usr/bin/env python
from distutils.core import *
import os,sys

def find_likely_directory():
    """Find (guess!) where Open Babel is installed.

    Order of precedence is:
      $OPENBABEL_INSTALL > /usr/local > /usr > ../../src
    """
    name = os.environ.get("OPENBABEL_INSTALL")
    if name: # OPENBABEL_INSTALL is set
        sys.stderr.write("INFO: Using the value of $OPENBABEL_INSTALL (%s)\n" % name)
        if not os.path.isdir(name):
            sys.stderr.write("ERROR: $OPENBABEL_INSTALL (%s) is not a directory\n" % name)
        else:
            return ([name+"/include/openbabel-2.0",name+"/include/openbabel-2.0/openbabel"],
                    [name+"/lib/openbabel"])

    else: # OPENBABEL_INSTALL is not set
        sys.stderr.write("WARNING: Environment variable OPENBABEL_INSTALL is not set\n")
        for dirname in ["/usr/local","/usr"]:
            # Look for each of these directories in turn for the directory include/openbabel-2.0
            # (This is version specific, so I may do as Andrew Dalke did for PyDaylight and use
            #  a regular expression to find the latest version of openbabel)
            if os.path.isdir(dirname+"/include/openbabel-2.0"):
                sys.stderr.write("INFO: Setting OPENBABEL_INSTALL to %s\n" % dirname)
                return ([dirname+"/include/openbabel-2.0",dirname+"/include/openbabel-2.0/openbabel"],
                        [dirname+"/lib/openbabel"])
            else:
                sys.stderr.write("WARNING: Open Babel does not appear to be globally installed\n" +
                                 "INFO: Looking for library and include files in ../../src\n")
                if os.path.isfile("../../src/atom.o"):
                    return ["../../src"],["../../src"]
                
    sys.stderr.write("ERROR: Cannot find Open Babel library directory\n")
    return (None,None)




OBinclude,OBlibrary = find_likely_directory()

obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=OBinclude,
                        library_dirs=OBlibrary,
                        libraries=['openbabel']
                        )

setup(name='openbabel',
      version='1.1.0',
      description='Python interface to Open Babel',
      author='Geoff Hutchison',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/',
      py_modules=['openbabel'],
      ext_modules=[obExtension])
