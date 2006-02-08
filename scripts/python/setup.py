###
## setup.py
from distutils.core import *

# lang = 'c++'

try:
      path = open("/usr/local/lib/libopenbabel.a")
      print " OK, found library "
      obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=['/usr/local/include/openbabel-2.0', '/usr/local/include/openbabel-2.0/openbabel'],
                        library_dirs=['/usr/local/lib'],
                        libraries=['openbabel']
                        )
except:
      print " Oops, can't open library"
      obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=['../../src'],
                        library_dirs=['../../src'],
                        libraries=['openbabel']
                        )

setup(name='openbabel',
      version='1.0.1',
      description='Chemistry interface to Open Babel',
      author='Geoff Hutchison',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/',
      py_modules=['openbabel'],
      ext_modules=[obExtension])
###
