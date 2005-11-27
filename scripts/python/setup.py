###
## setup.py
from distutils.core import *

# lang = 'c++'

obExtension = Extension('_openbabel',
                        ['openbabel_python.cpp'],
                        include_dirs=['../../src'],
                        library_dirs=['../../src'],
                        libraries=['openbabel']
                        )

setup(name='openbabel',
      version='1.0.0',
      description='Chemistry interface to Open Babel',
      author='Geoff Hutchison',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/',
      py_modules=['openbabel'],
      ext_modules=[obExtension])
###
