###
## setup.py
from distutils.core import *

lang = 'c++'

setup(name='openbabel',
      version='2.0.0b1',
      description='Chemistry interface to Open Babel',
      author='Geoff Hutchison',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/',
      py_modules=['openbabel'],
      ext_modules=[Extension('openbabel',
                             ['openbabel_python.cpp'],
                             language=lang,
                             include_dirs=['../../src'],
                             library_dirs=['../../src'],
                             libraries=['openbabel']
                            )]
      )
###
