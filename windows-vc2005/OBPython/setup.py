#!/usr/bin/env python
from distutils.core import *
import os, shutil

about = """The Open Babel package provides a Python wrapper
to the Open Babel C++ chemistry library. Open Babel is a project
designed to pick up where Babel left off, as a cross-platform program
and library designed to interconvert between many file formats used in
molecular modeling, computational chemistry, and many related
areas. It provides a broad base of chemical functionality for custom
development.
"""

# The following line is necessary because only one 'root package' location
# is possible: either "." or "../../scripts/python", but not both.
# The root package location is set by "package_dir" and defaults to "."
shutil.copy("../../scripts/python/pybel.py", ".")

setup(name='openbabel-python',
      version='1.0',
      author='The Open Babel development team',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/wiki/Python',
      license='http://www.gnu.org/copyleft/gpl.html',
      py_modules=['openbabel', 'pybel'],
      # package_dir = {'': '../../scripts/python'},
      data_files=[('Lib/site-packages',
                   ['_openbabel.pyd', '../libinchi.dll',
                    '../libxml2.dll', 'OpenBabelDLL.dll'])],
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