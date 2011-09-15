#!/usr/bin/env python
from distutils.core import *
import os, sys, shutil, glob

about = """Open Babel is a chemical toolbox designed to speak the
many languages of chemical data. It's an open, collaborative project
allowing anyone to search, convert, analyze, or store data from
molecular modeling, chemistry, solid-state materials, biochemistry,
or related areas.
"""

# The following line is necessary because only one 'root package' location
# is possible: either "." or "../../scripts/python", but not both.
# The root package location is set by "package_dir" and defaults to "."
shutil.copyfile("../../../scripts/python/pybel_py%dx.py" % sys.version_info[0],
                "pybel.py")
shutil.copyfile("../../build/bin/Release/_openbabel.pyd", "_openbabel.pyd")
shutil.copyfile("../../../scripts/python/openbabel.py", "openbabel.py")

setup(name='openbabel-python',
      version='1.7',
      author='The Open Babel development team',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/wiki/Python',
      license='http://www.gnu.org/copyleft/gpl.html',
      py_modules=['openbabel', 'pybel'],
      data_files=[('Lib/site-packages', ['_openbabel.pyd'])],
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
os.remove("openbabel.py")
os.remove("_openbabel.pyd")
