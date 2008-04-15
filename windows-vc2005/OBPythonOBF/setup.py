#!/usr/bin/env python
from distutils.core import *
import os, shutil, glob

about = """Open Babel is a chemical toolbox designed to speak the
many languages of chemical data. It's an open, collaborative project
allowing anyone to search, convert, analyze, or store data from
molecular modeling, chemistry, solid-state materials, biochemistry,
or related areas.
"""

# The following line is necessary because only one 'root package' location
# is possible: either "." or "../../scripts/python", but not both.
# The root package location is set by "package_dir" and defaults to "."
shutil.copy("../../scripts/python/pybel.py", ".")

setup(name='openbabel-python',
      version='1.2',
      author='The Open Babel development team',
      author_email='openbabel-scripting@lists.sourceforge.net',
      url='http://openbabel.sourceforge.net/wiki/Python',
      license='http://www.gnu.org/copyleft/gpl.html',
      scripts=["openbabel_postinstall.py"],
      py_modules=['openbabel', 'pybel'],
      # package_dir = {'': '../../scripts/python'},
      data_files=[('Lib/site-packages',
                   ['_openbabel.pyd', '../libinchi.dll',
                    '../libxml2.dll', 'OpenBabelDLL.dll',
                    '../zlib1.dll', '../OBFPRT.obf',
                    '../OBDESC.obf', '../obcommon.obf',
                    '../OBMore.obf', '../OBXML.obf',
                    '../OBUtil.obf', '../OBInchi.obf',
                    '../OBMCDL.obf']),
                  ('Lib/site-packages/openbabel_data',
                   glob.glob("../../data/*.txt"))
                 ],
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
