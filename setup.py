#!/usr/bin/env python
import os
import sys
from cmake_build_extension import BuildExtension, CMakeExtension
from setuptools import setup

# Path to scripts/python
base_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "scripts", "python")


setup(
    name="openbabel",
    version="3.1.1",
    author='Noel O\'Boyle',
    author_email='openbabel-discuss@lists.sourceforge.net',
    license='GPL-2.0',
    url='http://openbabel.org/',
    description='Python interface to the Open Babel chemistry library',
    long_description=open(os.path.join(base_dir, 'README.rst')).read(),
    ext_modules=[
        CMakeExtension(name="OpenBabel",
                       install_prefix="openbabel",
                       cmake_configure_options=[
                            "-DPYTHON_EXECUTABLE=" + sys.executable,
                            "-DCMAKE_BUILD_TYPE=Release",
                            "-DWITH_INCHI=ON",
                            "-DPYTHON_BINDINGS=ON",
                            "-DRUN_SWIG=ON",
                            "-DBUILD_BY_PIP=ON",
                       ]),
    ],
    cmdclass=dict(build_ext=BuildExtension),
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
        'Topic :: Software Development :: Libraries'
    ]
)
