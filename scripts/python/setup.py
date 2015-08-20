#!/usr/bin/env python
import os
import subprocess
import sys
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.errors import DistutilsExecError
from distutils.version import StrictVersion
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools import setup, Extension


__author__ = 'Noel O\'Boyle'
__email__ = 'openbabel-discuss@lists.sourceforge.net'
__version__ = '1.8.3'
__license__ = 'GPL'


if os.path.exists('README.rst'):
    long_description = open('README.rst').read()
else:
    long_description = '''
        The Open Babel package provides a Python wrapper to the Open Babel C++
        chemistry library. Open Babel is a chemical toolbox designed to speak
        the many languages of chemical data. It's an open, collaborative
        project allowing anyone to search, convert, analyze, or store data from
        molecular modeling, chemistry, solid-state materials, biochemistry, or
        related areas. It provides a broad base of chemical functionality for
        custom development.
    '''


class PkgConfigError(Exception):
    pass


def pkgconfig(package, option):
    """Wrapper around pkg-config command line tool."""
    try:
        p = subprocess.Popen(['pkg-config', option, package],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True)
        stdout, stderr = p.communicate()
        if stderr:
            raise PkgConfigError('package %s could not be found by pkg-config' % package)
        return stdout.strip()
    except OSError:
        raise PkgConfigError('pkg-config could not be found')


def locate_ob():
    """Try use pkgconfig to locate Open Babel, otherwise guess default location."""
    try:
        version = pkgconfig('openbabel-2.0', '--modversion')
        if not StrictVersion(version) >= StrictVersion('2.3.0'):
            print('Warning: Open Babel 2.3.0 or later is required. Your version (%s) may not be compatible.' % version)
        include_dirs = pkgconfig('openbabel-2.0', '--variable=pkgincludedir')
        library_dirs = pkgconfig('openbabel-2.0', '--variable=libdir')
        print('Open Babel location automatically determined by pkg-config:')
    except PkgConfigError as e:
        print('Warning: %s.\nGuessing Open Babel location:' % e)
        include_dirs = '/usr/local/include/openbabel-2.0'
        library_dirs = '/usr/local/lib'
    return include_dirs, library_dirs


class CustomBuild(build):
    """Ensure build_ext runs first in build command."""
    def run(self):
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    """Ensure build_ext runs first in install command."""
    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()


class CustomSdist(sdist):
    """Add swig interface files into distribution from parent directory."""
    def make_release_tree(self, base_dir, files):
        sdist.make_release_tree(self, base_dir, files)
        link = 'hard' if hasattr(os, 'link') else None
        self.copy_file('../stereo.i', base_dir, link=link)
        self.copy_file('../openbabel-python.i', base_dir, link=link)


class CustomBuildExt(build_ext):
    """Custom build_ext to set SWIG options and print a better error message."""
    def finalize_options(self):
        # Setting include_dirs, library_dirs, swig_opts here instead of in Extension constructor allows them to be
        # overridden using -I and -L command line options to python setup.py build_ext.
        build_ext.finalize_options(self)
        include_dirs, library_dirs = locate_ob()
        self.include_dirs.append(include_dirs)
        self.library_dirs.append(library_dirs)
        self.swig_opts = ['-c++', '-small', '-O', '-templatereduce', '-naturalvar']
        self.swig_opts += ['-I%s' % i for i in self.include_dirs]
        print('- include_dirs: %s\n- library_dirs: %s' % (self.include_dirs, self.library_dirs))

    def swig_sources(self, sources, extension):
        try:
            return build_ext.swig_sources(self, sources, extension)
        except DistutilsExecError:
            print('\nError: SWIG failed. Is Open Babel installed?\n'
                  'You may need to manually specify the location of Open Babel include and library directories. '
                  'For example:\n'
                  '  python setup.py build_ext -I/usr/local/include/openbabel-2.0 -L/usr/local/lib\n'
                  '  python setup.py install')
            sys.exit(1)


obextension = Extension('_openbabel', ['openbabel-python.i'], libraries=['openbabel'])


setup(name='openbabel',
      version=__version__,
      author=__author__,
      author_email=__email__,
      license=__license__,
      url='http://openbabel.org/',
      description='Python interface to the Open Babel chemistry library',
      long_description=long_description,
      zip_safe=True,
      cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall, 'sdist': CustomSdist},
      py_modules=['openbabel', 'pybel'],
      ext_modules=[obextension],
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
