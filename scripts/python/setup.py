#!/usr/bin/env python
import os
import re
import subprocess
import sys
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.errors import DistutilsExecError
from distutils.version import StrictVersion
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools import setup, Extension


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))


def find_version():
    """Extract the current version of these python bindings from the __init__.py file."""
    try:
        with open(os.path.join(base_dir, 'openbabel', '__init__.py')) as fp:
            for line in fp:
                version_match = re.match(r'^__version__ = "(.+?)"$', line)
                if version_match:
                    return version_match.group(1)
            raise Exception('Could not find version string in openbabel/__init__.py.')
    except IOError:
        raise Exception('Could not find openbabel/__init__.py.')


def pkgconfig(package, option):
    """Wrapper around pkg-config command line tool."""
    try:
        return subprocess.check_output(['pkg-config', option, package],
                                       universal_newlines=True).strip()
    except subprocess.CalledProcessError:
        raise Exception('Failed to run pkg-config')


def locate_ob():
    """Try use pkgconfig to locate Open Babel, otherwise guess default location."""
    try:
        # Warn if the (major, minor) version of the installed OB doesn't match these python bindings
        py_ver = StrictVersion(find_version())
        py_major_ver, py_minor_ver = py_ver.version[:2]
        pcfile = 'openbabel-{}'.format(py_major_ver)
        ob_ver = StrictVersion(pkgconfig(pcfile, '--modversion'))
        if not ob_ver.version[:2] == py_ver.version[:2]:
            print('Warning: Open Babel {}.{}.x is required. Your version ({}) may not be compatible.'
                  .format(py_major_ver, py_minor_ver, ob_ver))
        include_dirs = pkgconfig(pcfile, '--variable=pkgincludedir')
        library_dirs = pkgconfig(pcfile, '--variable=libdir')
        print('Open Babel location automatically determined by pkg-config:')
    except Exception as e:
        print('Warning: %s.\nGuessing Open Babel location:' % e)
        include_dirs = '/usr/local/include/openbabel3'
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
        install.run(self)


class CustomSdist(sdist):
    """Add swig interface files into distribution from parent directory."""
    def make_release_tree(self, base_dir, files):
        sdist.make_release_tree(self, base_dir, files)
        link = 'hard' if hasattr(os, 'link') else None
        pkg_dir = os.path.join(base_dir, 'openbabel')
        self.copy_file(os.path.join('..', 'stereo.i'), pkg_dir, link=link)
        self.copy_file(os.path.join('..', 'openbabel-python.i'), pkg_dir, link=link)


class CustomBuildExt(build_ext):
    """Custom build_ext to set SWIG options and print a better error message."""
    def finalize_options(self):
        # Setting include_dirs, library_dirs, swig_opts here instead of in Extension constructor allows them to be
        # overridden using -I and -L command line options to python setup.py build_ext.
        build_ext.finalize_options(self)
        self.ob_include_dir, self.ob_library_dir = locate_ob()
        self.include_dirs.append(self.ob_include_dir)
        self.library_dirs.append(self.ob_library_dir)
        self.swig_opts = ['-c++', '-small', '-O', '-templatereduce', '-naturalvar']
        self.swig_opts += ['-I%s' % i for i in self.include_dirs]
        print('- include_dirs: %s\n- library_dirs: %s' % (self.include_dirs, self.library_dirs))

    def swig_sources(self, sources, extension):
        try:
            return build_ext.swig_sources(self, sources, extension)
        except DistutilsExecError:
            print('\nError: SWIG failed. Is Open Babel installed?',
                  'You may need to manually specify the location of Open Babel include and library directories. '
                  'For example:',
                  '  python setup.py build_ext -I{} -L{}'.format(self.ob_include_dir, self.ob_library_dir),
                  '  python setup.py install',
                  sep='\n')
            sys.exit(1)


obextension = Extension(
    'openbabel._openbabel', [os.path.join('openbabel', 'openbabel-python.i')], libraries=['openbabel']
)


setup(
    name='openbabel',
    version=find_version(),
    author='Noel O\'Boyle',
    author_email='openbabel-discuss@lists.sourceforge.net',
    license='GPL-2.0',
    url='http://openbabel.org/',
    description='Python interface to the Open Babel chemistry library',
    long_description=open(os.path.join(base_dir, 'README.rst')).read(),
    zip_safe=False,
    cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall, 'sdist': CustomSdist},
    packages=['openbabel'],
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
