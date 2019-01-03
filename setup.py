#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import importlib

from glob import glob
from setuptools import find_packages, Command
from setuptools.command.test import test as TestCommand
from subprocess import check_call, CalledProcessError

from setuptools import setup
from setuptools import Command
from setuptools.command.install import install

from distutils.command.clean import clean
from configparser import ConfigParser

AXELIB_DIR = "cextern/aXe_c_code/"
CONF_H_NAME = os.path.join(AXELIB_DIR, "config.h")

# We only need to compile with these
AXE_SOURCES = glob(AXELIB_DIR+"/src/aXe*.c")
AXELIB_DEFINES = [("HAVE_CONFIG_H", "1")]

# Create the axe module extension
# no-strict-prototypes is inserted in order
# to cut back on the cfitsio lib warnings
# axe_module = Extension("axe",
#                        sources=AXE_SOURCES,
#                        include_dirs=[AXELIB_DIR],
#                        define_macros=AXELIB_DEFINES,
#                        depends=[CONF_H_NAME],
#                        language="C",
#                        extra_compile_args=['-Wno-strict-prototypes'],
#                        )

# hack building the sphinx docs with C source
try:
    import sphinx
    from sphinx.setup_command import BuildDoc

    class BuildSphinx(BuildDoc):
        """Build Sphinx documentation after compiling C source files"""

        description = 'Build Sphinx documentation with C extensions'

        def initialize_options(self):
            BuildDoc.initialize_options(self)

        def finalize_options(self):
            BuildDoc.finalize_options(self)

        def run(self):
            build_cmd = self.reinitialize_command('build_ext')
            build_cmd.inplace = 1
            self.run_command('build_ext')
            sphinx.build_main(['setup.py', '-b', 'html',
                               './docs', './docs/_build/html'])

except ImportError:
    class BuildSphinx(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print('!\n! Sphinx is not installed!\n!', file=sys.stderr)
            exit(1)


class PyTest(TestCommand):

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = [PACKAGENAME]

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        try:
            import pytest
        except ImportError:
            print('Unable to run tests...')
            print('To continue, please install "pytest":')
            print('    pip install pytest')
            exit(1)

        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


class MyClean(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("cleaning")
        try:
            check_call(["make", "clean"], cwd=AXELIB_DIR)
        except CalledProcessError as e:
            print(e)
            exit(1)

        if os.access(CONF_H_NAME, os.F_OK):
            os.remove(CONF_H_NAME)

        clean.run(self)


class BuildExtWithConfigure(install):
    """Build and install the aXe C code."""
    def run(self):
        CURRENT_ENV = sys.prefix
        try:
            check_call(["make", "clean"], cwd=AXELIB_DIR)
            check_call(["sh", "./configure",
                        "--with-cfitsio="+CURRENT_ENV,
                        "--with-wcstools="+CURRENT_ENV,
                        "--with-gsl="+CURRENT_ENV], cwd=AXELIB_DIR)
            check_call(["make"], cwd=AXELIB_DIR)
        except CalledProcessError as e:
            print(e)
            exit(1)
        install.run(self)


# Get some values from the setup.cfg
conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))
PACKAGENAME = metadata.get('package_name', 'pyaxe')
DESCRIPTION = metadata.get('description', 'aXe - spectral extraction for HST')
LONG_DESCRIPTION = metadata.get('long_description',
                                'aXe spectral extraction for HST minus the IRAF')
AUTHOR = metadata.get('author', 'STScI')
AUTHOR_EMAIL = metadata.get('author_email', 'https://hsthelp.stsci.edu')
LICENSE = metadata.get('license', '3-Clause BSD')
URL = metadata.get('url', 'http://axe-info.stsci.edu')


if os.path.exists('relic'):
    sys.path.insert(1, 'relic')
    import relic.release
else:
    try:
        import relic.release
    except ImportError:
        try:
            check_call(['git', 'clone',
                        'https://github.com/spacetelescope/relic.git'])
            sys.path.insert(1, 'relic')
            import relic.release
        except CalledProcessError as e:
            print(e)
            exit(1)

version = relic.release.get_info()
relic.release.write_template(version, PACKAGENAME)

# # Modifiy this if C/BLAS are not in /usr/lib.
# BLAS_LIB_DIR = '/usr/lib'

# # Default names of BLAS libraries
# BLAS_LIB = ['cblas']
# BLAS_EXTRA_LINK_ARGS = []

# # Set environment variable BLAS_NOUNDERSCORES=1 if your BLAS does
# # not use trailing underscores
# BLAS_NOUNDERSCORES = False


# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

# Define entry points for command-line scripts
# entry_points = {'console_scripts': []}

# entry_point_list = conf.items('entry_points')
# for entry_point in entry_point_list:
#     entry_points['console_scripts'].append('{0} = {1}'.format(entry_point[0],
#                                                               entry_point[1]))


# Note that requires and provides should not be included in the call to
# ``setup``, since these are now deprecated. See this link for more details:
# https://groups.google.com/forum/#!topic/astropy-dev/urYO8ckB2uM
package_info = {}
# package_info['ext_modules'] = [axe_module]

setup(
    name=PACKAGENAME,
    version=version.pep386,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    license=LICENSE,
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: C',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=['numpy',
                      'astropy',
                      'wcstools',
                      'cfitsio',
                      'gsl',
                      'stwcs',
                      'drizzlepac'],
    scripts=scripts,
    packages=find_packages(),
    tests_require=[
        'backports.tempfile',
        'pytest',
        'requests_mock',
        'pytest-catchlog'
    ],
    cmdclass={
        'test': PyTest,
        'build_sphinx': BuildSphinx,
        'build_ext': BuildExtWithConfigure,
        'install': BuildExtWithConfigure,
        'clean': MyClean,
    },
    #entry_points=entry_points,
    **package_info
)
