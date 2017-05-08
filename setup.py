#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import subprocess
from glob import glob
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser


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


# Get some values from the setup.cfg
conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))
PACKAGENAME = metadata.get('package_name', 'axenoiraf')
DESCRIPTION = metadata.get('description', 'aXe - spectral extraction for HST')
LONG_DESCRIPTION = metadata.get('long_description', 'aXe minus the IRAF')
AUTHOR = metadata.get('author', 'STScI')
AUTHOR_EMAIL = metadata.get('author_email', 'help@stsci.edu')
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
            subprocess.check_call(['git', 'clone',
                                   'https://github.com/jhunkeler/relic.git'])
            sys.path.insert(1, 'relic')
            import relic.release
        except subprocess.CalledProcessError as e:
            print(e)
            exit(1)


version = relic.release.get_info()
relic.release.write_template(version, PACKAGENAME)

# Modifiy this if C/BLAS are not in /usr/lib.
BLAS_LIB_DIR = '/usr/lib'

# Default names of BLAS libraries
BLAS_LIB = ['cblas']
BLAS_EXTRA_LINK_ARGS = []

# Set environment variable BLAS_NOUNDERSCORES=1 if your BLAS does
# not use trailing underscores
BLAS_NOUNDERSCORES = False

BUILD_GSL = 1

# Directory containing libgsl (used only when BUILD_GSL = 1).
GSL_LIB_DIR = '/usr/local/lib'

# Directory containing the GSL header files (used only when BUILD_GSL = 1).
GSL_INC_DIR = 'cextern/gsl'

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']


# Define entry points for command-line scripts
# entry_points = {'console_scripts': []}

# entry_point_list = conf.items('entry_points')
# for entry_point in entry_point_list:
#    entry_points['console_scripts'].append('{0} = {1}'.format(entry_point[0],
#                                                              entry_point[1]))


# Note that requires and provides should not be included in the call to
# ``setup``, since these are now deprecated. See this link for more details:
# https://groups.google.com/forum/#!topic/astropy-dev/urYO8ckB2uM

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
        'Programming Language :: Python :: C',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=['numpy', 'astropy'],
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
        'build_sphinx': BuildSphinx
    },
)
