# HSTaXe

`hstaxe` is a Python package to extract spectra from all
slitless spectroscopy modes provided by the [Hubble Space
Telescope](https://www.stsci.edu/hst) [Wide Field Camera
3 (WFC3)](https://www.stsci.edu/hst/instrumentation/wfc3)
and the [Advanced Camera for Surveys
(ACS)](https://www.stsci.edu/hst/instrumentation/acs). It supports the
features provided by the previous IRAF-based aXe package.


## Installation

### UPDATE (May 2022)

The WFC3 team has recently tested `hstaxe` installation and found
three distinct pathways to get the package installed successfully.
Each of these installations were tested with the [HSTaXe WFC3
cookbook](https://github.com/npirzkal/aXe_WFC3_Cookbook) and were found
to execute all the cookbook steps successfully. NOTE: These installation
pathways were tested on macOS Catalina.


### Option 1. Installing the PyPI release of HSTaXe (v.1.0.1):

Installing the latest released version of `hstaxe` with `pip` involves
the following steps (`conda` is used to install some non-Python
dependencies). Note that this installation method is currently
incompatible with python versions greater than 3.9:

    conda create --name hstaxe_test python=3.9 -y
    conda activate hstaxe_test
    conda install gsl cfitsio make automake autoconf libtool pkg-config -y
    conda install wcstools -c https://ssb.stsci.edu/astroconda -y
    pip install hstaxe

Optionally install Jupyter:

    pip install jupyter


### Option 2. Installing the Astroconda release of HSTaXe (v.1.0.0):

The following steps should be followed to install
`hstaxe` from using `conda` using the [Astroconda
channel](https://astroconda.readthedocs.io/en/latest).
This installation should start by downloading only the
[conda_environment.yml](https://raw.githubusercontent.com/spacetelescope/hstaxe/master/conda_environment.yml) file from this repository
followed by creating a conda environment using this file:

    conda create --name hstaxe_test --file conda_environment.yml

Activate this environment:

    conda activate hstaxe_test

Install HSTaXe from Astroconda:

    conda install hstaxe -c https://ssb.stsci.edu/astroconda


### Option 3. Installing the development version of HSTaXe from source distribution:

These instructions describe the installation of the latest development
version of the HSTaXe source code. This workflow involves the the
creation of a new conda environment, retrieval of the source code
from this repository (https://github.com/spacetelescope/hstaxe), and
the installation of the HSTaXe software and its dependencies in that
environment. This method has been confirmed to work with python 3.8 - 3.10.

Create and activate a conda environment where `hstaxe` and all its
dependencies will be installed in. We call it `hstaxe_test` in this
example:

    conda create --name hstaxe_test python=3 -y

    conda activate hstaxe_test

Install supporting packages before installing `hstaxe`:

    conda install gsl cfitsio make automake autoconf libtool pkg-config -y

    conda install wcstools -c conda-forge -y

Clone the `hstaxe` software from this repository:

    git clone https://github.com/spacetelescope/hstaxe.git

Install `hstaxe` (and its dependencies):

    cd hstaxe
    pip install .


## Notebooks and Examples

#### WFC3

The aXe WFC3 cookbook data can be downloaded by cloning this repository:
https://github.com/npirzkal/aXe_WFC3_Cookbook

1. git clone https://github.com/npirzkal/aXe_WFC3_Cookbook.git

2. Enter the `aXe_WFC3_Cookbook` directory and work through
   `aXe_WFC3_Cookbook.ipynb` within the conda environment created as
   part of the `hstaxe` installation.

#### ACS

The aXe ACS cookbook and associated data can be downloaded from the
following area: https://stsci.box.com/s/eo98zjtyccnoq7z73akfrx94jog3pg7j

#### Examples

There are also a couple simple test scripts in hstaxe/tests that can be
run against the example data in the aXe cookbook:

  * `run_acs_cookbook.py`: runs basic aXe against ACS data

  * `run_cookbook.py`: runs basic aXe against WFC3 data from the cookbook

  * `run_cookbook_part2.py`: runs axedrizzle against WFC3 data from the cookbook


## Software Development History

This software was originally developed by the ACS group of the Space
Telescope - European Coordinating Facility (ST-ECF). The ST-ECF is a
department jointly run by the European Space Agency and the European
Southern Observatory. It is located at the ESO headquarters at Garching
near Munich. The ST-ECF staff supports the European astronomical
community in exploiting the research opportunities provided by the
earth-orbiting Hubble Space Telescope.

The developers have included, roughly in order of who has worked on the
software:

    Norbert Pirzkal, ST-ECF/STScI
    Markus Demleitner, ST-ECF
    Martin Kuemmel, ST-ECF
    Richard Hook, ST-ECF/STScI
    Howard Bushouse, STScI
    Megan Sosey, STScI


### Installation Requirements (for reference)

`hstaxe` has the following Python requirements:

  * Python 3.7 or later
  * numpy
  * astropy
  * stwcs
  * stsci.imagestats
  * drizzlepac
  * photutils < 1.1.0
  * drizzle

In addition to the above Python packages, `hstaxe` depends on the
following C-based packages to be installed:

  * [cfitsio](https://heasarc.gsfc.nasa.gov/lheasoft/fitsio/fitsio.html)
  * [gsl](https://www.gnu.org/software/gsl/)
  * [libwcs](http://tdc-www.harvard.edu/software/wcstools/subroutines/libwcs.wcs.html) from [wcstools](http://tdc-www.harvard.edu/wcstools/)

These packages must be installed in your environment before you can
install hstaxe. They can be compiled and installed locally or installed
using pip and/or conda.

You will also need a C compiler suite (e.g., ``gcc`` or ``clang``) to
build `hstaxe` via [pip](https://pip.pypa.io/en/latest/) from either
PyPI or from the source distribution. In addition to the above required
packages, the following packages are required to build `hstaxe`:

  * make
  * automake
  * autoconf
  * libtool
  * pkg-config

Note that you may also need to install the
[c-blosc](https://github.com/Blosc/c-blosc) package as a dependency
of the [tables](https://pypi.org/project/tables/) package on certain
platforms (e.g., this is currently needed for macOS and Windows with
Python 3.9 due to missing `tables` wheels on PyPI).
