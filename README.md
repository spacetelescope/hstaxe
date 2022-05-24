# hstaxe

`hstaxe` is a Python package to extract spectra from all
slitless spectroscopy modes provided by the [Hubble Space
Telescope](https://www.stsci.edu/hst) [Wide Field Camera
3 (WFC3)](https://www.stsci.edu/hst/instrumentation/wfc3)
and the [Advanced Camera for Surveys
(ACS)](https://www.stsci.edu/hst/instrumentation/acs). It supports the
features provided by the previous IRAF-based aXe package.


## Installation

----
----
### UPDATE (May, 2022)!!!

The WFC3 team has recently tested `hstaxe` installation and found three distinct pathways to get the package installed successfully. Each of these installations were tested with The HSTaXe WFC3 cookbook (`aXe_WFC3_Cookbook`; available from this repository: https://github.com/npirzkal/aXe_WFC3_Cookbook) and were found to execute all the cookbook steps successfully. NOTE: These installation pathways were tested on macOS Catalina.

#### ----

#### A. Installing the development version of HSTaXe from source distribution:

These instructions describe the installation of the latest development version of the HSTaXe source code. This workflow involves the retrieval of the source code from this repository (https://github.com/spacetelescope/hstaxe), the creation of a new conda environment, and the installation of the HSTaXe software and its dependencies in that environment. Some of the dependancies need to be installed using pip owing to bugs in importing these packages if installed using conda.

Create a directory where `hstaxe` would be cloned. cd into that directory. We call this directory `test` in this example:

    mkdir test
    cd test

Create a conda environment where `hstaxe` and all its dependencies will be installed in. We call it `hstaxe_test` in this example:

    conda create --name hstaxe_test python=3

Activate this environment:

    conda activate hstaxe_test

Install supporting packages before installing `hstaxe`:

    conda install numpy astropy gsl cfitsio wcstools stwcs stsci.imagestats jupyter make automake autoconf libtool pkg-config

Due to recent issues with Astroconda, several other dependancies need to be installed in our environment via `pip`:

    pip install drizzle drizzlepac gwcs photutils tweakwcs spherical-geometry

Clone the `hstaxe` software from this repository:

    git clone https://github.com/spacetelescope/hstaxe.git

Install `hstaxe`:

    cd hstaxe
    python setup.py install
    ** OR (if installing with pip) **
    pip install .

#### ----

#### B. Installing the PyPI release of HSTaXe (v.1.0.1):

Installing the latest released version of `hstaxe` with `pip` (instead of retrieving the source code of the development version as described above) involves the following steps:

    conda create --name hstaxe_test python=3
    conda activate hstaxe_test
    conda install numpy astropy gsl cfitsio wcstools stwcs stsci.imagestats jupyter make automake autoconf libtool pkg-config
    pip install drizzle drizzlepac gwcs photutils tweakwcs spherical-geometry hstaxe

NOTE: The last step installs the latest released version of HSTaXe along with some of its dependencies.

#### ----

#### C. Installing the Astroconda release of HSTaXe (v.1.0.0):

The following steps should be followed to install `hstaxe` from Astroconda. This installation should start by downloading only the `conda_environment.yml` file from this repository (https://github.com/spacetelescope/hstaxe) followed by creating a conda environment using this file:

    conda create --name hstaxe_test --file conda_environment.yml

Activate this environment:

    conda activate hstaxe_test

Install HSTaXe from Astroconda:

    conda install hstaxe -c https://ssb.stsci.edu/astroconda

#### ----

Please follow the following steps to retrieve and run the aXe_WFC3_Cookbook:

1. git clone https://github.com/npirzkal/aXe_WFC3_Cookbook.git

2. Enter the 'aXe_WFC3_Cookbook' directory and work through aXe_WFC3_Cookbook.ipynb within the conda environment created as part of the `hstaxe` installation.

----
----

### Requirements

`hstaxe` has the following Python requirements:

  * Python 3.7 or later
  * Numpy
  * Astropy
  * stwcs
  * stsci.imagestats
  * drizzlepac
  * photutils < 1.1.0
  * drizzle

In addition to the above Python packages, `hstaxe` depends on the
following C-based packages to be installed:

  * cfitsio
      * [CFITSIO](https://heasarc.gsfc.nasa.gov/lheasoft/fitsio/fitsio.html)
  * gsl
      * [GSL](https://www.gnu.org/software/gsl/)
  * wcstools
      * [libwcs](http://tdc-www.harvard.edu/software/wcstools/subroutines/libwcs.wcs.html) from [wcstools](http://tdc-www.harvard.edu/wcstools/)

These packages must be installed in your environment before you can install hstaxe using pip or setup. They can be compiled and installed locally or installed using conda.


### Installing `hstaxe` using conda

`hstaxe` can be installed with
[conda](https://docs.conda.io/en/latest/) if you have installed
[Anaconda](https://www.anaconda.com/products/individual) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html)
`hstaxe` is available via the [Astroconda
channel](https://astroconda.readthedocs.io/en/latest) and can be
installed by running:

    conda install hstaxe -c https://ssb.stsci.edu/astroconda


### Installing using pip

You will need a C compiler suite (e.g., ``gcc`` or ``clang``) to build
`hstaxe` via [pip](https://pip.pypa.io/en/latest/) from either PyPI
or from the source distribution. In addition to the above required
packages, the following packages are required to build `hstaxe`:

  * make
  * automake
  * autoconf
  * libtool
  * pkg-config

Note that you may also need to install the
[c-blosc](https://github.com/Blosc/c-blosc) package as a dependency
of the [tables](https://pypi.org/project/tables/) package on certain
platforms (e.g., this is currently needed for OSX and Windows with
Python 3.9 due to missing `tables` wheels on PyPI).

After all of the required packages are installed, to install the latest
released version of `hstaxe` with `pip`, run (*NOTE this does not work
for versions <= 1.0.0*):

    pip install hstaxe

The latest development version of the `hstaxe` source code can be
retrieved and installed via:

    git clone https://github.com/spacetelescope/hstaxe.git
    cd hstaxe
    pip install .

Note: using pip to install will quiet all build messages by default. If you would like
to see the full output from all build commands, including the C compile statements, use `pip install . -v`


## Examples

There are a couple simple test scripts in hstaxe/tests that can be run
against the example data in the aXe cookbook:

  * `run_acs_cookbook.py`: runs basic aXe against ACS data

  * `run_cookbook.py`: runs basic aXe against WFC3 data from the cookbook

  * `run_cookbook_part2.py`: runs axedrizzle against WFC3 data from the cookbook

  * `aXe_WFC3_Cookbook.ipynb`: runs through a full WFC3 dataset cookbook style

The aXe WFC3 cookbook data can be downloaded by cloning this repository:
https://github.com/npirzkal/aXe_WFC3_Cookbook

The aXe ACS cookbook and associated data can be downloaded from the
following area: https://stsci.box.com/s/eo98zjtyccnoq7z73akfrx94jog3pg7j


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
