# hstaxe

`hstaxe` is a Python package to extract spectra from all
slitless spectroscopy modes provided by the [Hubble Space
Telescope](https://www.stsci.edu/hst) [Wide Field Camera
3 (WFC3)](https://www.stsci.edu/hst/instrumentation/wfc3)
and the [Advanced Camera for Surveys
(ACS)](https://www.stsci.edu/hst/instrumentation/acs). It supports the
features provided by the previous IRAF-based aXe package.


## Installation

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
  * gsl
  * [libwcs](http://tdc-www.harvard.edu/software/wcstools/subroutines/libwcs.wcs.html) from [wcstools](http://tdc-www.harvard.edu/wcstools/)


### Installing using conda

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
    pip install ".[all]"


## Examples

There are a couple simple test scripts in hstaxe/tests that can be run
against the example data in the aXe cookbook:

  * `run_acs_cookbook.py`: runs basic aXe against ACS data

  * `run_cookbook.py`: runs basic aXe against WFC3 data

  * `run_cookbook_part2.py`: runs axedrizzle against WFC3 data

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


Some versioning feature notes from the previous IRAF software package:

- Version 2.4.4:
    axeprep doesn't require direct image column for input list file

- Version 2.4.3:
    - #978 - added full support for astrodrizzle cases with dimen_info not zero
    - #779 - check for subarray input images in axeprep and return an error
    - #678 - dont do background subtraction in axeprep when < 10% pixels are deemed good
    - #1033 - I put a catch in so that if the C code reports 0 good sky pixels no background subtraction is done and it doesn't fatal error out

- Version 2.4.2:
    - #1017 resolved. User reported an error with the "SCI" extension when the segementation image was read in fcubeprep. The code was updated to be more generic in the way headers are read and constructed.

- Version 0.20: PyDrizzle will be automatically run to generate coeffs files
    - if not already present for input image.

- Version 0.21: Syntax for calling PyDrizzle updated for new 'bits' syntax.
    - Created wtran to use wtraxy and wtranback, and hence be MUCH faster in
      list mode.
    - Updated 'raise' statements to use Exception class, to conform to latest
      python syntax rules.

- Version 0.12 for STSDAS 3.3 release, October 2004
