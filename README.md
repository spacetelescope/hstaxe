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
