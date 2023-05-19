![aXe logo](docs/hstaxe/images/ACS_aXe02.png)

HSTaXe
======

HSTaXe is a Python package to provide a uniform process to perform spectral
extraction for the Hubble Space Telescope. HSTaXe supports all slitless
spectroscopy modes provided by the Wide Field Camera 3 (WFC3) and Advanced
Camera for Surveys (ACS).

The core ``hstaxe`` is written in ANSI C and is highly cross-platform by
leveraging ``CFITSIO``, ``GSL``, and ``WCSLIB``, which have been successfully
employed under Linux, Solaris, and MacOS X.

HSTaXe is the successor to aXe, a similar package written on PyRAF/IRAF.


Quickstart
----------
To install the latest release of ``hstaxe``, we recommend the following steps::

    conda create --name hstaxe-env "python>=3.8, <3.11"
    conda activate hstaxe-env
    conda install gsl cfitsio make automake autoconf libtool pkg-config -y
    conda install wcstools -c https://conda.anaconda.org/conda-forge/ --override-channels -y
    pip install hstaxe --no-cache-dir

For additional installation instructions, including instructions on installing older
or development versions of hstaxe, visit our full documentation:
https://hstaxe.readthedocs.io/en/latest/hstaxe/installing.html

Example notebooks can also be found on our full documentation:
https://hstaxe.readthedocs.io/en/latest/hstaxe/examples.html

To run the notebooks, you will need to install jupyter::

    pip install jupyter