.. image:: docs/hstaxe/images/ACS_aXe02.png
    :width: 400
    :alt: Jdaviz logo
    :align: center

# HSTaXe

`hstaxe` is a Python package to provide a uniform process to perform spectral
extraction for the Hubble Space Telescope. `hstaxe` supports all slitless
spectroscopy modes provided by the Wide Field Camera 3 (WFC3) and Advanced
Camera for Surveys (ACS).

The core `hstaxe` is written in ANSI C and is highly cross-platform by
leveraging `CFITSIO`, `GSL`, and `WCSLIB`, which have been successfully
employed under Linux, Solaris, and MacOS X.

`hstaxe` is the successor to aXe, a similar package written on PyRAF/IRAF.


## Quickstart
To install the latest release of `hstaxe`:

    conda create --name hstaxe-env -y
    conda activate hstaxe-env
    conda install gsl cfitsio make automake autoconf libtool pkg-config -y
    conda config --add channels conda-forge
    conda install wcstools -c conda-forge -y
    pip install hstaxe

For additional installation instructions, visit our full documentation:
https://hstaxe.readthedocs.io/en/latest/hstaxe/installing.html

Example notebooks can also be found on our full documentation:
https://hstaxe.readthedocs.io/en/latest/hstaxe/notebooks.html
