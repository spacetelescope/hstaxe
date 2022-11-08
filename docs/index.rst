.. hstaxe Manual documentation master file, created by
   sphinx-quickstart on Tue Feb 19 15:56:53 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: hstaxe/images/ACS_aXe02.png


Welcome to the hstaxe package documentation!
============================================

`hstaxe` is a Python package to provide a uniform process to perform spectral
extraction for the Hubble Space Telescope. `hstaxe` supports all slitless
spectroscopy modes provided by the Wide Field Camera 3 (WFC3) and Advanced
Camera for Surveys (ACS).

The core `hstaxe` is written in ANSI C and is highly cross-platform by
leveraging `CFITSIO`, `GSL`, and `WCSLIB`, which have been successfully
employed under Linux, Solaris, and MacOS X.

`hstaxe` is the successor to aXe, a similar package written on PyRAF/IRAF.

.. note::
   This documentation mostly overlaps with, but is not identical to, the original
   aXe manual. The last version of the aXe manual can be downloaded from
   `the hstaxe repository <https://github.com/spacetelescope/hstaxe/tree/main/docs/aXe_manual/axe-user-manual-v-2-3.pdf>`_.

Attribution
-----------
This software was originally developed by the ACS group of the
Space Telescope - European Coordinating Facility (ST-ECF). The ST-ECF is a
department jointly run by the European Space Agency and the European Southern
Observatory. It is located at the ESO headquarters at Garching near Munich.
The ST-ECF staff supports the European astronomical community in exploiting the
research opportunities provided by the earth-orbiting Hubble Space Telescope.

STECF supported the use of the aXe software and the slitless spectroscopic
modes of ACS and WFC3 until late 2010. Since the beginning of 2011 support has
been provided by STScI. Please send questions to: http://hsthelp.stsci.edu

The developers have included, roughly in order of who has worked on the
software:

    Norbert Pirzkal, ST-ECF/STScI
    Markus Demleitner, ST-ECF
    Martin Kuemmel, ST-ECF
    Richard Hook, ST-ECF/STScI
    Howard Bushouse, STScI
    Megan Sosey, STScI
    Richard O'Steen, STScI
    Duy Nguyen, STScI


Using HSTaXe
------------
.. toctree::
   :maxdepth: 3
   :numbered: 3

    Installation <hstaxe/installing.rst>
    Example Notebooks <hstaxe/examples.rst>
    Description <hstaxe/description.rst>
    Calibration Files <hstaxe/calibration_files.rst>
    Using Axe <hstaxe/using_axe.rst>
    aXe Tasks <hstaxe/axe_tasks.rst>
    Configuration Files <hstaxe/configuration_files.rst>
    File Formats <hstaxe/file_formats.rst>
    Appendix <hstaxe/appendix.rst>
    Bibliography <hstaxe/bibliography.rst>
