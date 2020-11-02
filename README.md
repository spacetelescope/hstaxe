
Current Requirements
====================
This version of the software no longer requires IRAF or PyRAF to be installed.
see `setup.py` and `setup.cfg` for additional install dependencies. The most notable here are
the GNU Science Library, cfitsio, and wcstools. 

There are some dependency packages which are only available through the Astroconda channel,
in order to ease setup of a conda environment with specific depenedencies, use the
`conda_environment.yml` file to create a working environment on your system. 


Examples
========
There are a couple simple test scripts in hstaxe/tests that can be run against the example data in the aXe cookbook:

`run_acs_cookbook.py`: runs basic aXe against ACS data
`run_cookbook.py`: runs basic aXe against WFC3 data
`run_cookbook_part2.py`: runs axedrizzle against WFC3 data
 

Software Development History
============================

This software was originally developed by the ACS group of the Space Telescope -
European Coordinating Facility (ST-ECF). The ST-ECF is a department jointly
run by the European Space Agency and the European Southern Observatory.
It is located at the ESO headquarters at Garching near Munich. The ST-ECF
staff supports the European astronomical community in exploiting the research
opportunities provided by the earth-orbiting Hubble Space Telescope.

The Developers have included, roughly in  order of who has worked on the software:

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



