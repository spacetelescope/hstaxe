
Current Requirements
====================
This version of the software no longer requires IRAF or PyRAF to be installed.
There are some deprecated packages which are still used:
    - pydrizzle


Software Development History
============================

This software was originally developed by the ACS group of the Space Telescope -
European Coordinating Facility (ST-ECF). The ST-ECF is a department jointly
run by the European Space Agency and the European Southern Observatory.
It is located at the ESO headquarters at Garching near Munich. The ST-ECF
staff supports the European astronomical community in exploiting the research
opportunities provided by the earth-orbiting Hubble Space Telescope.

The Developers have included:
Megan Sosey, STScI
Howard Bushouse, STScI
Richard Hook, ST-ECF/STScI
Martin Kuemmel, ST-ECF
Norbert Pirzkal, ST-ECF/STScI

- Version 0.12 for STSDAS 3.3 release, October 2004
- Version 0.20: PyDrizzle will be automatically run to generate coeffs files
if not already present for input image.
-Version 0.21: Syntax for calling PyDrizzle updated for new 'bits' syntax.
    - Created wtran to use wtraxy and wtranback, and hence be MUCH faster in
      list mode.
    - Updated 'raise' statements to use Exception class, to conform to latest
      python syntax rules.
