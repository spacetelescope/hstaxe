The source distribution for the aXe spectral extracion software
===============================================================

Authors of this package have included:
Martin Kuemmel, Markus Demleitner, Nor Pirzkal, Megan Sosey

Originally managed by STECF, the software was transferred to STScI
in 2012 for maintainence and updates.

1. The aXe software
-------------------
The aXe software was designed to extract spectra in a consistent
manner from all slitless spectroscopy modes provided by the
Advanced Camera for Surveys (ACS)  and the Wide Field Camera 3
which are installed on the Hubble Space Telescope. aXe was designed
to be quite general and can also be used for extraction of slitless
spectra from other instruments.

A full description of aXe and its capabilities can also be found
at: http://axe-info.stsci.edu/

The purpose of the current document is restricted to providing
an explanation on how to compile and install the C-code of aXe
(see also Chapter 2 of the aXe manual).


2. Requirements
---------------
The following are required to compile and install the aXe C code:

- GNU CC 2.95 or later compiler
- GNU Scientific Libraries 1.9 or later (http://www.gnu.org/software/gsl/)
- WCStools libraries 3.5 or later (http://tdc-www.harvard.edu/wcstools/)
- CFITSIO 3.x libraries (http://heasarc.nasa.gov/lheasoft/fitsio/fitsio.html)

Install the libraries according to the instructions given with the
respective source or binary distributions.

In case that you do not have the privileges to install the libraries
at their default location (e.g. "/usr/local/lib" or similar) it
may be wise to install them at a common location, e.g. under
"/yourpath/axelibs/", where "yourpath" is the location of a directory
accessible to you. In this case do not forget to:
a) set the location of the libraries on the library path LD_LIBRARY_PATH
   (e.g. by inserting the line
    "setenv  LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/yourpath/axelibs/lib"
    into your .cshrc);
b) make the GSL-script "gsl-config" accessible by either putting its
   location ("/yourpath/axelibs/bin") on your PATH-variable
   or defining the variable GSL_CONFIG by inserting the
   line "setenv GSL_CONFIG /yourpath/axelibs/bin/gsl-config"
   into your .cshrc.


3. Configuring aXe
------------------
The configure script generates a Makefile which can be used to compile
the aXe tasks. After all the above preparation it should be
straightforward to configure aXe. In the case that the libraries are
at the usual locations this is done by typing

./configure

in the directory where the aXe-source was extracted.

In the event that the libraries listed above are not installed in the
usual places, online parameters should be used to tell the configure
script where to find them. If you installed the libraries at the
locations recommended above, typing

>./configure --with-cfitsio-prefix=/yourpath/axelibs
--with-wcstools-prefix=/yourpath/axelibs/wcstools-3.x.x

should do the job. In case things do not go smoothly, follow the
instructions which go along with the error messages to fix problems.
Also the file "config.log" records continually the status information
of the configuration process. A detailed description of the errors
and possible solutions can be found there, too.

The configure script should create a proper Makefile on Linux, Solaris
and MacOSX platforms. Note that under MacOSX the option "-build=powerpc"
must be added to the  configure command for Power PC machines, and
the option "-build=i386-pc-macosx" for Intel Macs .

4. Compiling and installing aXe
-------------------------------
The next step consists of actually building the aXe library
and aXe tasks:

>make

to execute the Makefile and create the c-tasks. The tasks must be
installed in the bin directory within the iraf-directory. Simply execute

>make install

to move the binaries to their proper location.

5. Validating the aXe installation
----------------------------------
Test data with grism and prism images can can be obtained from the aXe
web site (http://axe-info.stsci.edu/).


