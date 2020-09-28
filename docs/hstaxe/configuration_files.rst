.. _configuration_files:

Configuration of aXe tasks
==========================

The aXe tasks are configured in three different ways:

-  environment Variables

-  configuration files

-  online parameters to aXe tasks


.. _environment_variables:

Environment Variables
---------------------

All aXe tasks use the following environment variables:

-  AXE\_IMAGE\_PATH: the path where the input data is located

-  AXE\_OUTPUT\_PATH: the path where all aXe outputs, except the drizzle
   related files, will be directed

-  AXE\_DRIZZLE\_PATH: the path where the drizzle outputs will be
   directed

-  AXE\_CONFIG\_PATH: the path where the aXe configuration files are
   located

These can be set before running the aXe tasks or by a PyRAF script which
runs all the aXe tasks in the desired order.

Using csh/tcsh:

::

    setenv AXE_IMAGE_PATH   /path/to/data/
    setenv AXE_OUTPUT_PATH  /path/to/output/directory/
    setenv AXE_DRIZZLE_PATH /path/to/drizzle/directory/
    setenv AXE_CONFIG_PATH  /path/to/the/axe/config/

Using bash:

::

    export AXE_IMAGE_PATH=/path/to/data/
    export AXE_OUTPUT_PATH=/path/to/output/directory/
    export AXE_DRIZZLE_PATH=/path/to/drizzle/directory/
    export AXE_CONFIG_PATH=/path/to/axe/config/


.. _configuration_files:

Configuration Files
-------------------


.. _main_configuration_file:

Main Configuration File
~~~~~~~~~~~~~~~~~~~~~~~

Many configuration parameters are read in by the aXe
tasks from a single text file which serves as the primary means to
configure the extraction process for a given mode of an instrument. A
separate Main Configuration File should be created for each of the
spectral modes of each of the instruments with which aXe tasks are to be
used. This configuration file contains a basic geometrical description
of where in the slitless image one would expect a given BEAM to be
located relative to the position of the source object in a direct image.

The character ”;” can be used to add comments to this file.
A general description of the format of the input data (location of the
science, error and data quality arrays) is also included in this file.

General configuration
^^^^^^^^^^^^^^^^^^^^^

The following keywords in the Main Configuration file are used to define
several parameters such as which extension of the input FITS images
contain the data, which keywords should be used to determine the
exposure time of the input data, etc...

-  INSTRUMENT [string] The name of the instrument to which this
   configuration file applies (optional).

-  CAMERA [string] The name of the camera (optional).

-  SCIENCE\_EXT [string or integer] The name of the FITS extension
   containing the data array (if a string) or the number of the
   extension containing the data array (if an integer).

-  ERRORS\_EXT [string or integer] The name of the FITS extension
   containing the error array (if a string) or the number of the
   extension containing the error array (if an integer). Set to ”-1” if
   no error array is to be read in.

-  DQ\_EXT [string or integer] The name of the FITS extension containing
   the data quality array (if a string) or the number of the extension
   containing the data quality array (if an integer). Set to ”-1” if no
   data quality array is to be read in.

-  DQMASK [integer] This integer value determines which bits in the data
   quality array must not be set in order that a given pixel is
   considered to be good. The integer value is logically AND’ed with the
   actual data quality value of each pixel. If the result is non-zero
   the pixel will be flagged as bad and ignored in aXe tasks. The data
   quality value assigned to each pixel in calacs and updated in axeprep
   has different flag values for the various pixel deficiencies (see the
   ACS and the WFC3 Data Handbooks for the exact codes). The flag values
   for ’new hot pixels’ and ’cosmic ray rejected pixels’ for example are
   16 and 8192, respectively. To flag both the new hot pixels and cosmic
   ray rejected pixels the value for DQMASK must be set to
   :math:`8192+16=8208`. To flag all non-zero values in the data quality
   array, DQMASK must be set to 16383.

-  EXPTIME [string or float] If set to a string, this keyword defines
   which FITS header keyword will be read in the data array FITS
   extension in order to define the exposure time of the data. If set to
   a float, then this value is used instead. The exposure time is used
   to compute pixel errors if an error extention in the flt image is
   missing and to compute the variance for optimal extraction.

-  POBJSIZE [float] The size (rms) of point-like objects. for the given
   configuration. Used for the smoothed flux conversion and as a minimum
   value for the object size (see Sect. [optpars] and box on page ).

-  SMFACTOR [float] Empirical correction factor to apply in the smoothed
   flux conversion (see Sect. [fluxconv]).

-  RDNOISE[float] The readnoise of the detector. This quantity is used
   to compute pixel errors in case that there is no explicit error
   extention in the flt-images and to compute the variance in optimal
   extraction.

-  PSFCOEFFS [float, float, ...] The numbers give the coefficients of
   the polynomial describing the variations of the PSF as a function of
   the wavelength in :math:`[nm]`

-  PSFRANGE [float, float] The two numbers give the lower and the upper
   wavelength range (in [nm]) for the application of the poynomial
   defined with the keyword PSFCOEFFS. Beyound that range the value at
   the border is used.

-  PSF\_OFFSET\_# [float] Inserting a grism or prism into the optical
   beam might degrade the PSF of the instrument, and the object widths
   as measured on a direct filter image must be corrected for this by
   applying an offset.

-  FFNAME [string] The name of the configuration file for the FITS data
   cube containing the flat-field model.

-  OPTKEY1 [string] The name of a keyword in the FITS headers to
   identify the proper extension to read (e.g. ”CCDTYPE”)

-  OPTVAL1 [string] The value that the FITS header keyword defined by
   OPTKEY1 must have in order to be selected (e.g. ”1”). The OPTKEY1,
   OPTVAL1 pair allow to select the proper chip from a multi-extension
   ACS WFC image for example.

-  REFX [int] The 2D field dependence contained in the configuration
   file is by default taken to be with respect to pixel (0,0). The
   parameter REFX and REFY can be set to different values. For example,
   these parameters can be used when a 2D field dependence with respect
   to the center of the image is required.

-  REFY [int] See REFX.

-  DRZRESOLA [real] The dispersion (in Å\ :math:`/pixel`) for the
   drizzled first order beams.

-  DRZSCALE [real] The pixelscale (in :math:`''` per pixel) in the
   cross-dispersion direction in the drizzled beams.

-  DRZLAMB0 [real] The reference wavelength (in Å) which is drizzled to
   the reference pixel in the drizzled beams.

-  DRZXINI [real] The x-value of the reference pixel in the drizzled
   images. The reference wavelength given in DRZLAM0 is drizzled to this
   reference pixel. The y-value of the reference pixel depends on the
   object width and the extraction width. For a given drizzled beam, the
   y-value of the reference pixel is at :math:`{\rm real}(ny/2)+1.0`
   where :math:`ny` is the number of rows in the drizzled beam.

-  DRZPFRAC [real] The pixfrac-value used in axedrizzle.

-  DRZKERNEL [string] The drizzle kernel to be used in axedrizzle. All
   kernels available in drizzle v2.92 are allowed. Those kernels are:
   :math:`[`\ *square,point,turbo,gaussian,tophat,lanczos2,lanczos3*\ :math:`]`.
   See the help for drizzle and astrodrizzle and note on page
   [drizz:sub:`k`\ ernel] for more details.

-  DRZROOT [string] The root name for the output files created in
   axedrizzle. The string ’hrcudf’ given as DRZROOT would result in the
   drizzled beams ’hrcudf\_ext\_ID1.fits’, ’hrcudf\_ext\_ID2.fits’, ...,
   the OAF/BAF ’hrcudf\_2.OAF/BAF’, the drizzle configuration file
   ’hrcudf.conf’, the list of drizzled images ’hrcudf\_2.lis’ and the
   dummy image ’hrcudf.fits’.

BEAM configuration
^^^^^^^^^^^^^^^^^^

There must be a description for each of the BEAMs (i.e. dispersion
orders) that are extracted. BEAMs are named using single letter
characters (’A’,’B’,’C’, etc.., for a maximum number of 26 BEAMs). All
pixel coordinates and offsets that appear in a BEAM description are in
fact offsets from the reference pixel in the BEAM (REFPIXEL## in
Aperture File). The following is defined for each BEAM:

-  Magnitude cutoffs

-  Trace description

-  Wavelength calibration description

-  Sensitivity

Magnitude cutoffs
^^^^^^^^^^^^^^^^^

[Magnitude cutoffs]

-  MMAG\_EXTRACT [float] The maximum magnitude listed in the input
   object catalog for this BEAM to be extracted during the extraction
   process. Objects fainter than this cutoff magnitude will not be
   extracted. They will however be avoided when computing the background
   estimate and will be used to flag extracted spectra for contamination
   (unless otherwise determined by the MMAG\_MARK parameter).

-  MMAG\_MARK [float] Objects that have an input catalog magnitude
   greater than this will be completely ignored and not accounted for.
   This BEAM will not be used at all for anything and will not be
   avoided when computing the background estimate.

Trace description
^^^^^^^^^^^^^^^^^

[Trace description] The following items apply to the BEAM ”#”. The
character ”A” through ”Z” should be substituted for ”#”.

-  BEAM# [int] [int] The extent of the spectrum in the row (X) direction
   with respect to the reference pixel of this BEAM. The location of the
   reference pixel of this beam with respect to the direct image
   position is defined by the parameters XOFF and YOFF listed below. The
   beam row extent is measured independently of the position angle and
   always along the column direction.

-  DYDX\_ORDER\_# [int] The order of the polynomial
   :math:`{\Delta y} = P({\Delta x}) = a_0+a_1*{\Delta x}+a_2*{\Delta x}^2+...`.
   :math:`\Delta x` and :math:`P({\Delta x})` which determines the
   actual location of the trace of the spectrum in this BEAM (See
   description of this process in Chapt.[fig\ :sub:`g`\ eometry1]).

-  DYDX\_#\_0 [int] [...] For each of the orders n as specified by the
   DISP\_ORDER\_A, an entry of the form DYDX\_#\_n must exist. This can
   be a field dependent representation (see note in Chapt. [fdepend]).

-  XOFF [float] A pixel row offset between the reference pixel of this
   BEAM and the position of the object in the Direct Image. This can be
   a field dependent representation (see note on page [fdepend]).

-  YOFF [float] A pixel column offset between the reference pixel of
   this BEAM and the position of the object in the Direct Image. This
   can be a field dependent representation (see note on page [fdepend]).

Wavelength calibration description for grisms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[Wavelength calibration description] The wavelength calibration is
handled using an :math:`n^{th}` order polynomial which, as is the case
for the Trace description, can be field dependent. The field dependence
format is the same as for the trace description.

-  DISP\_ORDER\_# [int] The order of the polynomial of the form
   :math:`\lambda(x_i) = a_0+a_1*x_i+a_2*{x_i}^2+...` which defines the
   wavelength at a distance :math:`x_i` along the spectral trace.

-  DLDP\_#\_0 [int] [..] Value of the parameter :math:`a_0`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLDP\_#\_1 [int] [..] Value of the parameter :math:`a_1`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLDP\_#\_2 [int] [..] Value of the parameter :math:`a_2`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLDP\_#\_\ :math:`n` Value of the parameter :math:`a_n`, which can be
   a field dependent representation as described for the Trace
   description.

Wavelength calibration description for prisms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PWavelength calibration description] The wavelength calibration is
handled using an :math:`n^{th}` order inverse polynomial which, as is
the case for the trace description, can be field dependent. The field
dependent format is the same as for the trace description.

-  DISP\_ORDER\_# [int] The order of the inverse polynomial of the form
   :math:`\lambda(x_i) = a_1+a_2/(x_i-a_0)+a_3/(x_i-a_0)^2+...`

-  DLD1P\_#\_0 [int] [..] Value of the parameter :math:`a_0`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLD1P\_#\_1 [int] [..] Value of the parameter :math:`a_1`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLD1P\_#\_2 [int] [..] Value of the parameter :math:`a_2`, which can
   be a field dependent representation as described for the Trace
   description.

-  DLD1P\_#\_\ :math:`n` Value of the parameter :math:`a_n`, which can
   be a field dependent representation as described for the trace
   description.

-  DLD1P\_#\_PRANGE [int] [int] In the form of the dispersion relation
   given above, the singularity at :math:`x_i=a_0` divides the inverse
   polynomial into the two branches :math:`x_i-a_0<0` and
   :math:`x_i-a_0>0`. The desired solution for the dispersion relation
   is on only one branch. The finite pointspread function and extended
   sources however require a beam definition which extends from the
   valid branch over the singularity at :math:`x_i=a_0` partly into the
   second, invalid branch. To avoid that pixels from the invalid branch
   enter the PET and the spectra, this keyword defines the minimum and
   maximum values for :math:`x_i-a_0` which are allowed in the PET. Thus
   pixels from the invalid branch can be excluded.

Sensitivity
^^^^^^^^^^^

[Sensitivity] The absolute sensitivity calibration is handled by
applying a sensitivity curve to the electron count rates at each
wavelength.

-  SENSITIVITY\_# [string] The name of a sensitivity FITS file. If no
   sensitivity is available this keyword can be set to "None" instead of
   a real filename.

Example of a Main Configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See Chapter [Main Configuration File].
