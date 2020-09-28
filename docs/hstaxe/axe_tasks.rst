.. _axe_tasks:

aXe tasks
=========

[aXe tasks]

This chapter gives detailed descriptions of each of the aXe tasks and of
their parameters. Examples of how to run each task are also included, as
well as descriptions of the files required and produced by each task. A
detailed description of all the output products can be found in
Chapt. [FIlesFormat].

The aXe tasks use environment variables (see Chapt.  [Env Var]) to
define the locations of the direct and slitless images. In addition, all
tasks are meant to work on specific extensions in the input FITS files
and the file names of the output products of each aXe task reflect this
by having a \_# (where # is an extension number) appended to the
original FITS file name (e.g. a grism\_1.SPC.fits file will be produced
if the input slitless FITS file was named grism.fits and if the science
data of interest was located in extension 1). Selection of the extension
to extract is defined in the configuration file (Chapt. [Main Config]) .

IOLPREP
-------

[IOLP] This task produces Input Object Lists for every input image of a
MultiDrizzle combined image. It projects the object positions from a
master catalogue, which contains all objects in the coordinate system of
the MultiDrizzled image, out into the coordinate system of each input
image. For each input image an Input Object List is generated, which
contains only the objects within the boundaries of the given input
image.

The names and drizzle parameters of the input images are retrieved from
the header of the MultiDrizzle combined image. The projection of the
object positions into the coordinate system of the input images is done
with the STSDAS task wtranback.

There is a parameter to influence the sensitive area to include objects
in the IOLs. This allows objects beyond the physical boundaries of the
input image to be included in the IOLs, to take into account partly
covered objects or to include bright objects outside of the FOV in
contamination estimates derived from those IOLs.

During the task execution, the drizzle coefficient files and the input
images must be available. For this reason it would be best practice to
run it in the directory that was used to combine the input images with
MultiDrizzle.

Usage
~~~~~

::

      iolprep comb_image master_cat dim_info

Parameters
~~~~~~~~~~

::

    mdrizzle_image: the name of the combined (MultiDrizzled) image

    input_cat:      the name of the SExtractor master catalogue

    dimension_info: four numbers to specify the modifications
                    [left, right, bottom, top] to the target area on the input
                    images. E.g. 100,500,10,0 would include in the Input Object
                    Lists all objects with -100 < x < x_size + 500 and
                    -10 < y < y_size.

    Example:
       axeprep mdrizzle_image='g800l_drz.fits' input_cat='bvri.cat'
               dimension_info=200,0,0,0

Output
~~~~~~

-  ./[grism filename]\_[science ext number].cat

FCUBEPREP
---------

[FPREP] The task produces fluxcubes for a set of grism images in a
standard scenario. The user should have prepared a MultiDrizzled grism
image and at least one MultiDrizzled direct image. Moreover there exists
a segmenation map that was produced together with the master catalogue
of objects in a SExtractor run. All the MultiDrizzled images must have
the identical coordinate system, which means every pixel :math:`(i,j)`
must represent the identical position :math:`(Ra, Dec)` on the sky

The task analyzes the header of the combined grism image and extracts
the name and the drizzle parameters of each grism input image. The
filter images are transformed from :math:`counts/sec` to flux units, and
the STSDAS task blot produces for each input grism image a set of cutout
images from the segmentation and the flux images. The cutout images for
each grism input image are finally combined to a fluxcube.

There exists the option to produce fluxcubes which have a different size
than their associated grism image. This allows to include objects
outside of the grism image FOV in the contamination estimate, or to
restrict the contamination analysis to only a restricted area.

As with iolprep the drizzle coefficient files and the grism images must
be available during task execution, and the best way to assure this
would be running the task in the directory used for the drizzle
combination of the grism images.

Usage
~~~~~

::

      fcubeprep grism segmentation filtinfo zeropoint dim_info

Parameters
~~~~~~~~~~

::

    grism_image:    the name of the combined (MultiDrizzled) grism image

    segm_image:     name of the segmentation image

    filter_info:    name, wavelength, zeropoint of the filter image(s).
                    If there are several filter images the comma separated
                    quantities are written in a file, and the name of this file
                    is given here.

    AB_zero:        boolean to indicate whether the zeropoints given in the
                    parameter 'filter_info' are in AB- or ST-magnitudes

    dimension_info: four numbers to specify the modifications
                    [left, right, bottom, top] to the target area on the grism
                    images. E.g. 100,500,10,0 would produce fluxcube images
                    which cover the area -100 < x < x_size + 500 and -10 < y < y_size
                    in the input grism images

    interpol:       the inpolation scheme used to compute flux values at
                    the interpolated wavelengths

    Example:
       fcubeprep grism_image='g800l_drz.fits' segm_image='f850_seg.fits'
                 filter_info='cubelisST.lis' ABzero='NO'

    with the file cubelisST.lis:
    f435w_drz.fits, 431.8, 25.157
    f606w_drz.fits, 591.8, 26.655
    f775w_drz.fits, 769.3, 26.393
    f850lp_drz.fits, 905.5, 25.954 

Output
~~~~~~

-  ./[grism filename]\_[ext number].FLX.fits

AXEPREP
-------

[PREP] This task prepares the science files (e.g. ACS and WFC3 \_flt
files produced by the on-the-fly pipeline or the calacs task) for
further processing within aXe. axeprep provides important keywords and
is mandatory if axedrizzle is to be used later on.

axeprep provides three different processing steps:

-  background subtraction:
   Provided that an Input Object List is given for the grism image,
   axeprep uses the tasks sex2gol, gol2af and backest to mark the beam
   areas on the grism image as well as on the master background image.
   Using the IRAF task imstat with 3 cycles of clipping pixels with
   values :math:`>3\sigma`, the median pixel values are derived for the
   unmarked pixels on both the grism image and on the master background
   image. The master background, scaled to the level of the grism image,
   is finally subtracted from the grism image.

-  exposure time normalization The input file is normalized by the
   exposure time to transform the images into counts per second.

-  gain correction The input file is multiplied by the gain conversion
   factor (electrons/ADU) to transform the images from units of detector
   counts to electrons. For HST data, this is usually only necessary for
   NICMOS images, because ACS and WFC3 images will normally already be
   in units of electrons.

Every processing step can be switched on/off independently by associated
boolean parameters.

For ACS/WFC and WFC3/UVIS images, the correspondence between the
configuration files, the background images and the IOL’s declared in the
Input Image List is explained in Fig. [inputext].

The file used by axeprep as inlist can be reused again in axecore,
drzprep and axedrizzle, perhaps extended with different dmag-values for
the grism images.

Usage
~~~~~

::

      axeprep inlist configs backgr backims mfwhm norm gaincorr

Parameters
~~~~~~~~~~

::

    inlist: Input Image List which gives on each line
            a) the name of the grism image to be processed (mandatory)
            b) the object catalog(s) (mandatory if back='yes',
               comma separated list if there is more than one catalogue)
            c) the direct image associated with the grism image (optional)

    configs: name of the aXe configuration file. If several image
             extensions are to be processed (e.g. for WFC images), one
             configuration file per extension must be given in a comma
             separated list.

    background: boolean to switch on/off background subtraction

    backims:name of the background image. If several image extensions
            are to be processed (e.g. for WFC images), one background 
            image per extension must be specified in a comma separated 
            list.

    mfwhm:  real number to specify the extent (as a multiple of A_IMAGE
            or B_IMAGE) of the area that is masked out perpendicular to
            the trace of each object before the background level is
            determined (see parameter mfwhm in gol2af).

    norm:   boolean to switch on/off the exposure time normalization

    gaincorr: boolean to switch on/off the gain conversion

    Example:
       axeprep inlist='imlist.lis', configs='conf1.conf,conf2.conf',
              back='YES', backims='back1.fits,back2.fits', fwhm=2.0,
              norm='YES' gain='NO'

If  back='YES' :

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].MSK.fits

AXECORE
-------

[AXECORE] This aXe task combines the Low Level Tasks sex2gol, gol2af,
af2pet, petff, petcont, pet2spc and stamps and offers the possibility to
make a complete aXe reduction based on the individual images in one
task. This also includes the reduction with background PETs (set
back=YES). The parameter list comprises all parameters for the
individual tasks, and as a consequence is rather long. For most of the
parameters the default value is appropriate, so the actual number of
parameters that will normally need to be edited by the user is quite
modest. In the listing below, the axecore parameters are organised
according to the low-level aXe tasks they affect.

The Input Image List used in the task axeprep as inlist can be reused
again in axecore, perhaps extended with individual dmag-values for the
grism images.

The sequence of configuration files must correspond to the sequence of
Input Object Lists and the sequence of background images in inlist (see
Fig. [inputext]).

If axedrizzle is not going to be used, the parameter drzfwhm should be
set to :math:`0.0`. Otherwise the values of parameters extrfwhm and
drzfwhm in the task axecore must correspond to the values to be used for
parameters infwhm and outfwhm in the task axedrizzle, respectively. The
extraction width used after drizzling, which is specified via drzfwhm
and outfwhm, must be larger than the extraction width for the PETs in
axecore, which is set by extrfwhm (see :ref:`axedrizzle_task`).

Usage
~~~~~

::

    axecore inlist configs extrfwhm drzfwhm back backfwhm orient slitless_geom exclude ...

Parameters
~~~~~~~~~~

::

    inlist:   Input Image List which gives on each line
              a) the name of the grism image to be processed (mandatory)
              b) the object catalog(s) (mandatory)
              c) the direct image associated with the grism image (optional)
              d) dmag value (see GOL2AF) for the grism image (optional)

    configs:  name of the axe configuration file. If several image
              extensions are to be processed (e.g. for WFC images), one
              configuration file per extension must be given in a comma
              separated list.

    back:     Boolean to switch on/off the creation of a background PET with
              mfwhm=backfwhm

    [The following parameters apply to GOL2AF:]

    extrfwhm: mfwhm value to specify the extraction width in gol2af

    drzfwhm:  mfwhm value to specify the extraction in axedrizzle

    backfwhm: mfwhm value to specify the width of the background PET

    orient:   enable tilted extraction

    slitless_geom: enable the best extraction for slitless spectroscopy

    exclude:  switch off the listing of faint objects

    lambda_mark: the wavelength at which to apply the cutoff magnitudes
                 MMAG_EXTRACT and MMAG_MARK

    [The following parameters apply to PETCONT:]

    cont_model:  name of the contamination model to be applied

    model_scale: scale factor for the gaussian contamination model

    interp_type: interpolation type for the flux values

    lambda_psf:  wavelength [nm] at which the object widths were measured

    [The following parameters apply to BACKEST:]

    np:       number of points for background estimation

    interp:   interpolation type for background determination
              (-1: GLOBAL median; 0: local median; 1: linear fit;
              2: quadratic fit)

    niter_med: number of kappa-sigma iterations around the median

    niter_fit: number of kappa-sigma iterations around the fit value

    kappa:     kappa value

    smooth_length: number of adjacent pixels on each side to use when
                   smoothing the background estimate in x-direction

    smooth_fwhm: FWHM of the Gaussian used in the background smoothing

    [The following parameters apply to PET2SPC and STAMPS:]

    spectr:   enable the creation of SPCs and STPs for each of the
              grism files individually

    weights:  compute and apply optimal weights

    adj_sens: adjust the sensitivity function for extended sources

    sampling: the sampling mode for the stamp images

    Example:
        axecore inlist='imlist.lis' configs='conf1,conf2' back='YES'
                extrfwhm=4.0 backfwhm=5.0 exclude='NO' cont_model='gauss'
                model_scale=4.0 interp_type='linear' lambda_psf=600.0
                slitless_geom='YES' np=10 interp=1 spectr='YES' adj_sens='YES'

Output
~~~~~~

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].cat

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].OAF

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].PET.fits

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].CONT.fits

If  back='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BAF

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BCK.fits

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BCK.PET.fits

If  drzfwhm>0  and cont_model='geometric':

-  $AXE\_OUTPUT\_PATH/[slitless filename][drzfwhm]\_[ext number].OAFthis
   file is used to recompute the contamination in the PET’s, using the
   value specified in drzfwhm as extraction width

If  spectr='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].STP.fits

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].SPC.fits

If  spectr='YES' and  weights='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number]\_opt.SPC.fits

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number]\_opt.WHT.fits


.. _drzprep_tasks:

DRZPREP
-------

This task produces a set of Drizzle PrePare (DPP) files for a
set of images given in an Input Image List. A DPP-file is a
multi-extension fits file with a pixel stamp image, an error stamp image
and a contamination stamp image for each first order beam in a grism
image. DRZPREP uses the PET file to derive the pixel/error/contamination
values for the stamp images and the background aperture file (BAFs) to
define a common geometry for the individual objects. The need for a
common geometry for all stamp images of a single object forces drzprep
to be run always on the set of images which later are also combined with
axedrizzle. If there is more than one set of PETs for each grism image
(as in the case of WFC data), the configuration files should be given as
a comma separated list in the parameter configs.

The task also derives and stores important keywords for axedrizzle. In
the Input Image List given with the parameter ’inlist’ the first item on
each line must be the name of the grism image. Further columns/items are
neglected by drzprep. Therefore the file used as inlist in axecore and
axeprep can be re-used in drzprep again.

Usage
~~~~~

::

      drzprep imagelist configs back

Parameters
~~~~~~~~~~

::

    inlist:  Input Image List which gives the name of the grism image to
             be processed as the first item on each line.

    configs:  name of the aXe configuration file. If several image
              extensions are to be processed (e.g. for WFC images), one
              configuration file per extension must be given in a comma
              separated list.

    opt_extr: boolean to generate also the necessary data  for optimal
              extraction in axedrizzle 
        
    back:     boolean to switch on the creation of background DPPs made
              by processing background PETs.     

    Example:
        drzprep inlist='axeprep.lis' configs='aXe_config1.conf,aXe_config2.conf'
                back='NO'


Output
~~~~~~

If  back='NO':

-  | $AXE\_DRIZZLE\_PATH/[slitless filename]\_[ext number].DPP.fits or

If  back='YES':

-  $AXE\_DRIZZLE\_PATH/[slitless filename]\_[ext number].BCK.DPP.fits


.. _axedrizzle_tasks:

AXEDRIZZLE
----------

[AXEDRIZZLE] This is the central task to the reduction method described
in Chapt. [drizzlingPETs]. This task takes the DPPs prepared by drzprep
as input. The extensions for the various objects are extracted from the
DPP, and the extracted stamp of each object are drizzled together to
form a deep, 2D drizzled grism image for each object. For a description
of the drizzle algorithm, see Fruchter & Hook (2002)

The drizzle coefficients computed by drzprep for each stamp image are
given as header keywords and are computed in such a way that the
combined 2D drizzled grism image resembles an ideal grism image with a
constant dispersion and a constant pixelscale in the cross-dispersion
direction. The trace of the drizzled spectra is parallel to the x-axis
of the image. The dispersion and the pixelscale (in cross-dispersion
direction) are set in the aXe configuration file with the keywords
DRZRESOLA and DRZSCALE , respectively (see Chapt. [Config Files]). At
present, only the first order beams of images can be drizzled.

Drizzling usually creates pixels with incomplete coverage at the borders
of the drizzle images. To avoid those pixels with their lower weight
entering the 1D extraction, the extraction width used in the 1D
extraction from the 2D drizzled grism images should be *smaller* than
the extraction width used to generate the PETs in axecore. The
extraction width (in multiples of the object fwhm) for the 1D extraction
must be specified with the parameter outfwhm, while the parameter infwhm
must be set to the value that was used in axecore to create the PETs
and therefore the DPPs. infwhm and outfwhm in the task axedrizzle
therefore directly correspond in axecore to the parameters extrfwhm and
drzfwhm, respectively. Then the task axedrizzle can recalculate the
extraction width for the 1D extraction. Typical value pairs for (infwhm,
outfwhm) and (extrfwhm, drzfwhm) are (4.0,3.0) or (3.0,2.0). A value
given in axecore cannot be corrected or changed in the task axedrizzle.

In addition to the 2D drizzled grism images, axedrizzle creates all the
necessary files to facilitate the extraction of the 1D spectra with the
tasks drz2pet and pet2spc. Usually, these additional steps are carried
out automatically within axedrizzle (if makespc=YES). To drizzle the
background DPPs, the task axedrizzle must be run with back=YES. If the
drizzling of the background is done *after* the drizzling of the object
DPPs, the background is correctly taken into account in the reduction of
the 1D spectra.

The Input Image List given with the parameter inlist must contain the
name of the grism image as the first item on each line. The name of the
corresponding DPP file(s) are then derived from the grism name and the
chip as specified in the configuration file(a). Further columns/items
are neglected. Therefore files used as inlist in axecore and axeprep
can be reused in axedrizzle again.

With driz\_separate=YES aXedrizzle identifies and excludes deviant
values coming from e.g. hot or cosmic ray affected pixels. In this mode,
axedrizzle works similar to MultiDrizzle () on direct images. As
indicated below, this method requires many parameters. Their name and
their meaning is identical as in Multidrizzle, and user are referred to
the MultiDrizzle documentation for details. axedrizzle with the
rejection of pixels works only if the background of the grism images had
been removed with global background subtraction (see task axeprep,
Chapt. [backsub]).

Usage
~~~~~

::

    axedrizzle inlist configs infwhm outfwhm back makespc

Parameters
~~~~~~~~~~

::

    inlist:   Input Image List with the input grism filename as the first item
              on each line.

    configs:  name of the aXe configuration file. If several image extensions
              (and therefore DPPs) are to be processed (e.g. for WFC images),
              one configuration file per extension must be given in a comma
              separated list.

    infwhm:   mfwhm for the input PETs and DPPs

    outfwhm:  mfwhm for the extraction of the objects in later steps

    back:     boolean to drizzle both the object and the background DPP


    clean:    boolean to remove temporary files

    makespc:  boolean to switch on/off whether SPCs shall be created directly
              from the coadded images

    adj_sens: adjust the sensitivity function for extended sources

    opt_extr: boolean to switch on the optimal extraction in addition to
              the regular, exposure time weighted extraction

    driz_separate: drizzling to separate image and CCR reject (="YES")
                   or "simple" aXedrizzle (="NO")

    **The following parameters apply only for driz_separate="YES":**

    **The following parameters apply for MEDIAN IMAGE parameters:**

    combine_type:    type of combine operation, (median|sum|minmed|minimum)

    combine_maskpt:  percentage of weight image below which it is flagged as bad

    combine_nsigma:  significance for accepting minimum instead of median
                     (for combine_type="minmed")

    combine_nlow:    number of low pixels to reject

    combine_nhigh:   number of high pixels to reject

    combine_lthresh: lower threshold for clipping input pixel values

    combine_hthresh: upper threshold for clipping input pixel values

    combine_grow:    radius (pixels) for neighbor rejection

    **The following parameters apply for BLOT BACK:**

    blot_interp:     interpolant (nearest,linear,poly3,poly5,sinc)

    blot_sinscl:     scale for sinc interpolation kernel

    **The following parameters apply for COSMIC RAYS REMOVAL:**

    driz_cr_snr:     driz_cr.SNR parameter

    driz_cr_grow:    driz_cr_grow parameter

    driz_cr_scale:   driz_cr.scale parameter

    Example:
        axedrizzle inlist="axegrism.lis" configs="HUDF.HRC.conf" infwhm=4.0
                   outfwhm=3.0 back="NO" makespc="YES" adj_sens="YES"
                   driz_separate="NO"

Output
~~~~~~

If for an input name ``./[drizzle root filename]_2.list``:

-  $AXE\_CONFIG\_PATH/[drizzle root filename].conf

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_2.OAF

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_mef\_ID[num].fits

-  BACK='YES': $AXE\_DRIZZLE\_PATH/[drizzle root
   filename]\_mef\_ID[num].BCK.fits

If makespc='YES':

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_2.SPC.fits

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_2.STP.fits

If  makespc='YES' and  opt_weight='YES':

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_2\_opt.STP.fits

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_2\_opt.WHT.fits

SEX2GOL
-------

[SEX2GOL] This task generates a Grism Object List file using an Input
Object List as input. There are three different kinds of Input Object
List that can be fed into aXe:

-  an Input Object List (in SExtractor format) of objects on a direct
   image covering (roughly) the same field as the grism image

-  an Input Object List in SExtractor format, which gives the objects on
   the grism image in world coordinates (RA, Dec and theta\_sky)

-  an Input Object List in SExtractor format, which lists the objects on
   the grism image in world coordinates and image coordinates (x\_image,
   y\_image and theta\_image)

A thorough description of the Input Object List is given in
Chapt. [SEX].

For the first two ways to specify object lists, the image coordinates of
the objects on the grism image will be recomputed using the WCS
information of the grism image and the direct image. This approach
therefore relies on the accuracy of the WCS information given in those
images. Refer to section [SEX] for a description of what values should
be in the the input catalog and which ones can be re-constructed by
SEX2GOL.

Usage
~~~~~

::

      sex2gol grism config in_sex use_direct direct dir_hdu
              spec_hdu out_sex

Parameters
~~~~~~~~~~

::

    grism:      name of the grism image to be processed.

    config:     name of the axe configuration file.

    in_sex:     name of the object file.

    use_direct: boolean to indicate that the Input Object List refers to a
                direct image

    direct:     name of the direct image

    dir_hdu:    direct image extension to be used

    spec_hdu:   grism/prism image extension to be used

    out_SEX:    overwrites the default output object catalog name

::

    Example:
         sex2gol grism='test_grismn.fits' config='SLIM.conf.test.0'
                 in_sex='test_0.cat' use_direct='NO'

Output
~~~~~~

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].cat

GOL2AF
------

[GOL2AF] This task generates an Aperture File using an input Grism
Object List and a valid configuration file which defines the length,
wavelength calibration and global offsets between direct and slitless
images. For positive numbers in mfwhm the extraction width of the BEAMs
is set up to be this number times the width listed in the Grism Object
List. Negative numbers specify the extraction width for all objects in
pixels directly (see Chapt. [ewidth] for a detailed discussion). Two
magnitude cutoffs are set in the Configuration File (Chapt.[Config
Files]). Sources which have magnitudes fainter than an extraction cutoff
magnitude are flagged so that they are not extracted, but will be
accounted for when computing spectral contamination and the background
estimates. Sources which have magnitudes fainter than another cutoff
magnitude are marked so that they will be completely ignored. The dmag
value can be used to globally adjust these cutoffs (to account for a
different signal-to-noise ratio in one dataset for example without
having to resort to editing of the configuration file).

This task can be used to generate both an Object Aperture File and a
Background Aperture File. While these files have a similar format, it is
often desirable to use different Aperture Files for the two cases. This
is because the former is used to extract counts from pixels which are
known to contain flux from the source, while the latter can be thought
to define a zone to avoid all source flux in the slitless image when
computing the background level (in the case that a master sky is not
used for background subtraction, see Chapt.[skyback]). In practice, a
larger extraction width multiplier should be used when generating the
Background Aperture File so that all the object flux is properly
isolated when generating a Background Estimate File (Chapt.[BEF]).

With orient=YES GOL2AF extracts the beams with an extraction angle
parallel to the semi-major-axis of the object. orient=YES forces a
vertical extraction perpendicular to the spectral trace of the beam. For
orient=YES and slitless\_geom=YES however GOL2AF adjusts the
extraction angle when the desired extraction angle forms too small an
angle with the spectral trace (:math:`|\alpha|<35^\circ`). Then the
extraction angle follows the semi-minor-axis instead of the
semi-major-axis of the object which results in more pixels being
extracted from the slitless image (see Chapt. [ewidth] and
Fig. [ext\ :sub:`w`\ idth] for more details).

Usage
~~~~~

::

      gol2af grism config mfwhm dmag back slitless_geom orient exclude
             sci_hdu out_af in_gol

Parameters
~~~~~~~~~~

::

    grism:       name of the grism image

    config:      name of the aXe configuration file

    mfwhm:       the extraction width multiplicative factor

    back:        to generate a BAF instead of an OAF file

    orient:      boolean to switch on/off tilted extraction

    slitless_geom: boolean to switch on/off automatic orientation
                 for the tilted extraction 

    exclude:     boolean to switch on the removal of  faint objects
                 in the result

    lambda_mark: the wavelength at which to apply the cutoff magnitudes
                 MMAG_EXTRACT and MMAG_MARK

    dmag:        a number to add to the MMAG_EXTRACT and MMAG_MARK 
                 values given in the configuration file

    out_af:      overwrites the default output OAF or BAF filename

    in_gol:      overwrites the default input catalog name

::

    Example: 
      gol2af grims='test_grismn.fits' config='SLIM.conf.test.0' mfwhm=4.0
             back='YES'

Output
~~~~~~

If  back='NO':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].OAF

If  back='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BAF

BACKEST
-------

[BE] This task uses the input slitless image and a Background Aperture
File to generate a Background Estimate File (Chapt.[BEF]) . This task is
applicable when a master sky is not used for background subtraction
(Chapt. [skyback]). The number of points to use and the order of the
interpolation to use to generate the Background Estimate File can be set
using online parameters. The values in the regions within each of the
BEAMs listed in the Background Estimate File are replaced by the median,
average, linear, or :math:`n^{th}` order polynomial interpolation of
pixels which are immediately above and below a BEAM (but not within any
other BEAM). The number of pixels to use for fitting is by default set
to 10 on each side below and above the BEAM (therefore 20 pixels in
total). The value given for the np option can be used to change this
default value. If the number of points is set to a value which is 0 or
less, then the entire column of an image will be used, ignoring any
pixels which are within any known BEAM. This option allows for a full
column background estimate to be created, instead of a local background
estimate. The type of interpolation is controlled by the parameter
interp:

-  interp= -1 ; Median

-  interp= 0 ; Average

-  interp= 1 ; Linear fit

-  interp= (:math:`n>1`) ; :math:`n^{th}` order polynomial fit ;

In case that bad pixels or cosmics are not marked in the dq-extention of
the flt-file, it is possible to execute a number of kappa-sigma klipping
iterations prior to the final fit. The kappa-sigma iterations exclude
extreme pixel values created by defects limit their impact on the fit.

To further suppress the noise in the background, it is possible to apply
a Gaussian smoothing in x-direction on the fitted background values.
This is controled by the parameters smooth\_length and smooth\_fwhm,
which give the number of adjacent pixels used for calculating the
smoothed value on each side of a background pixel and the FWHM of the
Gaussian, respectively.

Usage
~~~~~

::

      backest grism config np interp niter_med niter_fit kappa

Parameters
~~~~~~~~~~

::

    grism:     name of the grism image

    config:    name of the aXe configuration file

    np:        the number of pixels used on each side of a beam 
               to compute the median/average/fitted background

    interp:    the type of interpolation to perform

    niter_med: number of kappa-sigma iterations around the median

    niter_fit: number of kappa-sigma iterations around the fit value

    kappa:     kappa value

    smooth_length: number of adjacent pixels on each side to use when
                   smoothing the background estimate in x-direction

    smooth_fwhm: FWHM of the Gaussian used in the background smoothing

    mask:      create a mask image with the OAF file

    in_af:     overwrite the default input aperture filename

    out_back:  overwrite the default output background filename

::

    Example:   
      backest grism='test_grismn.fits' config='SLIM.conf.test.0' np=10 interp=1 

Output
~~~~~~

If  mask='NO':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BCK.fits

If  mask='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].MSK.fits

AF2PET
------

[AF2PET] This task uses the input slitless image together with an Object
Aperture File to generate a Pixel Extraction Table (PET) for the input
data. The same task should be used with the Background Estimate File and
the same Object Aperture File to generate a Background Pixel Extraction
Table containing information about the spectral background (BPET).

Usage
~~~~~

::

      af2pet grism config back out_pet

Parameters
~~~~~~~~~~

::

    grism:   name of the grism image

    config:  name of the aXe configuration file

    back:    generate a PET for a background image using
             a BAF file instead of a OAF file and using a
             background image generated by backest

    out_PET: overwrite the default output PET filename

::

    Example:
       af2pet grism='test_grismn.fits' config='SLIM.conf.test.0' back='YES'

Output
~~~~~~

If  back='NO':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].PET.fits

If  back='YES':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BCK.PET.fits

PETCONT
-------

[PETCONT] The task computes and stores the contamination information for
a given Pixel Extraction Table. There are two distinct ways to compute
the contamination:

-  The geometrical contamination records, for each PET pixel, how often
   it is a member of a different beam. If a pixel is a member of two
   separate beams, i.e. is in a region where two beams overlap, it is
   assigned a value of 1 in each of the two beam PET’s, thus indicating
   that this pixel is also part of another beam.

-  In quantitative contamination, the amount of contaminating flux from
   other beams is estimated for each PET pixel. This estimate is based
   on a model of the emitting sources. There are two different methods
   to establish an emission model , the **gaussian emission model** and
   the **fluxcube model**. Chapt. [quant\ :sub:`c`\ ont] gives a
   detailed discussion on the emission models and their implications for
   the contamination.

The right panel of Fig. [beams] is a geometrical contamination image,
which carries the basic information about contamination.

Usage
~~~~~

::

      petcont grism config cont_model model_scale inter_type cont_map

Parameters
~~~~~~~~~~

::

    grism:       name of the grism image

    config:      name of the aXe configuration file

    cont_model:  name of the contamination model to be applied

    model_scale: scale factor for the gaussian cont. model

    spec_models: name of the multi-extension fits table with model spectra

    object_models: name of the multi-extension fits image with object templates.

    interp_type: interpolation type for the flux values

    lambda_psf:  wavelength [nm] at which the object widths were measured

    cont_map:    write the contamination map into a FITS file

    in_af:       overwrites the input AF file name

::

    Example:
       petcont grism='test_grismn.fits' config='SLIM.conf.test.0' cont_map='YES'

Output
~~~~~~

-  Updates $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].PET.fits

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].CONT.fits

PETFF
-----

[PETFF] This task uses a flat-field calibration file to flat-field the
content of a Pixel Extraction Table (see Chapt.[PET]). The wavelength of
a pixel is used in conjunction with a flat-fielding data cube containing
the coefficients of a polynomial which can be used to compute at each
pixel (x,y):

:math:`FF(x,y,x) = a_0(x,y) + a_1(x,y)*x + .. +a_i * x^i`, where,

:math:`x = (\lambda - {\lambda}_{min})/(\lambda_{max} - \lambda_{min})`

The coefficients :math:`a_0(x,y)` are stored in the first data
extension of the flat-field cube, :math:`a_1(x,y)` in the second, etc...
The values for :math:`\lambda_{max}` and :math:`\lambda_{min}` are in
the FITS header keywords :math:`WMIN` and :math:`WMAX`. The name of the
flat-field cube is read from the aXe configuration file using the
parameter FFNAME . Chapter [Flat field] gives a detailed description of
the flatfield.

Usage
~~~~~

::

      petff grism config back ffname

Parameters
~~~~~~~~~~

::

    grism:  name of the grism image

    config: name of the aXe configuration file

    back:   apply FF to the background Pixel Extraction Table (BPET)

    ffname: overwrite the default input flat-field cube name

::

    Example:
       petff grism='test_grismn.fits' config='SLIM.conf.test.0' back='YES'

Output
~~~~~~

If  back='NO':

-  Updates ``$AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].PET.fits``

If  back='YES':

-  Updates ``$AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].BPET.fits``


.. _pet2spc_tasks:

PET2SPC
-------

This task is used to transform the content of an Object Pixel
Extraction Table into a set of 1D binned spectra in an Extracted Spectra
File (see Chapt.[SPC]). The binning process is explained in more detail
in Chapt. [binning] and allows the application of optimal weights (see
Chapt. [optimal weighting]).

The task can be used simultaneously with both an Object Pixel Extraction
Table and a Background Pixel Extraction Table, in which case a
background subtraction is performed. Care must be taken that both Object
and Background Pixel Extraction Tables were created with the same
Aperture File. Additionally, absolute flux calibration can be performed
if the proper information is included in the Main Configuration File.

Usage
~~~~~

::

      pet2spc grism config use_bpet adj_sens do_flux

Parameters
~~~~~~~~~~

::

    grism:     name of the grism image

    config:   name of the aXe configuration file

    use_bpet: use of a BPET file

    adj_sens: adjust the sensitivity function for extended sources

    weights:  compute and apply optimal weights

    do_flux:  do flux calibration

    drzpath:  use AXE_DRIZZLE_PATH for IN/Output?

    in_af:    overwrite the default input Aperture File name

    opet:     overwrite the default input Object PET file name

    bpet:     overwrite the default input Background PET file name

    out_spc:  overwrite the default output SPC file name

::

    Example:
       pet2spc grism='test_grismn.fits' config='SLIM.conf.test.0' use_bpet='YES'

Output
~~~~~~

If  drzpath='NO':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].SPC.fits

If  drzpath='YES':

-  $AXE\_DRIZZLE\_PATH/[slitless filename]\_[ext number].SPC.fits

If  weights='YES':

-  $AXE\_OUTPUT/DRIZZLE\_PATH/[slitless filename]\_[ext
   number]\_opt.STP.fits

-  $AXE\_OUTPUT/DRIZZLE\_PATH/[slitless filename]\_[ext
   number]\_opt.WHT.fits

STAMPS
------

[STAMP] This task uses the content of a Pixel Extraction Table (see
Chapt.[PET]) to generate a FITS Stamp Image File (see Chapt.[STP])
containing stamp images of the BEAMs that were extracted. This task can
output various types of stamp images. In addition to the usual
**trace**-stamp images which display the beams as they appear on the
grism images, the **rectified** stamp images order the PET pixels in a
rectangular grid using the XI and DIST values of the pixels. For
**drizzled** stamp images the pixels are resampled onto a rectangular
grid with wavelength and trace distance as the axes and with a constant
dispersion. The first order beams are resampled to the dispersion
specified in the keyword DRZRESOLA or, as all other beams, to the
average dispersion in the PET pixels.

The stamp images allow a quick visual check on the extraction process.
Moreover the drizzled stamp images can be used as an input for
alternative 1D extractions with other iraf or IDL tools.

Usage
~~~~~

::

      stamps grism config sampling

Parameters
~~~~~~~~~~

::

    grism:    name of the grism image

    config:   name of the aXe configuration file

    sampling: the sampling type

    drzpath:  use AXE_DRIZZLE_PATH for IN/Output?

    in_af:    non standard OAF name

    in_PET:   non standard PET name

    out_stp:  non standard STP name

::

    Example:
       stamps grism='test_grismn.fits' config='SLIM.conf.test.0' sampling='rectified'

Output
~~~~~~

If  drzpath='NO':

-  $AXE\_OUTPUT\_PATH/[slitless filename]\_[ext number].STP.fits

If  drzpath='YES':

-  $AXE\_DRIZZLE\_PATH/[slitless filename]\_[ext number].STP.fits

DRZ2PET
-------

[DRZ2PET]

This task produces one object PET (background BPET if back=’YES’) from a
set of images created with AXEDRIZZLE. On this PET the task pet2spc can
then perform the extraction of the 1D spectra for the drizzled grism
images.

All the necessary input files (OAF/BAF, image list, modified
configuration file) are automatically created by the AXEDRIZZLE task.
The sequence of the images in the image list must match the sequence of
the beams in the OAF. Interactive changes to the image list and/or the
OAF are not recommended.

The 1D extraction of the 2D drizzled grism spectra is usually done
within axedrizzle by calls to the tasks drz2pet and pet2spc.

The task drz2pet also sets the pixel weights to reflect the different
signal-to-noise (S/N) ratios in each pixel. The S/N variations are
caused by the masking of bad and cosmic ray affected pixels and by the
partial coverage of objects on the border of grism object. The pixels
that will be co-added into a single resolution element in the 1D spectra
are weighted according to their relative exposure times. Moreover it is
in addition possible to compute and store optimal weights to enhance the
signal-to-noise (S/N) ratio in the 1D extracted spectra.

Usage
~~~~~

::

    drz2pet inlist config opt_extr back

Parameters
~~~~~~~~~~

::

    inlist:   ascii list which gives the name of the grism image to be processed 
              as the first item on each line.

    config:   name of the aXe configuration file(s).

    opt_extr: boolean to set the computation and storage
              of optimal weights

    back:     boolean to switch on/off the creation of background PETs made
              from drizzled background images.

    in_af:    non standard OAF name

    out_pet:  non standard PET name

    Example:
        drz2pet inlist='aXedrizzle_2.lis' conifgs='axedrizzle.conf' back='NO'

Output
~~~~~~

If  back='NO':

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_[ext number].PET.fits

If  back='YES':

-  $AXE\_DRIZZLE\_PATH/[drizzle root filename]\_[ext
   number].BCK.PET.fits

AXEGPS
------

[AXEGPS] This task reports the spectral properties of a single pixel.
The spectral properties for individual pixels can only be assigned with
respect to a reference point or reference beam. axegps lists:

-  the wavelength at pixel center

-  the dispersion at pixel center

-  the trace distance of the section point

-  the distance of the pixel center to the section point

-  the data value of the pixel

The task axegps works on the .OAF file. The corresponding OAF file and
the reference beam therein must therefore exist before axegps can give a
result.

For numerical reasons a solution can only be guaranteed within the
bounding box of the specified beam. The extraction width as specified
with the parameter extrfwhmin axecore (or mfwhm in gol2af) has an
influence on the bounding box. In the case that the desired information
for the pixel of interest is not given, a repetition of axecore (or
gol2af) with a larger value of drzfwhm (mfwhm) may enlarge the
bounding box sufficiently to get a result from axegps.

Even in case of failure, the corner points which define the bounding box
of the beam are listed in the output such that the user can understand
why the pixel information could not be computed.

Usage
~~~~~

::

    axegps grism config  beam_ref xval yval

Parameters
~~~~~~~~~~

::

    grism:    name of the grism image

    config:   name of aXe configuration file used to create the OAF

    beam_ref: the beam to define the spectral solutions

    xval:     the x-coordinate of the pixel

    yval:     the y-coordinate of the pixel

    Example:
        axegps grism="j8m822qhq_flt.fits" config="HUDF.HRC.conf"
               beam_ref="3A" xval=102 yval=588

Output
~~~~~~

All output is directly printed to the standard output.
