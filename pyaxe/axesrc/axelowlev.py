import os
import subprocess
import logging
# from hstaxe.config import __AXE_BINDIR as AXE_BINDIR
from hstaxe.axeerror import aXeError
from hstaxe.config import getOUTPUT

# make sure there is a logger
_log = logging.getLogger(__name__)

# define the good return
# value for the binaries
GOOD_RETURN_VALUE = 0


class TaskWrapper(object):
    """General class to execute the C-tasks"""
    def __init__(self, taskname="", tshort=""):
        """
        Parameters
        ----------
        taskname: str
            name of the C-executable
        tshort: str
            shor name for task
        """
        self.taskname = taskname
        self.tshort = tshort

        # initialize the command list
        self.command_list = []

        # save a name for stdout
        self.stdout = getOUTPUT(tshort+'.stdout')

        # save a name for stderr
        self.stderr = getOUTPUT(tshort+'.stderr')

        # put the command into the list
        # self.command_list.append("/".join([AXE_BINDIR, taskname]))
        self.command_list.append(taskname)

    def _cleanup(self):
        """The method deletes the files created for stdout and stderr.

        This is a usual cleaning procedure in case nothing bad happened.
        """
        # delete stdout/stderr
        if os.path.isfile(self.stdout):
            os.unlink(self.stdout)
        if os.path.isfile(self.stderr):
            os.unlink(self.stderr)

    def _report_all(self, silent=True):
        """Print stdout and stderr on the screen

        The method gives a feedback in case of problems. stdout and stderr
        are both listed onto the screen for a further interactive analysis.

        Parameters
        ----------
        silent: bool
            indicates silent/noisy runs
        """
        # if silent, save error info to file
        if silent:
            # dump the files with
            # stdout and stderr onto the screen
            _log.info("\nThere was a problem in the task: {0:s}"
                  .format(self.taskname))
            _log.info("The output of the task ({0:s}) is:\n".format(self.stdout))
            _log.info("---------------------------------------------------------"
                  "-----------------------")
            for line in open(self.stdout):
                _log.info(line.strip())
            _log.info("\n\nThe error report is of the task ({0:s}) is\n:"
                  .format(self.stdout))
            _log.info("----------------------------------------------------------"
                  "----------------------")
            for line in open(self.stderr):
                _log.info(line.strip())

        # report an error
        raise aXeError(f"An error occurred in the aXe task: {self.taskname}")

    def run(self, silent=True):
        """Run the wrapped task

        The method executes the associated C-executable. The return code from
        the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent: bool
            for silent mode

        Returns
        -------
        retcode: int
            the return code of the C-executable
        """
        # print("command list: ",self.command_list)
        if silent:
            # open stdout/stderr
            sout = open(self.stdout, 'w+')
            serr = open(self.stderr, 'w+')

            # execute the task
            retcode = subprocess.call(self.command_list,
                                      stdout=sout,
                                      stderr=serr)

            # close stdout/stderr
            sout.close()
            serr.close()

        else:

            # execute the task with the default stdout and
            # stderr, which is the system one
            retcode = subprocess.call(self.command_list)


    def runall(self, silent=True):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent: bool
            for silent mode

        Returns
        -------
        recode: int
            he return code of the C-executable
        """
        # run the executable
        retcode = self.run(silent=silent)
        # check whether the run was good
        if retcode == GOOD_RETURN_VALUE:
            # do the cleaning
            self._cleanup()
        else:
            self._report_all(silent)


class aXe_AF2PET(TaskWrapper):
    """Wrapper around the aXe_AF2PET task"""
    def __init__(self, grism, config, **params):
        """Generate a pixel extraction table (PET) from a grism image and an aperture file.

        Parameters
        ----------
        grism : str
            Name of grism/prism image.
        config : str
            Filename of grism/prism extraction configuration file.
        back : bool
            Whether to generate a PET for a background image using a BAF file,
            instead of a OAF file, and using a background image generated by
            the "backest" aXe background extraction task.
        out_pet : str
            Name to use for the output PET file instead of the default.
            Filename conventions for aXe are described below.
        in_af : str or None
            Name to use for the input aperture file with the stamp images
            instead of the default.petcont

        Description
        -----------
        This task uses the input slitless image together with an Object
        Aperture File (OAF) to generate an Object Pixel Extraction Table
        (OPET) for the input data. The same task should be used with the
        Background Estimate File and the same Object Aperture File (OAF) to
        generate a Background PET containing information about the spectral
        background (BPET).

        FILE NAMING CONVENTION aXe tasks use default names for input and
        output files based on the given ame of the "grism" image. For this
        task the default input OAF would be called <grism-rootname>_<science
        extension>.OAF and the output PET file would be
        <grism-rootname>_<science extension>.PET.fits. If used to create
        background PET the output would be <grism-rootname>_<science
        extension>.bck.PET.

        Example
        -------
        This will take the aperture file corresponding to the grism image
        test_grism_1.fits, using the standard naming scheme, and create an
        output pixel extraction table. The aXe configuration parameters will
        be read from the file SLIM.conf.test.0.

        axelowlev.aXe_AF2PET('test_grism_1.fits',
                             'SLIM.conf.test.0',
                              back=back,
                              in_af=in_af,
                              out_pet=out_pet

        Notes
        -----
        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name
        """
        super().__init__('aXe_AF2PET', 'af2pet')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for background flag
        if (('back' in params) and (params['back'])):
            # put the bck-flag to the list
            self.command_list.append('-bck')

        # check for the in-AF name
        if (('in_af' in params) and (params['in_af'])):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the out-PET name
        if (('out_pet' in params) and (params['out_pet'] is not None)):
            # put the name to the list
            self.command_list.append("-out_PET={0:s}"
                                     .format(params['out_pet']))


class aXe_GPS(TaskWrapper):
    """Wrapper around the aXe_GPS task"""
    def __init__(self, grism, config, beam_ref, xval, yval):
        """Compute the spectral properties of a single pixel.

        Parameters
        ----------
        grism: str
            Name of the grism image.
        config: str
            Name of aXe configuration file used to create the OAF.
        beam_ref: str
            The beam to define the spectral solutions.
        xval: float
            The x-coordinate of the pixel.
        yval: float
            The y-coordinate of the pixel.

        Description
        -----------
        This task computes the spectral properties of a single pixel in the
        context of the beam to define the local spectral properties and
        the configuration file which defines the global spectral properties
        on the grism image.

        The tasks lists:
          1. the wavelength at pixel center
          2. the dispersion at pixel center
          3. the trace distance of the section point
          4. the distance of the pixel center to the section point
          5. the data value of the pixel.

        `axegps` works on the .OAF file. The corresponding OAF file and the
        reference beam therein must therefore exist(run `sex2gol` and
        `gol2af`) before axegps can give a result.

        For numerical reasons a solution can only be guaranteed within the
        bounding box of the specified beam. The extraction width as specified
        with the parameter `-mfwhm=` in `gol2af` has an influence on the
        bounding box. In case that you do not get the desired information for
        the pixel of your interest you may repeat `gol2af` with a larger value
        of `-mfwhm=` to make the bounding box sufficiently large. The corner
        points which define the bounding box of the beam are listed in the
        output such that the user can understand why the pixel information
        could not be computed.

        Example
        -------
        Print the spectral properties of the pixel 102,588 on the grism image
        "j8m822qhq_flt.fits", using the beam "3A" as the reference beam.

        axelowlev.aXe_GPS(grism="j8m822qhq_flt.fits",
                                 config="HUDF.HRC.conf",
                                 beam_ref="3A",
                                 xval=102,
                                 yval=588)

        """
        super().__init__('aXe_GPS', 'axegps')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # put the beam name to the list
        self.command_list.append(beam_ref)

        # put the x-value to the list
        self.command_list.append(str(xval))

        # put the y-value to the list
        self.command_list.append(str(yval))


class aXe_BE(TaskWrapper):
    """Wrapper around the aXe_BE task"""
    def __init__(self, grism, config, **params):
        """This method is a simple initializer for the class.

        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_BE', 'backest')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'np'
        if (('np' in params) and (params['np'] is not None)):
            self.command_list.append('-np={0:s}'.format(str(params['np'])))

        # append the parameter 'interp'
        if (('interp' in params) and (params['interp'] is not None)):
            self.command_list.append("-interp={0:s}"
                                     .format(str(params['interp'])))

        # append the parameter 'niter_med'
        if (('niter_med' in params) and (params['niter_med'] is not None)):
            self.command_list.append('-niter_med={0:s}'
                                     .format(str(params['niter_med'])))

        # append the parameter 'niter_fit'
        if (('niter_fit' in params) and (params['niter_fit'] is not None)):
            self.command_list.append('-niter_fit={0:s}'
                                     .format(str(params['niter_fit'])))

        # append the parameter 'kappa'
        if (('kappa' in params) and (params['kappa'] is not None)):
            self.command_list.append('-kappa={0:s}'
                                     .format(str(params['kappa'])))

        # append the parameter 'smooth_length'
        if (params['smooth_length'] is not None):
            self.command_list.append('-smooth_length={0:s}'
                                     .format(str(params['smooth_length'])))

        # append the parameter 'smooth_length'
        if (params['smooth_fwhm'] is not None):
            self.command_list.append('-fwhm={0:s}'
                                     .format(str(params['smooth_fwhm'])))

        # append the parameter 'af_file'
        if (params['in_af']):
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # append the parameter 'out_bck'
        if (params['out_bck'] is not None):
            self.command_list.append('-out_BCK={0:s}'
                                     .format(params['out_bck']))

        # append the flag 'old_bck'
        if (('old_bck' in params) and (params['old_bck'])):
            self.command_list.append('-nor_flag')

        # append the flag 'mask'
        if (('mask' in params) and (params['mask'])):
            self.command_list.append('-msk')


class aXe_DRZ2PET(TaskWrapper):
    """Wrapper around the aXe_DRZ2PET task"""
    def __init__(self, inlist, config, **params):
        """Generate a PET froma list of 2D drizzled grism images.

        Parameters
        ----------
        inlist: str
            List to give input images.
        config: str
            Name of aXe configuration file(s) to be used.
        opt_extr : bool
            Compute and store optimal weights.
        back : bool
            Create a background PET from 2D drizzled background images?
        in_af : str
            Name to use for the input Aperture file instead of the default.
        out_pet : str
            Name to use for the output PET file (instead of the default).

        Description
        -----------
        The task `drz2pet` combines the 2D drizzled grism images given in
        the `inlist` into a single PET. The task needs an OAF/BAF whose
        name is derived from the name of inlist. The sequence of beams in
        the OAF/BAF must coincide with sequence of images in "inlist".
        Also the configuration file must reflect the properties of the
        ideal 2D drizzled, ideal grism spectra (trace description).

        The task `axedrizzle` produces all the necessary input files for
        `drz2pet`, and usually `drz2pet` is implicitly run within `axedrizzle`
        (`makespc=True`).

        With `opt_extr=False`the task updates the weight image and sets the
        weight proportional to the relative exposure time within the set of
        pixels that will be contracted to an individual resolution element.
        This is necessary to properly take into account the uneven distributed
        exposure time introduced by the masking of bad pixels/cosmic ray hist
        and objects which sometimes are at the border of a chip.

        With `opt_extr=True` the task computes optimal weights from the
        extensions `[WEIGH]` and `[VAR]` in the 2D drizzled grism images.
        The optimal weights are stored in the PET and used in subsequent
        calls to `pet2spc` on that data set.

        FILE NAMING CONVENTION
        The name of the PET created by this task is based on the name of the
        input list `inlist`. The task creates a file <inlist>.PET.FITS
        for `back=False`and <inlist>.PET.fits for `back=True`, respectively.

        Example
        -------
        Create an object PET from the images given in "test_aXe_2.lis".
        The names for the image list and the configuration file are typical
        filenames given in `axedrizzle`. The PET name is then
        "test_aXe_2.PET.fits".

            axelowlev.aXe_DRZ2PET(inlist="test_aXe_2.lis",
                                  config="test_aXe.conf",
                                  back=False)

        Notes
        -----
        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name.
        """
        super().__init__('aXe_DRZ2PET', 'drz2pet')

        # put the grism name to the list
        self.command_list.append(inlist)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'out_pet'
        if (('out_pet' in params) and (params['out_pet'] is not None)):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'in_af'
        if (('in_af' in params) and (params['in_af'])):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # append the flag 'bck'
        if (('back' in params) and (params['back'])):
            # put the bck-flag to the list
            self.command_list.append('-bck')

        # append the flag 'bck'
        if (('opt_extr' in params) and (params['opt_extr'])):
            # put the bck-flag to the list
            self.command_list.append('-opt_extr')


class aXe_DRZPREP(TaskWrapper):
    """Wrapper around the aXe_DRZPREP task"""
    def __init__(self, inlist, configs, **params):
        """Generate DPP's for a list of images and object lists

        Parameters
        ----------
        inlist : str
            List to give input images and object lists.
        configs : str
            Name of aXe configuration file(s) to be used.
        opt_extr : bool
            Also create extension for optimal extraction.
        back : bool
            Create also background DPP's?.

        Description
        -----------
        The task produces a set of Dirzzle PrePare (DPP) files for the
        grism images specified in `inlist`. The PET's and the BAF's for
        those grism images must have been prepared, e.g. by a previous run
        of `axecore`.

        A DPP is a multi extension fits file. For each first order beam in the
        corresponding PET the DPP has three extensions, [beam_??A], [err_??a]
        and [cont_??a]. The three extensions are the stamp images for
        the data values, the error values and the contamination values of the
        beam with the ID number "??". All extensions have important keywords
        which later are extracted and used in `axedrizzle`. Among those
        keywords are the coefficients that allow `axedrizzle` to combine the
        stamp images with the same ID from different DPP's onto a deep 2D grism
        image with identical dispersion and pixel scale in cross dispersion direction.

        To do that `drzprep` defines a common frame for all identical objects
        in the all grism images given in "inlist". `axedrizzle` will give
        consistent results ONLY if its input DPP's were generated within the
        same `drzprep` run. It is also not possible to run `drzprep`
        separately on the two science extensions of the WFC images. The usual
        dithers done for ACS image sets let an object sometimes fall on chip 1
        and in some instances on chip 2. To establish the common frame needed
        in `axedrizzle` requires also in this case that the DPP's were
        generated in only one `drzprep` run. For WFC images with its two
        chips you will specify two configuration files in a comma separated
        list. Of course the two PET's, one for each extension number, must
        both exist.

        With 'opt_extr=True' also two additional extensions necessary
        for optimal weighting are created and saved. Those extensions are named
        [mod_??A] and [var_??A] and contain the emission model as calculated
        in quantitative contamination and the theoretical inverse variance,
        respectively.

        Setting "back=True" also generates background
        DPP's, which are then based on the background PET's. With `axedrizzle`
        those background DPP's can be drizzled to deep 2D backgrounds for
        the individual object ID's.

        The format of the image list given as 'inlist' is:

        grim_image1 object_cat11[,object_cat12] [direct_image1] [dmag1]
        grim_image2 object_cat21[,object_cat22] [direct_image2] [dmag2]
        grim_image3 object_cat31[,object_cat32] [direct_image3] [dmag3]


        For the task `drzprep` only the first column with the name of the
        grism image is extracted and used. All further columns are neglected.
        It is possible to use the same 'inlist' as in `axeprep` and `axecore`
        (also later in axedrizzle).

        FILE NAMING CONVENTION
        aXe tasks use default names for input and output files based on the
        given name of the "grism" image. If the input grism image would be
        called <grism-rootname>.fits, the task creates the file

        <grism-rootname>_<science extension>.DPP.fits, and if "back='YES'"
        <grism-rootname>_<science extension>.BCK.DPP.fits

        All files are stored in the directory pointed by $AXE_OUTPUT_ROOT.

        Example
        -------
        For ACS HRC images,  reduce both, object DPP's and background DPP's
        from the images and object lists given in "axeprep.lis".

            axelowlev.aXe_DRZPREP(inlist="axeprep.lis",
                                  configs="HUDF.HRC.conf",
                                  back=True)

        For ACS WFC images, generate only the object DPP's for the grism
        images in "axeprep.lis". Work on the PET's specified by the extensions
        specified in ACS.WFC.CHIP1.conf and ACS.WFC.CHIP2.conf (usin the usual
        file naming convention. Create also the extensions for optimal
        weighting.

            aXe_DRZPREP(nlist="axeprep.lis",
                        opt_extr=True,
                        configs="ACS.WFC.CHIP1.conf,ACS.WFC.CHIP2.conf",
                        back=False)
        """
        super().__init__('aXe_DRZPREP', 'drzprep')

        # put the grism name to the list
        self.command_list.append(inlist)

        # put the config file name to the list
        self.command_list.append(configs)


        # store the 'back'-flag
        if (('back' in params) and (params['back'])):
            # put the opt_extr-flag to the list
            self.bck = True
            self.command_list.append('-opt_extr')
        else:
            self.bck = False
            

    def runall(self, silent=True):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent : bool
            or silent mode

        Returns
        -------
        recode : int
            the return code of the C-executable
        """
        # run the method of the super class
        super().runall(silent=silent)

        # check for the background flag
        if self.bck:
            # put the bck-flag to the list
            self.command_list.append('-bck')

            # run the method of the super class
            super().runall(silent=silent)


class aXe_GOL2AF(TaskWrapper):
    """Wrapper around the aXe_GOL2AF task"""
    def __init__(self, grism, config,
                 slitless_geom=True,
                 orient=True,
                 back=False, **params):
        """Generates an aperture file from a grism object list.

        Parameters
        ----------
        grism : str
            Grism or Prism image filename.
        config : str
            aXe configuration file term

        mfwhm : int
            Extraction width multiplicative factor.

        back : bool
            Whether a BAF will be generated rather than an OAF.

        orient : bool
            Whether to use tilted extracted. This is the default. When set to
            False only vertical extraction (along columns) is performed.

        slitless_geom : bool
             Whether to use an extraction orientation which is optimized for
             slitless spectroscopy.

        exclude : bool
            Do not even list faint objects (mag>max_mag) in the output.

        lambda_mark : float
            the wavelength at which to apply the cutoff magnitudes MMAG_EXTRACT
            and MMAG_MARK

        dmag : float or None
            number to add to the MMAG_EXTRACT and MMAG_MARK values given in the
            configuration file

        out_af : str
            overwrites the default output OAF or BAF filename

        in_gol : str
            Use this name in preference to the default for the input object
            list catalog file name.

        Description
        -----------
        This task generates an Aperture File (AF) using an input Grism Object
        List (GOL) and a valid Main Configuration File which specifies the
        length, wavelength calibration and global offsets between direct and
        slitless images. Any offset between the two images produced by
        dithering of the telescope between the time the direct and the
        slitless images were taken is assumed to  be properly accounted for in
        the WCS of both images. The width of the BEAMs  is set up to be a
        given number times the full width at half maximum listed  in the GOL.
        If the exclude parameter is also set these faint objects will not
        appear in the output GOL at all.

        This task can be used to generate both an Object Aperture File and a
        Background Aperture File. These files have a similar format, but it is
        often  desirable to use different Aperture Files in both cases since
        the first is  used to extract counts from pixels which are known to
        contain flux  from the source while the second can be thought of as a
        definition of a zone  to avoid in the slitless image when computing
        the background level. In  practice, a larger extraction width
        multiplier should be used when generating  the Background Aperture
        File, so that the background estimation is not contaminated by flux
        from the object.

        The task defines how exactly the spectra will be extracted by the
        tasks pet2spc or axedrizzle in later stages of the aXe reduction. With
        'orient="yes"' the cross-dispersion direction follows the angle
        "THETA_IMAGE" specified for each image in the Input Object List. The
        parameter 'slitless_geom="yes"' automatically adjusts the
        cross-dispersion angle such that the extracted signal becomes maximal.

        With 'orient="no"' the cross-dispersion direction is, for every
        object, perpendicular to the trace.

        Positive values for 'mfwhm' are used to define the extraction width as
        'mfwhm'*'object_width' with 'object_with' either A_IMAGE/B_IMAGE
        (orient='yes') or their projection perpendicular to the trace
        (orient='no'). Negative values specify the extration width in pixels.
        A detailed description on the interplay between those three parameters
        is given in the aXe manual.

        When using quantitative contamination scheme such as the gaussian
        contamination there may exist information on the source brightness at
        different wavelengths. The parameter 'lambda_mark' specifies the
        wavelength where the cutoff magnitudes defined in the configuration
        file should be applied. The algorithm selects the wavlength closest to
        the wavelength given in 'lambda_mark' and applies the selection there.
        An interpolation between the source brightness at different
        wavelengths is not performed.

        FILE NAMING CONVENTION

        aXe tasks use default names for input and output files based on the
        given name of the "grism" image. For this task the default input GOL
        would be called <grism-rootname>_<science extension>.GOL and the
        output OAF file would be <grism-rootname>_<science extension>.OAF. If
        used to create a background aperture file the output would be
        <grism-rootname>_<science extension>.BAF.

        Example
        -------
        Create an object aperture file or the objects in the grism/prism image
        'test_grism_1.fits'. The cross-dispersion direction is along the angle
        "THETA_IMAGE" if this forms not a too small angle with the trace, and
        perpendicular to "THETA_IMAGE" otherwise. In the former case the
        extraction width is 2.0*A_IMAGE, in the latter case 2.0*B_IMAGE.

            axelowlev.aXe_GOL2AF("test_grism_1.fits",
                                 fwhm=2,
                                 orient=True,
                                 slitless_geom=False,
                                 )

        Create an object aperture file or the objects in the grism/prism
        image 'test_grism_1.fits'. The cross-dispersion direction is
        perpendicuar to the object trace, and the extraction width is fixed at
        3 pixels for all objects.

            axelowlev.aXe_GOL2AF("test_grism_1.fits",
                                 fwhm=-3,
                                 orient=False)
        """
        super().__init__('aXe_GOL2AF', 'gol2af')

        # add the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'mfwhm'
        if (('mfwhm' in params) and (params['mfwhm'] is not None)):
            # put the name to the list
            self.command_list.append('-mfwhm={0:s}'
                                     .format(str(params['mfwhm'])))

        # append the parameter 'dmag'
        if (('dmag' in params) and (params['dmag'] is not None)):
            # put the name to the list
            self.command_list.append('-dmag={0:s}'.format(str(params['dmag'])))

        # append the parameter 'lambda_mark'
        if (('lambda_mark' in params) and (params['lambda_mark'] is not None)):
            # put the name to the list
            self.command_list.append('-lambda_mark={0:s}'
                                     .format(str(params['lambda_mark'])))
        else:
            self.command_list.append('-lambda_mark=800.0')

        # append the parameter 'out_pet'
        if (('out_pet' in params) and (params['out_pet'] is not None)):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'in_af'
        if (('out_af' in params) and (params['out_af'])):
            # put the name to the list
            self.command_list.append('-out_AF={0:s}'.format(params['out_af']))

        # append the parameter 'in_gol'
        if (('in_gol' in params) and (params['in_gol'] is not None)):
            # put the name to the list
            self.command_list.append('-in_GOL={0:s}'.format(params['in_gol']))

        # append the flag 'slitless_geom'
        if (('slitless_geom' in params) and (params['slitless_geom'])):
            # put the slitless_Geom-flag to the list
            self.command_list.append('-slitless_geom=1')
        else:
            self.command_list.append('-slitless_geom=0')

        # append the flag 'orient'
        if (('orient' in params) and (params['orient'])):
            # put the orient-flag to the list
            self.command_list.append('-orient=1')
        else:
            self.command_list.append('-orient=0')

        # append the flag 'exclude'
        if (('exclude' in params) and (params['exclude'])):
            # put the exclude_faint-flag to the list
            self.command_list.append('-exclude_faint')

        #  append the flag 'bck'
        if (('back' in params) and (params['back'])):
            # put the bck-flag to the list
            self.command_list.append('-bck')
        _log.info("Command list: {}".format(self.command_list))


class aXe_INTPIXCORR(TaskWrapper):
    """Wrapper around the aXe_INTPIXCORR task"""
    def __init__(self, grism, config, **params):
        """
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_INTPIXCORR', 'ipixcorr')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_OAF' name
        if (params['in_OAF'] is not None):
            # put the name to the list
            self.command_list.append('-in_OAF={0:s}'
                                     .format(str(params['in_oaf'])))

        # check for the 'in_SPC' name
        if (params['in_SPC'] is not None):
            # put the name to the list
            self.command_list.append('-in_SPC={0:s}'
                                     .format(str(params['in_spc'])))

        # check for the 'out_SPC' name
        if (params['out_SPC'] is not None):
            # put the name to the list
            self.command_list.append('-out_SPC={0:s}'
                                     .format(params['out_spc']))

        # check for the 'max_ext' value
        if (params['max_ext'] is not None):
            # put the value to the list
            self.command_list.append('-max_ext={0:s}'
                                     .format(str(params['max_ext'])))
        else:
            # put the default to the list
            self.command_list.append('-max_ext=0.0')

        # append the flag 'ip_corr'
        if (params['ip_corr']):
            # put the ip_corr-flag to the list
            self.command_list.append('-ipixcorr')

        # append the flag 'nl_corr'
        if (params['nl_corr']):
            # put the nl_corr-flag to the list
            self.command_list.append('-nlincorr')


class aXe_NICBACK(TaskWrapper):
    """Wrapper"""
    def __init__(self, grism, config, master_bck, back_ped=None):
        """
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_NICBACK', 'nicback')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # put the master background
        self.command_list.append(master_bck)

        if ((back_ped is not None) and (len(back_ped) > 1)):
            # put the master background
            self.command_list.append(back_ped)


class aXe_SCALEBCK(TaskWrapper):
    """
    Wrapper around the aXe_PET2SPC task
    """
    def __init__(self, grism, mask, config, master_sky, to_master=False,
                 make_plis=False):
        """This method is a simple initializer for the class.

        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Inputs
        ------
        grism : string
            Name of the dispersed image
        mask : string
            Name of the mask image
        config : string
            Name of the aXe configuration file
        master_sky : string
            Name of the master sky image
        to_master : boolean
            Scale grism to master
        make_plis : boolean
            Generate the pixel list

        Outputs
        -------
        """
        super().__init__('aXe_SCALEBCK', 'scalebck')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the mask name to the list
        self.command_list.append(mask)

        # put the config file name to the list
        self.command_list.append(config)

        # put the master background
        self.command_list.append(master_sky)

        if to_master:
            # add the flag for scaling to the master image
            self.command_list.append('-toMaster')

        if make_plis:
            # add flag to generate pixel list
            self.command_list.append('-make_plis')


class aXe_PET2SPC(TaskWrapper):
    """Wrapper around the aXe_PET2SPC task"""
    def __init__(self, grism, config,
                 usr_bpet=True,
                 adj_sens=True,
                 weights=False,
                 do_flux=True,
                 drzpath="",
                 in_af="",
                 opet="",
                 bpet="",
                 use_bpet=False,
                 out_spc=""):
        """Bin contents of a Pixel Extraction Table into 1D spectra.

        Parameters
        ----------
        grism : str
            Name of grism/prism image.
        config : str
            Filename of grism/prism extraction configuration file.
        usr_bpet : bool
            Whether to use a background PET file.
        adj_sens : bool
            Adjust sensitivity for extended objects.
        weights : bool
            Compute and use optimal weights in 1D extraction.
        do_flux : bool
            Whether to perform flux calibration.
        out_spc : str
            Name to use for the output SPC file instead of the default.
        drz_path : str
            Use the directory indicated by the system variable
            $AXE_DRIZZLE_PATH to located the input PET and write out the
            output spectra (SPC).
        in_af : str
            Name to use for the input Aperture file instead of the default.
        opet : str
            Name to use for the input object pixel extraction table file
            instead of the default.
        bpet : str
            Name to use for the input background pixel extraction table file
            instead of the default.
        use_bpet : bool
            Put the ip_corr-flag to the list
        out_spc : str
            Name to use for the output file with the spectra instead of the
            default).

        Description
        -----------
        This task is used to transform the content of an Object Pixel
        Extraction Table (PET) into a set of 1D binned spectra in a Extracted
        Spectra File (SPC). The  binning process is explained in more detail
        in the aXe manual.

        In case that the object PET was created from a grism image with
        substantial sky background, "use_bpet='yes'" extracts a background
        spectrum from the corresponding background PET and subtracts the
        background from the object spectrum.

        The sensitivity file to convert from e/s to flux units can be adjusted
        for extended objects ("adj_sens='yes'").

        Also optimal weighting (weights='yes') can be used to enhance the
        signal-to-noise ratio. Using quantitative contamination is a
        prerequisite to be able to computing optimal weights.

        FILE NAMING CONVENTION

        aXe tasks use default names for input and output files based on the
        given name of the "grism" image. For this task the default input PET
        would be called <grism-rootname>_<science extension>.PET.fits and the
        output SPC file would be <grism-rootname>_<science
        extension>.SPC.fits.

        Example
        -------
        Extract 1D spectra from the PET
        which was derived from image 'test_grism_1.fits'. To give meaningful
        (without background) result, the sky background must have been removed
        from the grism/prism image before creating the PET.

            axelowlev.aXe_PET2SPC("test_grism_1.fits",
                                  "SLIM.conf.test.0")

        Same as above, but use the background PET to extract and subtract a
        background spectrum from the object spectrum. Apply optimal weights in
        the 1D extraction.

            axelowlev.aXe_PET2SPC("test_grism_1.fits",
                                  "SLIM.conf.test.0",
                                  use_bpet=True,
                                  weights=True)

        """
        super().__init__('aXe_PET2SPC', 'pet2spc')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_AF' name
        if in_af:
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(in_af))

        # check for the 'OPET' name
        if opet:
            # put the name to the list
            self.command_list.append('-OPET={0:s}'.format(opet))

        # check for the 'BPET' name
        if bpet:
            # put the name to the list
            self.command_list.append('-BPET={0:s}'.format(bpet))

        # check for the 'out_SPC' name
        if out_spc:
            # put the name to the list
            self.command_list.append('-out_SPC={0:s}'
                                     .format(out_spc))

        # append the flag 'drzpath'
        if drzpath:
            # put the ip_corr-flag to the list
            self.command_list.append('-drz')

        # append the flag 'noBPET'
        if not use_bpet:
            # put the ip_corr-flag to the list
            self.command_list.append('-noBPET')

        # append the flag 'opt_weights'
        if weights:
            # put the ip_corr-flag to the list
            self.command_list.append('-opt_weights')

        # append the flag 'noflux'
        if not do_flux:
            # put the ip_corr-flag to the list
            self.command_list.append('-noflux')

        # append the flag 'smooth_conv'
        if adj_sens:
            # put the ip_corr-flag to the list
            self.command_list.append('-smooth_conv')


class aXe_PETCONT(TaskWrapper):
    """Wrapper around the aXe_PETCONT task"""
    def __init__(self, grism, config,
                 cont_model="gauss",
                 model_scale=None,
                 spec_models="",
                 object_models="",
                 inter_type="linear",
                 lambda_psf=None,
                 cont_map=False,
                 no_pet=False,
                 silent=False,
                 in_af=""):
        """Check whether pixels in a PET are members of more than one beam

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            Name of grism/prism image.
        config: str
            Filename of grism/prism configuration file.
        cont_model : str {"gauss", "fluxcube", "geometric"}
            The contamination model to be applied.
        model_scale : float
            Scale factor for the gaussian emission model
        spec_models : str
            Spectral models to use.
        object_models : str
            Object models to use.
        inter_type : str {"linear","polynomial","spline"}
            Interpolation type for the flux values.
        lambda_psf : float
            Wavelength in [nm] A_IMAGE and B_IMAGE was measured at (gauss).
        cont_map : bool
            Whether to write the contamination map into a FITS file.
        in_af : str
            Name to use for the input aperture file (instead of the default).

        Description
        -----------
        This task determines for each pixel in a PET the degree of
        contamination experienced from other sources. There are two
        fundamentally different ways to describe the contamination:

        The quantitative contamination uses an emission model and computes the
        contaminating flux according to the emission model (cont_model='gauss'
        or 'fluxcube').

        The geometrical contamination (cont_model=geometric) just records the
        number of beams a PET pixel is used in.

        Depending on the contamination model used the requirements for the
        Input Object Lists and mandatory input files are different. Please
        read the aXe manual for details.

        FILE NAMING CONVENTION
        aXe tasks use default names for input and output files based on the
        given name of the "grism" image. For this task the default input PET
        would be called <grism-rootname>_<science extension>.PET.fits and the
        input  aperture file would be <grism-rootname>_<science
        extension>.OAF.

        Example
        -------
        Compute the contamination for the image 'test_grism_1.fits' using
        the gaussian emission model:

            axelowlev.aXe_PETCONT("test_grism_1.fits",
                                  "SLIM.conf.test.0",
                                cont_model='gauss')
        """
        super().__init__('aXe_PETCONT', 'petcont')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_AF' name
        if in_af:
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(in_af))

        # check for the 'specmodels' name
        if spec_models:
            # put the name to the list
            self.command_list.append('-model_spectra={0:s}'
                                     .format(spec_models))

        # check for the 'objectmodels' name
        if object_models:
            # put the name to the list
            self.command_list.append('-model_images={0:s}'
                                     .format(object_models))

        # append the flag 'cont_model'
        if (cont_model == 'gauss'):
            # put the gauss-flag to the list
            self.command_list.append('-cont_model=1')
        elif (cont_model == 'direct'):
            # put the according number to the list
            self.command_list.append('-cont_model=2')
        elif (cont_model == 'fluxcube'):
            # put the according number to the list
            self.command_list.append('-cont_model=3')
        elif (cont_model == 'geometric'):
            # put the according number to the list
            self.command_list.append('-cont_model=4')

        # append the flag 'inter_type'
        if (inter_type == 'linear'):
            # put the according number to the list
            self.command_list.append('-inter_type=1')
        elif (inter_type == 'polynomial'):
            # put the according number to the list
            self.command_list.append('-inter_type=2')
        elif (inter_type == 'spline'):
            # put the according number to the list
            self.command_list.append('-inter_type=3')

        # append the parameter 'model_scale'
        if (model_scale is not None):
            # put the flag to the list
            self.command_list.append('-model_scale={0:s}'
                                     .format(str(params['model_scale'])))

        # append the parameter 'lambda_psf'
        if (lambda_psf is not None):
            # put the number to the list
            self.command_list.append('-lambda_psf={0:s}'
                                     .format(str(lambda_psf)))
        else:
            self.command_list.append('-lambda_psf=800.0')

        # append the flag 'cont_map'
        if cont_map:
            # put the cont_map-flag to the list
            self.command_list.append('-cont_map')

        # append the flag 'no_pet' to indicate
        # that no PET exists
        if no_pet:
            # append the no-PET flagg
            self.command_list.append('-noPET')


class aXe_PETFF(TaskWrapper):
    """Wrapper around the aXe_PETFF task"""
    def __init__(self, grism, config,
                 back=False,
                 ffname=""):
        """Flat field the contents of a Pixel Extraction Table (PET).

        Parameters
        ----------
        grism: str
            Name of grism/prism image.
        config: str
            Filename of grism/prism extraction configuration file.
        back : bool
            Specify whether or not to apply the flat-field (FF) to the
            background Pixel Extraction Table (BPET).
        ffname : str
            Name to use for the output flat field cube file instead of the
            default

        Description
        -----------
        Apply a wavelength dependent flat-fielding to the COUNT and
        ERROR parts of a Pixel Extraction Table (PET)
        The wavelength of a pixel is used in conjunction with a flat-fielding
        data cube containing the coefficents of a polynomial which can be used
        to compute at each pixel (x,y):

            FF(x,y,lambda) = a_0(x,y) + a_1(x,y)*lambda + .. +a_i * lambda^i

        The coefficients a_0(x,y) are stored in the first data extension of the
        flat-field cube, a_1(x,y) in the second, etc... The name of the
        flat-field cube is read from the aXe configuration file.

        Note that that flat-fielding should be applied to both the OPET and
        the BPET  separately so that, when subtracted, the information
        contained in the BPET  has been flat-fielded with the same flat-field
        as the OPET (which contains  counts of objects plus counts from the
        background).

        FILE NAMING CONVENTION
        aXe tasks use default names for input and output files based on the
        given name of the "grism" image. For this task the default input PET
        would be called <grism-rootname>_<science extension>.PET.fits and will
        be flat-fielded "in place" using the flat field cube defined in the
        main configuration file.

        Example
        -------
        Perform in-place flat-fielding of the PET associated with the grism
        data file test_grism_1.fits. Use configuration parameters from the file
        SLIM.conf.test.0.

            axelowlev.aXe_PETFF("test_grism_1.fits",
                                "SLIM.conf.test.0")
        Notes
        -----
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        """
        super().__init__('aXe_PETFF', 'petff')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'FFNAME' name
        if ffname:
            # put the name to the list
            self.command_list.append('-FFNAME={0:s}'.format(params['ffname']))

        # append the flag 'bck'
        if back:
            # put the according flag to the list
            self.command_list.append('-bck')


class aXe_PETIPC(TaskWrapper):
    """Wrapper around the aXe_PETIPC task"""
    def __init__(self, grism, config, **params):
        """

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_PETIPC', 'petipc')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_OAF' name
        if (('in_OAF' in params) and (params['in_OAF'] is not None)):
            # put the name to the list
            self.command_list.append('-in_OAF={0:s}'.format(params['in_oaf']))

        # check for the 'in_PET' name
        if (('in_PET' in params) and (params['in_PET'] is not None)):
            # put the name to the list
            self.command_list.append('-in_PET={0:s}'.format(params['in_pet']))

        # check for the 'out_PET' name
        if (('out_PET' in params) and (params['out_PET'] is not None)):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'max_ext'
        if (('max_ext' in params) and (params['max_ext'] is not None)):
            # put the number to the list
            self.command_list.append('-max_ext={0:s}'
                                     .format(str(params['max_ext'])))

        # append the flag 'bck'
        if (params['back']):
            # put the according flag to the list
            self.command_list.append('-bck')

        # append the flag 'origname'
        if (('orig_name' in params) and (params['orig_name'])):
            # put the according flag to the list
            self.command_list.append('-origname')


class aXe_STAMPS(TaskWrapper):
    """Wrapper around the aXe_STAMPS task"""
    def __init__(self, grism, config, **params):
        """Generate stamp images from a pixel extraction table.

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism : str
            Name of grism/prism image
        config : str
            Name of grism/prism extraction configuration file
        sampling : str {"drizzle", "rectified", "trace"}
            Which sampling should be used generating the stamp?
            "drizzle" is the default
        drzpat : bool
            Use the directory indicated by the system variable $AXE_DRIZZLE_PATH
            to located the input PET and write the output stamp image (STP).
        in_af : str
            Name to use for the input AF file with the stamp images instead of the default.
        in_pet : str
            Name to use for the input PET file with the stamp images instead of the default.
        out_STP : str
            Name to use for the output file with the stamp images instead of the default.
        silent : bool
            Give feedback incase of problems


        Returns
        -------
        Each stamp is saved as a seperate extention inside a FITS file that is placed
        in the AXE_OUTPUT directory.


        Description
        -----------
        This task uses the content of a Pixel Extraction Table (PET) file to
        generate a Stamp Image File (STP) containing stamp images of the BEAMs that
        were extracted. Different kinds of stamp images can be created:

        sampling="drizzle"  : fully resampled stamp images in the
                              trace length - cross dispersion plane
                              with the trace axis parallel to the
                              x-axis
        sampling="rectified": not resampled stamp images in the
                              trace length - cross dispersion plane
                              with the trace axis parallel to the
                              x-axis
        sampling="trace":     trace axis and cross dispersion axis
                              as in the flt-image

        Users can quickly check that the extraction angle used was the
        appropriate one. Because this uses the content of a PET and not the
        original input slitless image, it offers a good check of what pixels were used
        during the extraction process.

        Using drizzle sampling the stamp images could even be used to combine
        spectra of the same object from different images (similar to axedrizzle).
        To aid this usage, each stamp image (extension '[BEAM_<object id><beam id>]')
        is accompanied by the corresponding weight image
        (extension '[WEIG_<object id><beam id>]').

        FILE NAMING CONVENTION
        aXe tasks use default names for input and output files based on the given
        name of the "grism" image. For this task the default input PET would be
        called <grism-rootname>_<science extension>.PET.fits and the output stamps
        image file would be <grism-rootname>_<science extension>.STP.fits.

        Example
        -------
        Create a stamp image for all beams in the corresponding PET file
        using the 'drizzle'-sampling:

        axelowlev.stamps(grism=fname,
                         config='G141.F140W.V4.31.conf',
                         sampling='drizzle',
                         drzpath=False,
                         in_af='',
                         in_pet=None,
                         out_stp=None,
                         silent=False)

            > stamps test_grism_1.fits SLIM.conf.test.0 sampling='drizzle'
        """
        super().__init__('aXe_STAMPS', 'stamps')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'in_af'
        if (params['in_af']):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the 'in_PET' name
        if (params['in_pet'] is not None):
            # put the name to the list
            self.command_list.append('-in_PET={0:s}'.format(params['in_pet']))

        # check for the 'out_STP' name
        if (params['out_stp'] is not None):
            # put the name to the list
            self.command_list.append('-out_STP={0:s}'
                                     .format(params['out_stp']))

        # append the flag 'drz'
        if (params['sampling'] is 'drizzle'):
            # put the flagg to the list
            self.command_list.append('-drzstamp')
        elif (params['sampling'] is 'rectified'):
            # put the flagg to the list
            self.command_list.append('-rectified')

        # append the flag 'drz'
        if (params['drzpath']):
            # put the according flag to the list
            self.command_list.append('-drz')
        _log.info(self.command_list)


class aXe_TRACEFIT(TaskWrapper):
    """Wrapper around the aXe_TRACEFIT task"""
    def __init__(self, grism, config, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_TRACEFIT', 'tracefit')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'in_af'
        if (('in_af' in params) and (params['in_af'])):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))


class aXe_FILET(TaskWrapper):
    """Wrapper around the aXe_STAMPS task"""
    def __init__(self, dppfile, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        @param dppfile: name of the dpp-file
        @type dppfile: string
        @param **params: all other parameters
        @type **params: dictionary
        """
        super().__init__('aXe_FILET', 'filet')
        _log.info("Called FILET")

        # put the grism name to the list
        self.command_list.append(dppfile)

        # append the 'opt_extr' flag
        if (('opt_extr' in params) and (params['opt_extr'])):
            self.command_list.append('-opt_extr')

            # append the parameter 'in_af'
        if (('drztmp' in params) and (params['drztmp'] is not None)):
            # put the name to the list
            self.command_list.append('-drztmp={0:s}'.format(params['drztmp']))


class aXe_DIRIMAGE(TaskWrapper):
    """Wrapper around the aXe_DIRIMAGE task"""
    def __init__(self, dirname, config, tpass_direct, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        dirname: str
            name of the direct image
        configfile: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        super().__init__('aXe_DIRIMAGE', 'dirimage')

        # put the direct image name to the list
        self.command_list.append(dirname)

        # put the aXe configuration file name to the list
        self.command_list.append(config)

        # put the total passband file name to the list
        self.command_list.append(tpass_direct)

        # check whether model spectra are given
        # append the file name to the list
        if (params['model_spectra'] is not None):
            # self.command_list.append('-model_spectra='+str(model_spectra))
            self.command_list.append('-model_spectra={0:s}'
                                     .format(params['model_spectra']))

        # check whether model images are given
        # append the file name to the list
        if (params['model_images'] is not None):
            self.command_list.append('-model_images={0:s}'
                                     .format(params['model_images']))

        # check whether model images are given
        # append the file name to the list
        if (params['tel_area'] is not None):
            self.command_list.append('-tel_area={0:f}'
                                     .format(params['tel_area']))

        # append the model scale, give default if necessary
        if (params['model_scale'] is not None):
            self.command_list.append('-model_scale={0:f}'
                                     .format(params['model_scale']))
        else:
            self.command_list.append('-model_scale=5.0')
