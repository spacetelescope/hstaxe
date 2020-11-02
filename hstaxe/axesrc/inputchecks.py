import os
import math
import logging
from astropy.io import fits
from astropy.table import Table
from astropy.io.registry import IORegistryError

from . import axeinputs
from . import configfile
# from . import axeiol

from hstaxe.axeerror import aXeError, aXeSIMError
from hstaxe import config as config_util

# make sure there is a logger
_log = logging.getLogger(__name__)

class InputChecker:
    def __init__(self, taskname, inlist=None, configs=None, backims=None):
        # store the parameters
        self.taskname = taskname

        # check whether an input image list exists
        if inlist is not None:
            # make sure the Input Image List does exist
            if not os.path.isfile(inlist):
                msg = ("{0:s}: The Input Image List {1:s} does not exist!"
                       .format(self.taskname, inlist))
                raise aXeError(msg)

            # create a table with the basic aXe inputs
            self.axe_inputs = axeinputs.aXeInput(inlist, configs, backims)

    def _is_prism_data(self):
        # define the default
        is_prism = False

        # make sure there are grism images
        for row in self.axe_inputs:
            one_grisim = config_util.getDATA(row['grisim'])
            # open the fits
            one_fits = fits.open(one_grisim, 'readonly')

            # read the keyword 'FILTER1'
            try:
                filter1 = one_fits[0].header['FILTER1']
            except KeyError:               
                filter1 = None

            # read the keyword 'FILTER2'
            try:
                filter2 = one_fits[0].header['FILTER2']
            except KeyError:
                filter2 = None

            # check whether it is prism data
            if ((filter1 and 'PR' in filter1) or
                    (filter2 and 'PR' in filter2)):
                is_prism = True

            # close the fits
            one_fits.close()

        # return the index
        return is_prism

    def _force_dirim(self):
        # go over all inputs
        for one_input in self.axe_inputs:

            # check whether there is a direct image
            if one_input['DIRIM'] is None:
                # error and out
                err_msg = ("{0:s}: The grism image: {1:s} does NOT have an "
                           "associated direct image!"
                           .format(self.taskname,
                                   config_util.getDATA(one_input['grisim'])))
                raise aXeError(err_msg)

    def _check_fluxcubes(self):
        # go over all inputs
        for one_input in self.axe_inputs:

            # load the config file and get the extension information
            conf = configfile.ConfigFile(config_util.getCONF(one_input['config']))
            ext_info = config_util.get_ext_info(config_util.getDATA(one_input['grisim']), conf)

            # derive the aXe names
            axe_names = config_util.get_axe_names(one_input['grisim'], ext_info)

            # check the fluxcube
            if not os.path.isfile(config_util.getDATA(axe_names['FLX'])):
                # error and out
                err_msg = ("{0:s}: The fluxcube file: {1:s} does not exist!"
                           .format(self.taskname,
                                   config_util.getDATA(axe_names['FLX'])))
                raise aXeError(err_msg)

    def _check_global_backsub(self):
        """Check for global background subtraction"""
        # go over all inputs
        for one_input in self.axe_inputs:

            # load the config file and get the extension information
            conf = configfile.ConfigFile(config_util.getCONF(one_input['config']))
            ext_info = config_util.get_ext_info(config_util.getDATA(one_input['grisim']), conf)

            # open the fits image
            gri_fits = fits.open(config_util.getDATA(one_input['grisim']), 'readonly')

            # go to the correct header
            act_header = gri_fits[ext_info['fits_ext']].header

            # make sure a sky background value is set
            if 'SKY_CPS' in act_header and act_header['SKY_CPS'] >= 0.0:
                # close the fits
                gri_fits.close()
            else:
                # close fits, complain and out
                gri_fits.close()
                err_msg = ("{0:s}: The grism image: \n{1:s}\nhas no keyword "
                           "SKY_CPS>=0.0 in the extension {2:d}. This means "
                           "it had NO global\nsky subtraction, which is "
                           "required for the CRR version of aXedrizzle!"
                           .format(self.taskname,
                                   config_util.getDATA(one_input['grisim']),
                                   ext_info['fits_ext']))
                raise aXeError(err_msg)

    def _check_dpps(self, back=False):
        # go over all inputs
        for one_input in self.axe_inputs:
            # load the config file and get the extension information
            conf = configfile.ConfigFile(config_util.getCONF(one_input['config']))
            ext_info = config_util.get_ext_info(config_util.getDATA(one_input['grisim']), conf)

            # derive the aXe names
            axe_names = config_util.get_axe_names(one_input['grisim'], ext_info)
            # check the DPP file
            if not os.path.isfile(config_util.getOUTPUT(axe_names['DPP'])):
                # error and out
                err_msg = (f"{self.taskname}: The DPP file: {config_util.getOUTPUT(axe_names['DPP'])}"
                           f" does not exist!")
                           
                raise aXeError(err_msg)

            # check for the background DPP file
            if back and not os.path.isfile(config_util.getOUTPUT(axe_names['BCK_DPP'])):
                # error and out
                err_msg = (f"{self.taskname}: The DPP file: {config_util.getOUTPUT(axe_names['BCK_DPP'])}"
                           f" does not exist!")
                raise aXeError(err_msg)

    def check_axeprep(self, backgr, backims):
        """Comprises all file and file format checks for AXEPREP"""

        # check for background subtraction
        if backgr:
            if not backims:
                err_msg = ("{0:s}: A background image must be given for the "
                           "background subtraction!".format(self.taskname))
                raise aXeError(err_msg)

    def _check_IOL(self):
        """Check the presence of all grism images."""
        IOL_list = []

        # go over all inputs
        for one_input in self.axe_inputs:
            if one_input['objcat'] in IOL_list:
                continue

            # check the prism image
            if not os.path.isfile(config_util.getDATA(one_input['objcat'])):
                raise aXeError("{0:s}: The direct image: '{1:s}' does not "
                               "exist!".format(self.taskname, config_util.getDATA(one_input['objcat'])))
            # load the IOL to check its format
            try:
                _log.info("reading {}".format(config_util.getDATA(one_input['objcat'])))
                Table.read(config_util.getDATA(one_input['objcat']), format='ascii.sextractor')
            except IORegistryError:
                raise aXeError("{0:s} is not in a known format".format(one_input['objcat']))
            # put the IOL to the list
            IOL_list.append(one_input['objcat'])

    def check_axecore(self, back, extrfwhm, drzfwhm, backfwhm, orient,
                      slitless_geom, np, interp, cont_model,
                      weights, sampling):
        """
        Comprises all file and file format checks for AXECORE
        """

        # check the IOL's
        self._check_IOL()

        # check the fluxcubes, if necessary
        if cont_model.lower() is 'fluxcube':
            self._check_fluxcubes()

        # check whether it is prism data
        if self._is_prism_data():
            #
            # NOTE: these checks are not exactly
            #       related to files.....
            #
            # make sure that there are
            # direct images
            self._force_dirim()

            # the fluxcube contamination does not work for prism data
            if cont_model.lower() is "fluxcube":
                err_msg = ("{0:s}: Fluxcube contamination is not possible for "
                           "prism data!".format(self.taskname))
                raise aXeError(err_msg)

            # drizzled stamp images are not supported for prism data
            if sampling.lower() is "drizzle":
                err_msg = ("{0:s}: Drizzle sampling for the stamp images is "
                           "not possible for prism data!".format(self.taskname))
                raise aXeError(err_msg)

        # the extraction width must be set!
        if not extrfwhm:
            err_msg = ("{0:s}: extrfwhm must be > 0.0 to create PETs, but "
                       "extrfwhm={1:0.1f}!".format(self.taskname, extrfwhm))
            raise aXeError(err_msg)

        # negative extraction width is significant ONLY
        # if orient="NO"
        if ((orient < 0.0) and (extrfwhm < 0.0)):
            err_msg = ("{0:s}: Negative width extrfwhm={1:0.1f} together with "
                       "extraction orient=yes does NOT make sense!"
                       .format(self.taskname, extrfwhm))
            raise aXeError(err_msg)

        # for background extraction the width must be set!
        if back and not backfwhm:
            err_msg = ("{0:s}: With back=yes, the parameter backfwhm must be "
                       "set to create background PETs!".format(self.taskname))
            raise aXeError(err_msg)

        # extraction width and drizzle extraction width
        # must have the same sign
        if (extrfwhm * drzfwhm < 0.0):
            err_msg = ("{0:s}: extrfwhm={1:0.1f} and drzfwhm={2:0.1f} must BOTH"
                       "be either positive or negative!".format(self.taskname,
                                                                extrfwhm,
                                                                drzfwhm))
            raise aXeError(err_msg)
        else:
            # the extractionwidth must be larger than the
            # drizzle extraction width
            if not math.fabs(extrfwhm) > math.fabs(drzfwhm):
                err_msg = ("{0:s}: fabs(extrfwhm) MUST be larger than "
                           "fabs(drzfwhm), but extrfwhm={1:0.1f} and "
                           "drzfwhm={2:0.1f}!".format(self.taskname,
                                                      extrfwhm,
                                                      drzfwhm))
                raise aXeError(err_msg)

        # extraction width and background extraction width
        # must have the same sign
        if back and extrfwhm*backfwhm < 0.0:
            err_msg = ("{0:s}: extrfwhm={1:0.1f} and backfwhm={2:0.1f} must "
                       "BOTH be either positive or negative!"
                       .format(self.taskname, extrfwhm, backfwhm))
            raise aXeError(err_msg)

        # the background extraction width must be larger than the
        # object extraction width
        elif back and math.fabs(extrfwhm) > math.fabs(backfwhm):
            err_msg = ("{0:s}: fabs(backfwhm) MUST be larger than fabs(extrfwhm"
                       "), but backfwhm={1:0.1f} and extrfwhm={2:0.1f}!"
                       .format(self.taskname, backfwhm, extrfwhm))
            raise aXeError(err_msg)

        # for background extraction the number of background
        # pixels must be set
        if back and not np:
            err_msg = ("{0:s}: The parameter 'np' must be set for the "
                       "background PETs!".format(self.taskname))
            raise aXeError(err_msg)

        # for background extraction the interpolation
        # type must be set
        if back and not interp:
            err_msg = ("{0:s}: The parameter 'interp' must be set for the "
                       "background PETs!".format(self.taskname))
            raise aXeError(err_msg)

        # check for proper contamination
        # to allow optimal extraction
        if ((cont_model is "geometric") and (weights)):
            err_msg = ("{0:s}: Optimal weigthing needs quantitative "
                       "contamination! Please change to either the 'gauss'"
                       " or 'fluxcube' contamination model or drop optimal "
                       "weighting!".format(self.taskname))
            raise aXeError(err_msg)

    def check_axedrizzle(self, infwhm, outfwhm, back=False):
        """Comprises all file and file format checks for AXEDRIZZLE"""

        # check the DPP files
        self._check_dpps(back)

        # make sure that fabs(infwhm) and fabs(outfwhm) > 0.0
        if ((math.fabs(infwhm) == 0.0) or (math.fabs(outfwhm) == 0.0)):
            err_msg = ("{0:s}: fabs(infwhm) AND fabs(outfwhm) must be larger "
                       "than 0.0, but infwhm={1:0.1f} and outfwhm={2:0.1f}!"
                       .format(self.taskname, infwhm, outfwhm))
            raise aXeError(err_msg)

        # make sure that fabs(infwhm) > fabs(outfwhm)
        if (math.fabs(infwhm) < math.fabs(outfwhm)):
            err_msg = ("{0:s}: fabs(infwhm) MUST be larger than fabs(outfwhm),"
                       " but infwhm={1:0.1f} and outfwhm={2:0.1f}!"
                       .format(self.taskname, infwhm, outfwhm))
            raise aXeError(err_msg)

        # make sure that infwhm and outfwhm
        # have consistent sign
        if ((infwhm * outfwhm) < 0.0):
            err_msg = ("{0:s}: infwhm={1:0.1f} and outfwhm={2:0.1f} must BOTH"
                       "be either positive or negative!".format(self.taskname,
                                                                infwhm,
                                                                outfwhm))
            raise aXeError(err_msg)

    def check_axecrr(self, back):
        """
        Comprises all checks for the CRR version of AXEDRIZZLE
        """
        # make sure that background drizzling is off
        if back:
            err_msg = ("{0:s}: Background drizzling is NOT possible in the CRR"
                       "version of aXedrizzle!".format(self.taskname))
            raise aXeError(err_msg)

        # check for global background subtraction
        self._check_global_backsub()

    def check_simdispim_input(self, incat, config, lambda_psf,
                              model_spectra, model_images,
                              nx, ny, exptime, bck_flux, extraction,
                              extrfwhm, orient, slitless_geom, adj_sens):
        """Does basic checks on the parameters

        The method checks whether all input values are reasonable, e.g.
        the exposure time and background flux >= 0.0 and similar.
        Input files are checked for existence. Also the input type is
        checked for the numbers.

        Parameters
        ----------
        incat: str
            name of model object table
        config: str
            aXe configuration file name
        lambda_psf: float
            wavelength the object shapes were determined at
        model_spectra: str
            name of model spectra
        model_images: str
            name of model images
        nx: int
            number of pixels in x
        ny: int
            number of pixels in y
        exptime: float
            exposure time
        bck_flux: float
            flux in background
        extraction: bool
            flag for default extraction
        extrfwhm: float
            multiplier for extraction width
        orient: bool
            flag for tilted extraction
        slitless_geom: bool
            flag for slitless optimized extraction
        adj_sens: bool
            flag for adjusted flux conversion
        """

        # do the setup
        config_util.axe_setup(axesim=True)

        # check the existence of the
        # model object table
        if not os.path.isfile(config_util.getDATA(incat)):
            msg = ("The Model Object Table does not exist: {}"
                   .format(config_util.getDATA(incat)))
            raise aXeSIMError(msg)

        # check the existence of the
        # axe configuration file
        if not os.path.isfile(config_util.getCONF(config)):
            msg = ("The aXe configuration file does not exist: {}"
                   .format(config_util.getCONF(config)))
            raise aXeSIMError(msg)

        else:
            # load the aXe configuration file
            conf = configfile.ConfigFile(config_util.getCONF(config))

            # make the internal checks
            n_sens = conf.check_files(check_glob=False)

            # make sure there is
            # at least one sens. file
            if n_sens < 1:
                msg = ("There must be at least one sensitivity file in: {}"
                       .format(config_util.getCONF(config)))
                raise aXeSIMError(msg)

        # check whether the configuration files
        # allows the requested extraction
        if extraction and (slitless_geom or adj_sens):
            extr_ready = conf.confirm_extrkeys()

            # error and out
            if not extr_ready:
                msg = ("It is not possible to perform the requested"
                       "extraction. The likely cause is that the configuration"
                       "file does NOT contain the keywords 'POBJSIZE' or "
                       "'SMFACTOR' or their values are NOT reasonable "
                       "(e.g. <0.0)!")
                raise aXeSIMError(msg)

        # check the lambda_psf-value
        if ((lambda_psf is not None) and (lambda_psf <= 0.0)):
            msg = ("Value for 'lambda_psf' must be positive: {0:s}"
                   .format(str(lambda_psf)))
            raise aXeSIMError(msg)

        if (model_spectra is not None):
            # check the existence of the
            # model spectra file
            if not os.path.isfile(config_util.getDATA(model_spectra)):
                msg = ("The model spectra file does not exist: {}"
                       .format(config_util.getDATA(model_spectra)))
                raise aXeSIMError(msg)

        if model_images is not None:
            # check the existence of the
            # model images file
            if not os.path.isfile(config_util.getDATA(model_images)):
                msg = ("The model images file does not exist: "
                       .format(config_util.getDATA(model_images)))
                raise aXeSIMError(msg)

        # check the nx-value
        if ((nx is not None) and (nx <= 0.0)):
            msg = ("Value for 'nx' or 'nx_disp' must be positive: {0:g}"
                   .format(nx))
            raise aXeSIMError(msg)

        # check the ny-value
        if ((ny is not None) and (ny <= 0)):
            error_message = ("Value for 'ny' or 'ny_disp' must be "
                             "positive: {0:g}".format(ny))
            raise aXeSIMError(error_message)

        # check the exptime-value
        if ((exptime is not None) and (exptime < 0)):
            error_message = ("Value for 'exptime' or 'exptime_disp' must be "
                             "positive: {0:g}".format(exptime))
            raise aXeSIMError(error_message)

        # the extraction width must be set!
        if not extrfwhm:
            error_message = ("Value for 'extrfwhm' must not be 0.0 to create"
                             "PETs, but extrfwhm={0:0.1f}!".format(extrfwhm))
            raise aXeSIMError(error_message)

        # negative extraction width is significant ONLY
        # if orient="NO"
        if orient and extrfwhm < 0.0:
            error_message = ("Negative width extrfwhm={0:0.1f} together with "
                             "extraction orient=yes does NOT make sense!"
                             .format(extrfwhm))
            raise aXeSIMError(error_message)

        try:
            # convert to float
            bck = float(bck_flux)

            # check for positive value
            if bck < 0:
                error_message = ("Value for 'bck_flux' or 'bck_flux_disp'"
                                 " most be positive: {0:g}".format(bck_flux))
                raise aXeSIMError(error_message)

        # catch a string
        except ValueError:
            # check for existence of file
            if not os.path.isfile(config_util.getCONF(bck_flux)):
                error_message = ("The background file does not exist: {0}"
                                 .format(config_util.getCONF(bck_flux)))
                raise aXeSIMError(error_message)

    def check_simdirim_input(self, incat, config, tpass_direct,
                             model_spectra, model_images,
                             nx, ny, exptime, bck_flux):
        """Does basic checks on the parameters

        The method checks whether all input values are reasonable, e.g.
        the exposure time and background flux >= 0.0 and similar.
        Input files are checked for existence. Also the input type is
        checked for the numbers.

        Parameters
        ----------
        incat: str
            name of model object table
        config: str
            aXe configuration file name
        tpass_direct: str
            total passband file
        model_spectra: str
            name of model spectra
        model_images: str
            name of model images
        nx: int
            number of pixels in x
        ny: int
            number of pixels in y
        exptime: float
            exposure time
        bck_flux: float
            flux in background
        """

        # do the setup
        config_util.axe_setup(axesim=True)

        # check the existence of the
        # model object table
        if not os.path.isfile(config_util.getDATA(incat)):
            error_message = ("The Model Object Table does not exist: {0}"
                             .format(config_util.getDATA(incat)))
            raise aXeSIMError(error_message)

        # check the existence of the
        # axe configuration file
        if not os.path.isfile(config_util.getCONF(config)):
            error_message = ("The aXe configuration file does not exist: {0}"
                             .format(config_util.getCONF(config)))
            raise aXeSIMError(error_message)

        else:
            # load the aXe configuration file
            conf = configfile.ConfigFile(config_util.getCONF(config))

            # make the internal checks
            n_sens = conf.check_files(check_glob=False)

            # make sure there is
            # at least one sens. file
            if n_sens < 1:
                error_message = ("There must be at least one sensitivity "
                                 "file in: {0}".format(config_util.getCONF(config)))
                raise aXeSIMError(error_message)

        # check the existence of the
        # total passband file
        if not os.path.isfile(config_util.getSIMDATA(tpass_direct)):
            error_message = ("The total passband file does not exist: {0}"
                             .format(config_util.getSIMDATA(tpass_direct)))
            raise aXeSIMError(error_message)

        if model_spectra is not None:
            # check the existence of the
            # model spectra file
            if not os.path.isfile(config_util.getDATA(model_spectra)):
                error_message = ("The model spectra file does not exist: {0}"
                                 .format(config_util.getDATA(config)))
                raise aXeSIMError(error_message)

        if model_images is not None:
            # check the existence of the
            # model images file
            if not os.path.isfile(config_util.getDATA(model_images)):
                error_message = ("The model images file does not exist: {0}"
                                 .format(config_util.getDATA(config)))
                raise aXeSIMError(error_message)

        # check the nx-value
        if ((nx is not None) and (nx <= 0.0)):
            error_message = ("Value for 'nx' or 'nx_dir' must be positive: "
                             "{0:s}".format(str(nx)))
            raise aXeSIMError(error_message)

        # check the ny-value
        if ((ny is not None) and (ny <= 0)):
            error_message = ("Value for 'ny' or 'ny_dir' must be positive: "
                             "{0:s}".format(str(ny)))
            raise aXeSIMError(error_message)

        # check the exptime-value
        if ((exptime is not None) and (exptime < 0)):
            error_message = ("Value for 'exptime' or 'exptime_dir' must be "
                             "positive: {0:s}".format(str(exptime)))
            raise aXeSIMError(error_message)

        if bck_flux is not None:
            # check the bck_flux-value
            try:
                # convert to float
                bck = float(bck_flux)

                # check for positive value
                if bck < 0:
                    error_message = ("Value for 'bck_flux' or 'bck_flux_dir'"
                                     " must be positive: {0:s}"
                                     .format(str(bck_flux)))
                    raise aXeSIMError(error_message)

                # catch a string
            except ValueError:
                # check for existence of file
                if not os.path.isfile(config_util.getCONF(bck_flux)):
                    error_message = ("The background file does not exist: {0}"
                                     .format(config_util.getCONF(bck_flux)))
                    raise aXeSIMError(error_message)
