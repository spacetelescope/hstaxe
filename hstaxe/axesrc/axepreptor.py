import os
import logging
from astropy.io import fits
import stsci.imagestats as imagestats

from hstaxe import config as config_util
from hstaxe.axeerror import aXeError


from . import axelowlev
from . import configfile
from . import axetasks

# make sure there is a logger
_log = logging.getLogger(__name__)


class aXePrepArator:
    """This task prepares the science files for further processing within aXe.

    `axeprep` provides important keywords and is mandatory if axedrizzle is to be
    used later on.

    Inputs
    ------
    grism: string
        The input grism image
    objcat: string
        The input object catalog
    dirim: string
        The input direct image
    config: string
        The name of the aXe configuration file
    dmag: float or None
        a number to add to the MMAG_EXTRACT and MMAG_MARK values given in the
        configuration file.

        Sources which have magnitudes fainter than another cutoff magnitude are
        marked so that they will be completely ignored. The dmag value can be
        used to globally adjust these cutoffs (to account for a different
        signal-to-noise ratio in one dataset for example without having to
        resort to editing of the configuration file). (see GOL2AF) for the
        grism image (optional)
    master_bck: string
        Name of the master background image

    Notes
    -----
    axeprep provides three different processing steps:

    background subtraction:  Provided that an Input Object List is given for
        the grism image, axeprep uses the tasks sex2gol, gol2af and backest
        to mark the beam areas on the grism image as well as on the master
        background image. The median pixel values are derived for the unmarked
        pixels on both the grism image and on the master background image. 
        The master background, scaled to the level of the grism image, is
        finally subtracted from the grism image.

    exposure time normalization:
        The input file is normalized by the exposure time to transform the
        images into counts per second.

    gain correction:
        The input file is multiplied by the gain conversion factor
        (electrons/ADU) to transform the images from units of detector
        counts to electrons. For HST data, this is usually only necessary
        for NICMOS images, because ACS and WFC3 images will normally
        already be in units of electrons.

    """

    def __init__(self, grisim="",
                       objcat="",
                       dirim="",
                       config="",
                       dmag=None,
                       **params):
        if grisim:
            self.grisim = config_util.getDATA(grisim)
        else:
            raise ValueError("No grisim image specified for axeprep")
        self.objcat = config_util.getDATA(objcat)
        _log.info("\n**Using object catalog: {0}\n".format(self.objcat))
        self.dirim = config_util.getDATA(dirim)
        # make the config into a list of configs
        self.config_name = config
        self.config_obj = configfile.ConfigFile(config_util.getCONF(config))
        self.dmag = dmag
        self.params = params

        # store the master background
        if 'master_bck' in params:
            self.master_bck = config_util.getCONF(params['master_bck'])

    def _is_nicmos_data(self):
        """Check whether the data comes from NICMOS

        Note: the package hasn't been recently tested
        with NICMOS data

        """
        try:
            header = fits.getval(self.grisim, 'INSTRUME')
            return "NICMOS" in header
        except KeyError:
            _log.info("\nData is not from NICMOS\n")
            return False

    def _is_wfc3ir_data(self):
        """Check whether the data comes from WFC3 IR channel"""
        # open the image
        header = fits.getheader(self.grisim)
        try:
            if header['INSTRUME'] == 'WFC3':
                if header['DETECTOR'] == 'IR':
                    return True
            return False
        except KeyError:
            _log.info("\nData is not from HST/WFC3")
            return False

    def _make_mask(self):
        """Make the background mask file"""

        # set the use_direct flag
        # NOTE: the flag is only useful
        #       for the C-version, the python
        #       version does not need it!
        if self.dirim is not None:
            use_direct = True
        else:
            use_direct = False

        # run SEX2GOL
        axetasks.sex2gol(grism=self.grisim,
                         config=self.config_name,
                         in_sex=self.objcat,
                         use_direct=use_direct,
                         direct=self.dirim,
                         dir_hdu=None,
                         spec_hdu=None,
                         out_sex=None)

        # run GOL2AF
        axetasks.gol2af(grism=os.path.split(self.grisim)[-1],
                        config=self.config_name,
                        mfwhm=self.params['mfwhm'],
                        back=False,
                        orient=True,
                        slitless_geom=True,
                        exclude=False,
                        lambda_mark=None,
                        dmag=self.dmag,
                        out_af=None,
                        in_gol=None)

        #  run BACKEST
        axetasks.backest(grism=os.path.split(self.grisim)[-1],
                         config=self.config_name,
                         np=0,
                         interp=-1,
                         niter_med=None,
                         niter_fit=None,
                         kappa=None,
                         smooth_length=None,
                         smooth_fwhm=None,
                         old_bck=False,
                         mask=True,
                         in_af="",
                         out_bck=None)

    def _subtract_sky(self, ext_info, flag=-1.0e10):
        """Perform a classical background subtraction."""

        # Derive the name of all aXe products for a given image
        axe_names = config_util.get_axe_names(self.grisim, ext_info)
        msk_image_sc = axe_names['MSK'] + '[SCI]'

        # check for a previous background subtraction
        with fits.open(self.grisim, mode='update') as grism_file:
            if 'AXEPRBCK' in grism_file[ext_info['fits_ext']].header:
                # warn that this is the second time
                _log.info("WARNING: Image %25s seems to be already background "
                      "subtracted!".format(self.grisim))

            # Compute the ratio of the grism SCI image to the background image
            sci_data = grism_file['SCI', ext_info['ext_version']].data
            sci_header = grism_file['SCI', ext_info['ext_version']].header
            npix = int(sci_header["NAXIS1"]) * int(sci_header["NAXIS2"]) 

            bck_data = fits.getdata(self.master_bck)
            ratio_data = sci_data / bck_data

            # Flag pixels in the ratio image based on the grism image DQ array
            grism_dq_data = grism_file['DQ', ext_info['ext_version']].data
            ratio_data[grism_dq_data > 0.5] = flag

            # Flag pixels in the ratio image based on the grism image MSK file
            msk_file = fits.open(config_util.getOUTPUT(msk_image_sc.split("[")[0]),  'readonly')
            msk_data = msk_file['SCI'].data
            msk_file.close()

            ratio_data[msk_data < -900000] = flag

            # Flag pixels in the background image based on the grism image DQ
            # and MSK file
            bck_data[grism_dq_data > 0.5] = flag
            bck_data[msk_data < -900000] = flag

            # Compute stats for the ratio image
            stats = imagestats.ImageStats(ratio_data[ratio_data > flag],
                                          fields='midpt,stddev,npix', lower=None,
                                          upper=None, nclip=3, lsig=3.0, usig=3.0,
                                          binwidth=0.01)

            # Compute stats for the background image
            bstats = imagestats.ImageStats(bck_data[bck_data > flag],
                                          fields='midpt,stddev,npix', lower=None,
                                          upper=None, nclip=3, lsig=3.0, usig=3.0,
                                          binwidth=0.01)

            # Subtract the scaled background from the grism image
            # Reload a clean version of background
            bck_data = fits.getdata(self.master_bck)

            grism_file['SCI', ext_info['ext_version']].data -= bck_data * stats.midpt 
            grism_header = grism_file['SCI', ext_info['ext_version']].header

            # write some header iformation
            grism_header['SKY_SCAL'] = (float(stats.midpt),  'scaling value for the master background')
            grism_header['SKY_MAST'] = (float(bstats.midpt),  'average value of the master background')
            grism_header['SKY_IMG'] = (self.master_bck, 'name of the master background image')
            grism_header['F_SKYPIX'] = (float(stats.npix)/float(npix), 'fraction of pixels used for scaling')
            grism_header['AXEPRBCK'] = ('Done',          'flag that background subtraction was done')
        return 0

    def _subtract_nicsky(self, ext_info):
        """Special sky subtraction for NICMOS images"""
        # get the axe names
        axe_names = config_util.get_axe_names(self.grisim, ext_info)

        # check for a previous background subtraction
        fits_img = fits.open(self.grisim, 'readonly')
        fits_head = fits_img[ext_info['fits_ext']].header
        # npix = int(fits_head['NAXIS1']) * int(fits_head['NAXIS1'])

        if 'AXEPRBCK' in fits_head:
            # warn that this is the second time
            _log.info(f"WARNING: Image {self.grisim} seems to be already background "
                  "subtracted!")

        # close the fits
        fits_img.close()

        # do the special background fitting for NICMOS
        if self.params['backped'] is not None:
            nicback = axelowlev.aXe_NICBACK(self.grisim,
                                            self.config_name,
                                            self.master_bck,
                                            self.params['backped'])
        else:
            nicback = axelowlev.aXe_NICBACK(self.grisim,
                                            self.config_name,
                                            self.master_bck)
        nicback.runall()
        del nicback

        # check whether the background image exists
        if not os.path.isfile(config_util.getOUTPUT(axe_names['NBCK'])):
            err_msg = ("The background image: {0:s} does NOT exist!"
                       .format(config_util.getOUTPUT(axe_names['NBCK'])))
            raise aXeError(err_msg)

        # Subtract the scaled background image from the grism image
        # copy the image to the output directory first
        sci_file = fits.open(self.grisim, mode='update')
        bck_file = fits.open(config_util.getOUTPUT(axe_names['NBCK']),'readonly')
        sci_file['SCI', ext_info['ext_version']].data -= bck_file[1].data
        sci_file.close()
        bck_file.close()

        # open the background image
        fits_img = fits.open(config_util.getOUTPUT(axe_names['NBCK']),
                                                'readonly')
        fits_head = fits_img['BCK'].header

        # open the grism image and isolate the correct extension header
        grism_img = fits.open(self.grisim, mode='update')
        grism_header = grism_img[ext_info['fits_ext']].header

        if 'SKY_SCAL' in fits_head and 'F_SKYPIX' in fits_head:

            # transfer important keywords
            # to the grism image
            grism_header['SKY_SCAL'] = (float(fits_head['SKY_SCAL']),
                                        'scaling value of background')
            grism_header['F_SKYPIX'] = (float(fits_head['F_SKYPIX']),
                                        'fraction of pixels used for scaling')

        # close the fits again
        fits_img.close()

        # write some keywords
        grism_header['AXEPRBCK'] = ('Done', 'flag that background subtraction was done')
        grism_header['SKY_IMG'] = (self.master_bck, 'name of the 1st master background image')
        if self.params['backped'] is not None:
            grism_header['SKY_IMG2'] = (self.params['backped'],
                                        'name of the 2nd master background image')

        # close the image again
        grism_img.close()

        return True

    def _subtract_wfc3irsky(self, ext_info):
        """Special sky subtraction for WFC3 IR images"""
        # get the axe names
        axe_names = config_util.get_axe_names(self.grisim, ext_info)

        # check for a previous background subtraction
        try:
            fits.getval(self.grisim,
                        'AXEPRBCK', exten=ext_info['fits_ext'])
            _log.info(f"WARNING: Image {self.grisim} seems to be already background "
                  "subtracted! Continuing anyways...")
        except KeyError:
            _log.info("Previous subtraction not recorded, proceeding with background subtraction.")

        scalebck = axelowlev.aXe_SCALEBCK(os.path.split(self.grisim)[-1],
                                          os.path.split(axe_names['MSK'])[-1],
                                          os.path.split(self.config_name)[-1],
                                          os.path.split(self.master_bck)[-1])
        try:
            scalebck.run()
        except aXeError:
            _log.info("There was a problem with the background subtraction, "
                  "continuing without it")
            return False
        # check whether the background image exists
        bckfilename = config_util.getOUTPUT(axe_names['SGRI'])

        if not os.path.isfile(bckfilename):
            err_msg = ("The background image: {0:s} does NOT exist!"
                       .format(bckfilename))
            raise aXeError(err_msg)

        fits_image = fits.open(bckfilename, ext=0, mode='readonly')
        sky_frac = fits_image[0].header["FRACFIN"]
        scal_val = fits_image[0].header["SCALVAL"]
        bck_data = fits_image[1].data
        fits_image.close()

        if sky_frac < 0.1:
            _log.info("Low fraction of sky pixels found (<10%) continuing WITHOUT"
                  " sky subtraction")
            return False

        # Subtract the scaled background image from the grism image
        grism_file = fits.open(self.grisim, mode='update')
        grism_file[ext_info['ext_version']].data -= bck_data
        grism_header = grism_file[ext_info['fits_ext']].header

        # write some information into the
        # grism image header
        grism_header['AXEPRBCK'] = ('Done', 'flag that background subtraction was done')
        grism_header['SKY_IMG'] = (self.master_bck, 'name of the 1st master background image')

        # write some scaling information into the header
        grism_header['F_SKYPIX'] = (float(sky_frac),
                                    'fraction of pixels used for scaling')
        grism_header['SKY_CPS'] = (float(scal_val),
                                   'scale used for master sky == sky value [cps]')
        # close the grism image and save and the scaled image
        grism_file.close()
        return True

    def _subtract_background(self, ext_info):
        """Determine and subtract the background"""
        goodreturn = True

        # generate the mask image
        self._make_mask()

        # check whether we have NICMOS
        if self._is_nicmos_data():
            # make normal background subtraction
            goodreturn = self._subtract_nicsky(ext_info)

        elif self._is_wfc3ir_data():
            # make normal background subtraction
            goodreturn = self._subtract_wfc3irsky(ext_info)

        else:
            # make normal background subtraction
            goodreturn = self._subtract_sky(ext_info)

        if goodreturn:
            pstring = (f"AXEPREP: Image {self.grisim}[SCI,{str(ext_info['ext_version'])}] sky-subtracted.")
                       
            _log.info(pstring)

    def _check_low_skyfrac(self, frac):
        """Check for a low fraction of background pixels"""

        msg = (f"\nAXEPREP Image {self.grisim}: Only {frac*100.0} percent of the pixels "
               "were used in the background scaling!")
                                                             
        raise aXeError(msg)

    def _check_second_normalization(self):
        """Check whether the data is already normalized.
        """

        msg = (f"AXEPREP: Image {self.grisim} has already been normalized! Will not renormalize")
        raise aXeError(msg)

    def _check_second_gaincorr(self):
        """Check whether the gain correction had already been applied.
        """

        msg = (f"AXEPREP: Image: {self.grisim} has already been gain corrected! Will not reapply.")
        raise aXeError(msg)

    def _check_gain_correction(self):
        """Check whether the gain correction had already been applied
        """

        msg = ("AXEPREP: Non-NICMOS images such as: {self.grisim} usually are "
               "already gain corrected!")
        raise aXeError(msg)

    def _transform_to_cps(self, ext_info, conf):
        """Transform the image from [e] to [e/s]"""

        # check for a previous normalization
        with fits.open(self.grisim, mode='update') as grism_image:
            grism_header = grism_image[ext_info['fits_ext']].header

            if 'AXEPRNOR' in grism_header:
                # check whether a second normalization
                # is really desired
                self._check_second_normalization()

            # check for the exposure time keyword
            if not conf.get_gkey('EXPTIME'):
                exptime_kword = 'EXPTIME'
            else:
                exptime_kword = conf.get_gvalue('EXPTIME')

            # get the exposure time
            exptime = grism_image[0].header[exptime_kword]

            pstring = ("AXEPREP: Image {0:25s}[SCI,{1:s}] exposure time "
                       "{2:7.1f}.".format(self.grisim,
                                          str(ext_info['ext_version']),
                                          exptime))
            _log.info(pstring)

            # Divide the grism image SCI and ERR arrays by exptime
            grism_image['SCI', ext_info['ext_version']].data /= exptime
            grism_image['ERR', ext_info['ext_version']].data /= exptime

            # write a header entry
            grism_header['AXEPRNOR'] = ('Done', 'flag that exposure time normalization was done')

            # extract the value of the constant
            # sky subtracted in multidrizzle
            if 'MDRIZSKY' in grism_image:
                mdrizsky = grism_header['MDRIZSKY']
            else:
                mdrizsky = None

            # check whether a global bias subtraction was done
            # in axeprep
            # if backgr == 'YES':
            if ('backgr' in self.params) and (self.params['backgr']):
                # read the desctiptors written
                # in the global bias subtraction
                sky_scal = grism_header['SKY_SCAL']
                sky_mast = grism_header['SKY_MAST']

                # compute the total subtracted background level
                # in multidrizzle and axeprep (in cps)
                if mdrizsky is not None:
                    sky_cps = (sky_scal * sky_mast + mdrizsky) / exptime
                else:
                    sky_cps = (sky_scal * sky_mast) / exptime
            else:

                # compute the total subtracted bias level (in cps)
                if mdrizsky is not None:
                    sky_cps = mdrizsky / exptime
                else:
                    sky_cps = 0.0

            # write the subtracted background level
            # into a defined descriptor for later use
            grism_header['SKY_CPS'] = (sky_cps,  'sky level in cps')

        return 0

    def _apply_gain_correction(self, ext_info):
        """Apply the gain correction"""
        # open the fits image
        fits_img = fits.open(self.grisim, 'readonly')

        # get the gain
        # for NICMOS:
        gain = float(fits_img[0].header['ADCGAIN'])

        # get the header of the target extension
        fits_head = fits_img[ext_info['fits_ext']].header
        # if 'AXEGAINC' in fits_head:
        #    dec = self._check_second_gaincorr()

        # close the fits image again
        fits_img.close()

        # multiply both the sci and the error array by the gain
        file_a = fits.open(self.grisim, mode='update')
        file_a["SCI"].data *= gain
        file_a["ERR"].data *= gain
        file_a.close()

        # open the fits image
        fits_img = fits.open(self.grisim, mode='update')

        # write the flag into the science extension
        fits_img[ext_info['fits_ext']].header['AXEGAINC'] = ('Done', 'flag that gain correction was done')

        # write the flag into the error extension, assuming it to be next to the science extension
        fits_img[ext_info['fits_ext']+1].header['AXEGAINC'] = ('Done', 'flag that gain correction was done')

        # close the image
        fits_img.close()

    def run(self):
        """Run AXEPREP on one slitless image"""

        # get the extension info
        ext_info = config_util.get_ext_info(self.grisim, self.config_obj)

        # make a background PET if necessary
        if 'backgr' in self.params and self.params['backgr']:
            self._subtract_background(ext_info)

        # make a background PET if necessary
        if 'norm' in self.params and self.params['norm']:
            self._transform_to_cps(ext_info, self.config_obj)

        # check wheter the gain correction is desired
        if 'gcorr' in self.params and self.params['gcorr']:

            # check whether we have NICMOS data
            if not self._is_nicmos_data():

                # check whether gain correction IS desired
                if self._check_gain_correction():

                    # make the gain correction
                    self._apply_gain_correction(ext_info)

            else:

                # make the gain correction
                self._apply_gain_correction(ext_info)

        # return something
        return 1
