import os
from astropy.io import fits
import stsci.imagestats as imagestats

from . import axelowlev
from . import configfile
from . import axetasks

from .. import config
from pyaxe.axeerror import aXeError


class aXePrepArator(object):
    def __init__(self, grisim, objcat, dirim, config, dmag, **params):
        # store the input
        self.grisim = grisim
        self.objcat = objcat
        self.dirim = dirim
        self.config = config
        self.dmag = dmag
        print(params)

        print(self.grisim, self.objcat)
        self.params = params
        # store the master background
        if 'master_bck' in params:
            self.master_bck = params['master_bck']

    def _is_nicmos_data(self):
        """Check whether the data comes from NICMOS"""
        # open the image
        header = fits.getvalue(config.getIMAGE(self.grisim), 'INSTRUME')
        return "NICMOS" in header

    def _is_wfc3ir_data(self):
        """Check whether the data comes from WFC3 IR channel"""
        # open the image
        header = fits.getvalue(config.getIMAGE(self.grisim), 'DETECTOR')
        return "IR" in header

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
        print("config = {}".format(self.config))
        axetasks.sex2gol(grism=self.grisim,
                         config=self.config,
                         in_sex=self.objcat,
                         use_direct=use_direct,
                         direct=self.dirim,
                         dir_hdu=None,
                         spec_hdu=None,
                         out_sex=None)

        # run GOL2AF
        axetasks.gol2af(grism=self.grisim,
                        config=self.config,
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
        axetasks.backest(grism=self.grisim,
                         config=self.config,
                         np=0,
                         interp=-1,
                         niter_med=None,
                         niter_fit=None,
                         kappa=None,
                         smooth_length=None,
                         smooth_fwhm=None,
                         old_bck=False,
                         mask=True,
                         in_af=None,
                         out_bck=None)

    def _subtract_sky(self, ext_info):
        """Make a classical background subtraction"""
        # get the axe names
        axe_names = config.get_axe_names(self.grisim, ext_info)

        msk_image_sc = axe_names['MSK'] + '[SCI]'

        # check for a previous background subtraction
        fits_img = fits.open(config.getIMAGE(self.grisim), 'readonly')
        fits_head = fits_img[ext_info['fits_ext']].header
        npix = int(fits_head['NAXIS1']) * int(fits_head['NAXIS1'])

        if 'AXEPRBCK' in fits_head:
            # warn that this is the second time
            print("WARNING: Image %25s seems to be already background "
                  "subtracted!".format(config.getIMAGE(self.grisim)))

        # close the fits
        fits_img.close()

        # Compute the ratio of the grism SCI image to the background image
        sci_file = fits.open(config.getIMAGE(self.grisim), 'readonly')
        sci_data = sci_file['SCI', ext_info['ext_version']].data
        bck_file = fits.open(config.getCONF(self.master_bck), 'readonly')
        bck_data = bck_file[0].data
        sci_data /= bck_data

        # Flag pixels in the ratio image based on the grism image DQ array
        dq_data = sci_file['DQ', ext_info['ext_version']].data
        sci_data[dq_data > 0.5] = -1.0e10

        # Flag pixels in the ratio image based on the grism image MSK file
        msk_file = fits.open(config.getOUTPUT(msk_image_sc.split("[")[0]), 'readonly')
        msk_data = msk_file['SCI'].data
        sci_data[msk_data < -900000] = -1.0e10

        # Flag pixels in the background image based on the grism image DQ
        # and MSK file
        bck_data[dq_data > 0.5] = -1.0e10
        bck_data[msk_data < -900000] = -1.0e10

        # Compute stats for the ratio image
        stats = imagestats.ImageStats(sci_data[sci_data > -1.0e9],
                                      fields='midpt,stddev,npix', lower=None,
                                      upper=None, nclip=3, lsig=3.0, usig=3.0,
                                      binwidth=0.01)
        flt_ave = stats.midpt
        # flt_std = stats.stddev
        flt_npx = stats.npix
        frac_pix = float(flt_npx)/float(npix)

        # Compute stats for the background image
        stats = imagestats.ImageStats(bck_data[bck_data > -1.0e9],
                                      fields='midpt,stddev,npix', lower=None,
                                      upper=None, nclip=3, lsig=3.0, usig=3.0,
                                      binwidth=0.01)
        mst_ave = stats.midpt
        # mst_std = stats.stddev
        # mst_npx = stats.npix

        sci_file.close()
        bck_file.close()
        msk_file.close()

        # Subtract the scaled background from the grism image
        sci_file = fits.open(config.getIMAGE(self.grisim), 'update')
        bck_file = fits.open(config.getCONF(self.master_bck), 'readonly')
        sci_file['SCI', ext_info['ext_version']].data -= flt_ave*bck_file[0].data
        sci_file.close()
        bck_file.close()

        # open the fits image ands isolate the correct extension
        grism_img = fits.open(config.getIMAGE(self.grisim), 'update')
        grism_header = grism_img[ext_info['fits_ext']].header

        # write some header iformation
        grism_header['SKY_SCAL'] = (float(flt_ave),  'scaling value for the master background')
        grism_header['SKY_MAST'] = (float(mst_ave),  'average value of the master background')
        grism_header['SKY_IMG'] = (self.master_bck, 'name of the master background image')
        grism_header['F_SKYPIX'] = (frac_pix,        'fraction of pixels used for scaling')
        grism_header['AXEPRBCK'] = ('Done',          'flag that background subtraction was done')

        # save the image
        grism_img.close()

        return 0

    def _subtract_nicsky(self, ext_info):
        """Special sky subtraction for NICMOS images"""
        # get the axe names
        axe_names = config.get_axe_names(self.grisim, ext_info)

        # check for a previous background subtraction
        fits_img = fits.open(config.getIMAGE(self.grisim), 'readonly')
        fits_head = fits_img[ext_info['fits_ext']].header
        # npix = int(fits_head['NAXIS1']) * int(fits_head['NAXIS1'])

        if 'AXEPRBCK' in fits_head:
            # warn that this is the second time
            print("WARNING: Image {0:25s} seems to be already background "
                  "subtracted!".format(config.getIMAGE(self.grisim)))

        # close the fits
        fits_img.close()

        # do the special background fitting for NICMOS
        if self.params['backped'] is not None:
            nicback = axelowlev.aXe_NICBACK(self.grisim,
                                            self.config,
                                            self.master_bck,
                                            self.params['backped'])
        else:
            nicback = axelowlev.aXe_NICBACK(self.grisim,
                                            self.config,
                                            self.master_bck)
        nicback.runall()
        del nicback

        # check whether the background image exists
        if not os.path.isfile(config.getOUTPUT(axe_names['NBCK'])):
            err_msg = ("The background image: {0:s} does NOT exist!"
                       .format(config.getOUTPUT(axe_names['NBCK'])))
            raise aXeError(err_msg)

        # Subtract the scaled background image from the grism image
        sci_file = fits.open(config.getIMAGE(self.grisim), 'update')
        bck_file = fits.open(config.getOUTPUT(axe_names['NBCK']),'readonly')
        sci_file['SCI', ext_info['ext_version']].data -= bck_file[1].data
        sci_file.close()
        bck_file.close()

        # open the background image
        fits_img = fits.open(config.getOUTPUT(axe_names['NBCK']),
                                                'readonly')
        fits_head = fits_img['BCK'].header

        # open the grism image and isolate the correct extension header
        grism_img = fits.open(config.getIMAGE(self.grisim), 'update')
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
        """Special sky subtraction for NICMOS images"""
        # get the axe names
        axe_names = config.get_axe_names(self.grisim, ext_info)

        # check for a previous background subtraction
        fits_img = fits.open(config.getIMAGE(self.grisim), 'readonly')
        fits_head = fits_img[ext_info['fits_ext']].header
        # npix = int(fits_head['NAXIS1']) * int(fits_head['NAXIS1'])

        if 'AXEPRBCK' in fits_head:
            # warn that this is the second time
            print("WARNING: Image {0:25s} seems to be already background "
                  "subtracted!".format(config.getIMAGE(self.grisim)))

        # close the fits
        fits_img.close()

        scalebck = axelowlev.aXe_SCALEBCK(self.grisim, axe_names['MSK'],
                                          self.config, self.master_bck)
        try:
            scalebck.runall()
        except aXeError:
            print("There was a problem with the background subtraction, "
                  "continuing without it")
            return False

        # check whether the background image exists
        bckfilename = config.getOUTPUT(axe_names['SGRI'])
        if not os.path.isfile(bckfilename):
            err_msg = ("The background image: {0:s} does NOT exist!"
                       .format(bckfilename))
            raise aXeError(err_msg)

        # in case of a low  value, dont do the subtraction if less than 10%
        sky_frac = fits.getval(bckfilename, "FRACFIN", ext=0)

        if sky_frac < 0.1:
            print("Low fraction of sky pixels found (<10%) continuing WITHOUT"
                  " sky subtraction")
            return False

        # Subtract the scaled background image from the grism image
        sci_file = fits.open(config.getIMAGE(self.grisim), 'update')
        bck_file = fits.open(bckfilename)
        sci_file['SCI', ext_info['ext_version']].data -= bck_file[1].data
        sci_file.close()
        bck_file.close()

        # open the background image
        fits_img = fits.open(bckfilename, 'readonly')
        fits_head = fits_img[0].header

        # open the grism image and isolate the correct extension header
        grism_img = fits.open(config.getIMAGE(self.grisim), 'update')
        grism_header = grism_img[ext_info['fits_ext']].header

        # write some information into the
        # grism image header
        grism_header['AXEPRBCK'] = ('Done', 'flag that background subtraction was done')
        grism_header['SKY_IMG'] = (self.master_bck, 'name of the 1st master background image')

        # write some scaling information into the header
        grism_header['F_SKYPIX'] = (float(fits_head['FRACFIN']),
                                    'fraction of pixels used for scaling')
        grism_header['SKY_CPS'] = (float(fits_head['SCALVAL']),
                                   'scale used for master sky == sky value [cps]')

        # close the grism image
        # and the scaled image
        fits_img.close()
        grism_img.close()

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
            pstring = ("AXEPREP: Image {0:25s}[SCI,{1:s}] sky-subtracted."
                       .format(self.grisim, str(ext_info['ext_version'])))
            print(pstring)

    def _check_low_skyfrac(self, frac):
        """Check for a low fraction of background pixels"""

        msg = ("\nAXEPREP Image {0:s}: Only {1:0.1f} percent of the pixels "
               "were used in the background scaling!".format(self.grisim,
                                                             frac*100.0))
        raise aXeError(msg)

    def _check_second_normalization(self):
        """Check whether the data is already normalized"""

        msg = ("AXEPREP: Image %25s has already been normalized!"
               .format(self.grisim))
        raise aXeError(msg)

    def _check_second_gaincorr(self):
        """Check whether the gain correction had already been applied"""

        msg = ("AXEPREP: Image: {0:s} has already been gain corrected!"
               .format(self.grisim))
        raise aXeError(msg)

    def _check_gain_correction(self):
        """Check whether the gain correction had already been applied"""

        msg = ("AXEPREP: Non-NICMOS images such as: {0:s} usually are "
               "already gain corrected!".format(self.grisim))
        raise aXeError(msg)

    def _transform_to_cps(self, ext_info, conf):
        """Transform the image from [e] to [e/s]"""
        dec = 1

        # check for a previous normalization
        fits_img = fits.open(config.getIMAGE(self.grisim), 'readonly')
        fits_head = fits_img[ext_info['fits_ext']].header

        if 'AXEPRNOR' in fits_head:
            # check whether a second normalization
            # is really desired
            dec = self._check_second_normalization()

        # check for the exposure time keyword
        if not conf.get_gkey('EXPTIME'):
            exptime_kword = 'EXPTIME'
        else:
            exptime_kword = conf.get_gvalue('EXPTIME')

        # get the exposure time
        exptime = fits_img[0].header[exptime_kword]

        if dec:
            pstring = ("AXEPREP: Image {0:25s}[SCI,{1:s}] exposure time "
                       "{2:7.1f}.".format(self.grisim,
                                          str(ext_info['ext_version']),
                                          exptime))
            print(pstring)

            # Divide the grism image SCI and ERR arrays by exptime
            file_a = fits.open(config.getIMAGE(self.grisim), 'update')
            file_a['SCI', ext_info['ext_version']].data /= exptime
            file_a['ERR', ext_info['ext_version']].data /= exptime
            file_a.close()

            # open the grism image and isolate the correct header
            grism_img = fits.open(config.getIMAGE(self.grisim), 'update')
            grism_header = grism_img[ext_info['fits_ext']].header

            # write a header entry
            grism_header['AXEPRNOR'] = ('Done', 'flag that exposure time normal was done')

            # extract the value of the constant
            # sky subtracted in multidrizzle
            if 'MDRIZSKY' in fits_head:
                mdrizsky = fits_head['MDRIZSKY']
            else:
                mdrizsky = None

            # check whether a global bias subtraction was done
            # in axeprep
            # if backgr == 'YES':
            if ('backgr' in self.params) and (self.params['backgr']):
                # read the desctiptors written
                # in the global bias subtraction
                sky_scal = fits_head['SKY_SCAL']
                sky_mast = fits_head['SKY_MAST']

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

            # close the fits image
            grism_img.close()

        else:
            return -1

        fits_img.close()

        return 0

    def _apply_gain_correction(self, ext_info):
        """Apply the gain correction"""
        # open the fits image
        fits_img = fits.open(config.getIMAGE(self.grisim), 'readonly')

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
        file_a = fits.open(config.getIMAGE(self.grisim), 'update')
        file_a["SCI"].data *= gain
        file_a["ERR"].data *= gain
        file_a.close()

        # open the fits image
        fits_img = fits.open(config.getIMAGE(self.grisim), 'update')

        # write the flag into the science extension
        fits_img[ext_info['fits_ext']].header['AXEGAINC'] = ('Done', 'flag that gain correction was done')

        # write the flag into the error extension, assuming it to be next to the science extension
        fits_img[ext_info['fits_ext']+1].header['AXEGAINC'] = ('Done', 'flag that gain correction was done')

        # close the image
        fits_img.close()

    def run(self):
        """Run AXEPREP on one slitless image"""

        # load the configuration files;
        # get the extension info
        conf = configfile.ConfigFile(config.getCONF(self.config))
        ext_info = config.get_ext_info(config.getIMAGE(self.grisim), conf)

        # make a background PET if necessary
        if 'backgr' in self.params and self.params['backgr']:
            self._subtract_background(ext_info)

        # make a background PET if necessary
        if 'norm' in self.params and self.params['norm']:
            self._transform_to_cps(ext_info, conf)

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

        del conf

        # return something
        return 1
