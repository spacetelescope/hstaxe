import os
import numpy as np
import logging

from astropy.io import fits

from stsci.imagestats import ImageStats
import stsci.convolve as convolve
from stsci.image.numcombine import numCombine


from drizzlepac.drizCR import quickDeriv
from drizzlepac import minmed
from drizzlepac import adrizzle
from drizzlepac.astrodrizzle import ablot

from hstaxe.axeerror import aXeError

# make sure there is a logger
_log = logging.getLogger(__name__)

class Drizzle:
    """Class to wrap drizzle command"""
    def __init__(self):

        _log.info("No init in class")

    def run(self, data,
                  outdata,
                  in_mask,
                  outweig,
                  coeffs,
                  wt_scl,
                  drizzle_params,
                  img_nx,
                  img_ny):
        """
        drizzle using drizzlepac.astrodrizzle

        Parameters
        ----------
        data: str
            The name of the input image file or association table which is to be "drizzled"

        in_mask: str
            input mask for blocking pixels

        outdata: str
            The name for the output data image desired by the user

        outweig: str
            The name for the output weight image

        coeffs: str
            distortion coefficients

        wt_scl: str
            Weighting factor for input image. If wt_scl=exptime then wt_scl
            will be set equal to the exposure time found in the image header.
            This is the standard operation and is recommended. It is also
            possible to give wt_scl=expsq for weighting by the square of
            exposure time. The latter is optimal for read-noise dominated
            images.


        drizzle_params:

        img_nx: int
            output size in x

        img_ny: int
            output size in y


        """

        # Check for file names that are too long for the drizzle task
        # this is a fits limitation
        if len(data) > 80:
            err_msg = ("File name '{0:s}' is too long (>80 chars)"
                       "for drizzle task".format(data))
            raise aXeError(err_msg)
        if len(outdata) > 80:
            err_msg = ("File name '{0:s}' is too long (>80 chars) "
                       "for drizzle task".format(outdata))
            raise aXeError(err_msg)

        ret = adrizzle.drizzle(data, outdata,
                               outweig=outweig,
                               in_mask=in_mask,
                               wt_scl=wt_scl,
                               coeffs=coeffs,
                               outnx=img_nx,
                               outny=img_ny,
                               in_un=drizzle_params['IN_UN'],
                               out_un=drizzle_params['OUT_UN'],
                               pixfrac=drizzle_params['PFRAC'],
                               scale=drizzle_params['PSCALE'],
                               kernel=drizzle_params['KERNEL'], Stdout=1)

        for i in range(len(ret)):
            _log.info(ret[i])


class MedianCombine:
    """Class to median-combine individual drizzles"""
    def __init__(self,
                 contributors,
                 drizzle_params,
                 ext_names):
        """Initialize the class"""

        # store the parameters
        self.combine_maskpt = drizzle_params['combine_maskpt']
        self.combine_type = drizzle_params['combine_type']
        self.combine_nsigma1 = drizzle_params['combine_nsigma1']
        self.combine_nsigma2 = drizzle_params['combine_nsigma2']
        self.combine_nlow = drizzle_params['combine_nlow']
        self.combine_nhigh = drizzle_params['combine_nhigh']
        self.combine_lthresh = drizzle_params['combine_lthresh']
        self.combine_hthresh = drizzle_params['combine_hthresh']
        self.combine_grow = drizzle_params['combine_grow']
        self.rdnoise = drizzle_params['RDNOISE']

        self.ext_names = ext_names

        # store the name of the median image
        self.median_image = ext_names['MED']

        self.input_data = self._get_inputs(contributors)

    def _get_inputs(self, contributors):
        """
        Extract the inputs for the median combine
        """
        # generate an empty dictionary
        input_data = {}

        # go over all contributing objects
        sci_imgs = []
        wht_imgs = []
        exp_vals = []
        rdn_vals = []
        sky_vals = []
        exp_tot = 0.0
        for one_contrib in contributors:

            # put the image names to the list
            sci_imgs.append(one_contrib.ext_names['SING_SCI'])
            wht_imgs.append(one_contrib.ext_names['SING_WHT'])

            # put image properties to the list
            exp_vals.append(one_contrib.info['EXPTIME'])
            rdn_vals.append(self.rdnoise)
            if 'SKY_CPS' in one_contrib.info:
                sky_vals.append(one_contrib.info['SKY_CPS'])
            else:
                err_msg = ("Sky value missing for image: {0:s}!"
                           .format(input_data['sci_imgs']))
                raise aXeError(err_msg)

            # compose the total exposure time
            exp_tot += one_contrib.info['EXPTIME']

        # put the images to the dictionary
        input_data['sci_imgs'] = sci_imgs
        input_data['wht_imgs'] = wht_imgs

        # put the values to the dictionary
        input_data['exp_vals'] = exp_vals
        input_data['rdn_vals'] = rdn_vals
        input_data['sky_vals'] = sky_vals
        input_data['exp_tot'] = exp_tot

        # return the dictionary
        return input_data

    def run(self):
        """
        Run the median combine step

        The code was either directly stolen from the corresponding
        pydrizzle version or done after this version. Necessary
        adjustments to the slitless data were applied.
        """
        sci_data = []

        for one_image in self.input_data['sci_imgs']:
            if os.access(one_image, os.F_OK):
                in_fits = fits.open(one_image, 'readonly')
                sci_data.append(in_fits[0].data)
                in_fits.close()

        wht_data = []
        for one_image in self.input_data['wht_imgs']:
            if os.access(one_image, os.F_OK):
                in_fits = fits.open(one_image, 'readonly')
                wht_data.append(in_fits[0].data)
                in_fits.close()
            else:
                _log.info("{0:s} not found/created by drizzle"
                      "...skipping it.".format(one_image))

        if len(sci_data) != len(wht_data):
            _log.info("The number of single_sci images created by "
                  "drizzle does not match the number of single_wht"
                  " files created!")
            raise aXeError("drizzle error")

        weight_mask_list = []

        # added the except so that if the image area contains only
        # zeros then the zero value is returned which is better for later
        # processing
        # we dont understand why the original lower=1e-8 value was
        # supplied unless it was for the case of spectral in the normal
        # field of view see #1110
        for wht_arr in wht_data:
            try:
                tmp_mean_value = self.combine_maskpt * ImageStats(wht_arr,lower=1e-8,lsig=None,usig=None,fields="mean",nclip=0).mean
            except (ValueError, AttributeError):
                tmp_mean_value = 0.
                _log.info("tmp_mean_value set to 0 because no good "
                      "pixels found; {0:s}".format(self.ext_names["MEF"]))
            except:
                tmp_mean_value = 0.
                _log.info("tmp_mean_value set to 0; possible uncaught "
                      "exception in dither.py; {0:s}"
                      .format(self.ext_names["MEF"]))

            weight_mask = np.zeros(wht_arr.shape, dtype=np.uint8)
            np.putmask(weight_mask, np.less(wht_arr, tmp_mean_value), 1)

            weight_mask_list.append(weight_mask)

        if len(sci_data) < 2:
            _log.info('\nNumber of images to flatten: %i!' % len(sci_data))
            _log.info('Set combine type to "minimum"!')
            self.combine_type = 'minimum'

        if (self.combine_type == "minmed"):
            # Create the combined array object using the minmed algorithm
            result = minmed(sci_data,  # list of input data to be combined.
                            wht_data,# list of input data weight images to be combined.
                            self.input_data['rdn_vals'],  # list of readnoise values to use for the input images.
                            self.input_data['exp_vals'],  # list of exposure times to use for the input images.
                            self.input_data['sky_vals'],  # list of image background values to use for the input images
                            weightMaskList = weight_mask_list,  # list of imput data weight masks to use for pixel rejection.
                            combine_grow = self.combine_grow,  # Radius (pixels) for neighbor rejection
                            combine_nsigma1 = self.combine_nsigma1,  # Significance for accepting minimum instead of median
                            combine_nsigma2 = self.combine_nsigma2  # Significance for accepting minimum instead of median
                            )
        else:
            # _log.info 'going to other', combine_type
            # Create the combined array object using the numcombine task
            result = numCombine(sci_data,
                                numarrayMaskList=weight_mask_list,
                                combinationType=self.combine_type,
                                nlow=self.combine_nlow,
                                nhigh=self.combine_nhigh,
                                upper=self.combine_hthresh,
                                lower=self.combine_lthresh
                                )

        # _log.info result.combArrObj
        hdu = fits.PrimaryHDU(result.combArrObj)
        hdulist = fits.HDUList([hdu])
        hdulist[0].header['EXPTIME'] = (self.input_data['exp_tot'],
                                        'total exposure time')
        hdulist.writeto(self.median_image)

        # delete the various arrays
        for one_item in sci_data:
            del one_item
        del sci_data
        for one_item in wht_data:
            del one_item
        del wht_data
        for one_item in weight_mask_list:
            del one_item
        del weight_mask_list


class Blot:
    """
    Class to wrap the blot command
    """
    def __init__(self):
        pass


    def run(self, in_data, out_data, out_nx, out_ny,
            drizzle_params):
    #     """
    #     Do the actual blot
    #     """
    #     drizzle.blot(data=in_data,
    #                outdata=out_data,
    #                scale=drizzle_params['PSCALE'],
    #                coeffs=coeffs,
    #                outnx=out_nx,
    #                outny=out_ny,
    #                interpol=mult_drizzle_par['blot_interp'],
    #                sinscl=mult_drizzle_par['blot_sinscl'],
    #                in_un=drizzle_params['IN_UN'],
    #                out_un=drizzle_params['OUT_UN'],
    #                expkey='exptime',
    #                expout='input')
        self._a_blot_image(in_data, out_data, 
                           sinscl=drizzle_params['blot_sinscl'],
                           out_nx=out_nx,
                           out_ny=out_ny,
                           interp=interpol)

    def _a_blot_image(self,
                      image_to_blot,
                      flt_image,
                      blotted_output,
                      sinscl=sinscl,
                      interp=interp):
            """
            Blot one image.

            Thats just a simple wrapper around the task blot in astrodrizzle

            Parameters
            ----------
            image_to_blot: str
                the input image name, either the grism or direct drizzled image

            blotted_output: str
                the name of the output blotted image

            """
            try:
                blot_header = fits.getheader(image_to_blot)
                blot_wcs = HSTWCS(image_to_blot)  # assume simple
                image_data = fits.getdata(image_to_blot)
                flt_header = fits.getheader(flt_image)
                flt_wcs = HSTWCS(flt_image)
            except:
                return IOError("File type of fits image is not "
                               "supported {0:s}".format(image_to_blot))

            # outimage is just the data array
            outimage = ablot.do_blot(image_data.astype(np.float32),
                                                  blot_wcs,
                                                  flt_wcs,
                                                  1.,
                                                  interp=interp,
                                                  sinscl=1.,
                                                  coeffs=True,
                                                  wcsmap=None,
                                                  stepsize=10)


            try:
                newimage = fits.PrimaryHDU()
                newimage.data = outimage
                newimage.header = flt_header
                newimage.header.update(flt_wcs.to_header())
                newimage.verify('silentfix')
                newimage.writeto(blotted_output)
            except:
                raise IOError("Problem writing fits image {0:s}".format(blotted_output))


class Deriv:
    """
    Class for the deriv-command
    """
    def __init__(self):
        """
        Initializes the class
        """
        pass

    def _absoluteSubtract(self, array, tmpArray, outArray):
        """
        Subtract the absolute value of two images
        """
        # subtract shifted image from imput image
        tmpArray = array - tmpArray
        # take the absolute value of tmpArray
        tmpArray = numpy.fabs(tmpArray)
        # save maximum value of outArray or tmpArray and save in outArray
        outArray = numpy.maximum(tmpArray, outArray)
        # zero out tmpArray before reuse
        tmpArray = tmpArray * 0.

        return (tmpArray, outArray)

    def _qderiv(self, array):
        """
        Take the absolute derivate of an image in memory
        """

        # Create 2 empty arrays in memory of the same dimensions as 'array'
        tmpArray = numpy.zeros(array.shape, dtype=numpy.float64)
        outArray = numpy.zeros(array.shape, dtype=numpy.float64)

        # Get the length of an array side
        (naxis1, naxis2) = array.shape

        # Main derivate loop:
        # Shift images +/- 1 in Y.
        for y in range(-1, 2, 2):
            if y == -1:
                # shift input image 1 pixel right
                tmpArray[0:(naxis1-1), 1:(naxis2-1)] = array[0:(naxis1-1),
                                                             0:(naxis2-2)]

            else:
                # shift input image 1 pixel left
                tmpArray[0:(naxis1-1), 0:(naxis2-2)] = array[0:(naxis1-1),
                                                             1:(naxis2-1)]

            # subtract the arrays
            (tmpArray, outArray) = self._absoluteSubtract(array,
                                                          tmpArray,
                                                          outArray)

        # Shift images +/- 1 in X.
        for x in range(-1, 2, 2):
            if x == -1:
                # shift input image 1 pixel right
                tmpArray[1:(naxis1-1), 0:(naxis2-1)] = array[0:(naxis1-2),
                                                             0:(naxis2-1)]

            else:
                # shift input image 1 pixel left
                tmpArray[0:(naxis1-2), 0:(naxis2-1)] = array[1:(naxis1-1),
                                                             0:(naxis2-1)]

            # subtract the arrays
            (tmpArray, outArray) = self._absoluteSubtract(array,
                                                          tmpArray,
                                                          outArray)

        # delete the tmp-array
        del tmpArray

        # return the result
        return outArray.astype(numpy.float32)

    def run(self, in_name, out_name):
        """Code stolen from Multidrizzle.deriv()"""
        # store the names
        self.in_name = in_name
        self.out_name = out_name

        # make sure the input image exists
        if not os.path.isfile(self.in_name):

            # complain and out if not
            err_msg = "Image missing: %s!" % self.in_name
            raise aXeError(err_msg)

        # delete output name if existing
        if os.path.isfile(self.out_name):
            os.unlink(self.out_name)

        _log.info("Running quickDeriv on ", self.in_name)
        # OPEN THE INPUT IMAGE IN READ ONLY MODE
        img = fits.open(self.in_name, mode='readonly', memmap=0)

        # calling qderiv with the assumption that the
        # input file is a simple FITS file.
        absderiv = quickDeriv.qderiv(img["PRIMARY"].data)
        # absderiv = self._qderiv(img["PRIMARY"].data)

        # WRITE THE OUTPUT IMAGE TO A FITS FILE
        outfile = fits.open(self.out_name, 'append')
        outhdu = fits.PrimaryHDU(data=absderiv)
        outfile.append(outhdu)

        # CLOSE THE IMAGE FILES
        outfile.close()
        img.close()
        del outfile
        del img


class CRIdent:
    def __init__(self, drizzle_params):
        """Initializes the class. """
        self.driz_cr_scale = (float(drizzle_params['driz_cr_scale'].split()[0]),
                              float(drizzle_params['driz_cr_scale'].split()[1]))
        self.driz_cr_snr = (float(drizzle_params['driz_cr_snr'].split()[0]),
                            float(drizzle_params['driz_cr_snr'].split()[1]))
        self.driz_cr_grow = int(drizzle_params['driz_cr_grow'])
        self.driz_cr_ctegrow = 0

        # store the readout noise
        self.rdnoise = drizzle_params['RDNOISE']

    def _identify_crr(self, in_img, blot_img, blotder_img, exptime, sky_val):
        """Identify cosmic rays and other deviant pixels.

        The code was taken from muldidrizzle.DrizCR. Small adjustments and
        re-factoring was done.
        """

        # create an empty file
        __crMask = numpy.zeros(in_img.shape, dtype=numpy.uint8)

        # Part 1 of computation:
        # flag the central pixels
        # Create a temp array mask
        __t1 = numpy.absolute(in_img - blot_img)
        __ta = numpy.sqrt(numpy.absolute(blot_img * exptime
                                         + sky_val * exptime) +
                          self.rdnoise*self.rdnoise)
        __t2 = self.driz_cr_scale[0] * blotder_img + self.driz_cr_snr[0] * __ta / exptime
        __tmp1 = numpy.logical_not(numpy.greater(__t1, __t2))

        # mop up
        del __ta
        del __t1
        del __t2

        # Create a convolution kernel that is 3 x 3 of 1's
        __kernel = numpy.ones((3, 3), dtype=numpy.uint8)
        # Create an output tmp file the same size as the input temp mask array
        __tmp2 = numpy.zeros(__tmp1.shape, dtype=numpy.int16)
        # Convolve the mask with the kernel
        convolve.convolve2d(__tmp1,
                            __kernel,
                            output=__tmp2,
                            fft=0,
                            mode='nearest',
                            cval=0)
        del __kernel
        del __tmp1

        # Part 2 of computation
        # flag the neighboring pixels
        # Create the CR Mask
        __xt1 = numpy.absolute(in_img - blot_img)
        __xta = numpy.sqrt(numpy.absolute(blot_img * exptime +
                                          sky_val * exptime) +
                           self.rdnoise*self.rdnoise)
        __xt2 = self.driz_cr_scale[1] * blotder_img + self.driz_cr_snr[1] * __xta / exptime

        # It is necessary to use a bitwise 'and' to create the mask with numarray objects.
        __crMask = numpy.logical_not(numpy.greater(__xt1, __xt2) & numpy.less(__tmp2,9) )

        del __xta
        del __xt1
        del __xt2
        del __tmp2

        # Part 3 of computation - flag additional cte 'radial'
        # and 'tail' pixels surrounding CR pixels as CRs
        # In both the 'radial' and 'length' kernels below, 0->good and
        # 1->bad, so that upon
        # convolving the kernels with __crMask, the convolution
        # output will have low->bad and high->good
        # from which 2 new arrays are created having 0->bad and 1->good.
        # These 2 new arrays are then 'anded'
        # to create a new __crMask.

        # recast __crMask to int for manipulations below;
        # will recast to Bool at end
        __crMask_orig_bool = __crMask.copy()
        __crMask = __crMask_orig_bool.astype(numpy.int8)

        # make radial convolution kernel and convolve it with original __crMask
        # kernel for radial masking of CR pixel
        cr_grow_kernel = numpy.ones((self.driz_cr_grow, self.driz_cr_grow))
        cr_grow_kernel_conv = __crMask.copy()   # for output of convolution
        convolve.convolve2d(__crMask,
                            cr_grow_kernel,
                            output=cr_grow_kernel_conv)

        # make tail convolution kernel and convolve it with original __crMask
        cr_ctegrow_kernel = numpy.zeros((2*self.driz_cr_ctegrow+1,
                                         2*self.driz_cr_ctegrow+1))  # kernel for tail masking of CR pixel
        cr_ctegrow_kernel_conv = __crMask.copy()  # for output convolution

        # which pixels are masked by tail kernel depends on sign of
        # ctedir (i.e., readout direction):
        ctedir = 0
        if (ctedir == 1):  # HRC: amp C or D ; WFC: chip = sci,1 ; WFPC2
            cr_ctegrow_kernel[0:ctegrow, ctegrow] = 1  # 'positive' direction
        if (ctedir == -1):  # HRC: amp A or B ; WFC: chip = sci,2
            cr_ctegrow_kernel[ctegrow+1:2*ctegrow+1, ctegrow ] = 1    #'negative' direction
        if (ctedir == 0):  # NICMOS: no cte tail correction
            pass

        # do the convolution
        convolve.convolve2d(__crMask, cr_ctegrow_kernel, output = cr_ctegrow_kernel_conv)

        # select high pixels from both convolution outputs; then 'and' them to create new __crMask
        where_cr_grow_kernel_conv = numpy.where(cr_grow_kernel_conv < self.driz_cr_grow*self.driz_cr_grow,0,1 )        # radial
        where_cr_ctegrow_kernel_conv = numpy.where(cr_ctegrow_kernel_conv < self.driz_cr_ctegrow, 0, 1 )     # length
        __crMask = numpy.logical_and(where_cr_ctegrow_kernel_conv, where_cr_grow_kernel_conv) # combine masks

        __crMask = __crMask.astype(numpy.uint8)  # cast back to Bool

        del __crMask_orig_bool
        del cr_grow_kernel
        del cr_grow_kernel_conv
        del cr_ctegrow_kernel
        del cr_ctegrow_kernel_conv
        del where_cr_grow_kernel_conv
        del where_cr_ctegrow_kernel_conv

        # get back the result
        return __crMask

    def _createcrmaskfile(self, crName = None, crmask = None, header = None, in_imag=None):
        """
        Create a fits file containing the generated cosmic ray mask.
        """

        # migrate the data over
        _cr_file = numpy.zeros(in_imag.shape,numpy.uint8)
        _cr_file = numpy.where(crmask,1,0).astype(numpy.uint8)

        # rmove file if it exists
        if os.path.isfile(crName):
            os.unlink(crName)

        # Create the output file
        fitsobj = fits.HDUList()

        if (header is not None):
            del(header['NAXIS1'])
            del(header['NAXIS2'])
            if 'XTENSION' in header:
                del(header['XTENSION'])
            if 'EXTNAME' in header:
                del(header['EXTNAME'])
            if 'EXTVER' in header:
                del(header['EXTVER'])
            if 'NEXTEND' in header:
                header['NEXTEND'] = 0

            hdu = fits.PrimaryHDU(data=_cr_file, header=header)
            del hdu.header['PCOUNT']
            del hdu.header['GCOUNT']

        else:
            hdu = fits.PrimaryHDU(data=_cr_file)

        fitsobj.append(hdu)
        fitsobj.writeto(crName)

        # close the fits image
        fitsobj.close()

        # mop up
        del fitsobj
        del _cr_file

    def run(self, in_image, blot_image, blotder_image, exptime, sky_val, crr_image):
        """
        Do the identification
        """

        # open the input image
        inImage = fits.open(in_image, 'readonly')

        # open the blot image
        blotImage = fits.open(blot_image, 'readonly')

        # open the blot image
        blotDerImage = fits.open(blotder_image, 'readonly')

        # identify the CR's
        crr_data = self._identify_crr(inImage[0].data,
                                      blotImage[0].data,
                                      blotDerImage[0].data,
                                      exptime, sky_val)

        # save the image
        self._createcrmaskfile(crr_image,
                               crr_data,
                               inImage[0].header,
                               inImage[0].data)

        # delete the array
        del crr_data

        # close the images
        inImage.close()
        blotImage.close()
        blotDerImage.close()
