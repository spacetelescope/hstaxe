import os
import numpy as np
import math
import logging

from astropy.io import fits
from astropy.table import Table

from drizzlepac import astrodrizzle

from stsci.tools import fileutil
from stwcs.wcsutil import HSTWCS

from hstaxe.axeerror import aXeError
from hstaxe import config as config_util

# make sure there is a logger
_log = logging.getLogger(__name__)

class DrizzleImage:
    """Class for working with drizzled images"""
    def __init__(self, drizzle_image, ext=0):
        # store the image name
        self.image = drizzle_image

        # store the extension number
        self.ext = ext

        # determine and store the
        # number of drizzled images
        self._get_nfcube()

    def _get_nfcube(self):
        """Get the number of requested fluxcubes

        The method looks in the header of the grism image
        for the name of the input images used in the drizzled
        process. This number is saved to the object.

        """

        # open the fits file and get the header
        img = fits.open(self.image, 'readonly')
        head = img[self.ext].header

        # create the keyname for the first input image
        ID = 1
        keyname = 'D{0:03d}DATA'.format(ID)

        # create the keyname for subsequent input images
        # and continue until the keynames do not exist
        while keyname in head:
            ID = ID+1
            keyname = 'D{0:03d}DATA'.format(ID)

        # close the grism image
        img.close()

        # correct the number
        ID = ID-1

        # return the number
        self.ndrizzle = ID


class FluxCube:
    """The class for the fluxcube images.

    The creation of the fluxcube images is the main purpose of the module.
    """
    def __init__(self, grism_image, index):
        """Extracts the name of an input file from the header of a drizzled image

        The name of the corresponding fluxcube is derived and stored.
        Also all other drizzle information for the input image
        image is extracted and stored

        Parameters
        ----------
        grism_image: string
            the name of the drizzled grism image
        index: int
            the index number in the header of the mult. gr. im.

        the index is the input data file in the header of the grism image,
        a flux cube will eventually be made for each input image that went into
        the drizzled grism image which was specified.

        """
        self.grism_image_name = grism_image
        # open the image and get the header
        grism_head = fits.getheader(grism_image)

        # subsequently extract all information
        # on the particular input
        keyname = 'D{0:03d}DATA'.format(index)
        self.data = grism_head[keyname]

        keyname = 'D{0:03d}DEXP'.format(index)
        self.dexp = grism_head[keyname]

        keyname = 'D{0:03d}OUDA'.format(index)
        self.ouda = grism_head[keyname]

        keyname = 'D{0:03d}OUWE'.format(index)
        self.ouwe = grism_head[keyname]

        keyname = 'D{0:03d}OUCO'.format(index)
        self.ouco = grism_head[keyname]

        keyname = 'D{0:03d}MASK'.format(index)
        self.mask = grism_head[keyname]

        keyname = 'D{0:03d}WTSC'.format(index)
        self.wtsc = grism_head[keyname]

        keyname = 'D{0:03d}KERN'.format(index)
        self.kern = grism_head[keyname]

        keyname = 'D{0:03d}PIXF'.format(index)
        self.pixf = grism_head[keyname]

        keyname = 'D{0:03d}COEF'.format(index)
        self.coef = grism_head[keyname]

        keyname = 'D{0:03d}SCAL'.format(index)
        self.scal = grism_head[keyname]

        keyname = 'D{0:03d}FVAL'.format(index)
        self.fval = grism_head[keyname]

        # get information on the flux cube
        self.fcube_info = self._get_fcube_info(self.data)

        # get the fluxcube name
        self.fcube_name = self._get_fcubename(self.fcube_info)

        # get the name of the grism input image
        # self.inima_name = string.split(self.data, '.fits')[0] + '.fits'
        if not os.path.isfile(self.fcube_info['fits']):
            msg = ("FCUBEPREP: File {0:s} does not exist!"
                   .format(self.fcube_info['fits']))
            raise aXeError(msg)

        self.inima_dims = self._get_indims(self.fcube_info)

    def _get_fcube_info(self, header_name):
        """Get all info on the extension

        Parameters
        ----------
        header_name: str
            the header name

        Returns
        -------
        cube_info: dict
            information on the fluxcube
        """
        # create an empty dict
        fcube_info = {}

        # get the bracket information
        exte_str = header_name.split('.fits')[1]

        # get the inside of the brackets
        exte_data = exte_str[1:len(exte_str)-1]

        # collect fits name, extension name and extension
        # version in the dictionary
        fcube_info['root'] = header_name.split('.fits')[0]
        fcube_info['fits'] = header_name.split('.fits')[0] + '.fits'
        fcube_info['ext_nam'] = exte_data.split(',')[0]
        fcube_info['ext_ver'] = int(exte_data.split(',')[1])

        # return the extension information
        return fcube_info

    def _a_blot_image(self,
                      image_to_blot,
                      tempname,
                      x_excess,
                      y_excess,
                      interp):
        """
        Blot one image.

        Thats just a simple wrapper around the task blot in astrodrizzle

        Parameters
        ----------
        image_to_blot: str
            the input image name, either the grism or direct drizzled image

        tempname: str
            the name of the output blotted image

        """
        # the drizzle coeff information for adriz is taken
        # from the self.data image,
        excess_x = 0
        excess_y = 0

        # use the current data as reference image for output
        # self.data comes from the header of the input grism or flux image
        # and is one of the input images used to make the drizzled
        # image_to_blot
        input_image = (self.data).split("[")[0]
        bunit = fits.getval(image_to_blot, 'BUNIT')

        flt_header = fits.getheader(input_image)
        flt_wcs = HSTWCS(self.data)

        # now look at the image to blot, this is a drizzled image
        ftype = fileutil.isFits(image_to_blot)[1]

        if (ftype == 'mef'):
            blot_wcs = HSTWCS(image_to_blot,
                              ext=(str(self.fcube_info["ext_nam"]),
                                   str(self.fcube_info["ext_ver"])))
            image_data = fits.getdata(image_to_blot,
                                      extname=str(self.fcube_info["ext_nam"]),
                                      extver=int(self.fcube_info["ext_nam"]))

        elif (ftype is 'simple'):
            blot_wcs = HSTWCS(image_to_blot)  # assume simple
            image_data = fits.getdata(image_to_blot)

        else:
            return IOError("File type of fits image is not "
                           "supported {0:s}".format(image_to_blot))

        # edit the wcs header information to add any dim_info shifts that
        # we need, expanding the size of the output image
        # make sure this gets saved to the output extension header.
        # The lambda image will be bigger than the segment image
        if x_excess > 0:
            excess_x = int(flt_wcs.naxis1+x_excess*2.)
            flt_wcs.naxis1 = excess_x
            crpix = flt_wcs.wcs.crpix
            newx = int(crpix[0]) + x_excess
            flt_wcs.wcs.crpix = np.array([newx, int(crpix[1])])
            flt_wcs.sip.crpix[0] = newx

        if y_excess > 0:
            excess_y = int(flt_wcs.naxis2+y_excess*2.)
            flt_wcs.naxis2 = excess_y
            crpix = flt_wcs.wcs.crpix
            newy = int(crpix[1]) + y_excess
            flt_wcs.wcs.crpix = np.array([int(crpix[0]), newy])
            flt_wcs.sip.crpix[1] = newy

        # outimage is just the data array
        outimage = astrodrizzle.ablot.do_blot(image_data.astype(np.float32),
                                              blot_wcs,
                                              flt_wcs,
                                              1.,
                                              interp=interp,
                                              sinscl=1.,
                                              coeffs=True,
                                              wcsmap=None,
                                              stepsize=10)

        # update the flt_header with the flt_wcs information I created
        flt_header['CRPIX1'] = flt_wcs.wcs.crpix[0]
        flt_header['CRPIX2'] = flt_wcs.wcs.crpix[1]

        try:
            newimage = fits.PrimaryHDU()
            newimage.data = outimage
            newimage.header = flt_header
            newimage.header['BUNIT'] = bunit
            newimage.header.update(flt_wcs.to_header())
            newimage.verify('silentfix')
            newimage.writeto(tempname)
        except:
            raise IOError("Problem writing fits image {0:s}".format(tempname))

    def _a_blot_segment_image(self, image_to_blot, tempname, x_excess,
                              y_excess, interp):

        """Blot the segmentation or other nondrizzled image as if it were
        assume self.grism_image is always used as the source wcs reference
        Thats just a simple wrapper around the task blot in astrodrizzle

        Parameters
        ----------
        image_to_blot: str
            the input image name, either the grism or direct drizzled image
        tempname: str
            the name of the output blotted image

        Notes
        -----
        exposure time is hard coded to 1 since it was made from the
        drizzled image and we dont want the id numbers rescaled by the
        exposure time that was used for blotting

        """
        excess_x = 0
        excess_y = 0

        # the drizzle coeff information for adriz is taken
        # from the self.data image
        # use the current data as reference image
        input_image = (self.data).split("[")[0]
        flt_header = fits.getheader(input_image)
        flt_wcs = HSTWCS(self.data)

        # check to see if this is a simple fits or MEF and grab
        # the science information.
        ftype = fileutil.isFits(self.grism_image_name)[1]

        if (ftype is 'mef'):
            grism_wcs = HSTWCS(self.grism_image_name,
                               ext=(str(self.fcube_info["ext_nam"]),
                                    self.fcube_info["ext_ver"]))
        elif (ftype is 'simple'):
            grism_wcs = HSTWCS(self.grism_image_name)
        else:
            return IOError("File type of fits image is not "
                           "supported {0:s}".format(image_to_blot))

        ftype = fileutil.isFits(image_to_blot)[1]
        if (ftype is 'mef'):
            image_data = fits.getdata(image_to_blot,
                                      ext=(str(self.fcube_info["ext_nam"]),
                                           self.fcube_info["ext_ver"]))
        elif (ftype is 'simple'):
            image_data = fits.getdata(image_to_blot)
        else:
            return IOError("Input image is not a supported FITS "
                           "type: {0:s}".format(image_to_blot))

        # edit the wcs header information to add any dim_info shifts that we
        # need the segment image needs to be the same sky area cut without
        # the added pixels
        type(x_excess)
        if x_excess > 0:
            excess_x = int(flt_wcs.naxis1 + x_excess * 2.)
            flt_wcs.naxis1 = excess_x
            crpix = flt_wcs.wcs.crpix
            newx = int(crpix[0]) + x_excess
            flt_wcs.wcs.crpix = np.array([newx, crpix[1]])
            flt_wcs.sip.crpix[0] = newx

        if y_excess > 0:
            excess_y = int(flt_wcs.naxis2 + y_excess * 2.)
            flt_wcs.naxis2 = excess_y
            crpix = flt_wcs.wcs.crpix
            newy = int(crpix[1]) + y_excess
            flt_wcs.wcs.crpix = np.array([int(crpix[0]), newy])
            flt_wcs.sip.crpix[1] = newy

        # returns a numpy.ndarray which is just the data
        outimage = astrodrizzle.ablot.do_blot(image_data.astype(np.float32),
                                              grism_wcs,
                                              flt_wcs,
                                              1.,
                                              interp=interp,
                                              sinscl=1.,
                                              coeffs=True,
                                              wcsmap=None,
                                              stepsize=10)

        # update the flt_header with the flt_wcs information I created
        flt_header['CRPIX1'] = flt_wcs.wcs.crpix[0]
        flt_header['CRPIX2'] = flt_wcs.wcs.crpix[1]

        # if the input flt was an MEF we need to write an MEF out
        try:
            newimage = fits.PrimaryHDU()
            newimage.data = outimage
            newimage.header = flt_header
            newimage.header.update(flt_wcs.to_header())
            newimage.verify('silentfix')
            newimage.writeto(tempname)
        except:
            raise IOError("Problem writing fits image {0:s}".format(tempname))

    def _get_fcubename(self, fcube_info):
        """Get the name of the fluxcube

        Parameters
        ----------
        fcube_info: dict
            the fluxcube info

        Returns
        -------
        fcube_name: str
            the name of the fluxcube
        """
        fcube_name = "{0:s}_{1:d}.FLX.fits".format(fcube_info['root'],
                                                   3*fcube_info['ext_ver']-1)
        return fcube_name

    def _get_indims(self, fcube_info):
        """Get the dimensions of the image

        @param fcube_info: the fluxcube info
        @type fcube_info: dict

        @return: the image dimension
        @rtype: [int,int]
        """
        # make sure the fits image exists
        if not os.path.isfile(fcube_info['fits']):
            msg = "Image: {0:s} does not exist!".format(fcube_info['fits'])
            raise aXeError(msg)

        # open the image and get the header
        in_img = fits.open(fcube_info['fits'], 'readonly')
        in_head = in_img[fcube_info['ext_nam'], fcube_info['ext_ver']].header

        # extract the keywords for the image size from
        # the header
        dims = [in_head['NAXIS1'], in_head['NAXIS2']]

        # close the image
        in_img.close()

        # return the list with the image dimension
        return dims

    def _get_indims_old(self, image_name, extnum):
        """
        Parameters
        ----------
        image_name: str
            the image name
        extnum: int
            the extension number

        Returns
        -------
        dims: list
            a two-entry list with the x/y dimensions of the imag[sci,extnum]

        Notes
        -----
        The method determines the size of the science extension of an image.

        Returns
        -------
        dims: int
            the number of input images
        """

        if not os.path.isfile(image_name):
            err_msg = "Image: {0:s} does not exist!".format(self.grism_image)
            raise aXeError(err_msg)

        # open the image and get the header
        in_img = fits.open(image_name, 'readonly')
        in_head = in_img['SCI', extnum].header

        # extract the keywords for the image size from
        # the header
        dims = [in_head['NAXIS1'], in_head['NAXIS2']]

        # close the image
        in_img.close()

        # return the list with the image dimension
        return dims

    def create_fitscube(self, segm_image, filter_images, dim_info, interpol):
        """Creates one fitscube

        This method creates a fluxcube fits-image. The input is evaluated
        and then the various calls to the blot routine - one for the
        segmentation image and one for every flux image - are executed.
        Also the headers are filled according to the specifications.

        Parameters
        ----------
        segm_image: str
            the name of the segmentation image
        filter_images: list
            list of the triples fluximage,wavelength,ST_zero
        dim_info: list
            the list with the excess pixels
        interpol: str
            the interpolation method for the flux images
        """
        # look for the maximum excess in x and y
        x_excess = max([dim_info[0], dim_info[1]])
        y_excess = max([dim_info[2], dim_info[3]])

        # given that in the blot you add the maximum
        # requested pixels on each side, derive the
        # image section area string to get back to the requested
        # pixels on either side
        x_start = x_excess - dim_info[0] + 1
        y_start = y_excess - dim_info[2] + 1

        self.x_offs = -1 * dim_info[0]
        self.y_offs = -1 * dim_info[2]

        x_end = self.inima_dims[0] + 2*x_excess - (x_excess - dim_info[1])
        y_end = self.inima_dims[1] + 2*y_excess - (y_excess - dim_info[3])

        _log.info(f'Creating {self.fcube_name}')

        # delete a just existing, previous fluxcube image
        if os.path.isfile(self.fcube_name):
            os.unlink(self.fcube_name)

        # set up the fluxcube image, and store the primary
        # header; also store the keywords XOFFS/YOFFS in
        # the primary header
        mex_hdu = fits.HDUList()
        hdrpr = fits.PrimaryHDU()
        mex_hdu.append(hdrpr)
        hdr = mex_hdu[0].header
        hdr['XOFFS'] = (self.x_offs, 'X-OFFSET between flt and fluxcube')
        hdr['YOFFS'] = (self.y_offs, 'Y-OFFSET between flt and fluxcube')
        mex_hdu.writeto(self.fcube_name)

        # get a tmp-filename
        tmpname = config_util.get_random_filename('', '.fits')

        # blot the segmenation image
        self._a_blot_segment_image(segm_image,
                                   tmpname,
                                   x_excess,
                                   y_excess,
                                   'nearest')

        # copy the appropriate image section to the fluxcube
        tmp_fits = fits.open(tmpname, 'readonly')
        tmp_fits[0].header['EXTNAME'] = 'SEGM'
        tmp_fits[0].header['EXTVER'] = 1
        fits.append(self.fcube_name,
                    tmp_fits[0].data[y_start-1:y_end, x_start-1:x_end],
                    tmp_fits[0].header)
        tmp_fits.close()

        # delete the tmp-file
        os.unlink(tmpname)

        # go over all filter images
        for fimage in filter_images:

            # get the name of the fluximage and the
            # wavelength; prepare the keword entry for
            # the wavelength
            fluximg = fimage.get_fluxname()
            wavelength = fimage.get_wavelength()

            # get a tmp-filename
            tmpname = config_util.get_random_filename('', '.fits')

            _log.info(f"Using excess pixels of {x_excess}, {y_excess} ")
            self._a_blot_image(fluximg,
                               tmpname,
                               x_excess,
                               y_excess,
                               interpol)

            # store the wavelength in the header
            tmp_fits = fits.open(tmpname, 'update')
            tmp_fits[0].header['WAVELENG'] = (wavelength,
                                              'wavelength for the image')
            tmp_fits.close()

            # copy the appropriate section to te fluxcube
            tmp_fits = fits.open(tmpname, 'readonly')
            tmp_fits[0].header['EXTNAME'] = "LAMBDA"+str(int(wavelength))
            tmp_fits[0].header['EXTVER'] = 1
            fits.append(self.fcube_name, tmp_fits[0].data, tmp_fits[0].header)
            tmp_fits.close()

            # delete the tmp-file again
            os.unlink(tmpname)

        _log.info(' Done')


class FluxImage:
    """The class for the flux images.

    The flux images are the drizzled filter images, however transformed
    from counts per second to flux. An according segment of each
    flux image is stored in the fluxcubes.
    """
    def __init__(self, image_name='', wavelength=0., mag_zero=0., AB_input=0):
        """
        Parameters
        ----------
        image_name: str
            the name of the direct image
        wavelength: int
            the wavelength of the direct image
        mag_zero: float
            the zero point for the filter
        AB_zero: bool
            1 if the zeropoint is in AB, 0 if in ST


        Notes
        -----
        The input data is stored as class data.
        """
        self.image_name = image_name + '[SCI]'
        self.wavelength = wavelength
        if AB_input:
            self.st_zero = self._get_stmag_from_magab(mag_zero, wavelength)
        else:
            self.st_zero = mag_zero

        self.flux_name = self._get_fluxname(self.image_name)

    def _get_fluxname(self, image_name):
        """Get the name of the flux image

        The method derives the name of the flux image which
        is derived for each direct image.
        Example: gaga_drz.fits[SCI]  -->  gaga_drz.FLX.fits
        """
        return image_name.split('.fits')[0] + '.FLX.fits'

    def _get_flambda_from_magab(self, mag, wlength):
        """Compute f_lambda from mag_AB and lambda [nm]"""
        fnu = math.pow(10.0, -0.4 * (mag+48.6))
        if wlength != 0:
            flambda = 2.99792458e+16 * fnu / (wlength * wlength)
        else:
            flambda = 0.

        return flambda

    def _get_stmag_from_magab(self, mag, wlength):
        """Compute STmagzero from ABmagzero"""
        if wlength != 0:
            flambda = self._get_flambda_from_magab(mag, wlength)
            stmag = -2.5 * math.log10(flambda) - 21.10
        else:
            stmag = 0.

        return stmag

    def get_fluxname(self):
        """Get the flux image name

        The method delivers the name of
        the fluximage to the outside world.
        """
        return self.flux_name

    def get_wavelength(self):
        """Get the proper wavelength

        The method delivers the wavelength of
        the fluximage to the outside world.
        """
        return self.wavelength

    def transform_toflux(self, segm_image):
        """Transform an image from cps to flux

        The method creates a flux image for a direct drizzled
        image. The direct image comes directly from drizzled,
        its unit is [cts/s], however the fluxcubes need the information
        in ergs/cm^2/s/AA.

        Parameters
        ----------
        segm_image: str
            the segmentation for the data set

        Notes
        -----
        Formulas:
        STmag = -2.5*log10(image) + self.st_zero
        STmag = -2.5*log10(F_lam) -21.10
        F_lam = 10**(-0.4*(st_zero+21.10)) * image
        """
        # get the name of the flux image
        self.flux_name = self._get_fluxname(self.image_name)

        # delete the flux image if it already exists
        if os.path.isfile(self.flux_name):
            os.unlink(self.flux_name)

        # the expression to apply to the image data
        expon = 10.0 ** (-0.4 * (21.10 + self.st_zero))

        _log.info(f" image_name: {self.image_name.split('[')[0]}")
        _log.info(f" segm_image: {segm_image}")
        _log.info(f" flux_name: {self.flux_name}")

        # open the input files and apply the conversion
        file_a = fits.open(self.image_name.split('[')[0], 'readonly')
        file_b = fits.open(segm_image, 'readonly')
        ftype = fileutil.isFits(file_b.filename())[1]
        if (ftype == 'mef'):
            bgt1 = file_b[1].data >= 1  # assume in first
            blt1 = file_b[1].data < 1  # assume in first
        else:
            bgt1 = file_b[0].data >= 1
            blt1 = file_b[0].data < 1

        ftype = fileutil.isFits(file_a.filename())[1]

        if (ftype == 'mef'):
            file_a['SCI'].data[bgt1] *= expon
            file_a['SCI'].data[blt1] = 0.0
        else:
            file_a[0].data[bgt1] *= expon
            file_a[0].data[blt1] = 0.0

        # concatenate the primary and extension headers for output
        # to a simple FITS file
        prihdr = file_a[0].header

        del prihdr['nextend']

        if (ftype == 'mef'):
            scihdr = file_a['SCI'].header
            scihdr._strip()
            newhdr = fits.Header(prihdr + scihdr)
            data = file_a['SCI'].data
        else:
            newhdr = fits.Header(prihdr)
            data = file_a[0].data

        # write the new flux image in simple fits format
        # regardless of input format
        fits.writeto(self.flux_name, data, newhdr)
        file_a.close()
        file_b.close()

        _log.info(f"\n{self.image_name} --> {self.flux_name}")


class FluxCubeMaker:
    """Central class to take the input and to create the fluxcubes
    for the list of images extracted from the header of the
    drizzled grism image.
    """
    def __init__(self, grism_image, segm_image, filter_info, AB_zero,
                 dim_term, interpol):
        """
        The class data is set. While doing that, basic checks
        on the input is done. The existence of the images is
        checked, also the data type of the various real
        or integer numbers.

        Parameters
        ----------
        grism_image: str
            the name of the drizzled grism image
        segm_image: str
            the name of the segmentation image
        filter_info: list or string
            information on drizzled direct images
        dim_term: str
            description of the additional rows/column for the fluxcubes
        interpol: str
            the interpolation method for the flux images
        """
        self.filter_images = []
        self.dim_info = []
        self.fcube_list = []

        # check whether the grism image exist,
        # store the name if it exists
        if not os.path.isfile(grism_image):
            raise aXeError(f"File {grism_image} does not exist!")
        else:
            self.grism_image = grism_image

        # check whether the segmentation image exist,
        # store the name if it exists
        if not os.path.isfile(segm_image):
            raise aXeError("File: {0:s} does not exist!".format(segm_image))
        else:
            self.segm_image = segm_image

        # check whether the input for the flux images is a file
        if not os.path.isfile(filter_info):
            # split the input into its items
            finfo = filter_info.split(',')

            if len(finfo) != 3:
                msg = ("There must be 3 items in {0:s} not {1:d}!"
                       .format(filter_info, str(len(finfo))))
                raise aXeError(msg)

            # in case of a string, create a fluximage
            self.filter_images.append(self._make_fluxim(finfo[0],
                                      finfo[1],
                                      finfo[2],
                                      AB_zero))
        else:
            # in case of a file, let the information be extraced
            # and create the fluximages
            self.filter_images.extend(self._evaluate_finfolist(filter_info,
                                                               AB_zero))

        # resolve and get the dimension information
        self.dim_info = self._get_dimension_info(dim_term)

        # store the interpolation method
        self.interpol = interpol

    def _get_dimension_info(self, dimension_term):
        """Get the dimension information"""
        # initialize the array
        dim_info = []

        # check the dimension input
        dim_entries = dimension_term.split(',')

        # is the number of items correct?
        if len(dim_entries) != 4:
            msg = ("There must be 4 entries in the term: {0:s}, "
                   "not {1:d}!".format(dimension_term, len(dim_entries)))
            raise aXeError(msg)

        # check whether each item is an integer
        for item in dim_entries:
            if self._toInt(item.strip()) is None:
                err_msg = 'Item: ' + item + ' must be integer!'
                raise aXeError(err_msg)
            # store the item
            dim_info.append(self._toInt(item.strip()))

        # return the array
        return dim_info

    def _fill_fcubelist(self):
        """Makes a list of fluxcubes to be created

        The method determines the number of input images
        in the header of the drizzled grism images
        and then creates for each input image a fluxcube
        object
        """
        fcubes = []

        # determine the number of iput images
        n_fcubes = self._get_nfcube()

        # create a fluxcube object for each input image
        for index in range(1, n_fcubes+1):
            fcubes.append(FluxCube(self.grism_image, index))

        # return the fluxcube list
        return fcubes

    def _get_nfcube(self):
        """Get the number of fluxcubes to be created

        The method looks in the header of the grism image
        for the name of the input images used in the drizzle
        process. This number is returned
        """
        # open the fits file and get the header
        grism_img = fits.open(self.grism_image, 'readonly')
        grism_head = grism_img[0].header

        # create the keyname for the first input image
        ID = 1
        keyname = 'D{0:03d}DATA'.format(ID)

        # create the keyname for subsequent input images
        # and continue until the keynames do not exist
        while keyname in grism_head:
            ID = ID+1
            keyname = 'D{0:03d}DATA'.format(ID)

        # close the grism image
        grism_img.close()

        # correct the number
        ID = ID-1

        # return the number
        return ID

    def _toFloat(self, fcheck):
        """Check for float

        The module checks whether an expression is a
        float or not. The float representation of the
        expresion or None is returned
        """
        if isinstance(fcheck, float):
            return fcheck
        else:
            try:
                fconv = float(fcheck)
            except (TypeError, NameError, ValueError):
                fconv = None
            return fconv

    def _toInt(self, icheck):
        """Check for integer

        The module checks whether an expression is an
        integer or not. The integer representation of the
        expression or None is returned
        """
        if isinstance(icheck, int):
            return icheck
        else:
            try:
                iret = int(icheck)
            except (TypeError, NameError, ValueError):
                iret = None
            return iret

    def _make_fluxim(self, img_name, wav_pivot, zeropoint, AB_zero):
        """Generate a flux image"""
        # check whether the first item is the name of an existing file
        if not os.path.isfile(img_name):
            raise aXeError("File {0:s} does not exist!".format(img_name))

        # check whether the second item is a float
        if self._toFloat(wav_pivot) is None:
            raise aXeError("Expression {} must be a float!".format(wav_pivot))

        # check whether the third item is a float
        if self._toFloat(zeropoint) is None:
            err_msg = 'Expression: ' + zeropoint + ' must be float!'
            raise aXeError(err_msg)

        # create a new fluximage and return it
        return FluxImage(img_name,
                         float(wav_pivot),
                         float(zeropoint),
                         AB_zero)

    def _evaluate_finfolist(self, filter_file='', AB_zero=0.):
        """Evaluate the filter list.

        The method opens and extracts direct image information
        from the file given in the input. For each direct image
        specified there a flux image object is created, and
        the list of flux image objects is returned.

        Parameters
        ----------
        filter_file: str
            the name of the file with the direct image information
        AB_zero: str
            the zeropoint in AB

        Returns
        -------
        filter_images: list
            the list of filter image objects created
        """

        # make an empty list
        filter_images = []

        # open the file
        f_list = Table.read(filter_file, format='ascii.no_header')

        # go over all rows
        for row in f_list:

            # get image name, pivot wavelength and zeropoint
            img_name = row[0].strip()
            wav_pivot = row[1]
            zeropoint = row[2]

            filter_images.append(self._make_fluxim(img_name,
                                                   wav_pivot,
                                                   zeropoint,
                                                   AB_zero))

        # return the object list
        return filter_images

    def run(self):
        """Make all fluxcubes

        This method is responsible to actually create the
        fluxcube image. Other internal methods as well as
        methods of aother classes are successively
        called to create the fluxcubes associated to the
        grism images listed in the header of the
        drizzled grism image.
        """
        # create the list of fluxcube instances that
        # will be created
        self.fcube_list.extend(self._fill_fcubelist())

        # prepare the direct images:
        for fimage in self.filter_images:
            fimage.transform_toflux(self.segm_image)

        # adjust the image headers
        for fcube in self.fcube_list:
            res = fcube.create_fitscube(self.segm_image, self.filter_images,
                                        self.dim_info, self.interpol)
            if res:
                # something went wrong, delete the fcube
                if os.path.isfile(fcube.fcube_name):
                    os.unlink(fcube.fcube_name)
