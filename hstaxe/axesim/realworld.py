from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import sys
from copy import deepcopy
import numpy as np
from astropy.io import fits

from ..axeerror import aXeSIMError
from ..axeutils import get_random_filename


class RealWorld(object):
    """Class to add poisson noise"""
    def __init__(self, image_name, extname='0', exptime=1.0, bck_flux=0.0,
                 rdnoise=0.0, instrument=None):
        """Initializes the RealWorld class

        The various input for the noise model is passed as
        parameters. Reasonable defaults are defined as well.
        'None' as input is converted to the corresponding default.

        Parameters
        ----------
        image_name: str
            the name of the input image
        extname: str
            the extension to use
        exptime: float
            he exposure time in seconds
        bck_flux: float
            the background flux in electrons
        rdnoise: float
            the readout noise in electrons
        """
        # check whether the image exists
        if not os.path.isfile(image_name):
            err_msg = ("\nImage: {0:s} does not exist!".format(image_name))
            raise aXeSIMError(err_msg)

        # save the parameters
        self.image_name = image_name
        self.extname = extname
        self.exptime = exptime

        # determine and store the image dimensions
        self.dimension = self._get_dimension(self.image_name, self.extname)

        # store the parameter for
        # the default background;
        # convert 'None' to default
        if (bck_flux is not None):
            self.bck_flux = bck_flux
        else:
            self.bck_flux = 0.0

        # store the parameter for
        # the CR-boolean
        # convert 'None' to default
        if (rdnoise is not None):
            self.rdnoise = rdnoise
        else:
            self.rdnoise = 0.0

        # store the parameter for
        # instrument keyword
        if (instrument is not None):
            self.instrument = instrument
        else:
            self.instrument = 'aXeSIM'

    def _set_keywords(self):
        """Set header kewords in output image

        The method sets header keywords in the zero extension
        header of the output image.
        """

        # open the fits image
        inima = fits.open(self.image_name, 'update')

        # write the instrument name in the header
        inima[0].header['INSTRUME'] = (self.instrument, 'instrument name')

        # write the exposure time in the header;
        # 'None' is converted to 1.0
        if (self.exptime is None):
            inima[0].header['EXPTIME'] = (1.0, 'default exposure time')
        else:
            inima[0].header['EXPTIME'] = (self.exptime, 'exposure time')

        # close the image
        inima.close()

    def _compose_multiext_image(self, sci_ext, err_ext, dq_ext, del_input=1):
        """
        Compose the multi-extension image from individual layers

        Parameters
        ----------
        sci_ext: str
            the science extension image
        err_ext: str
            the error extension image
        dq_ext: str
            the dq extension image
        """
        # set the various keywords
        self._set_keywords()

        # put the noisy stuff on the input image,
        # replacing the original first extension
        current_image = fits.open(self.image_name, mode='update')
        current_image['SCI'].data = sci_ext
        current_image['ERR'].data = err_ext
        current_image['DQ'].data = dq_ext
        current_image.close()

        # delete the input images,
        # if reuqested
        if del_input:
            os.unlink(sci_ext)
            os.unlink(err_ext)
            os.unlink(dq_ext)

    def _get_dimension(self, image_name, extname):
        """Get the image dimension

        Parameters
        ----------
        image_name: str
            the name of the image
        extname: str
            the extension name

        Returns
        -------
        dimension: (int, int)
            the image-dimension (yaxis, xaxis)
        """
        # open the fits
        f_img = fits.open(image_name)

        # extract the image dimension
        dimension = f_img[extname].data.shape

        # close the image
        f_img.close()

        # return the dimension
        return dimension

    def _make_real_sciimage(self, in_image, bck_flux, exptime):
        """Create the science extension

        Starting from a simulated image in [e/s], the module adds
        background and scales to [e]. The output image is returned.

        Parameters
        ----------
        in_image: str
            the simulated image
        bck_flux: float
            background flux [e/s]
        exptime: float
            the exposure time

        Returns
        -------
        tmpfile1: str
            The name of the output file

        """
        # get a random filename
        tmpfile1 = get_random_filename('t', '.fits')

        # add background; scale by exptime
        image = fits.open(in_image)
        image.data = (image.data + bck_flux) * exptime
        image.writeto(tmpfile1)
        image.close()

        # return the image name
        return tmpfile1

    def _add_noise(self, in_image):
        """Add noise to an image

        The module adds noise to an image in [e].

        Parameters
        ----------
        in_image: str
            the input image, assumed to be an MEF fits image with data
            in self.extname

        Returns
        -------
        in_image is updated inplace

        Notes
        -----

        readnoise is in electrons
        gain is used for scaling the readnoise parameter

        Translated from the iraf.mknoise task, Poisson noise is added as:
            out = P((in+background)*gain) / gain

        where P(x) is a poisson deviate with mean x,  in  and  out  are  the
        input  and  final  pixel  values,  and  background  and gain are the
        parameter values of the same name

        the readnoise values is used as the Gaussian sigma
        The  sigma is divided by the specified gain to convert to  image
        data  units.   Gaussian  random numbers  of mean zero are then
        generated for each pixel and added to the image, or background value
        for  new  images,  after  the  photon noise is computed.


        First,  the  gaussian approximation is  used  for  data  values
        greater  than  20  (after applying  the  background  and  gain).

        The square root of the data value is used as the gaussian  sigma
        about  the  data  value. For values  less  than  20  a  true
        poisson  deviate is generated.  The second speed up is to allow
        storing a number of normalized  gaussian values   given   by   the
        package  parameter  ranbuf  as  they  are generated.  If more values
        than this  are  desired  then  a  uniform random  number  is  used
        to select one of these stored values.  This applies  to  both  the
        read  noise  and  poisson   noise   gaussian approximation  though
        not  the  true  poisson evaluation.

        For most purposes this approximation is good and one would need to
        look  very hard  to  detect  the  nonrandomness  in the noise.
        However, if one wants to take the extra  computational  time  then
        by  setting  the ranbuf  parameter  to  zero  each  gaussian  random
        number  will be generated independently.


        This function duplicates the gaussian approximation, depending on
        runtime tests we could remove that
        """
        with fits.open(in_image, mode='update') as image:
            # background = 0.0
            # gain = 1.0
            new_image = deepcopy(image[self.extname].data)

            # select large values and use simple shot noise
            lv = np.where(new_image >= 20.)
            new_image[lv] += np.sqrt(new_image[lv])

            # create normal, where readnoise is the sigma and the pixel
            # value is the mean
            y, x = new_image.shape
            noise = self.rdnoise * np.random.randn(y, x)
            sv = np.where(new_image < 20.)
            new_image[sv] += noise[sv]
            image[self.extname].data = new_image
            image.close()

    def _compute_err_ext(self, sci_ext):
        """Compute the error image for a science image

        For a, image in [extension], the module computes the associated
        error image assuming a simple noise model with photon shot noise
        and readout noise.

        Parameters
        ----------
        sci_ext: str
            the input image

        Returns
        -------
        The name of the tempfile that's created
        """
        # get a random filename
        tmpfile1 = get_random_filename('t', '.fits')

        try:
            image = fits.open(self.image)
        except IOError:
            print("Problem opening image: {0:s}".format(self.image))

        image[sci_ext].data = np.sqrt(image[sci_ext].data + self.rdnoise**2)
        image.writeto(tmpfile1)
        image.close()

        return tmpfile1

    def _scale_image(self, in_image):
        """Converts from [e] to [e/s]

        @param in_image: the input image
        @type in_image: string
        """

        try:
            image = fits.open(in_image, mode='update')
        except IOError:
            print("Problem opening image: {0:s}".format(in_image))

        # divide by exposure time
        image.data /= self.exptime
        image.close()

    def _make_dq_extension(self):
        """Creates an empty dq-extension image

        Returns
        -------
        Saves the output to a new simple fits image and returns the name

        Notes
        -----
        iraf.mkpattern(input=tmpfile1, output="", pattern="constant",
                       option="replace",nlines=self.dimension[0],
                       v1=0.0, v2=0.0, size=1, title="", pixtype="integer",
                       ndim=2, ncols=self.dimension[1],
                       n3=1, n4=1, n5=1, n6=1, n7=1, header="")
        """
        # get a random filename
        tmpfile1 = get_random_filename('t', '.fits')
        image = np.zeros((self.dimension[1], self.dimension[0]), dtype=np.int)
        hdu = fits.PrimaryHDU(image)
        hdu.writeto(tmpfile1)
        return tmpfile1

    def _make_const_image(self, value):
        """Creates a constant image

        Notes
        -----
        iraf.mkpattern(input=tmpfile1, output="", pattern="constant",
                   option="replace",v1=value, v2=0.0, size=1, title="",
                   pixtype="real",ndim=2, ncols=self.dimension[1],
                   nlines=self.dimension[0],n3=1, n4=1, n5=1, n6=1, n7=1,
                   header="")
        """
        # get a random filename
        tmpfile1 = get_random_filename('t', '.fits')
        image = float(value) * np.ones((self.dimension[1], self.dimension[0]),
                                       dtype=np.float64)
        hdu = fits.PrimaryHDU(image)
        hdu.writeto(tmpfile1)
        return tmpfile1

    def make_real(self):
        """Create a 'natural' image

        Depending on the class data, method adds background and noise
        in order to create a 'natural' image from the plain,
        simulated image.
        """
        print("\nCompleting image ...")
        sys.stdout.flush()
        if (self.exptime is not None):
            self.make_real_exptime()
        else:
            self.make_real_noexp()
        print(" Done\n")

    def make_real_noexp(self):
        """Create a 'natural' image

        The method creates a natural image without noise
        """
        # get the science image plus background
        # in electrons
        imname = "{}[{}]".format(self.image_name, self.extname)

        sci_scaled = self._make_real_sciimage(imname, self.bck_flux, 1.0)
        # compute the error extension
        err_scaled = self._make_const_image(self.rdnoise)

        # make a dq-extension
        dq_ext = self._make_dq_extension()

        # compose the final multi-extension image
        self._compose_multiext_image(sci_scaled, err_scaled, dq_ext, 1)

    def make_real_exptime(self):
        """Create a 'natural' image

        The method creates a natural image with noise
        """
        # get the science image plus background
        # in electrons
        imname = "{}[{}]".format(self.image_name, self.extname)
        sci_scaled = self._make_real_sciimage(imname,
                                              self.bck_flux,
                                              self.exptime)

        # add readout and poisson errors
        self._add_noise(sci_scaled)

        # compute the error extension
        err_scaled = self._compute_err_ext(sci_scaled)

        # check whether the exposure time is non-zero
        if self.exptime != 0.0:
            # scale both extensions
            # by the exposure time
            self._scale_image(sci_scaled)
            self._scale_image(err_scaled)

        # make a dq-extension
        dq_ext = self._make_dq_extension()

        # compose the final multi-extension image
        self._compose_multiext_image(sci_scaled, err_scaled, dq_ext, 1)
