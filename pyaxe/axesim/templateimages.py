import os
import os.path
from astropy.io import fits
from hstaxe import axe_asciidata

from hstaxe.axeerror import aXeSIMError
from hstaxe.config import getDATA, getOUTSIM, getSIMDATA


class ArtImaList:
    """Class for the image list"""
    def __init__(self, inlist):
        """Initializer for the class

        Parameters
        ----------
        inlist: str
            name of list with the fits images
        """
        # initialize the list of spectra
        self._imalist = []

        # load the list with input images
        ilist = axe_asciidata.open(inlist)

        # go over the image names
        for item in ilist[0]:

            # check whether the file is there
            if not os.path.isfile(getSIMDATA(item)):
                error_message = ("\nDid not find image: {0:s} !!"
                                 .format(str(getSIMDATA(item))))
                raise aXeSIMError(error_message)

            # apend the fits to the liast
            self._imalist.append(ArtImage(getSIMDATA(item.strip())))

    def tofits(self, fitsname, indata_copy=0):
        """
        Converts and stores the spectra in a fits file

        Converts all images stored in the class instance to a
        multi-extension fits image with a name gvien as input.
        A flagg indicates whether, besides the normal output to
        AXE_OUTSIM_PATH, a copy to AXE_IMAGE_PATH is desired.

        Parameters
        ----------
        fitsname: str
            name for the MEX fits file
        indata_copy: int
            flag to save also a copy
        """
        # create a HDU list
        hdulist = fits.HDUList()

        # create an empty primary HDU
        phdu = fits.PrimaryHDU()

        # put the primary to the list
        hdulist.append(phdu)

        # give a linefeed
        print()

        # check whther the fitsname
        # ends with '.fits'
        if ('.fits' not in fitsname[-5:]):
            fitsname += '.fits'

        # initialize an index
        index = 0

        # go over all spectra
        for ima in self._imalist:
            # enhance the counter
            index += 1

            # append the the image HDU
            # to the output fits
            hdulist.append(ima.imgHDU)

            # print what you do on the screen
            print("Adding: {0:s} to {1:s}, ext: {2:s}"
                  .format(os.path.basename(ima.filename),
                          fitsname,
                          str(index)))

        # delete older versions
        # of the fits name
        if os.path.isfile(getOUTSIM(fitsname)):
            os.unlink(getOUTSIM(fitsname))

        print("\nWriting images to file: {0:s} ..."
              .format(getOUTSIM(fitsname)))
        # write it to fits
        hdulist.writeto(getOUTSIM(fitsname))
        # give an end notice
        print('Done')

        # check whether a copy
        # is needed at AXE_IMAGE_PATH directory
        if indata_copy:
            # delete older versions
            # of the fits name
            if os.path.isfile(getDATA(fitsname)):
                os.unlink(getDATA(fitsname))

            print("Writing images to file: {0:s} ..."
                  .format(getDATA(fitsname)))
            # write it to fits
            hdulist.writeto(getDATA(fitsname))
            # give an end notice
            print('Done')

        # add an extra linefeed
        print('')


class ArtImage:
    """Class for one image template"""
    def __init__(self, filename):
        """Initializer for the ArtImage class

        Parameters
        ----------
        filename: str
            name of the spectrum
        """
        # store the filename
        self.filename = filename

        # define the image name,
        # which is the file name until
        # the last '.', e.g. 'psf.fits' --> 'psf'
        self.imaname = os.path.basename(filename[:filename.rfind('.')])

        # create an image HDU from the fits
        self.imgHDU = self._makeImgHDU(self.filename, self.imaname)

    def _makeImgHDU(self, filename, imaname):
        """Extract and return the first non-empty image HDU from the fits

        The method opens a fits image and goes along its extension.
        The first extension with a non-zero data part is returned
        after updating the header with a given image name.

        Parameters
        ----------
        filename: str
            name of the fits file
        imaname: str
            name of the image

        Returns
        -------
        hdr: fits.HDU
            the image hdu
        """
        # intialize the image HDU
        imgHDU = None

        # open the fits
        fitsHDU = fits.open(filename, 'update')

        # initialize an index
        index = 0

        # go over all HDU's in the fits
        for HDU in fitsHDU:

            # check for an existing data HDU
            if (HDU.data is not None):

                # write the image name to a keyword
                HDU.header['IMANAME'] = imaname

                # normalize the data
                normData = self._norm_data(HDU.data)

                # create an image HDU
                imgHDU = fits.ImageHDU(normData,
                                       header=HDU.header,
                                       name=imaname)
                print("File {0:s}, ext: {1:s} loaded".format(filename,
                                                             str(index)))

                # exit after taking the first HDU with data
                break

        # close the fits file
        fitsHDU.close()

        # return the image HDU
        return imgHDU

    def _norm_data(self, imgdata):
        """Normalize the image data

        The method normalizes the image data. It uses
        methods of the data class (numpy or numarray).

        Parameters
        ----------
        imgdata: ImgData
            the image data

        Returns
        -------
        imgdata: ImgData
            the normalized image data
        """
        # get the sum of all pixels
        cts_sum = imgdata.sum()

        # create the  normalized image data
        normData = imgdata / cts_sum

        # return the normalized image data
        return normData
