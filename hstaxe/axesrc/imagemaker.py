import logging

from astropy.io import fits

from hstaxe import axeerror

from . import configfile

# make sure there is a logger
_log = logging.getLogger(__name__)

class DummyImages:
    """
    Small class to create dummy images for the simulations

    This class creates empty dispersed and direct images as dummies.
    Since all aXe C-executables need the dispersed image as an input,
    such an image must exist for making the simualtions.
    """
    def __init__(self, confname, griname=None, dirname=None, nx=None, ny=None):
        """
        Initializes the class

        @param confname: name of the aXe configuration file
        @type confname: string
        @param griname: name of the dispersed image
        @type griname: string
        @param dirname: name of the direct image
        @type dirname: string
        @param nx: image dimension in x
        @type nx: int
        @param ny: image dimension in y
        @type ny: int
        """

        # set the variable
        self.WCSimage = None
        self.WCSext = None

        # load the aXe configuration file
        self.conf = configfile.ConfigFile(confname)

        # load the pre-defined image information
        image_data = self._get_image_data(self.conf)

        # check whether a grism image
        # shall be created
        if image_data['grism'] != None and griname != None:
            self.griname  = griname
            self.gridata  = image_data['grism']
            self.WCSimage = griname
            self.WCSext   = '[SCI]'
        else:
            self.griname = None
            self.gridata = None

        # check whether a direct image
        # shall be created
        if dirname is not None:
            self.dirname = dirname
            self.WCSimage = dirname
            self.WCSext = '[SCI]'

            if image_data['direct'] is not None:
                self.dirdata = image_data['direct']
            else:
                self.dirdata = image_data['grism']
        else:
            self.dirname = None
            self.dirdata = None

        # store the drizzle metadata
        if 'drizzle' in image_data:
            self.drzdata = image_data['drizzle']
        else:
            self.drzdata = None

        # store the x-dimension
        if nx is not None:
            self.nx = nx
        else:
            self.nx = image_data['dimension'][0]

        # store the y-dimension
        if ny is not None:
            self.ny = ny
        else:
            self.ny = image_data['dimension'][1]

    def _get_image_data(self, conf):
        """
        Determines the pre-defined image information

        @param conf: aXe configuration object
        @type conf: <ConfigFile>

        @return: image information
        @rtype: {}
        """
        from . import WCSdata

        # load the WCS for the various grism modes
        if self.conf['CAMERA'] == 'HRC':
            # in case that several orders exist OR
            # the XOFFSET value is larger than -100, its HRG/G800L
            if self.conf['B'] != None or float(self.conf['A']['XOFF_'].split()[0]) > -100.0:
                image_data = WCSdata.get_HRC_G800L_WCS()
            else:
                image_data = WCSdata.get_HRC_PR200L_WCS()

        elif self.conf['CAMERA'] == 'SBC':
            # in case that the singularity is more than
            # -90 away from the reference point, its PR110L
            if float(self.conf['A'].disp[0][0]) < -90.0:
                image_data = WCSdata.get_SBC_PR110L_WCS()
            else:
                # otherwise, it is PR130L
                image_data = WCSdata.get_SBC_PR130L_WCS()

        # the WFC is fortunately unique
        elif self.conf['CAMERA'] == 'WFC':
            image_data = WCSdata.get_WFC_G800L_WCS()

        # check whether the mode is WFC3/IR
        elif self.conf['INSTRUMENT'] == 'WFC3' and self.conf['CAMERA'] == 'IR':
            image_data = WCSdata.get_WFC3_IR_G102_WCS()

        # check whether the mode is WFC3/IR
        elif self.conf['INSTRUMENT'] == 'WFC3' and self.conf['CAMERA'] == 'UV':
            image_data = WCSdata.get_WFC3_UV_G280_WCS()

        # check whether the mode is NICMOS/G141
        elif self.conf['INSTRUMENT'] == 'NICMOS' and self.conf['CAMERA'] == 'NIC3':
            image_data = WCSdata.get_NICMOS3_G141_WCS()

        else:
            # HRC/G800L is the dummy
            image_data = WCSdata.get_HRC_G800L_WCS()

        return image_data

    def deleteImages(self):
        """
        Deletes all images

        The method deletes the dummy images of the class instance.
        """
        # check whether there should exist a grism image
        if self.griname is not None:
            # check whether it exists and delete it
            if os.path.isfile(self.griname):
                os.unlink(self.griname)

        # check whether there should exist a direct image
        if self.dirname is not None:
            # check whether it exists and delete it
            if os.path.isfile(self.dirname):
                os.unlink(self.dirname)

    def makeImages(self):
        """
        Makes all images

        The method lets all images be generated.
        """
        # check whether a grism image
        # shall be created
        if self.griname is not None:
            # make the grism image
            self.makeOneImage(self.griname, self.nx, self.ny, self.gridata, self.drzdata)

        # check whether a direct image
        # shall be created
        if self.dirname is not None:
            # make the direct image
            self.makeOneImage(self.dirname, self.nx, self.ny, self.dirdata, self.drzdata)

    def makeOneImage(self, imgname, nx, ny, metadata, drzmeta=None):
        """
        Creates one dummy image

        The method creates a dummy image with the given name and dimension.
        A list of general metadata and a list of drizzle metadata is added
        to the zero extension header of the image.

        @param imgname: name of the image to be created
        @type imgname: string
        @param nx: image dimension in x
        @type nx: int
        @param ny: image dimension in y
        @type ny: int
        @param metadata: the list of general image metadata
        @type metadata: []
        @param drzmeta: the list of drizzle metadata
        @type drzmeta: []
        """
        from pyraf import iraf
        from iraf import noao, artdata

        # delete a previous version
        if os.path.isfile(imgname):
            os.unlink(imgname)

        # open a HDU-list
        mex_hdu = fits.HDUList()

        # create a primary HDU,
        # append it to the list
        hdrpr = fits.PrimaryHDU()
        mex_hdu.append(hdrpr)

        # go the the header and put
        # the exposure time
        hdr = mex_hdu[0].header
        hdr['EXPTIME'] = (1.0, 'dummy exposure time')

        if drzmeta != None:
            # update the header
            for item in drzmeta:
                hdr[item[0]] = (item[1], item[2])

        # write the image and close it
        mex_hdu.writeto(imgname)
        mex_hdu.close()

        # get a random filename
        tmpfile = get_random_filename('t', '.fits')

        # create a tmp-image with the right dimension
        iraf.mkpattern(input=tmpfile, output="", pattern="constant",
                       option="replace", v1=0.0, v2=0.0, size=1,
                       title="aXeSIM simulation", pixtype="real", ndim=2, ncols=nx,
                       nlines=ny, n3=1, n4=1, n5=1, n6=1, n7=1, header="")

        # copy the tmp-image to the empty dummy image
        # as science extension
        iraf.imcopy(input=tmpfile,
                    output=(imgname+'[SCI,1,append]'),
                    verbose='YES', Stdout=1)

        # open the dummy image, go to the header
        img = fits.open(imgname, 'update')
        hdr = img[1].header

        # update the header
        for item in metadata:
            hdr[item[0]] = (item[1], item[2])

        # write to disk and close
        img.flush()
        img.close

        # delete the tmp-image
        if os.path.isfile(tmpfile):
            os.unlink(tmpfile)
