import sys
import shutil
import logging

from astropy.io import fits

from hstaxe import config as config_util

from . import axetasks
from .axeerror import aXeSIMError


# make sure there is a logger
_log = logging.getLogger(__name__)

"""
The following deal with axe simulations
"""

class DispImator:
    """Class to create a dispersed image"""
    def __init__(self, dummyImages, configfile, simobjects, lambda_psf=None,
                 model_spectra=None, model_images=None):
        """

        Parameters
        ----------
        dummyImages: DummyImages()
            dummy image structure
        configfile: str
            name of the aXe configuration file
        simobjects: str
            name of the model object table
        lambda_psf: float
            reference wavelength for the psf
        model_spectra: str
            name of the model spectra file
        model_images: str
            name of the model image file
        """
        # save the naked name of the grism image
        self.grismname = os.path.basename(dummyImages.griname)

        # check whether there is a direct image
        if dummyImages.dirname is not None:
            # save the direct image name
            self.dirname = os.path.basename(dummyImages.dirname)
        else:
            # set the direct image name to 'None'
            self.dirname = None

        # save all other input to class variables
        self.configfile = configfile
        self.iolname = simobjects
        self.model_spectra = model_spectra
        self.model_images = model_images
        self.lambda_psf = lambda_psf

        # check whether model images are given
        # append the file name to the list
        if self.model_images is not None:
            self.cont_model = 'direct'
        else:
            self.cont_model = 'gauss'

    def run(self, silent=True):
        """Generates a simulated dispersed image

        The method executes the series of aXe tasks necessary to generate
        a simulated dispersed image. The input from the class data is
        supplemented with default values.

        Parameters
        ----------
        silent: boolean for silent mode
        """

        # run SEX2GOL
        _log.info('Running task "sex2gol" ...', end=' ')
        sys.stdout.flush()
        axetasks.sex2gol(grism=self.grismname,
                         config=self.configfile,
                         in_sex=self.iolname,
                         use_direct=True,
                         direct=self.dirname,
                         silent=silent)
        _log.info(' Done')

        # run GOL2AF
        _log.info('Running task "gol2af"  ...')
        axetasks.gol2af(grism=self.grismname,
                        config=self.configfile,
                        orient=1,
                        slitless_geom=1,
                        silent=silent)
        _log.info(' Done')

        # run PETCONT
        _log.info('Running task "petcont" ...', end=' ')
        sys.stdout.flush()
        axetasks.petcont(grism=self.grismname,
                         config=self.configfile,
                         cont_model=self.cont_model,
                         spec_models=self.model_spectra,
                         object_models=self.model_images,
                         lambda_psf=self.lambda_psf,
                         no_pet=True,
                         silent=silent)
        _log.info(' Done')

    def mopup(self):
        """Deleting GOL and OAF files"""

        # get the root name of the dispersed image
        pos = self.grismname.rfind('.fits')
        root_name = self.grismname[:pos]

        # delete the GOL, the OAF and the PET
        result_cat = config_util.getOUTPUT(root_name + '_2.cat')
        if os.path.isfile(result_cat):
            os.unlink(result_cat)
        result_oaf = config_util.getOUTPUT(root_name + '_2.OAF')
        if os.path.isfile(result_oaf):
            os.unlink(result_oaf)

class DirImator:
    """Class to create a direct image"""
    def __init__(self, dummyImages, configfile, simobjects, tpass_direct,
                 model_spectra=None, model_images=None, tel_area=None):
        """

        Parameters
        ----------
        dummyImages: dDummyImages()
            dummy image structure
        configfile: str
            name of the aXe configuration file
        simobjects: str
            name of the model object table
        tpass_direct: str
            name of the total passband file
        model_spectra: str
            name of the model spectra file
        model_images: str
            name of the model image file
        tel_area: float
            the collecting area of the telescope
        """
        # save the naked name of the direct image
        self.dirname = os.path.basename(dummyImages.dirname)

        # save all other input to local variables
        self.configfile = configfile
        self.iolname = simobjects
        self.tpass_direct = tpass_direct
        self.model_spectra = model_spectra
        self.model_images = model_images
        self.tel_area = tel_area

    def run(self, silent=False):
        """Generates a simulated direct image

        The method executes the series of aXe tasks necessary to generate
        a simulated direct image.

        Parameters
        ----------
        silent: bool
            for silent mode
        """
        # run SEX2GOL
        _log.info('Running task "sex2gol"  ...', end=' ')
        sys.stdout.flush()
        axetasks.sex2gol(grism=self.dirname,
                         config=self.configfile,
                         in_sex=self.iolname,
                         use_direct=False,
                         silent=silent)
        _log.info(' Done')

        # run GOL2AF
        _log.info('Running task "gol2af"   ...')
        sys.stdout.flush()
        axetasks.gol2af(grism=self.dirname,
                        config=self.configfile,
                        silent=silent)
        _log.info(' Done')

        # run DIRIMAGE
        _log.info('Running task "dirimage" ...')
        sys.stdout.flush()
        axetasks.axedirim(dirname=self.dirname,
                          config=self.configfile,
                          tpass_direct=self.tpass_direct,
                          model_spectra=self.model_spectra,
                          model_images=self.model_images,
                          tel_area=self.tel_area,
                          silent=silent)
        _log.info(' Done')

    def mopup(self):
        """Deleting GOL and OAF files"""

        # get the root name of the dispersed image
        pos = self.dirname.rfind('.fits')
        root_name = self.dirname[:pos]

        # delete the GOL, the OAF and the PET
        result_cat = config_util.getOUTPUT(root_name + '_2.cat')
        if os.path.isfile(result_cat):
            os.unlink(result_cat)
        result_oaf = config_util.getOUTPUT(root_name + '_2.OAF')
        if os.path.isfile(result_oaf):
            os.unlink(result_oaf)


class DummyExtractor:
    """Class to make a dummy extraction"""
    def __init__(self, dummyImages, grism_image, configfile, simobjects,
                 bck_flux, extrfwhm=3.0, orient=True, slitless_geom=True,
                 adj_sens=True, lambda_mark=None):
        """

        Parameters
        ----------
        dummyImages: DummyImages()
            dummy image structure
        grism_image: str
            the simulated grism image name
        configfile: str
            name of the aXe configuration file
        simobjects: str
            name of the model object table
        bck_flux: float or string
            backround-flux or image
        extrfwhm: float
            multiplier for extraction width
        orient: bool
            flag for tilted extraction
        slitless_geom: bool
            flag for slitless optimized extraction
        adj_sens: bool
            flag for adjusted flux conversion
        lambda_mark: float
            wavelength to apply cutoff magnitudes
        """
        # save the direct image name
        self.direct_image = os.path.basename(dummyImages.dirname)
        self.simul_grisim = os.path.basename(grism_image)
        self.dispersed_image = None

        # save all other input to local variables
        self.configfile = configfile
        self.iolname = simobjects
        self.bck_flux = bck_flux
        self.extrfwhm = extrfwhm
        self.orient = orient
        self.slitless_geom = slitless_geom
        self.adj_sens = adj_sens
        self.lambda_mark = lambda_mark

        # check whether everything
        # is where it is supposed to be
        self._check_files()

    def _check_files(self):
        """Checks the existence of the input files"""
        # check the direct image
        if not os.path.isfile(config_util.getDATA(self.direct_image)):
            err_msg = ("\nThe direct image is not available: {0:s}"
                       .format(config_util.getDATA(self.direct_image)))
            raise aXeSIMError(err_msg)

        # check the configuration file
        if not os.path.isfile(config_util.getCONF(self.configfile)):
            err_msg = ("\nThe configuration file is not available: {0:s}"
                       .format(config_util.getCONF(self.configfile)))
            raise aXeSIMError(err_msg)

        # check the simulated grism image
        if not os.path.isfile(config_util.getOUTSIM(self.simul_grisim)):
            err_msg = ("\nThe grism image is not available: {0:s}"
                       .format(config_util.getOUTSIM(self.simul_grisim)))
            raise aXeSIMError(err_msg)

        # check the IOL
        if not os.path.isfile(self.iolname):
            err_msg = ("\nThe Input Object List is not available: {0:s}"
                       .format(self.iolname))
            raise aXeSIMError(err_msg)

        try:
            float(self.bck_flux)
        except ValueError:
            # check the background image
            if not os.path.isfile(config_util.getCONF(self.bck_flux)):
                err_msg = ("\nThe background imagage is not available: {0:s}"
                           .format(config_util.getCONF(self.bck_flux)))
                raise aXeSIMError(err_msg)

    def _decorate_PET(self):
        """Write the 'CONTAM'-keyword into the PET

        The method determines the name of the PET and
        sets the contamination keyword in the zero-extension header
        """
        # get the root name of the dispersed image
        root_name = self.dispersed_image.split('.fits')[0]

        # compose the name of the PET
        result_pet = config_util.getOUTPUT(root_name + '_2.PET.fits')

        # open the PET
        pet_fits = fits.open(result_pet, mode='update')

        # update the PET header
        comment_str = 'dummy flag - no quantitative contamination'
        pet_fits[0].header['CONTAM'] = ('GEOM', comment_str)

        # close and out
        pet_fits.close()

    def prepare_extraction(self):
        """Prepares the aXe extraction

        The module does some preparatory stuff before the extraction
        can start. This includes copying the simulated dispersed image
        to AXE_IMAGE_PATH and subtracting the background on this copy.
        """
        # give brief feedback
        _log.info('Dummy extraction on the dispersed image:')
        sys.stdout.flush()

        # get a random filenames
        tmpfile1 = config_util.get_random_filename('t', '.fits')

        # copy the grism image to AXE_IMAGE_PATH
        shutil.copy(config_util.getOUTSIM(self.simul_grisim),
                    config_util.getDATA(tmpfile1))

        # subtract the background from
        # the grism image
        # expression = "(a - b)"
        # iraf.imexpr(expr=expression, output=tmpfile2,
        #    a=config_util.getDATA(tmpfile1)+'[SCI]', b=self.bck_flux, Stdout=1)

        in_image = fits.open(config_util.getDATA(tmpfile1))
        in_image['sci'].data -= self.bck_flux
        in_image.close()

        # store the name of the background
        # subtracted grism image - this was tmpfile
        self.dispersed_image = tmpfile1

    def mopup(self):
        """Deleting tmp-files, copying SPC's, STP's"""

        # get the root name of the dispersed image
        root_name = self.dispersed_image.split('.fits')[0]

        #  get the root name of the simulated image
        result_root = self.simul_grisim.split('.fits')[0]

        # move and rename the SPC-file
        out_spc = config_util.getOUTPUT(root_name + '_2.SPC.fits')
        result_spc = config_util.getOUTSIM(result_root + '_2.SPC.fits')
        shutil.move(out_spc, result_spc)

        # move and rename the STP-file
        out_stp = config_util.getOUTPUT(root_name + '_2.STP.fits')
        result_stp = config_util.getOUTSIM(result_root + '_2.STP.fits')
        shutil.move(out_stp, result_stp)

        # delete the background subtracted
        # grism image
        os.unlink(config_util.getDATA(self.dispersed_image))

        # delete the GOL, the OAF and the PET
        result_cat = config_util.getOUTPUT(root_name + '_2.cat')
        if os.path.isfile(result_cat):
            os.unlink(result_cat)
        result_oaf = config_util.getOUTPUT(root_name + '_2.OAF')
        if os.path.isfile(result_oaf):
            os.unlink(result_oaf)
        result_pet = config_util.getOUTPUT(root_name + '_2.PET.fits')
        if os.path.isfile(result_pet):
            os.unlink(result_pet)

    def run(self, silent=True):
        """Generates a simulated dispersed image

        The method executes the series of aXe tasks necessary to generate
        a simulated dispersed image. The input from the class data is
        supplemented with default values.

        Parameters
        ----------
        silent: bool
            boolean for silent mode
        """
        # run SEX2GOL
        _log.info('Running task "sex2gol" ...', end=' ')
        sys.stdout.flush()
        axetasks.sex2gol(grism=self.dispersed_image,
                         config=self.configfile,
                         in_sex=self.iolname,
                         use_direct=True,
                         direct=self.direct_image,
                         silent=silent)
        _log.info(' Done')

        # run GOL2AF
        _log.info('Running task "gol2af"  ...', end=' ')
        sys.stdout.flush()
        axetasks.gol2af(grism=self.dispersed_image,
                        config=self.configfile,
                        mfwhm=self.extrfwhm,
                        orient=self.orient,
                        slitless_geom=self.slitless_geom,
                        lambda_mark=self.lambda_mark,
                        ilent=silent)
        _log.info(' Done')

        # run AF2PET
        _log.info('Running task "af2pet"  ...', end=' ')
        sys.stdout.flush()
        axetasks.af2pet(grism=self.dispersed_image,
                        config=self.configfile,
                        silent=silent)
        _log.info(' Done')

        # -----------------------------------------------
        # set the contamination keyword
        #
        # Running PECONT, is, since we are doing
        # a simulation, not very reasonable.
        # However PET2SPC complains if the contmaintion
        # keyword is not set. A solution is just to set
        # the geometrical contamination keyword to make
        # the warning in PET2SPC disappear.
        self._decorate_PET()
        # -----------------------------------------------

        # run PET2SPC
        _log.info('Running task "pet2spc" ...', end=' ')
        sys.stdout.flush()
        axetasks.pet2spc(grism=self.dispersed_image,
                         config=self.configfile,
                         adj_sens=self.adj_sens,
                         silent=silent)
        _log.info(' Done')

        # run STAMPS
        _log.info('Running task "stamps"  ...', end=' ')
        sys.stdout.flush()
        axetasks.stamps(grism=self.dispersed_image,
                        config=self.configfile,
                        sampling='rectified',
                        silent=silent)
        _log.info(' Done')
