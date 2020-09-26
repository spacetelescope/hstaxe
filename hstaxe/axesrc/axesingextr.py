import logging
from hstaxe import config as config_util

from . import axetasks
from . import configfile
from . import nlincoeffs


# make sure there is a logger
_log = logging.getLogger(__name__)

class aXeSpcExtr:
    def __init__(self, grisim, objcat, dirim, config, dmag, **params):
        # store the input
        self.grisim = grisim
        self.objcat = objcat
        self.dirim = dirim
        self.config = config
        self.dmag = dmag

        self.params = params

    def _make_bckPET(self):
        """Generate the object PET."""
        # run GOL2AF
        axetasks.gol2af(grism=self.grisim,
                        config=self.config,
                        mfwhm=self.params['backfwhm'],
                        back=True,
                        orient=self.params['orient'],
                        slitless_geom=self.params['slitless_geom'],
                        exclude=self.params['exclude'],
                        lambda_mark=self.params['lambda_mark'],
                        dmag=self.dmag,
                        out_af=None,
                        in_gol=None)

        #  run BACKEST
        axetasks.backest(grism=self.grisim,
                         config=self.config,
                         np=self.params['np'],
                         interp=self.params['interp'],
                         niter_med=self.params['niter_med'],
                         niter_fit=self.params['niter_fit'],
                         kappa=self.params['kappa'],
                         smooth_length=self.params['smooth_length'],
                         smooth_fwhm=self.params['smooth_fwhm'],
                         old_bck=False,
                         mask=False,
                         in_af="",
                         out_bck=None)

        # run AF2PET
        axetasks.af2pet(grism=self.grisim,
                        config=self.config,
                        back=True,
                        out_pet=None)

        # run PETFF
        axetasks.petff(grism=self.grisim,
                       config=self.config,
                       back=True,
                       ffname=None)

    def _make_objPET(self):
        """Generate the object PET."""

        # set the use_direct flag
        # NOTE: the flag is only usefull
        #       for the C-version, the python
        #       version does not need it!

        use_direct = False
        if self.dirim is not None:
            use_direct = True

        # run SEX2GOL
        axetasks.sex2gol(grism=self.grisim,
                         config=config_util.getCONF(self.config),
                         in_sex=config_util.getDATA(self.objcat),
                         use_direct=use_direct,
                         direct=config_util.getDATA(self.dirim),
                         dir_hdu=None,
                         spec_hdu=None,
                         out_sex=None)

        # run GOL2AF
        axetasks.gol2af(grism=self.grisim,
                        config=self.config,
                        mfwhm=self.params['extrfwhm'],
                        back=False,
                        orient=self.params['orient'],
                        slitless_geom=self.params['slitless_geom'],
                        exclude=self.params['exclude'],
                        lambda_mark=self.params['lambda_mark'],
                        dmag=self.dmag,
                        out_af=None, in_gol=None)

        # run AF2PET
        axetasks.af2pet(grism=self.grisim,
                        config=self.config,
                        back=False,
                        out_pet=None)

        # run PETCONT
        axetasks.petcont(grism=self.grisim,
                         config=self.config,
                         cont_model=self.params['cont_model'],
                         model_scale=self.params['model_scale'],
                         spec_models=None,
                         object_models=None,
                         inter_type=self.params['inter_type'],
                         lambda_psf=self.params['lambda_psf'],
                         cont_map=True,
                         in_af="")

        # run PETFF
        axetasks.petff(grism=self.grisim,
                       config=self.config,
                       back=False,
                       ffname=None)

    def _make_spectra(self):
        """Extract the spectra."""
        # set the switch for using
        # a background PET
        if (('back' in self.params) and (self.params['back'])):
            use_bpet = True
        else:
            use_bpet = False

        # run PET2SPC
        axetasks.pet2spc(grism=self.grisim,
                         config=self.config,
                         use_bpet=use_bpet,
                         adj_sens=self.params['adj_sens'],
                         weights=self.params['weights'],
                         do_flux=True,
                         drzpath=False,
                         in_af="",
                         opet=None,
                         bpet=None,
                         out_spc=None)

        # run STAMPS
        axetasks.stamps(grism=self.grisim,
                        config=self.config,
                        sampling=self.params['sampling'],
                        drzpath=False,
                        in_af="",
                        in_pet=None,
                        out_stp=None)

    def _make_drzgeocont(self, ext_info):
        # get the aXe names
        axe_names = config_util.get_axe_names(self.grisim, ext_info)

        # for the name of a special contamination OAF
        cont_oaf = config_util.getOUTPUT(axe_names['OAF'].replace('.OAF', '_{0:s}.OAF'.format(int(self.params['drzfwhm']*10.0))))

        # run GOL2AF,
        # getting the special OAF as output
        axetasks.gol2af(grism=self.grisim,
                        config=self.config,
                        mfwhm=self.params['drzfwhm'],
                        back=False,
                        orient=self.params['orient'],
                        slitless_geom=self.params['slitless_geom'],
                        exclude=self.params['exclude'],
                        lambda_mark=self.params['lambda_mark'],
                        dmag=self.dmag,
                        out_af=cont_oaf, in_gol=None)

        # run PETCONT,
        # using the special OAF as input
        axetasks.petcont(grism=self.grisim,
                         config=self.config,
                         cont_model=self.params['cont_model'],
                         model_scale=self.params['model_scale'],
                         spec_models=None,
                         object_models=None,
                         inter_type=self.params['inter_type'],
                         lambda_psf=self.params['lambda_psf'],
                         cont_map=True, in_af=cont_oaf)

    def run(self):
        # load the configuration files;
        # get the extension info
        conf = configfile.ConfigFile(config_util.getCONF(self.config))
        ext_info = config_util.get_ext_info(config_util.getDATA(self.grisim), conf)
        del conf

        # Does this harm data that was astrodrizzled?
        if (('drzfwhm' in self.params) and
            (self.params['drzfwhm']) or
            (('cont_model' in self.params) and
            (config_util.is_quant_contam(self.params['cont_model'])))):

            # generate the non-linear distortions from the IDCTAB;
            # and store them in the fits-file header
            this_data = config_util.getDATA(self.grisim)
            _log.info("Generating and storing nonlinear distortions in {0}".format(this_data))
            nlins = nlincoeffs.NonLinCoeffs(this_data, ext_info)
            nlins.write_file()
            nlins.store_coeffs()
            del nlins

        # make the object PET's
        self._make_objPET()

        # make a background PET if necessary
        if 'back' in self.params and self.params['back']:
            _log.info("\nMaking backpet\n")
            self._make_bckPET()

        # extract the spectra
        if 'spectr' in self.params and self.params['spectr']:
            _log.info("\nMaking spectra\n")
            self._make_spectra()

        # make the proper non-quantitative contamination
        if ('drzfwhm' in self.params and self.params['drzfwhm']) and \
           ('cont_model' in self.params and not config_util.is_quant_contam(self.params['cont_model'])):
            _log.info("\nmaking non quant contam\n")
            self._make_drzgeocont(ext_info)
