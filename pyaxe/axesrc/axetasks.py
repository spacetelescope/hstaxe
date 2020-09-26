import logging

from . import axeinputs
from . import axelowlev
from . import axepreptor
from . import axesingextr
from . import dppdumps
from . import drizzleobjects
from . import fcubeobjs
from . import inputchecks
from . import iolmaking
from . import mefobjects
from . import pysex2gol

from hstaxe.config import axe_setup
from hstaxe.axeerror import aXeError

# make sure there is a logger
_log = logging.getLogger(__name__)

def iolprep(drizzle_image='',
            input_cat='',
            dimension_in='0,0,0,0'):
    """Convenience function for the aXe task IOLPREP.

    This task produces Input Object Lists for every input image of a
    Drizzle combined image. It projects the object positions from a
    master catalogue, which contains all objects in the coordinate system
    of the Drizzled image, out into the coordinate system of each
    input image. For each input image an Input Object List is generated,
    which contains only the objects within the boundaries of the given input image.

    Inputs
    ------
    drizzle_image : string
      The name of the drizzled image
    input_cat : string
      The name of the input SExtractor master catalog
    dimension_in : string
      The numbers of border pixels to extend the output image.
      [left, right, bottom, top] to the target area on the input
      images. E.g. 100,500,10,0 would include in the Input Object
      Lists all objects with -100 < x < x_size + 500 and
      -10 < y < y_size.

    Outputs
    -------
    Creates the catalog for each input grism image.

    Notes
    -----
    The names and drizzle parameters of the input images are retrieved from
    the header of the Astrodrizzle combined image. The projection of the object
    positions into the coordinate system of the input images is done by going
    through the WCS information. 

    There is a parameter to influence the sensitive area to include objects
    in the IOL's. This allows objects beyond the physical boundaries of the
    input image to be included in the IOL's, to take into account partly
    covered objects or to include bright objects outside of the FOV in
    contamination estimates derived from those IOL's.

    During the task execution, the drizzle coefficient files and the input
    images must be available. For this reason it would be best practice to
    run it in the directory that was used to combine the input images with
    Drizzle.

    """
    iol_maker = iolmaking.IOLMaker(drizzle_image,
                                   input_cat,
                                   dimension_in)
    iol_maker.run()


def fcubeprep(grism_image='',
              segm_image='',
              filter_info=None,
              AB_zero=True,
              dim_info='0,0,0,0',
              interpol='nearest'):
    """Convenience function for the aXe task FCUBEPREP"""

    # run the main command
    fcmaker = fcubeobjs.FluxCubeMaker(grism_image, segm_image, filter_info,
                                      AB_zero, dim_info, interpol)
    fcmaker.run()


def axeprep(inlist='',
            configs='',
            backgr=True,
            backims='',
            backped=None,
            mfwhm=None,
            norm=True,
            gcorr=False):
    """Convenience function for the aXe task AXEPREP.

    Inputs
    ------
    inlist: string
      Input Image List which gives on each line
        a) the name of the grism image to be processed (mandatory)
        b) the object catalog(s) (mandatory if back='yes',
           comma separated list if there is more than one catalogue)
        c) the direct image associated with the grism image (optional)

        The file used by axeprep as inlist can be reused again in axecore,
        drzprep and axedrizzle, perhaps extended with different dmag-values
        for the grism images.

    configs: string
      name of the aXe configuration file. If several image
      extensions are to be processed (e.g. for WFC images), one
      configuration file per extension must be given in a comma
      separated list.

    background: boolean
      Switch on/off background subtraction

    backims: string
      name of the background image. If several image extensions
      are to be processed (e.g. for WFC images), one background 
      image per extension must be specified in a comma separated 
      list.
    backped: None
      UNKNOWN

    mfwhm: float
      real number to specify the extent (as a multiple of A_IMAGE
      or B_IMAGE) of the area that is masked out perpendicular to
      the trace of each object before the background level is
      determined (see parameter mfwhm in gol2af).

    norm: boolean
      switch on/off the exposure time normalization

    gcorr: boolean
      switch on/off the gain conversion

    Notes
    -----
    AXEPREP changes the SCI (and potentially ERR) extensions
    of the grism images! It is highly recommended to work only
    on copies of the original files in order to be able to repeat
    the reduction with different parameters.
    """

    # do all the input checks
    inchecks = inputchecks.InputChecker('AXEPREP', inlist, configs, backims)
    inchecks.check_axeprep(backgr, backims)

    # create a list with the basic aXe inputs
    axe_inputs = axeinputs.aXeInput(inlist, configs, backims)

    # go over all the input
    for row in axe_inputs:
        # make a prepare-object; run the prepare
        aXePrep = axepreptor.aXePrepArator(row['grisim'],
                                           row['objcat'],
                                           row['dirim'],
                                           row['config'],
                                           row['dmag'],
                                           backgr=backgr,
                                           master_bck=row['fringe'],
                                           backped=backped,
                                           mfwhm=mfwhm,
                                           norm=norm,
                                           gcorr=gcorr)

        aXePrep.run()
        del aXePrep


def axecore(inlist='',
            configs='',
            fconfterm='',
            back=False,
            extrfwhm=None,
            drzfwhm=None,
            backfwhm=None,
            lambda_mark=None,
            slitless_geom=True,
            orient=True,
            exclude=False,
            cont_model='gauss',
            model_scale=None,
            inter_type='linear',
            lambda_psf=None,
            np=None,
            interp=None,
            niter_med=None,
            niter_fit=None,
            kappa=None,
            smooth_length=None,
            smooth_fwhm=None,
            spectr=True,
            adj_sens=True,
            weights=False,
            sampling='drizzle'):
    """Convenience function for the aXe task AXECORE"""
    axe_setup()

    # do all the file checks
    inchecks = inputchecks.InputChecker('AXECORE', inlist, configs)
    inchecks.check_axecore(back, extrfwhm, drzfwhm, backfwhm, orient,
                           slitless_geom, np, interp, cont_model, weights,
                           sampling)

    # create a list with the basic aXe inputs
    axe_inputs = axeinputs.aXeInput(inlist, configs, fconfterm)

    # go over all the input
    for row in axe_inputs:

        # make an extraction object
        _log.info("image is located: {0}".format(row['grisim']))
        aXeNator = axesingextr.aXeSpcExtr(row['grisim'],
                                          row['objcat'],
                                          row['dirim'],
                                          row['config'],
                                          row['dmag'],
                                          back=back,
                                          extrfwhm=extrfwhm,
                                          drzfwhm=drzfwhm,
                                          backfwhm=backfwhm,
                                          lambda_mark=lambda_mark,
                                          slitless_geom=slitless_geom,
                                          orient=orient,
                                          exclude=exclude,
                                          cont_model=cont_model,
                                          model_scale=model_scale,
                                          inter_type=inter_type,
                                          lambda_psf=lambda_psf,
                                          np=np,
                                          interp=interp,
                                          niter_med=niter_med,
                                          niter_fit=niter_fit,
                                          kappa=kappa,
                                          smooth_length=smooth_length,
                                          smooth_fwhm=smooth_fwhm,
                                          spectr=spectr,
                                          adj_sens=adj_sens,
                                          weights=weights,
                                          sampling=sampling)
        aXeNator.run()
        del aXeNator


def drzprep(inlist='',
            configs='',
            back=False,
            opt_extr=True):
    """Convenience function for the aXe task DRZPREP"""
    axe_setup()

    # create a list with the basic aXe inputs
    configlist=list(configs.split(','))
    if len(configlist) < 1:
      raise aXeError("No configuration file input to drzprep")

    for conf in configlist:
    # make the objects task object, run it an do the cleaning
        prepArator = axelowlev.aXe_DRZPREP(inlist, 
                                           conf,
                                           back=back,
                                           opt_extr=opt_extr)
    prepArator.run()

    del prepArator


def axecrr(inlist='',
           configs='',
           infwhm=0.0,
           outfwhm=0.0,
           back=False,
           clean=True,
           makespc=True,
           adj_sens=True,
           opt_extr=True,
           driz_separate=False):

    """Function for aXedrizzle with CosmicRay-rejection"""
    axe_setup(tmpdir=True)

    # do all the input checks
    inchecks = inputchecks.InputChecker('AXEDRIZZLE', inlist, configs)
    inchecks.check_axedrizzle(infwhm, outfwhm, back)

    # unload the DPP's
    dpps = dppdumps.DPPdumps(inlist, configs, False)
    dpps.filet_dpp(opt_extr)

    # get the contamination information
    cont_info = dpps.is_quant_contam()

    # delete the object
    del dpps

    # assemble the drizzle parameters
    drizzle_params = drizzleobjects.DrizzleParams(configs)

    # make a list of drizzle objects
    dols = drizzleobjects.DrizzleObjectList(drizzle_params,
                                            cont_info,
                                            opt_extr, back=False)

    _log.info(f"checking files {dols}")
    dols.check_files()

    # prepare and perform the drizzling
    dols.prepare_drizzle()
    dols.drizzle()

    # if there are no background
    # files, immediately extract the spectra
    if not back and makespc:
        # extract spectra from the deep 2D stamps
        mefs = mefobjects.MEFExtractor(drizzle_params,
                                       dols, 
                                       opt_extr=opt_extr)
        mefs.extract(infwhm, outfwhm, adj_sens)
        del mefs

        # delete files
        if clean:
            dols.delete_files()

        del dols

    if back:
        # do all the input checks
        inchecks = inputchecks.InputChecker('AXEDRIZZLE', inlist, configs)
        inchecks.check_axedrizzle(infwhm, outfwhm, back)

        # unload the DPP's
        dpps = dppdumps.DPPdumps(inlist, configs, back=back)
        dpps.filet_dpp(opt_extr)

        # get the contamination information
        # cont_info = dpps.is_quant_contam()

        del dpps

        # make a list of drizzle objects
        back_dols = drizzleobjects.DrizzleObjectList(drizzle_params, None,
                                                     opt_extr, back=back)

        # check all files
        back_dols.check_files()

        # prepare and do the drizzling
        back_dols.prepare_drizzle()
        back_dols.drizzle()

        # extract the spectra,
        if makespc:
            mefs = mefobjects.MEFExtractor(drizzle_params, dols, back_dols,
                                           opt_extr=opt_extr)
            mefs.extract(infwhm, outfwhm, adj_sens)

        if clean:
            dols.delete_files()
            back_dols.delete_files()


def axeddd(inlist='',
           configs='',
           drizzle_par=None,
           infwhm=0.0,
           outfwhm=0.0,
           back=False,
           clean=True,
           makespc=True,
           adj_sens=True,
           opt_extr=True,
           driz_separate=False):
    """Function for aXedrizzle"""
    # make the general setup
    axe_setup(tmpdir=True)

    # do all the input checks
    inchecks = inputchecks.InputChecker('AXEDRIZZLE', inlist, configs)
    inchecks.check_axedrizzle(infwhm, outfwhm, back)
    inchecks.check_axecrr(back)

    # unload the DPP's
    dpps = dppdumps.DPPdumps(inlist, configs, False)
    dpps.filet_dpp(opt_extr)

    # get the contamination information
    cont_info = dpps.is_quant_contam()

    # assemble the drizzle parameters
    drizzle_params = drizzleobjects.DrizzleParams(configs)

    # make a list of drizzle objects
    dols = drizzleobjects.DrizzleObjectList(drizzle_params, drizzle_par,
                                     cont_info, opt_extr, back)

    # check all files
    dols.check_files()

    # prepare and do the drizzling
    dols.prepare_drizzle()
    dols.drizzle()

    # if there are no background files, immediately extract
    # the spectra
    if makespc:
        # extract spectra from the deep 2D stamps
        mefs = mefobjects.MEFExtractor(drizzle_params, dols, opt_extr=opt_extr)
        mefs.extract(infwhm, outfwhm, adj_sens)

        # delete tmp objects if
        # requested
        if clean:
            dols.delete_files()


def sex2gol(grism='',
            config='',
            in_sex='',
            use_direct=True,
            direct=None,
            dir_hdu=None,
            spec_hdu=None,
            out_sex=None,
            silent=False):
    """Function for the aXe task SEX2GOL"""
    # make the general setup
    axe_setup()
    sex2gol = pysex2gol.Sex2GolPy(grism, config,
                                  in_sex=in_sex,
                                  dirname=direct,
                                  out_sex=out_sex,
                                  spec_hdu=spec_hdu,
                                  dir_hdu=dir_hdu)
    sex2gol.runall(silent)


def gol2af(grism='',
           config='',
           mfwhm=None,
           back=False,
           orient=True,
           slitless_geom=True,
           exclude=False,
           lambda_mark=None,
           dmag=None,
           out_af="",
           in_gol=None):
    """Function for the aXe task GOL2AF"""
    # check for required environment variables
    axe_setup()

    # run GOL2AF
    gol2af = axelowlev.aXe_GOL2AF(grism, config,
                                  back=back,
                                  mfwhm=mfwhm,
                                  orient=orient,
                                  slitless_geom=slitless_geom,
                                  exclude=exclude,
                                  lambda_mark=lambda_mark,
                                  dmag=dmag,
                                  in_gol=in_gol,
                                  out_af=out_af)
    gol2af.run()


def af2pet(grism='',
           config='',
           back=False,
           in_af="",
           out_pet=None):
    """Function for the aXe task AF2PET"""
    # check for required environment variables
    axe_setup()

    # run AF2PET
    af2pet = axelowlev.aXe_AF2PET(grism, config,
                                  back=back,
                                  in_af=in_af,
                                  out_pet=out_pet)
    af2pet.run()


def petcont(grism='',
            config='',
            cont_model="",
            model_scale=None,
            spec_models="",
            object_models="",
            inter_type='linear',
            lambda_psf=None,
            cont_map=True,
            in_af="",
            no_pet=False,
            silent=False):
    """Function for the aXe task PETCONT"""
    # check for required environment variables
    axe_setup()

    # run PETCONT
    petcont = axelowlev.aXe_PETCONT(grism, config,
                                    cont_model=cont_model,
                                    model_scale=model_scale,
                                    spec_models=spec_models,
                                    object_models=object_models,
                                    inter_type=inter_type,
                                    lambda_psf=lambda_psf,
                                    cont_map=cont_map,
                                    in_af=in_af,
                                    no_pet=no_pet,
                                    silent=silent)

    petcont.run()


def petff(grism='',
          config='',
          back=False,
          ffname=None):
    """Function for the aXe task PETFF"""
    # check for required environment variables
    axe_setup()

    # run PETFF
    petff = axelowlev.aXe_PETFF(grism, config, back=back, ffname=ffname)
    petff.run()


def backest(grism='',
            config='',
            np=None,
            interp=None,
            niter_med=None,
            niter_fit=None,
            kappa=None,
            smooth_length=None,
            smooth_fwhm=None,
            old_bck=False,
            mask=False,
            in_af="",
            out_bck=None):
    """Function for the aXe task BACKEST"""
    # check for required environment variables
    axe_setup()

    #  run BACKEST
    backest = axelowlev.aXe_BE(grism, config,
                               np=np,
                               interp=interp,
                               niter_med=niter_med,
                               niter_fit=niter_fit,
                               kappa=kappa,
                               smooth_length=smooth_length,
                               smooth_fwhm=smooth_fwhm,
                               old_bck=old_bck,
                               mask=mask,
                               in_af=in_af,
                               out_bck=out_bck)
    backest.run()


def pet2spc(grism='',
            config='',
            use_bpet=False,
            adj_sens=True,
            weights=False,
            do_flux=True,
            drzpath=False,
            in_af="",
            opet=None,
            bpet=None,
            out_spc=None):
    """Function for the aXe task PET2SPC"""
    # check for required environment variables
    axe_setup()

    pet2spc = axelowlev.aXe_PET2SPC(grism, config,
                                    use_bpet=use_bpet,
                                    adj_sens=adj_sens,
                                    weights=weights,
                                    do_flux=do_flux,
                                    drzpath=drzpath,
                                    in_af=in_af,
                                    opet=opet,
                                    bpet=bpet,
                                    out_spc=out_spc)
    pet2spc.run()


def stamps(grism='',
           config='',
           sampling='trace',
           drzpath=False,
           in_af="",
           in_pet=None,
           out_stp=None):
    """Function for the aXe task STAMPS"""
    # check for required environment variables
    axe_setup()

    # run STAMPS
    stamps = axelowlev.aXe_STAMPS(grism, config,
                                  sampling=sampling,
                                  drzpath=drzpath,
                                  in_af=in_af,
                                  in_pet=in_pet,
                                  out_stp=out_stp)

    stamps.run()


def drz2pet(inlist='',
            config='',
            opt_extr=False,
            back=False,
            in_af="",
            out_pet=None):
    """Function for the aXe task DRZ2PET"""
    # check for required environment variables
    axe_setup()

    # run the DRZ2PET task
    drz2pet = axelowlev.aXe_DRZ2PET(inlist=inlist,
                                    config=config,
                                    opt_extr=opt_extr,
                                    back=back,
                                    in_af=in_af,
                                    out_pet=out_pet)
    drz2pet.run()


def axegps(grism='',
           config='',
           beam_ref='',
           xval=None,
           yval=None):
    """Function for the aXe task AXEGPS"""
    # check for required environment variables
    axe_setup()
    axegps = axelowlev.aXe_GPS(grism, config, beam_ref, xval, yval)
    axegps.runall()


def axedirim(dirname='',
             config='',
             tpass_direct='',
             model_spectra=None,
             model_images=None,
             model_scale=None,
             tel_area=None,
             silent=False):
    """Function for the aXe task AXEDIRIM"""
    # check for required environment variables
    axe_setup()

    # run the command and delete what's left
    axedirim = axelowlev.aXe_DIRIMAGE(dirname, config, tpass_direct,
                                      model_spectra=model_spectra,
                                      model_images=model_images,
                                      model_scale=model_scale,
                                      tel_area=tel_area)
    return axedirim.runall(silent)
