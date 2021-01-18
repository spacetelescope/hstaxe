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
    """Convenience function for the aXe task FCUBEPREP.

    Parameters
    ----------
    grism_image:    the name of the combined grism image

    segm_image:     name of the segmentation image

    filter_info:    name, wavelength, zeropoint of the filter image(s).
                    If there are several filter images the comma separated
                    quantities are written in a file, and the name of this file
                    is given here.

    AB_zero:        boolean to indicate whether the zeropoints given in the
                    parameter 'filter_info' are in AB- or ST-magnitudes

    dimension_info: four numbers to specify the modifications
                    [left, right, bottom, top] to the target area on the grism
                    images. E.g. 100,500,10,0 would produce fluxcube images
                    which cover the area -100 < x < x_size + 500 and -10 < y < y_size
                    in the input grism images

    interpol:       the inpolation scheme used to compute flux values at
                    the interpolated wavelengths
    """

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
    """Convenience function for the aXe task AXECORE.

    Parameters
    ----------
    inlist: str
      Input Image List which gives on each line
          a) the name of the grism image to be processed (mandatory)
          b) the object catalog(s) (mandatory)
          c) the direct image associated with the grism image (optional)
          d) dmag value (see GOL2AF) for the grism image (optional)

    configs: str
      name of the axe configuration file. If several image
      extensions are to be processed (e.g. for WFC images), one
      configuration file per extension must be given in a comma
      separated list.

    back: Bool
      to switch on/off the creation of a background PET with
      mfwhm=backfwhm

    extrfwhm: float
      mfwhm value to specify the extraction width in gol2af

    drzfwhm: float
      mfwhm value to specify the extraction in axedrizzle

    backfwhm: float
      mfwhm value to specify the width of the background PET

    orient: bool
      enable tilted extraction

    slitless_geom: bool
      enable the best extraction for slitless spectroscopy

    exclude: bool
     switch off the listing of faint objects

    lambda_mark: float
      the wavelength at which to apply the cutoff magnitudes
      MMAG_EXTRACT and MMAG_MARK

    cont_model: str
      name of the contamination model to be applied

    model_scale: float
      scale factor for the gaussian contamination model

    interp_type: str
      interpolation type for the flux values

    lambda_psf: float
      wavelength [nm] at which the object widths were measured

    np: int
      number of points for background estimation

    interp: int
      interpolation type for background determination
      (-1: GLOBAL median; 0: local median; 1: linear fit;
      2: quadratic fit)

    niter_med: int
      number of kappa-sigma iterations around the median

    niter_fit: int
      number of kappa-sigma iterations around the fit value

    kappa: float
      kappa value

    smooth_length: int
      number of adjacent pixels on each side to use when
      smoothing the background estimate in x-direction

    smooth_fwhm: float
      FWHM of the Gaussian used in the background smoothing

    spectr: bool
      enable the creation of SPCs and STPs for each of the
      grism files individually

    weights: bool
      compute and apply optimal weights

    adj_sens: bool
       adjust the sensitivity function for extended sources

    sampling: str
      the sampling mode for the stamp images

    """
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
            opt_extr=False):
    """Convenience function for the aXe task DRZPREP.

    Parameters
    ----------
    inlist: str
      Input Image List which gives the name of the grism image to
      be processed as the first item on each line.

    configs: str
      name of the aXe configuration file. If several image
      extensions are to be processed (e.g. for WFC images), one
      configuration file per extension must be given in a comma
      separated list.

    opt_extr: bool 
      to generate also the necessary data  for optimal
      extraction in axedrizzle

    back: bool
      to switch on the creation of background DPPs made
      by processing background PETs.

    Output
    ------
    if back = False
      $AXE_DRIZZLE_PATH/[slitless filename]_[ext number].DPP.fits
    if back = True
      $AXE_DRIZZLE_PATH/[slitless filename]_[ext number].BCK.DPP.fits

    """
    axe_setup()

    # create a list with the basic aXe inputs
    configlist=list(configs.split(','))
    if len(configlist) < 1:
      raise aXeError("No configuration file input to drzprep")

    prepArator = axelowlev.aXe_DRZPREP(inlist, 
                                       configs,
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
           opt_extr=False,
           driz_separate=False):

    """Function for aXedrizzle with CosmicRay-rejection.

    Parameters
    ----------


    """
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
                                            opt_extr, back=back)

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
    """Function for the aXe task SEX2GOL.

    Parameters
    ----------
    grism: str
      name of the grism image to be processed.

    config: str
      name of the axe configuration file.

    in_sex: str
      name of the object file.

    use_direct: bool
      indicate that the Input Object List refers to a
      direct image

    direct: str
      name of the direct image

    dir_hdu: int
      direct image extension to be used

    spec_hdu: int
      grism/prism image extension to be used

    out_SEX: str
      overwrites the default output object catalog name

    silent: bool
      print messages

    """
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
    """Function for the aXe task GOL2AF.

    Parameters
    ----------
    grism: str
      name of the grism image

    config: str
      name of the aXe configuration file

    mfwhm: float
      the extraction width multiplicative factor

    back: bool
      to generate a BAF instead of an OAF file

    orient: bool
      switch on/off tilted extraction

    slitless_geom: bool
      switch on/off automatic orientation for the tilted extraction

    exclude: bool to switch on the removal of faint
      objects in the result

    lambda_mark: float
      the wavelength at which to apply the cutoff magnitudes
      MMAG_EXTRACT and MMAG_MARK

    dmag: float
      a number to add to the MMAG_EXTRACT and MMAG_MARK
      values given in the configuration file

    out_af: str
      overwrites the default output OAF or BAF filename

    in_gol: str
      overwrites the default input catalog name


    """
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
    """Function for the aXe task AF2PET.

    This task uses the input slitless image together with an Object
    Aperture File (OAF) to generate an Object Pixel Extraction Table
    (OPET) for the input data.

    Parameters
    ----------
    grism: str
      name of the grism image

    config: str
      name of the aXe configuration file

    back:  bool
      generate a PET for a background image using
      a BAF file instead of a OAF file and using a
      background image generated by backest

    in_af : str
      Name to use for the input aperture file with the stamp images
      instead of the default.petcont

    out_pet : str
      Name to use for the output PET file instead of the default.

    """
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
            cont_model='',
            model_scale=None,
            spec_models='',
            object_models='',
            inter_type='linear',
            lambda_psf=None,
            cont_map=True,
            in_af='',
            no_pet=False,
            silent=False):
    """Function for the aXe task PETCONT.

    The task computes and stores the contamination information for a 
    given Pixel Extraction Table. There are two distinct ways to
    compute the contamination:

    The geometrical contamination records, for each PET pixel, how often
    it is a member of a different beam. If a pixel is a member of two 
    separate beams, i.e. is in a region where two beams overlap, it is 
    assigned a value of 1 in each of the two beam PET’s, thus indicating 
    that this pixel is also part of another beam. In quantitative contamination,
    the amount of contaminating flux from other beams is estimated for each 
    PET pixel. This estimate is based on a model of the emitting sources.
    There are two different methods to establish an emission model, 
    the gaussian emission model and the fluxcube model.
    
    Parameters
    ----------
    grism: str
      name of the grism image

    config: str
      name of the aXe configuration file

    cont_model: str
      name of the contamination model to be applied

    model_scale: float
      scale factor for the gaussian cont. model

    spec_models: str
      name of the multi-extension fits table with model spectra

    object_models: str
      name of the multi-extension fits image with object templates.

    interp_type: str
      interpolation type for the flux values

    lambda_psf: float
      wavelength [nm] at which the object widths were measured

    cont_map: bool
      write the contamination map into a FITS file

    in_af: str
      overwrites the input AF file name

    no_pet: bool
      whether a PET exists

    """
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
