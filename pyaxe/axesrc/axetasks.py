from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from . import iolmaking
from . import fcubeobjs
from . import axeutils
from . import axeinputs
from . import axepreptor
from . import inputchecks
from . import axesingextr
from . import axelowlev
from . import dppdumps
from . import drizzleobjects
from . import mefobjects
from . import mdrzobjects
from . import pysex2gol


def iolprep(drizzle_image='',
            input_cat='',
            dim_info='0,0,0,0'):
    """Function for the aXe task IOLPREP"""
    iol_maker = iolmaking.IOL_Maker(drizzle_image,
                                    input_cat,
                                    dim_info)
    iol_maker.run()


def fcubeprep(grism_image='',
              segm_image='',
              filter_info=None,
              AB_zero=True,
              dim_info='0,0,0,0',
              interpol='nearest'):
    """Function for the aXe task FCUBEPREP"""

    # run the main command
    fcmaker = fcubeobjs.FluxCube_Maker(grism_image, segm_image, filter_info,
                                       AB_zero, dim_info, interpol)
    fcmaker.run()


def axeprep(inlist='',
            configs='',
            backgr=True,
            backims='',
            backped=None,
            mfwhm=None,
            norm=True,
            gcorr=False,
            histogram=False):
    """Function for the aXe task AXEPREP"""
    # only temporarily here
    axeutils.axe_setup()

    # do all the input checks
    inchecks = inputchecks.InputChecker('AXEPREP', inlist, configs, backims)
    inchecks.check_axeprep(backgr, backims)

    # create a list with the basic aXe inputs
    axe_inputs = axeinputs.aXeInputList(inlist, configs, backims)

    # go over all the input
    for item in axe_inputs:

        # make a prepare-object; run the prepare
        aXePrep = axepreptor.aXePrepArator(item['GRISIM'],
                                           axeutils.getIMAGE(item['OBJCAT']),
                                           item['DIRIM'],
                                           item['CONFIG'],
                                           item['DMAG'],
                                           backgr=backgr,
                                           master_bck=item['FRINGE'],
                                           backped=backped,
                                           mfwhm=mfwhm,
                                           norm=norm,
                                           gcorr=gcorr)

        # delete the object
        aXePrep.run()


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
            sampling='drzizzle'):
    """Function for the aXe task AXECORE"""
    # only temporarily here
    axeutils.axe_setup()

    # do all the file checks
    inchecks = inputchecks.InputChecker('AXECORE', inlist, configs)
    inchecks.check_axecore(back, extrfwhm, drzfwhm, backfwhm, orient,
                           slitless_geom, np, interp, cont_model, weights,
                           sampling)

    # create a list with the basic aXe inputs
    axe_inputs = axeinputs.aXeInputList(inlist, configs, fconfterm)

    # go over all the input
    for item in axe_inputs:

        # make an extraction object
        aXeNator = axesingextr.aXeSpcExtr(item['GRISIM'],
                                          axeutils.getIMAGE(item['OBJCAT']),
                                          item['DIRIM'],
                                          item['CONFIG'],
                                          item['DMAG'],
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
                                          np=np, interp=interp,
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


def drzprep(inlist='',
            configs='',
            back=False,
            opt_extr=False):
    """Function for the aXe task DRZPREP"""
    # make the general setup;
    # needed to define
    # the BIN-directory
    axeutils.axe_setup()

    # make the objects task object, run it an do the cleaning
    prepArator = axelowlev.aXe_DRZPREP(inlist, configs, back=back,
                                       opt_extr=opt_extr)
    prepArator.runall()


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
    """Function for aXedrizzle with CR-rejection"""
    # make the general setup
    axeutils.axe_setup(tmpdir=True)

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
    dols = drizzleobjects.DrizzleObjectList(drizzle_params, cont_info,
                                            opt_extr, back=False)

    # check all files
    dols.check_files()

    # prepare and do the drizzling
    dols.prepare_drizzle()
    dols.drizzle()

    # if there are no background
    # files, immediately extract
    # the spectra
    if not back and makespc:
        # extract spectra from the deep 2D stamps
        mefs = mefobjects.MEFExtractor(drizzle_params, dols, opt_extr=opt_extr)
        mefs.extract(infwhm, outfwhm, adj_sens)
        del mefs

        # delete files
        if clean:
            dols.delete_files()

        # delete the object
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

        # delete the object
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
        # if desired
        if makespc:
            mefs = mefobjects.MEFExtractor(drizzle_params, dols, back_dols,
                                           opt_extr=opt_extr)
            mefs.extract(infwhm, outfwhm, adj_sens)

        # delete the files
        if clean:
            dols.delete_files()
            back_dols.delete_files()


def axeddd(inlist='',
           configs='',
           mult_drizzle_par=None,
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
    axeutils.axe_setup(tmpdir=True)

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
    dols = mdrzobjects.MulDrzObjList(drizzle_params, mult_drizzle_par,
                                     cont_info, opt_extr, back)

    # check all files
    dols.check_files()

    # prepare and do the drizzling
    dols.prepare_drizzle()
    dols.multidrizzle()

    # if there are no background
    # files, immediately extract
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
    axeutils.axe_setup()

    sex2gol = pysex2gol.Sex2GolPy(grism, config, in_sex=in_sex,
                                  dirname=direct, out_sex=out_sex,
                                  spec_hdu=spec_hdu, dir_hdu=dir_hdu)
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
           out_af=None,
           in_gol=None,
           silent=False):
    """Function for the aXe task GOL2AF"""
    # make the general setup
    axeutils.axe_setup()

    # run GOL2AF
    gol2af = axelowlev.aXe_GOL2AF(grism, config, back=back, mfwhm=mfwhm,
                                  orient=orient, slitless_geom=slitless_geom,
                                  exclude=exclude, lambda_mark=lambda_mark,
                                  dmag=dmag, in_gol=in_gol, out_af=out_af)
    gol2af.runall(silent)


def af2pet(grism='',
           config='',
           back=False,
           in_af=None,
           out_pet=None,
           silent=False):
    """Function for the aXe task AF2PET"""
    # make the general setup
    axeutils.axe_setup()

    # run AF2PET
    af2pet = axelowlev.aXe_AF2PET(grism, config, back=back, in_af=in_af,
                                  out_pet=out_pet)
    af2pet.runall(silent)


def petcont(grism='',
            config='',
            cont_model=None,
            model_scale=None,
            spec_models=None,
            object_models=None,
            inter_type=None,
            lambda_psf=None,
            cont_map=True,
            in_af=None,
            no_pet=False,
            silent=False):
    """Function for the aXe task PETCONT"""
    # make the general setup
    axeutils.axe_setup()

    # run PETCONT
    petcont = axelowlev.aXe_PETCONT(grism, config, cont_model=cont_model,
                                    model_scale=model_scale,
                                    spec_models=spec_models,
                                    object_models=object_models,
                                    inter_type=inter_type,
                                    lambda_psf=lambda_psf, cont_map=cont_map,
                                    in_af=in_af, no_pet=no_pet)

    petcont.runall(silent)


def petff(grism='',
          config='',
          back=False,
          ffname=None):
    """Function for the aXe task PETFF"""
    # make the general setup
    axeutils.axe_setup()

    # run PETFF
    petff = axelowlev.aXe_PETFF(grism, config, back=back, ffname=ffname)
    petff.runall()


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
            in_af=None,
            out_bck=None):
    """Function for the aXe task BACKEST"""
    # make the general setup
    axeutils.axe_setup()

    #  run BACKEST
    backest = axelowlev.aXe_BE(grism, config, np=np, interp=interp,
                               niter_med=niter_med, niter_fit=niter_fit,
                               kappa=kappa, smooth_length=smooth_length,
                               smooth_fwhm=smooth_fwhm,
                               old_bck=old_bck, mask=mask, in_af=in_af,
                               out_bck=out_bck)
    backest.runall()


def pet2spc(grism='',
            config='',
            use_bpet=False,
            adj_sens=True,
            weights=False,
            do_flux=True,
            drzpath=False,
            in_af=None,
            opet=None,
            bpet=None,
            out_spc=None,
            silent=False):
    """Function for the aXe task PET2SPC"""
    # make the general setup
    axeutils.axe_setup()

    pet2spc = axelowlev.aXe_PET2SPC(grism, config, use_bpet=use_bpet,
                                    adj_sens=adj_sens, weights=weights,
                                    do_flux=do_flux, drzpath=drzpath,
                                    in_af=in_af, opet=opet, bpet=bpet,
                                    out_spc=out_spc)
    pet2spc.runall(silent)


def stamps(grism='',
           config='',
           sampling='rectified',
           drzpath=False,
           in_af=None,
           in_pet=None,
           out_stp=None,
           silent=False):
    """Function for the aXe task STAMPS"""
    # make the general setup
    axeutils.axe_setup()

    # run STAMPS
    stamps = axelowlev.aXe_STAMPS(grism, config, sampling=sampling,
                                  drzpath=drzpath, in_af=in_af,
                                  in_pet=in_pet, out_stp=out_stp)

    stamps.runall(silent)


def drz2pet(inlist='',
            config='',
            opt_extr=False,
            back=False,
            in_af=None,
            out_pet=None):
    """Function for the aXe task DRZ2PET"""
    # make the general setup
    axeutils.axe_setup()

    # run the DRZ2PET task
    drz2pet = axelowlev.aXe_DRZ2PET(inlist=inlist, config=config,
                                    opt_extr=opt_extr, back=back,
                                    in_af=in_af, out_pet=out_pet)
    drz2pet.runall()


def axegps(grism='',
           config='',
           beam_ref='',
           xval=None,
           yval=None):
    """Function for the aXe task AXEGPS"""
    # make the general setup
    axeutils.axe_setup()
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
    # make the general setup
    axeutils.axe_setup()

    # run the command and delete what's left
    axedirim = axelowlev.aXe_DIRIMAGE(dirname, config, tpass_direct,
                                      model_spectra=model_spectra,
                                      model_images=model_images,
                                      model_scale=model_scale,
                                      tel_area=tel_area)
    return axedirim.runall(silent)
