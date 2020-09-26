from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import shutil
from .axeerror import aXeSIMError
from ..axeutils import (getDATA, getOUTSIM, getCONF, get_random_filename,
                        getSIMDATA, getOUTPUT, axe_setup)


def simdata(incat=None, config=None, output_root=None, silent=True,
            inlist_spec=None, tpass_flux=None, inlist_ima=None,
            lambda_psf=None, nx_disp=None, ny_disp=None, exptime_disp=None,
            bck_flux_disp=0.0, extraction=True, extrfwhm=3.0, orient=True,
            slitless_geom=True, adj_sens=True, tpass_direct=None, nx_dir=None,
            ny_dir=None, exptime_dir=None, bck_flux_dir=0.0):
    """Main function for the task SIMDATA

    This module is the high level wrapper function for the
    task SIMDATA. All necessary actions are done, feedback
    is given to the user

    Parameters
    ----------
    incat: str
        name of model object table
    config: str
        aXe configuration file name
    output_root: str
        root for naming output file
    silent: boolean
        flag for silent run
    inlist_spec: str
        spectrum template list
    tpass_flux: str
        total passband file for flux normalization
    inlist_ima: str
        image template list
    lambda_psf: float
        wavelength the object shapes were determined at
    nx_disp: int
        number of pixels in x for dispersed image
    ny_disp: int
        number of pixels in y for dispersed image
    exptime_disp: float
        exposure time for dispersed image
    bck_flux_disp: float
        flux in background for dispersed image
    extraction: bool
        flag for default extraction
    extrfwhm: float
        multiplier for extraction width
    orient: bool
        flag for tilted extraction
    slitless_geom: bool
        flag for slitless optimized extraction
    adj_sens: bool
        flag for adjusted flux conversion
    tpass_direct: str
        total passband file for direct image
    nx_dir: int
        number of pixels in x for direct image
    ny_dir: int
        number of pixels in y for direct image
    exptime_dir: float
        exposure time for direct image
    bck_flux_dir: float
        flux in background for direct image
    """
    from .inputchecks import InputChecker

    # give brief feedback
    print('\nSIMDATA: Starting ...')

    # check for minimal input,
    # print 'doc' if not given
    if ((incat is None) or (config is None)):
        print(__doc__)
        return 1

    # check the input parameters
    in_check = InputChecker(taskname='simdata')

    # for the 'simdisp'-task
    in_check.check_simdispim_input(incat, config,
                                   lambda_psf,
                                   model_spectra=None,
                                   model_images=None,
                                   nx=nx_disp,
                                   ny=ny_disp,
                                   exptime=exptime_disp,
                                   bck_flux=bck_flux_disp,
                                   extraction=extraction,
                                   extrfwhm=extrfwhm,
                                   orient=orient,
                                   slitless_geom=slitless_geom,
                                   adj_sens=adj_sens)

    # check whether a direct image
    # is requested
    if tpass_direct is not None:
        # for the 'simdirim'-task
        in_check.check_simdirim_input(incat, config,
                                      tpass_direct,
                                      model_spectra=None,
                                      model_images=None,
                                      nx=nx_dir,
                                      ny=ny_dir,
                                      exptime=exptime_dir,
                                      bck_flux=bck_flux_dir)

    # check whether there is all input for creating
    # a model spectra file
    if (((inlist_spec is not None) and (tpass_flux is None))
       or ((inlist_spec is None) and (tpass_flux is not None))):
        err_msg = ("\nTo compile model spectra you HAVE to specify a list\n"
                   "of model spectra AND a total passband file for "
                   "normalization.\nOne of them is missing!")
        raise aXeSIMError(err_msg)

    # check whether model spectra are computed
    if ((inlist_spec is not None) and (tpass_flux is not None)):
        # check whether an output root was given
        if (output_root is not None):
            # define a root for the model spectra
            model_spectra = ("{0:s}_spectra.fits".format(output_root))
        else:
            # derive the output name
            # naming convention replicated from module 'prepspectra.py'
            pos = inlist_spec.rfind('.')
            if pos < 0:
                model_spectra = inlist_spec + '.fits'
            else:
                model_spectra = inlist_spec[:pos] + '.fits'

        # compute the model spectra
        prepspectra(inlist=inlist_spec, incat=incat,
                    tpass_flux=tpass_flux,
                    model_spectra=model_spectra, indata_copy=1)

    else:
        # set the name to None
        model_spectra = None

    # check whether model images are computed
    if (inlist_ima is not None):
        # check whether a output root was given
        if (output_root is not None):
            # define a root for the model images
            model_images = "{0:s}_images.fits".format(output_root)
        else:
            # derive the output name;
            # naming convention replicated from module 'prepimages.py'
            pos = inlist_ima.rfind('.')

            if pos < 0:
                model_images = "{0:s}.fits".format(inlist_ima)
            else:
                model_images = "{0:s}.fits".format(inlist_ima[:pos])

        # compute the model images
        prepimages(inlist=inlist_ima, model_images=model_images, indata_copy=1)

    else:
        # set the name to None
        model_images = None

    # check whether a root
    # for the output was given
    if (output_root is not None):
        # define a root for the dispersed image
        dispim_name = "{0:s}_slitless.fits".format(output_root)

        # define a root for the direct image
        dirim_name = "{0:s}_direct.fits".format(output_root)

    else:
        dispim_name = None
        dirim_name = None

    # compute the dispersed image
    simdispim(incat=incat, config=config, dispim_name=dispim_name,
              lambda_psf=lambda_psf, model_spectra=model_spectra,
              model_images=model_images, nx=nx_disp, ny=ny_disp,
              exptime=exptime_disp, bck_flux=bck_flux_disp,
              extraction=extraction, extrfwhm=extrfwhm, orient=orient,
              slitless_geom=slitless_geom, adj_sens=adj_sens,
              silent=silent)

    # check whether a direct images
    # is requested
    if (tpass_direct is not None):
        simdirim(incat=incat, config=config,
                 tpass_direct=tpass_direct, dirim_name=dirim_name,
                 model_spectra=model_spectra, model_images=model_images,
                 nx=nx_dir, ny=ny_dir, exptime=exptime_dir,
                 bck_flux=bck_flux_dir, silent=silent)

    # give brief feedback
    print('SIMDATA: Done ...')

    # return 'success'
    return 0


def simdirim(incat=None, config=None, tpass_direct=None, dirim_name=None,
             model_spectra=None, model_images=None, nx=None, ny=None,
             exptime=None, bck_flux=0.0, silent=True):
    """
    Main function for the task SIMDIRIM

    This module is the high level wrapper function for the
    task SIMDIRIM. All necessary actions are done, feedback
    is given to the user

    @param incat: name of model object table
    @type incat: string
    @param config: aXe configuration file name
    @type config: string
    @param tpass_direct: total passband file
    @type tpass_direct: string
    @param dirim_name: name of direct image
    @type dirim_name: string
    @param model_spectra: name of model spectra file
    @type model_spectra: string
    @param model_images: name of model images
    @type model_image: string
    @param nx: number of pixels in x
    @type nx: int
    @param ny: number of pixels in y
    @type ny: int
    @param exptime: exposure time
    @type exptime: dloat
    @param bck_flux: flux in background
    @type bck_flux: float
    @param silent: flag for silent run
    @type silen: boolean
    """
    from . import interpolator
    from . import modspeclist
    from . import axecommands
    from . import realworld
    from . import configfile
    from . import imagemaker
    from .inputchecks import InputChecker

    # give brief feedback
    print('\nSIMDIRIM: Starting ...')

    # only temporarily here
    axe_setup(axesim=True)

    if ((incat is None) or (config is None) or (tpass_direct is None)):
        print(__doc__)
        return 1

    # check the input parameters
    in_check = InputChecker(taskname='simdirim')
    # for the 'simdisp'-task
    in_check.check_simdirim_input(incat, config, tpass_direct,
                                  model_spectra, model_images, nx,
                                  ny, exptime, bck_flux)

    if (dirim_name is None):
        # derive the output name
        pos = incat.rfind('.')
        if pos < 0:
            dirima_name = "{0:s}_direct.fits".format(incat)
        else:
            dirima_name = "{0:s}_direct.fits".format(incat[:pos])
    else:
        dirima_name = dirim_name

    # make a full path to the
    # direct image as dummy and as final output
    dummy_dirima_path = getDATA(dirima_name)
    final_dirima_path = getOUTSIM(dirima_name)

    try:
        # to convert the background value
        # to a float
        bck_flux = float(bck_flux)
    except ValueError:
        # now it must be a file;
        # check for its existence
        if not os.path.isfile(getCONF(bck_flux)):
            err_msg = "Missing background image: {0}".format(getCONF(bck_flux))
            raise aXeSIMError(err_msg)

        # store the path to the
        # background image
        bck_flux = getCONF(bck_flux)

    # load the aXe configuration file
    conf = configfile.ConfigFile(getCONF(config))

    # make the simulation configuration
    # file pointing the correct extensions
    config_simul = conf.axesim_prep()

    # delete the object
    # explicitly
    del conf

    # load the simulation configuration file
    conf_simul = configfile.ConfigFile(getCONF(config_simul))

    print("SIMDIRIM: Input Model Object List:\t{0:s}".format(getDATA(incat)))
    print("SIMDIRIM: Input aXe configuration file:\t{0:s}"
          .format(getCONF(config)))
    print("SIMDIRIM: Input Total Passband file:\t{0:s}"
          .format(getSIMDATA(tpass_direct)))
    if (model_spectra is not None):
        print("SIMDIRIM: Input Model Spectra:\t{0:s}"
              .format(getDATA(model_spectra)))
    if (model_images is not None):
        print("SIMDIRIM: Input Model Spectra:\t{0:s}"
              .format(getDATA(model_images)))
    print("SIMDIRIM: Background flux/image:\t{0:s}".format(str(bck_flux)))
    if (exptime is not None):
        print("SIMDIRIM: Input exposure time:\t{0:s}".format(str(exptime)))
    if ((nx is None) and (ny is None)):
        print("SIMDIRIM: Input image dimensions:\tAUTO")
    else:
        print("SIMDIRIM: Input image dimensions:\t({0:s},{1:s})"
              .format(str(nx), str(ny)))
    print("SIMDIRIM: Output dispersed image:\t{0:s}\n"
          .format(final_dirima_path))

    # check whether the name ends with '.fits'
    if ('.fits' in tpass_direct[-5:]):
        # load the ascii list
        new_interp = interpolator.Interpolator(getSIMDATA(tpass_direct))

        # save it as fits
        new_name = new_interp.writetofits(colname1='WAVELENGTH',
                                          colname2='THROUGHPUT')

        # give feedback
        print("Ascii list: {0:s} converted to fits: {1:s}".format(tpass_direct,
                                                                  new_name))

        # overwrite the old name
        tpass_direct = os.path.basename(new_name)

    # create the dummy image maker
    i_maker = imagemaker.DummyImages(getCONF(config_simul),
                                     dirname=dummy_dirima_path,
                                     nx=nx,
                                     ny=ny)
    # nake the dummy images
    i_maker.makeImages()

    # load the model object table
    inobjects = modspeclist.ModelObjectTable(getDATA(incat))
    # fill the model object table
    inobjects.fill_columns(i_maker.WCSimage, i_maker.WCSext)

    # load the object to make the grism simulations
    dirmator = axecommands.DirImator(i_maker,
                                     config_simul,
                                     getDATA(incat),
                                     tpass_direct,
                                     model_spectra,
                                     model_images,
                                     float(conf_simul['TELAREA']))
    dirmator.run(silent=silent)
    dirmator.mopup()

    # delete the dummy images
    i_maker.deleteImages()

    # get the name of the result image, which is the contamination image
    result_image = getOUTPUT(dirima_name.replace('.fits', '_2.CONT.fits'))

    # convert the 'contamination' image into
    # a full output image with three extensions
    # and noise (if desired)
    rworld = realworld.RealWorld(result_image,
                                 extname='SCI',
                                 exptime=exptime,
                                 bck_flux=bck_flux,
                                 rdnoise=conf_simul['RDNOISE'],
                                 instrument=conf_simul['INSTRUMENT'])
    rworld.make_real()

    # move the resulting image to the correct
    # name and place
    shutil.move(result_image, final_dirima_path)

    # give brief feedback
    print('SIMDIRIM: Done ...')
    return 0


def simdispim(incat=None, config=None, lambda_psf=None, dispim_name=None,
              model_spectra=None, model_images=None, nx=None, ny=None,
              exptime=None, bck_flux=0.0, extraction=True, extrfwhm=3.0,
              orient=True, slitless_geom=True, adj_sens=True, silent=True):
    """
    Main function for the task SIMDISPIM

    This module is the high level wrapper function for the
    task SIMDISPIM. All necessary actions are done, feedback
    is given to the user

    @param incat: name of model object table
    @type incat: string
    @param config: aXe configuration file name
    @type config: string
    @param lambda_psf: wavelength the object shapes were determined at
    @type lambda_psf: float
    @param dispim_name: name of dispersed image
    @type dispim_name: string
    @param model_spectra: name of model spectra file
    @type model_spectra: string
    @param model_images: name of model images
    @type model_image: string
    @param nx: number of pixels in x
    @type nx: int
    @param ny: number of pixels in y
    @type ny: int
    @param exptime: exposure time
    @type exptime: dloat
    @param bck_flux: flux in background
    @type bck_flux: float
    @param extraction: flag for default extraction
    @type extraction: boolean
    @param extrfwhm: multiplier for extraction width
    @type extrfwhm: float
    @param orient: flag for tilted extraction
    @type orient: boolean
    @param slitless_geom: flag for slitless optimized extraction
    @type slitless_geom: boolean
    @param adj_sens: flag for adjusted flux conversion
    @type adj_sens: boolean
    @param silent: flag for silent run
    @type silen: boolean
    """
    from . import imagemaker
    from . import modspeclist
    from . import axecommands
    from . import realworld
    from . import configfile
    from .inputchecks import InputChecker

    # give brief feedback
    print('\nSIMDISPIM: Starting ...')

    # just set the environments
    axe_setup(axesim=True)

    if ((incat is None) or (config is None)):
        print(__doc__)
        return 1

    # check the input parameters
    in_check = InputChecker(taskname='simdispim')
    # for the 'simdisp'-task
    in_check.check_simdispim_input(incat, config, lambda_psf,
                                   model_spectra, model_images, nx,
                                   ny, exptime, bck_flux, extraction,
                                   extrfwhm, orient, slitless_geom,
                                   adj_sens=adj_sens)

    if (dispim_name is None):
        # derive the output name
        pos = incat.rfind('.')
        if pos < 0:
            dirima_name = incat + '_direct.fits'
            grisima_name = incat + '_slitless.fits'
        else:
            dirima_name = incat[:pos] + '_direct.fits'
            grisima_name = incat[:pos] + '_slitless.fits'
    else:
        dirima_name = dispim_name.replace('.fits', '_direct.fits')
        grisima_name = dispim_name

    # make a full path to the
    # direct image as dummy and as final output
    dummy_dirima_path = getDATA(get_random_filename('t', '.fits'))
    dummy_grisima_path = getDATA(get_random_filename('t', '.fits'))

    final_dirima_path = getOUTSIM(dirima_name)
    final_grisima_path = getOUTSIM(grisima_name)

    try:
        # to convert the background value
        # to a float
        bck_flux = float(bck_flux)
    except ValueError:
        # now it must be a file;
        # check for its existence
        if not os.path.isfile(getCONF(bck_flux)):
            err_msg = 'Missing background image: ' + getCONF(bck_flux)
            raise aXeSIMError(err_msg)

        # store the path to the
        # background image
        bck_flux = getCONF(bck_flux)

    # load the aXe configuration file
    conf = configfile.ConfigFile(getCONF(config))

    # make the simulation configuration
    # file pointing the correct extensions
    config_simul = conf.axesim_prep()

    # delete the object
    # explicitly
    del conf

    # load the simulation configuration file
    conf_simul = configfile.ConfigFile(getCONF(config_simul))

    # make sure a reasonable default
    # for lambda_psf is given if needed
    if (lambda_psf is None):
        lambda_psf = conf_simul.confirm_lambda_psf()

    print("SIMDISPIM: Input Model Object List:\t{0:s}".format(getDATA(incat)))
    print("SIMDISPIM: Input aXe configuration file:\t{0:s}"
          .format(getCONF(config)))
    if (model_spectra is not None):
        print("SIMDISPIM: Input Model Spectra:\t{0:s}"
              .format(getDATA(model_spectra)))
    if (model_images is not None):
        print("SIMDISPIM: Input Model Spectra:\t{0:s}"
              .format(getDATA(model_images)))
    print("SIMDISPIM: Fixed wavlength for PSF:\t{0:s}".format(str(lambda_psf)))
    print("SIMDISPIM: Background flux/image:\t{0:s}".format(str(bck_flux)))
    if (exptime is not None):
        print("SIMDISPIM: Input exposure time:\t{0:s}".format(str(exptime)))
    if ((nx is None) and (ny is None)):
        print("SIMDISPIM: Input image dimensions: AUTO")
    else:
        print("SIMDISPIM: Input image dimensions:\t(%s,%s)".format(str(nx),
                                                                   str(ny)))
    print("SIMDISPIM: Output dispersed image:\t(%s,%s)"
          .format(final_grisima_path))
    if extraction:
        print("SIMDISPIM: Extraction width scaling:\t{0:.2f}".format(extrfwhm))
        print("SIMDISPIM: Extraction tilted: \t({})".format(str(orient)))
        print("SIMDISPIM: Extraction slitless optimized:\t({})"
              .format(str(slitless_geom)))
        print("SIMDISPIM: Size-adjusted flux conversion:\t({})"
              .format(str(adj_sens)))
        print("SIMDISPIM: Output extracted spectra:\t({})"
              .format(final_grisima_path.replace(".fits", "_2.SPC.fits")))
        print("SIMDISPIM: Output stamp images:\t({})"
              .format(final_grisima_path.replace(".fits", "_2.STP.fits")))

    # create the dummy image maker
    i_maker = imagemaker.DummyImages(getCONF(config_simul), dummy_grisima_path,
                                     dummy_dirima_path, nx, ny)
    # nake the dummy images
    i_maker.makeImages()

    # load the model object table
    inobjects = modspeclist.ModelObjectTable(getDATA(incat))
    # fill the model object table
    inobjects.fill_columns(i_maker.WCSimage, i_maker.WCSext)

    # load the object to make the grism simulations
    grismator = axecommands.DispImator(i_maker, config_simul, getDATA(incat),
                                       lambda_psf, model_spectra, model_images)
    grismator.run(silent=silent)
    grismator.mopup()

    # get the name of the result image, which is the contamination image
    result_image = getOUTPUT(os.path.basename(dummy_grisima_path
                                              ).replace(".fits",
                                                        "_2.CONT.fits"))

    # convert the 'contamination' image into
    # a full output image with three extensions
    # and noise (if desired)
    rworld = realworld.RealWorld(result_image,
                                 extname="SCI",
                                 exptime=exptime,
                                 bck_flux=bck_flux,
                                 rdnoise=conf_simul["RDNOISE"],
                                 instrument=conf_simul["INSTRUMENT"])
    rworld.make_real()

    # move the resulting image to the correct
    # name and place
    shutil.move(result_image, final_grisima_path)

    # check whether an extraction
    # is desired
    if extraction:
        # create and extractor
        extractor = axecommands.DummyExtractor(i_maker,
                                               final_grisima_path,
                                               config_simul,
                                               getDATA(incat),
                                               bck_flux,
                                               extrfwhm,
                                               orient,
                                               slitless_geom,
                                               adj_sens,
                                               lambda_mark=lambda_psf)
        # make the extraction
        extractor.prepare_extraction()
        extractor.run(silent=silent)
        extractor.mopup()

    # delete the dummy images
    i_maker.deleteImages()

    # give brief feedback
    print('SIMDISPIM: Done ...')
    return 0


def prepimages(inlist=None, model_images=None, indata_copy=0):
    """
    Main function for the task PREPIMAGES

    This module is the high level wrapper function for the
    task PREPIMAGES. All necessary actions are done, feedback
    is given to the user

    @param incat: name of the image template list
    @type incat: string
    @param model_images: name of model image file
    @type model_images: string
    @param indata_copy: flag to save also a copy
    @type indata_copy:int
    """
    from . import templateimages

    # give brief feedback
    print("\nPREPIMAGES: Starting ...")

    # just set the environments
    axe_setup(axesim=True)

    # check whether the model images filename
    # is given; if not, derive it
    if (model_images is None):
        pos = inlist.rfind('.')
        if pos < 0:
            model_images = "{0:s}.fits".format(inlist)
        else:
            model_images = "{0:s}.fits".format(inlist[:pos])

    # give feedback onto screen
    print("PREPIMAGES: Input Image Template List:\t{0:s}".format(inlist))
    print("PREPIMAGES: Output Model Image:\t{0:s}"
          .format(getOUTSIM(model_images)))
    if indata_copy:
        print("PREPIMAGES: Output Model Image:\t{0:s}\n"
              .format(getDATA(model_images)))

    # load the images defined in the image template list
    imglist = templateimages.ArtImaList(inlist)

    # save and write to disk
    imglist.tofits(model_images, indata_copy)

    # give brief feedback
    print('PREPIMAGES: Done ...')

    # return 'success'
    return 0


def prepspectra(inlist=None, incat=None, tpass_flux=None,
                model_spectra=None, indata_copy=0):
    """
    Main function for the task PREPIMAGES

    This module is the high level wrapper function for the
    task PREPSPECTRA. All necessary actions are done, feedback
    is given to the user

    @param incat: name of the model object table
    @type incat: string
    @param tpass_flux: total passband for flux normalization
    @type tpass_flux: string
    @param model_spectra: name of model spectra file
    @type model_spectra: string
    @param indata_copy: flag to save also a copy
    @type indata_copy:int
    """
    from . import templatespectra

    # give brief feedback
    print("\nPREPSPECTRA: Starting ...")

    # just set the environments
    axe_setup(axesim=True)

    # check whether a output name was given
    if (model_spectra is None):
        # derive the output name
        pos = inlist.rfind('.')
        if pos < 0:
            model_spectra = "{0:s}.fits".format(inlist)
        else:
            model_spectra = "{0:s}.fits".format(inlist[:pos])

    # give user feedback onto the screen
    print("PREPSPECTRA: Input Spectrum Template List:\t{0:s}".format(inlist))
    print("PREPSPECTRA: Input Model Object Table:\t{0:s}"
          .format(getDATA(incat)))
    if (len(tpass_flux.split(',')) > 1):
        print("PREPSPECTRA: Input Total Passband expression:\t{0:s}"
              .format(tpass_flux))
    else:
        print("PREPSPECTRA: Input Total Passband:\t{0:s}"
              .format(getSIMDATA(tpass_flux)))
    print("PREPSPECTRA: Output Model Spectra:\t{0:s}"
          .format(getOUTSIM(model_spectra)))
    if indata_copy:
        print("PREPSPECTRA: Output Model Spectra:\t{0:s}\n"
              .format(getDATA(model_spectra)))

    # load the spectra given in the spectrum template list
    speclist = templatespectra.TemplateSpectrumList(inlist, incat, tpass_flux)

    # check for wavelength
    # coverage and similar
    speclist.check()

    # move in redshit- and flux- space
    speclist.beamto()

    # write the model spectra file to disk
    speclist.tofits(model_spectra, indata_copy)
    speclist.save_modspec_list(incat)

    # give brief feedback
    print('PREPSPECTRA: Done ...')

    # return 'success'
    return 0
