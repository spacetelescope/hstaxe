import iraf

no  = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile  = 'axe$axedrizzle.par'
_taskname = 'axedrizzle'

######
# Set up Python IRAF interface here
######
def axedrzcrr_iraf(inlist,
                   configs,
                   infwhm,
                   outfwhm,
                   back,
                   driz_separate,
                   clean,
                   combine_type,
                   combine_maskpt,
                   combine_nsigmas,
                   combine_nlow,
                   combine_nhigh,
                   combine_lthresh,
                   combine_hthresh,
                   combine_grow,
                   blot_interp,
                   blot_sinscl,
                   driz_cr_snr,
                   driz_cr_grow,
                   driz_cr_scale,
                   makespc,
                   adj_sens,
                   opt_extr):

    # properly format the strings
    inlist  = axesrc.straighten_string(inlist)
    configs = axesrc.straighten_string(configs)


    # transform the IF booleans to python
    if back == yes:
        back = True
    else:
        back = False
    if clean == yes:
        clean = True
    else:
        clean = False
    if driz_separate == yes:
        driz_separate = True
    else:
        driz_separate = False
    if makespc == yes:
        makespc = True
    else:
        makespc = False
    if adj_sens == yes:
        adj_sens = True
    else:
        adj_sens = False
    if opt_extr == yes:
        opt_extr = True
    else:
        opt_extr = False

    # check whether something should be done
    if inlist != None and configs != None:

        # check for FULL multidrizzle
        if driz_separate:

            # make an empty dict
            mult_drizzle_par={}
            
            # fill with input parameters
            mult_drizzle_par['combine_maskpt']   = combine_maskpt
            mult_drizzle_par['combine_type']     = combine_type
            mult_drizzle_par['combine_nsigma1']  = float(combine_nsigmas.split()[0])
            mult_drizzle_par['combine_nsigma2']  = float(combine_nsigmas.split()[1])
            mult_drizzle_par['combine_nlow']     = combine_nlow
            mult_drizzle_par['combine_nhigh']    = combine_nhigh
            mult_drizzle_par['combine_lthresh']  = combine_lthresh
            mult_drizzle_par['combine_hthresh']  = combine_hthresh
            mult_drizzle_par['combine_grow']     = combine_grow
            mult_drizzle_par['blot_interp']      = blot_interp
            mult_drizzle_par['blot_sinscl']      = blot_sinscl
            mult_drizzle_par['driz_cr_snr']      = driz_cr_snr
            mult_drizzle_par['driz_cr_grow']     = driz_cr_grow
            mult_drizzle_par['driz_cr_scale']    = driz_cr_scale
            
            # call the main function
            axesrc.axeddd(inlist=inlist,
                          configs=configs,
                          mult_drizzle_par=mult_drizzle_par,
                          infwhm=infwhm,
                          outfwhm=outfwhm,
                          back=back,
                          clean=clean,
                          makespc=makespc,
                          opt_extr=opt_extr,
                          adj_sens=adj_sens)
        else:
            # do a simple, traditional aXedrizzle run
            axesrc.axecrr(inlist=inlist,
                          configs=configs,
                          infwhm=infwhm,
                          outfwhm=outfwhm,
                          back=back,
                          clean=clean,
                          makespc=makespc,
                          opt_extr=opt_extr,
                          adj_sens=adj_sens)
    else:
        # print the help
        iraf.help(_taskname)

parfile = iraf.osfn(_parfile)
multid = iraf.IrafTaskFactory(taskname=_taskname, value=parfile,
        pkgname=PkgName, pkgbinary=PkgBinary, function=axedrzcrr_iraf)
