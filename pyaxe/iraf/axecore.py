from axe import axesrc

# Point to default parameter file for task
_parfile  = 'axe$axecore.par'
_taskname = 'axecore'

def axecore_iraf(inlist,
                 configs,
                 fconfigs,
                 back,
                 extrfwhm,
                 drzfwhm,
                 backfwhm,
                 orient,
                 slitless_geom,
                 exclude,
                 lambda_mark,
                 cont_model,
                 model_scale,
                 inter_type,
                 lamdbda_psf,
                 np,
                 interp,
                 niter_med,
                 niter_fit,
                 kappa,
                 smooth_length,
                 smooth_fwhm,
                 spectr,
                 adj_sens,
                 weights,
                 sampling):

    # properly format the strings
    inlist   = axesrc.straighten_string(inlist)
    configs  = axesrc.straighten_string(configs)
    fconfigs = axesrc.straighten_string(fconfigs)

    # check whether something should be done
    if inlist is not  None and configs is not  None:
        # run the main code
        axesrc.axecore(inlist,
                       configs,
                       fconfigs,
                       back,
                       extrfwhm,
                       drzfwhm,
                       backfwhm,
                       lambda_mark,
                       slitless_geom,
                       orient,
                       exclude,
                       cont_model,
                       model_scale,
                       inter_type,
                       lamdbda_psf,
                       np,
                       interp,
                       niter_med,
                       niter_fit,
                       kappa,
                       smooth_length,
                       smooth_fwhm,
                       spectr,
                       adj_sens,
                       weights,
                       sampling)
    else:
        raise ValueError("Nothing to run!")
