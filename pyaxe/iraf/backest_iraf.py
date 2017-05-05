import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$backest.par'
_taskname = 'backest'

######
# Set up Python IRAF interface here
######
def backest_iraf(grism,
                 config,
                 np,
                 interp,
                 niter_med,
                 niter_fit,
                 kappa,
                 smooth_length,
                 smooth_fwhm,
                 old_bck,
                 mask,
                 in_af,
                 out_back):

    # properly format the strings
    grism    = axesrc.straighten_string(grism)
    config   = axesrc.straighten_string(config)
    in_af    = axesrc.straighten_string(in_af)
    out_back = axesrc.straighten_string(out_back)


    # transform the IF booleans to python
    if old_bck == True:
        old_bck = True
    else:
        old_bck = False
    if mask == True:
        mask = True
    else:
        mask = False

    # check whether there is something to start
    if grism != None and config != None:
        axesrc.backest(grism=grism,
                       config=config,
                       np=np,
                       interp=interp,
                       niter_med=niter_med,
                       niter_fit=niter_fit,
                       kappa=kappa,
                       smooth_length=smooth_length ,
                       smooth_fwhm=smooth_fwhm,
                       old_bck=old_bck,
                       mask=mask,
                       in_af=in_af,
                       out_bck=out_back)

    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=backest_iraf)
