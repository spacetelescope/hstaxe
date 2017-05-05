import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$stamps.par'
_taskname = 'stamps'

######
# Set up Python IRAF interface here
######
def stamps_iraf(grism,
                config,
                sampling,
                drzpath,
                in_af,
                in_pet,
                out_stp):

    # properly format the strings
    grism   = axesrc.straighten_string(grism)
    config  = axesrc.straighten_string(config)
    in_af   = axesrc.straighten_string(in_af)
    in_pet  = axesrc.straighten_string(in_pet)
    out_stp = axesrc.straighten_string(out_stp)

    # transform the IF booleans to python
    if drzpath == True:
        drzpath = True
    else:
        drzpath = False

    # check for minimal input
    if grism != None and config != None:
        axesrc.stamps(grism=grism,
                      config=config,
                      sampling=sampling,
                      drzpath=drzpath,
                      in_af=in_af,
                      in_pet=in_pet,
                      out_stp=out_stp)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=stamps_iraf)
