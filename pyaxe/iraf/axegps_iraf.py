import iraf

no  = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$axegps.par'
_taskname = 'axegps'

######
# Set up Python IRAF interface here
######
def axegps_iraf(grism,
                config,
                beam_ref,
                xval,
                yval):

    # properly format the strings
    grism    = axesrc.straighten_string(grism)
    config   = axesrc.straighten_string(config)
    beam_ref = axesrc.straighten_string(beam_ref)



    # check whether something should be done
    if grism != None and config != None and beam_ref != None:
        # call the python task
        axesrc.axegps(grism, config, beam_ref, xval, yval)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=axegps_iraf)
