import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$petff.par'
_taskname = 'petff'

######
# Set up Python IRAF interface here
######
def petff_iraf(grism,
               config,
               back,
               ffname):

    # properly format the strings
    grism  = axesrc.straighten_string(grism)
    config = axesrc.straighten_string(config)
    ffname = axesrc.straighten_string(ffname)

    # transform the IF booleans to python
    if back == True:
        back = True
    else:
        back = False

    # check whether something should be done
    if grism != None and config!= None:
        axesrc.petff(grism=grism,
                     config=config,
                     back=back,
                     ffname=ffname)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=petff_iraf)
