import iraf

no  = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile  = 'axe$drzprep.par'
_taskname = 'drzprep'

######
# Set up Python IRAF interface here
######
def drzprep_iraf(inlist,
                 configs,
                 opt_extr,
                 back):

    # properly format the strings
    inlist  = axesrc.straighten_string(inlist)
    configs = axesrc.straighten_string(configs)

    # transform the IF booleans to python
    if opt_extr == yes:
        opt_extr = True
    else:
        opt_extr = False
    if back == yes:
        back = True
    else:
        back = False

    # check whether something should be done
    if inlist != None and configs != None:
        # call the main function
        axesrc.drzprep(inlist=inlist,
                       configs=configs,
                       opt_extr=opt_extr,
                       back=back)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=drzprep_iraf)
