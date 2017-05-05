import iraf

no  = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile  = 'axe$axeprep.par'
_taskname = 'axeprep'

######
# Set up Python IRAF interface here
######
def axeprep_iraf(inlist,
                 configs,
                 backgr,
                 backims,
                 backped,
                 mfwhm,
                 norm,
                 gcorr):

    # properly format the strings
    inlist  = axesrc.straighten_string(inlist)
    configs = axesrc.straighten_string(configs)
    backims = axesrc.straighten_string(backims)
    backped = axesrc.straighten_string(backped)


    # transform the IF booleans to python
    if backgr == yes:
        backgr = True
    else:
        backgr = False
    if norm == yes:
        norm = True
    else:
        norm = False
    if gcorr == yes:
        gcorr = True
    else:
        gcorr = False

    # check whether something should be done
    if inlist != None and configs != None:
        # call the main function
        axesrc.axeprep(inlist=inlist,
                       configs=configs,
                       backgr=backgr,
                       backims=backims,
                       backped=backped,
                       mfwhm=mfwhm,
                       norm=norm,
                       gcorr=gcorr)
    else:
        # print the help
        iraf.help(_taskname)

parfile = iraf.osfn(_parfile)
multid = iraf.IrafTaskFactory(taskname=_taskname, value=parfile,
        pkgname=PkgName, pkgbinary=PkgBinary, function=axeprep_iraf)
