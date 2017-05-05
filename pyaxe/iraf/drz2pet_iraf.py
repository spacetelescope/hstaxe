import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$drz2pet.par'
_taskname = 'drz2pet'

######
# Set up Python IRAF interface here
######
def drz2pet_iraf(inlist,
                 config,
                 opt_extr,
                 back,
                 in_af,
                 out_pet):

    # properly format the strings
    inlist  = axesrc.straighten_string(inlist)
    config  = axesrc.straighten_string(config)
    in_af   = axesrc.straighten_string(in_af)
    out_pet = axesrc.straighten_string(out_pet)

    # transform the IF booleans to python
    if back == True:
        back = True
    else:
        back = False
    if opt_extr == True:
        opt_extr = True
    else:
        opt_extr = False

    # check for minimal input
    if inlist!= None and config != None:
        axesrc.drz2pet(inlist=inlist,
                       config=config,
                       opt_extr=opt_extr,
                       back=back,
                       in_af=in_af,
                       out_pet=out_pet)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=drz2pet_iraf)
