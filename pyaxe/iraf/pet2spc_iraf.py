import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$pet2spc.par'
_taskname = 'pet2spc'

######
# Set up Python IRAF interface here
######
def pet2spc_iraf(grism,
                 config,
                 use_bpet,
                 adj_sens,
                 weights,
                 do_flux,
                 drzpath,
                 in_af,
                 opet,
                 bpet,
                 out_spc):

    # properly format the strings
    grism   = axesrc.straighten_string(grism)
    config  = axesrc.straighten_string(config)
    in_af   = axesrc.straighten_string(in_af)
    opet    = axesrc.straighten_string(opet)
    bpet    = axesrc.straighten_string(bpet)
    out_spc = axesrc.straighten_string(out_spc)


    # transform the IF booleans to python
    if use_bpet == True:
        use_bpet = True
    else:
        use_bpet = False
    if adj_sens == True:
        adj_sens = True
    else:
        adj_sens = False
    if weights == True:
        weights = True
    else:
        weights = False
    if do_flux == True:
        do_flux = True
    else:
        do_flux = False
    if drzpath == True:
        drzpath = True
    else:
        drzpath = False

    # check whether something should be done
    if grism != None and config != None:
        axesrc.pet2spc(grism=grism,
                       config=config,
                       use_bpet=use_bpet,
                       adj_sens=adj_sens,
                       weights=weights,
                       do_flux=do_flux,
                       drzpath=drzpath,
                       in_af=in_af,
                       opet=opet,
                       bpet=bpet,
                       out_spc=out_spc)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=pet2spc_iraf)
