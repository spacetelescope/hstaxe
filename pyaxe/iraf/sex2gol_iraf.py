import iraf

no  = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$sex2gol.par'
_taskname = 'sex2gol'

######
# Set up Python IRAF interface here
######
def sex2gol_iraf(grism,
                 config,
                 in_sex,
                 use_direct,
                 direct,
                 dir_hdu,
                 spec_hdu,
                 out_sex):

    # properly format the strings
    grism   = axesrc.straighten_string(grism)
    config  = axesrc.straighten_string(config)
    in_sex  = axesrc.straighten_string(in_sex)
    direct  = axesrc.straighten_string(direct)
    out_sex = axesrc.straighten_string(out_sex)


    # transform the IF booleans to python
    if use_direct == True:
        use_direct = True
    else:
        use_direct = False
        direct = None

    # check whether something should be done
    if grism != None and config != None:
        axesrc.sex2gol(grism=grism,
                       config=config,
                       in_sex=in_sex,
                       use_direct=use_direct,
                       direct=direct,
                       dir_hdu=dir_hdu,
                       spec_hdu=spec_hdu,
                       out_sex=out_sex)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=sex2gol_iraf)
