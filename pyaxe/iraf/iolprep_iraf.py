import iraf

no  = iraf.no
yes = True


from axe import axesrc

# Point to default parameter file for task
_parfile = 'axe$iolprep.par'
_taskname = 'iolprep'

######
# Set up Python IRAF interface here
######
def iolprep_iraf(mdrizzle_image,
                 input_cat,
                 dim_info,
                 useMdriz):

    # properly format the strings
    mdrizzle_image = axesrc.straighten_string(mdrizzle_image)
    input_cat  = axesrc.straighten_string(input_cat)
    dim_info    = axesrc.straighten_string(dim_info)

    if useMdriz == iraf.no:
        useMdriz = False
    else:
        useMdriz = True

    # check whether something should be done
    if mdrizzle_image != None and input_cat != None and dim_info != None:
        # call the main function
        axesrc.iolprep(mdrizzle_image=mdrizzle_image,
                       input_cat=input_cat,
                       dim_info=dim_info,
                       useMdriz=useMdriz)
    else:
        # print the help
        iraf.help(_taskname)

parfile = iraf.osfn(_parfile)
multid = iraf.IrafTaskFactory(taskname=_taskname, value=parfile,
                              pkgname=PkgName, pkgbinary=PkgBinary,
                              function=iolprep_iraf)
