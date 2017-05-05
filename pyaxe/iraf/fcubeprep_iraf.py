import iraf

no  = iraf.no
yes = True

from axe import axesrc

# Point to default parameter file for task
_parfile = 'axe$fcubeprep.par'
_taskname = 'fcubeprep'

######
# Set up Python IRAF interface here
######
def fcubeprep_iraf(grism_image,
                   segm_image,
                   filter_info,
                   AB_zero,
                   dim_info,
                   interpol,
                   useMdriz):

    # properly format the strings
    grism_image = axesrc.straighten_string(grism_image)
    segm_image  = axesrc.straighten_string(segm_image)
    filter_info = axesrc.straighten_string(filter_info)
    dim_info    = axesrc.straighten_string(dim_info)

    # transform the IF booleans to python
    if AB_zero == yes:
        AB_zero = True
    else:
        AB_zero = False

    if useMdriz == no:
        useMdriz = False
    else:
        useMdriz = True
        
    # check whether something should be done
    if grism_image != None and segm_image != None and filter_info != None:
        # call the main function
        axesrc.fcubeprep(grism_image=grism_image,
                         segm_image=segm_image,
                         filter_info=filter_info,
                         AB_zero=AB_zero,
                         dim_info=dim_info,
                         interpol=interpol, 
                         useMdriz=useMdriz)
    else:
        # print the help
        iraf.help(_taskname)


parfile = iraf.osfn(_parfile)
multid = iraf.IrafTaskFactory(taskname=_taskname, value=parfile,
                              pkgname=PkgName, pkgbinary=PkgBinary,
                              function=fcubeprep_iraf)
