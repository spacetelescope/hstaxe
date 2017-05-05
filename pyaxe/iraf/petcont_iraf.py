import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$petcont.par'
_taskname = 'petcont'

######
# Set up Python IRAF interface here
######
def petcont_iraf(grism,
                 config,
                 cont_model,
                 model_scale,
                 spec_models,
                 object_models,
                 inter_type,
                 lambda_psf,
                 cont_map,
                 in_af):

    # properly format the strings
    grism         = axesrc.straighten_string(grism)
    config        = axesrc.straighten_string(config)
    spec_models   = axesrc.straighten_string(spec_models)
    object_models = axesrc.straighten_string(object_models)
    in_af         = axesrc.straighten_string(in_af)



    # transform the IF booleans to python
    if cont_map == True:
        cont_map = True
    else:
        cont_map = False

    # check whether something should be done
    if grism != None and config != None:
        axesrc.petcont(grism=grism,
                       config=config,
                       cont_model=cont_model,
                       model_scale=model_scale,
                       spec_models=spec_models,
                       object_models=object_models,
                       inter_type=inter_type,
                       lambda_psf=lambda_psf,
                       cont_map=cont_map,
                       in_af=in_af)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=petcont_iraf)
