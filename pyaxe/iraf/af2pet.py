from axe import axesrc

# Point to default parameter file for task
_parfile = 'axe$af2pet.par'
_taskname = 'af2pet'

def af2pet(grism,
           config,
           back,
           in_af,
           out_pet):

    # properly format the strings
    grism   = axesrc.straighten_string(grism)
    config  = axesrc.straighten_string(config)
    in_af   = axesrc.straighten_string(in_af)
    out_pet = axesrc.straighten_string(out_pet)

    # check whether there is something to start
    if ( (grism is None) and (config is None) ):
        raise ValueError("Nothing to start")
    else:
        axesrc.af2pet(grism=grism,
                      config=config,
                      back=back,
                      in_af=in_af,
                      out_pet=out_pet)
