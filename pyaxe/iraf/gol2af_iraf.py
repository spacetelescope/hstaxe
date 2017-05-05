import iraf

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$gol2af.par'
_taskname = 'gol2af'

######
# Set up Python IRAF interface here
######
def gol2af_iraf(grism,
                config,
                mfwhm,
                back,
                orient,
                slitless_geom,
                exclude,
                lambda_mark,
                dmag,
                out_af,
                in_gol):

    # properly format the strings
    grism  = axesrc.straighten_string(grism)
    config = axesrc.straighten_string(config)
    out_af = axesrc.straighten_string(out_af)
    in_gol = axesrc.straighten_string(in_gol)

    # transform the IF booleans to python
    if back == True:
        back = True
    else:
        back = False
    if orient == True:
        orient = True
    else:
        orient = False
    if slitless_geom == True:
        slitless_geom = True
    else:
        slitless_geom = False
    if exclude == True:
        exclude = True
    else:
        exclude = False

    # check whether something should be done
    if grism != None and config != None:
        axesrc.gol2af(grism=grism,
                      config=config,
                      mfwhm=mfwhm,
                      back=back,
                      orient=orient,
                      slitless_geom=slitless_geom,
                      exclude=exclude,
                      lambda_mark=lambda_mark,
                      dmag=dmag,
                      out_af=out_af,
                      in_gol=in_gol)
    else:
        # print the help
        iraf.help(_taskname)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=gol2af_iraf)
