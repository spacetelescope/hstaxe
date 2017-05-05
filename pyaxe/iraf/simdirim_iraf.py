"""
$Revision: 1.2 $ $Date: 2010/05/18 07:56:00 $
Author: Martin Kuemmel (mkuemmel@stecf.org)
Affiliation: Space Telescope - European Coordinating Facility
WWW: http://www.stecf.org/software/slitless_software/axesim/
"""
import os
import iraf
import sys

no = iraf.no
yes = True

from axe import axesrc
#import axesrc

# Point to default parameter file for task
_parfile = 'axe$simdirim.par'
_taskname = 'simdirim'


######
# Set up Python IRAF interface here
######
def simdirim_iraf(incat, config, tpass_direct, dirim_name=None,
                  model_spectra=None, model_images=None, nx=None, ny=None,
                  exptime=None, bck_flux=0.0, silent=None):

    # properly format the strings
    dirim_name    = axesrc.straighten_string(dirim_name)
    model_spectra = axesrc.straighten_string(model_spectra)
    model_images  = axesrc.straighten_string(model_images)
    tpass_direct  = axesrc.straighten_string(tpass_direct)

    # stratify the input
    # for the background level
    if bck_flux == None:
        bck_flux = 0.0

    # convert the iraf-value
    # to a true boolean
    if silent != no:
        silent = True
    else:
        silent = False

    if incat==None or config==None or incat==None:
        # print the help
        iraf.help(_taskname)

    else:
        axesrc.simdirim(incat=incat.strip(), config=config.strip(),
                        tpass_direct=tpass_direct.strip(), dirim_name=dirim_name,
                        model_spectra=model_spectra, model_images=model_images,
                        nx=nx, ny=ny, exptime=exptime, bck_flux=bck_flux,
                        silent=silent)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=simdirim_iraf)
