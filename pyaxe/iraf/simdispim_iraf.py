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
_parfile = 'axe$simdispim.par'
_taskname = 'simdispim'


######
# Set up Python IRAF interface here
######
def simdispim_iraf(incat, config, dispim_name=None, lambda_psf=None,
                   model_spectra=None, model_images=None, nx=None, ny=None,
                   exptime=None, bck_flux=None, extraction=None, extrfwhm=None,
                   orient=None, slitless_geom=None, adj_sens=None, silent=None):

    # properly format the strings
    dispim_name   = axesrc.straighten_string(dispim_name)
    model_spectra = axesrc.straighten_string(model_spectra)
    model_images  = axesrc.straighten_string(model_images)


    # stratify the input
    # for the real numbers
    if bck_flux == None:
        bck_flux = 0.0
    if extrfwhm == None:
        extrfwhm = 0.0

    # convert the iraf-value
    # to a true boolean
    if extraction != no:
        extraction = True
    else:
        extraction = False
    if orient != no:
        orient = True
    else:
        orient = False
    if slitless_geom != no:
        slitless_geom = True
    else:
        slitless_geom = False
    if adj_sens != no:
        adj_sens = True
    else:
        adj_sens = False
    if silent != no:
        silent = True
    else:
        silent = False

    if (incat == None or len(incat) < 1) or (config==None or len(config) < 1):
        # print the help
        iraf.help(_taskname)

    else:
        axesrc.simdispim(incat=incat.strip(), config=config.strip(),
                         lambda_psf=lambda_psf, dispim_name=dispim_name,
                         model_spectra=model_spectra, model_images=model_images,
                         nx=nx, ny=ny, exptime=exptime, bck_flux=bck_flux,
                         extraction=extraction, extrfwhm=extrfwhm, orient=orient,
                         slitless_geom=slitless_geom, adj_sens=adj_sens,
                         silent=silent)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=simdispim_iraf)
