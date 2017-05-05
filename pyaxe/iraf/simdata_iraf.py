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
_parfile = 'axe$simdata.par'
_taskname = 'simdata'


######
# Set up Python IRAF interface here
######
def simdata_iraf(incat, config, output_root=None, silent=None, inlist_spec=None,
                 tpass_flux=None, inlist_ima=None, lambda_psf=None,
                 nx_disp=None, ny_disp=None, exptime_disp=None,
                 bck_flux_disp=0.0, extraction=None, extrfwhm=None, orient=None,
                 slitless_geom=None, adj_sens=None, tpass_direct=None,
                 nx_dir=None, ny_dir=None, exptime_dir=None,
                 bck_flux_dir=0.0, version=None):

    # properly format the strings
    output_root   = axesrc.straighten_string(output_root)
    inlist_spec   = axesrc.straighten_string(inlist_spec)
    tpass_flux    = axesrc.straighten_string(tpass_flux)
    inlist_ima    = axesrc.straighten_string(inlist_ima)
    tpass_direct  = axesrc.straighten_string(tpass_direct)

    # stratify the input
    # for the background level
    if bck_flux_disp == None:
        bck_flux_disp = 0.0
    if bck_flux_dir == None:
        bck_flux_dir = 0.0
    if extrfwhm == None:
        extrfwhm = 0.0

    # convert the iraf-value
    # to true booleans
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

    if incat == None or config==None:
        # print the help
        iraf.help(_taskname)

    else:
        axesrc.simdata(incat=incat.strip(), config=config.strip(),
                       output_root=output_root, silent=silent,
                       inlist_spec=inlist_spec,
                       tpass_flux=tpass_flux, inlist_ima=inlist_ima,
                       lambda_psf=lambda_psf, nx_disp=nx_disp,
                       ny_disp=ny_disp, exptime_disp=exptime_disp,
                       bck_flux_disp=bck_flux_disp, extraction=extraction,
                       extrfwhm=extrfwhm, orient=orient,
                       slitless_geom=slitless_geom, adj_sens=adj_sens,
                       tpass_direct=tpass_direct, nx_dir=nx_dir, ny_dir=ny_dir,
                       exptime_dir=exptime_dir, bck_flux_dir=bck_flux_dir)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=simdata_iraf)
