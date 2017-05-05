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
_parfile = 'axe$prepspectra.par'
_taskname = 'prepspectra'


######
# Set up Python IRAF interface here
######
def prepspectra_iraf(inlist, incat, tpass_flux, model_spectra=None):

    # properly format the strings
    inlist        = axesrc.straighten_string(inlist)
    incat         = axesrc.straighten_string(incat)
    model_spectra = axesrc.straighten_string(model_spectra)

    # check whether there is enough input
    if inlist == None or incat == None or tpass_flux == None:

        # print the help if not
        iraf.help(_taskname)

    else:

        # execute the python function
        axesrc.prepspectra(inlist=inlist.strip(), incat=incat.strip(),
                           tpass_flux=tpass_flux.strip(),
                           model_spectra=model_spectra)

# Initialize IRAF Task definition now...
parfile = iraf.osfn(_parfile)
a = iraf.IrafTaskFactory(taskname=_taskname,value=parfile,
            pkgname=PkgName, pkgbinary=PkgBinary, function=prepspectra_iraf)
