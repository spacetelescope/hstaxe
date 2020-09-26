"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

Notes
-----
aXe_BE: Interpolation order 
        (-1=median, 0= cst, 1= linear etc..)

"""
import os
from hstaxe import axetasks


def test_be():
    """test the background estimation task.

    The data file and the corresponding OAF file
    are used to produce a mask file. 
    From the C task:
         aXe task to compute an estimate of the p/grism image background
         or to create a mask file for background subtraction (option -msk).
         This tasks uses the existing p/grism image and an existing
         aperture file defining the position, orientation, and width
         of all the beams (orders) that are believed to be in the image.
         The values in the regions within each of these beams (orders) are
         replaced by the median, average, linear, or n^th order
         interpolation of pixels which are immediately above and below a
         beam (but not within any other beam). The number of pixels to use
         for this is by default se to 10 both below and above each aperture.
         The -np option can be used to change this default value.
         If the number of points is set to a value which is 0 or less, then
         the entire column of an image will be used, ignoring any pixel
         which are within any known beam. This option allows for a global
         background estimate to be created instead of a local background
         estimate.
         The type of interpolation is controlled by the -interp option:
                       -interp= -1    ; Median
                       -interp= 0     ; Average
                       -interp= 1     ; Linear fit
                       -interp= (n>1) ; n^th order polynomial fit
        
         The output of this task is a FITS image containing two extensions:
         A SCI extension containing the actual image of the background
         A ERR extension containing an estimate of the error in the fit
        
        Input FITS mages are looked for in $AXE_IMAGE_PATH
        aXe config file is looked for in $AXE_CONFIG_PATH
        All outputs are writen to $AXE_OUTPUT_PATH


    """

    infiles = ['ib6o23rsq_flt.fits',
               'ib6o23ruq_flt.fits',
               'ib6o23ryq_flt.fits',
               'ib6o23s0q_flt.fits']

    for filename in infiles:
        axetasks.backest(grism=filename,
                         config='G141.F140W.V4.31.conf',
                         np=0,
                         interp=-1,
                         niter_med=None,
                         niter_fit=None,
                         kappa=None,
                         smooth_length=None,
                         smooth_fwhm=None,
                         old_bck=False,
                         mask=True,
                         in_af=None,
                         out_bck=None)

        # validate the output file
        maskname = "OUTPUT/"+filename.split(".")[0] + "_2.MSK.fits"
        assert os.path.isfile(maskname)
