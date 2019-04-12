"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

Notes
-----
aXe_BE: Interpolation order 
        (-1=median, 0= cst, 1= linear etc..)
"""
import os
from pyaxe import axetasks


def test_be():
    """test the axeprep task"""

    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for filename in outfiles:
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

