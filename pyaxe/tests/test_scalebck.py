"""Licensed under a 3-clause BSD style license - see LICENSE.rst.
"""
from hstaxe.axesrc import axelowlev


def test_scalebck():
    """test the back ground scaling task"""

    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for filename in outfiles:  
        maskname = filename.split(".fits")[0] + "_2.MSK.fits"
        axelowlev.aXe_SCALEBCK(filename,
                               maskname,
                               "G141.F140W.V4.31.conf",
                               "WFC3.IR.G141.sky.V1.0.fits")

