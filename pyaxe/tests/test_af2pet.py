"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from hstaxe import axetasks


def test_af2pet():
    """test the af2pet task"""

    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in outfiles:
        axetasks.af2pet(grism=fname,
                        config="G141.F140W.V4.31.conf",
                        back=False,
                        out_pet=None)