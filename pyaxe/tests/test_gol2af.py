"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks


def test_gol2af():
    """test the sex2gol task.

    This should turn the translated catalog files
    into grism object list (GOL) files
    """

    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in outfiles:
        axetasks.gol2af(grism=fname,
                        config='G141.F140W.V4.31.conf',
                        mfwhm=3.,
                        back=False,
                        orient=True,
                        slitless_geom=True,
                        exclude=False,
                        lambda_mark=4000,
                        dmag=None,
                        out_af=None,
                        in_gol=None,
                        silent=False)
