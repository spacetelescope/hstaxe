"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from pyaxe import axetasks


def test_stamps():
    """test the aXe_stamps task"""

    # These are the files in the OUTPUT dir
    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in outfiles:
        axetasks.stamps(grism=fname,
                        config='G141.F140W.V4.31.conf',
                        sampling='rectified',
                        drzpath=False,
                        in_af=None,
                        in_pet=None,
                        out_stp=None,
                        silent=False)
