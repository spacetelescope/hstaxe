"""Licensed under a 3-clause BSD style license - see LICENSE.rst.
pure test of the run
"""
import os
from hstaxe import axetasks


def test_af2pet():
    """test the af2pet task"""

    infiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']
    outfiles = ['OUTPUT/ib6o23rsq_flt_2.PET.fits',
                'OUTPUT/ib6o23ruq_flt_2.PET.fits',
                'OUTPUT/ib6o23ryq_flt_2.PET.fits',
                'OUTPUT/ib6o23s0q_flt_2.PET.fits']


    for fname in infiles:
        axetasks.af2pet(grism=fname,
                        config="G141.F140W.V4.31.conf",
                        back=False,
                        out_pet=None)

    # make some basic existence checks
    for image in infiles:
        for out in outfiles:
            assert os.path.isfile(out)
            stats = os.stat(out)
            assert stats.st_size > 0

