"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks

def test_iolprep():
    """test the iolprep task"""
    outfiles = ['ib6o23rtq_flt_1.cat',
                'ib6o23rwq_flt_1.cat',
                'ib6o23rzq_flt_1.cat',
                'ib6o23s2q_flt_1.cat']

    os.chdir('F140W')
    axetasks.iolprep(drizzle_image='F140W_drz.fits',
                     input_cat='cookbook.cat',
                     dimension_in='183,85,50,50')

    #  basic check for existance and non-zero size
    for f in outfiles:
        assert os.path.isfile(f)
        stats = os.stat(f)
        assert stats.st_size > 0

