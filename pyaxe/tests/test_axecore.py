"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks

axelist = """ib6o23rsq_flt.fits ib6o23rtq_flt_1.cat ib6o23rtq_flt.fits
ib6o23ruq_flt.fits ib6o23rwq_flt_1.cat ib6o23rwq_flt.fits
ib6o23ryq_flt.fits ib6o23rzq_flt_1.cat ib6o23rzq_flt.fits
ib6o23s0q_flt.fits ib6o23s2q_flt_1.cat ib6o23s2q_flt.fits

"""

def test_axecore():
    """test the axecore task"""

    filename = 'aXe.lis'
    with open(filename, 'w') as fd:
        fd.write(axelist)
        fd.close()
    axetasks.axecore(filename,
                     "G141.F140W.V4.31.conf",
                     extrfwhm=4.,
                     drzfwhm=3.,
                     backfwhm=0.,
                     orient=False,
                     weights=True,
                     slitless_geom=False,
                     exclude=True)
    os.remove(filename)
