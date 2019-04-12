"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks

axelist = """ib6o23rsq_flt.fits ib6o23rtq_flt_1.cat ib6o23rtq_flt.fits
ib6o23ruq_flt.fits ib6o23rwq_flt_1.cat ib6o23rwq_flt.fits
ib6o23ryq_flt.fits ib6o23rzq_flt_1.cat ib6o23rzq_flt.fits
ib6o23s0q_flt.fits ib6o23s2q_flt_1.cat ib6o23s2q_flt.fits

"""

def test_axeprep():
    """test the axeprep task"""

    filename = 'aXe.lis'
    with open(filename, 'w') as fd:
        fd.write(axelist)
        fd.close()

    axetasks.axeprep(inlist="aXe.lis",
                     configs="G141.F140W.V4.31.conf",
                     backgr=True,
                     backims="WFC3.IR.G141.sky.V1.0.fits",
                     norm=False,
                     mfwhm=3.0)

    os.remove(filename)
