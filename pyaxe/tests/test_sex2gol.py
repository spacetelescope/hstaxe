"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks

axelist = """ib6o23rsq_flt.fits ib6o23rtq_flt_1.cat ib6o23rtq_flt.fits
ib6o23ruq_flt.fits ib6o23rwq_flt_1.cat ib6o23rwq_flt.fits
ib6o23ryq_flt.fits ib6o23rzq_flt_1.cat ib6o23rzq_flt.fits
ib6o23s0q_flt.fits ib6o23s2q_flt_1.cat ib6o23s2q_flt.fits

"""


def test_sex2gol():
    """test the sex2gol task"""

    filename = "aXe.lis"
    configs = "G141.F140W.V4.31.conf"
    backims = "WFC3.IR.G141.sky.V1.0.fits"

    with open(filename, 'w') as fd:
        fd.write(axelist)
        fd.close()

    axe_inputs = axetasks.axeinputs.aXeInput(filename, configs, backims)

    for row in axe_inputs:
        axetasks.sex2gol(grism='DATA/'+row['grisim'],
                         config=configs,
                         in_sex='DATA/'+row['objcat'],
                         use_direct=True,
                         direct='DATA/'+row['dirim'],
                         dir_hdu=None,
                         spec_hdu=None,
                         out_sex=None,
                         silent=False)
        assert os.path.isfile('OUTPUT/' + row['grisim'].split('.fits')[0] + '_2.cat')

    os.remove(filename)


