"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
import pytest
from hstaxe import axetasks


@pytest.fixture(scope='module')
def global_data():
    yield {"config": "G141.F140W.V4.31.conf",
           "background_image": "WFC3.IR.G141.sky.V1.0.fits"}


@pytest.fixture
def axe_inputs(global_data, scope='module'):
    axe_lis_filename = "aXe.lis"
    axelist = """ib6o23rsq_flt.fits ib6o23rtq_flt_1.cat ib6o23rtq_flt.fits
ib6o23ruq_flt.fits ib6o23rwq_flt_1.cat ib6o23rwq_flt.fits
ib6o23ryq_flt.fits ib6o23rzq_flt_1.cat ib6o23rzq_flt.fits
ib6o23s0q_flt.fits ib6o23s2q_flt_1.cat ib6o23s2q_flt.fits

"""

    with open(axe_lis_filename, 'w') as fd:
            fd.write(axelist)
            fd.close()
    yield axetasks.axeinputs.aXeInput(axe_lis_filename,
                                      global_data['config'],
                                      global_data['background_image'])
    print("\nTearing down axe_inputs")
    os.remove(axe_lis_filename)


def test_sex2gol(axe_inputs, global_data):
    """test the sex2gol creates basic catalog files"""
    for row in axe_inputs:
        axetasks.sex2gol(grism='DATA/'+row['grisim'],
                         config=global_data['config'],
                         in_sex='DATA/'+row['objcat'],
                         use_direct=True,
                         direct='DATA/'+row['dirim'],
                         dir_hdu=None,
                         spec_hdu=None,
                         out_sex=None,
                         silent=False)
        assert os.path.isfile('OUTPUT/' + row['grisim'].split('.fits')[0] + '_2.cat')


