"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
import pytest
from hstaxe import axetasks
from hstaxe.axesrc import axesingextr


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


def test_axecore(axe_inputs, global_data):
    """test the axecore task"""

    for row in axe_inputs:
        runner = axesingextr.aXeSpcExtr(row['grisim'],
                                        row['objcat'],
                                        row['dirim'],
                                        global_data['config'],
                                        row['dmag'],
                                        extrfwhm=4.,
                                        drzfwhm=3.,
                                        backfwhm=0.,
                                        orient=False,
                                        weights=True,
                                        slitless_geom=False,
                                        exclude=True,
                                        lambda_mark=None,
                                        cont_model='gauss',
                                        model_scale=None,
                                        inter_type='linear',
                                        lambda_psf=None,
                                        np=None,
                                        interp=None,
                                        niter_med=None,
                                        niter_fit=None,
                                        kappa=None,
                                        smooth_length=None,
                                        smooth_fwhm=None,
                                        spectr=True,
                                        adj_sens=True,
                                        sampling='drizzle')
        runner.run()
