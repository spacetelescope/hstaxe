"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from hstaxe import axetasks
import os

# available contamination models
cont_models = ['gauss', 'direct', 'fluxcube', 'geometric']
inter_types = ['linear', 'polynomial', 'spline']


def test_petcont():
    """test the petcont task."""

    # These are the files in the OUTPUT dir
    infiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']
    outfiles = ['OUTPUT/ib6o23rsq_flt_2.CONT.fits',
                'OUTPUT/ib6o23ruq_flt_2.CONT.fits',
                'OUTPUT/ib6o23ryq_flt_2.CONT.fits',
                'OUTPUT/ib6o23s0q_flt_2.CONT.fits']

    for fname in infiles:
        axetasks.petcont(grism=fname,
                         config='G141.F140W.V4.31.conf',
                         cont_model='gauss',
                         model_scale=3.,
                         spec_models=None,
                         object_models=None,
                         inter_type='linear',
                         lambda_psf=800.,
                         cont_map=True,
                         in_af=None,
                         no_pet=False,
                         silent=False)

    # make some basic existence checks
    for image in outfiles:
        assert os.path.isfile(image)
        stats = os.stat(image)
        assert stats.st_size > 0

