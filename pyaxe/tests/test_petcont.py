"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from hstaxe import axetasks

# available contamination models
cont_models = ['gauss', 'direct', 'fluxcube', 'geometric']
inter_types = ['linear', 'polynomial', 'spline']


def test_petcont():
    """test the petcont task.

    This should turn the translated catalog files
    into grism object list (GOL) files that still
    look like source extractor output catalogs.
    """

    # These are the files in the OUTPUT dir
    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in outfiles:
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
