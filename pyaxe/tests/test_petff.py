"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from hstaxe import axetasks


def test_petff():
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
        axetasks.petff(grism=fname,
                       config='G141.F140W.V4.31.conf',
                       back=True,
                       ffname=None)
