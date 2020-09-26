"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from hstaxe import axetasks


def test_pet2spc():
    """test the pet2spc task which bins the contents of
    the Pixle Extraction Table into a 1D spectrum"""

    # These are the files in the OUTPUT dir
    outfiles = ['ib6o23rsq_flt.fits']
                # 'ib6o23ruq_flt.fits',
                # 'ib6o23ryq_flt.fits',
                # 'ib6o23s0q_flt.fits']

    for fname in outfiles:
        axetasks.pet2spc(grism=fname,
                         config='G141.F140W.V4.31.conf',
                         use_bpet=False,
                         adj_sens=True,
                         weights=False,
                         do_flux=True,
                         drzpath="",
                         in_af="",
                         opet="",
                         bpet=False,
                         out_spc="",
                         silent=False)
