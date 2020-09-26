"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
from hstaxe import axetasks


def test_stamps():
    """test the aXe_stamps task

    Notes
    -----
    SPC files contained 1D spectra, opt.SPC files contained optimally extracted
    spectra using gaussian profiles STP files contain 2D stamps.
    CONT files contain the gaussian based contamination estimate

    This test should produce STP FITS files where each stamp is contained in
    its own extension. All output files should be in the directory pointed to
    by the AXE_OUTPUT_PATH environment variable.
    """

    # These are the files in the OUTPUT dir
    outfiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in outfiles:
        axetasks.stamps(grism=fname,
                        config='G141.F140W.V4.31.conf',
                        sampling='rectified',
                        drzpath=False,
                        in_af="",
                        in_pet=None,
                        out_stp=None,
                        silent=False)
