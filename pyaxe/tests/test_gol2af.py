"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from hstaxe import axetasks


def test_gol2af():
    """test the gol2af task.

    This should turn the grism object catalog files
    into aperture files. The C function is called
    without additions from python.

    This task can be used to generate both an Object Aperture File and a
    Background Aperture File. These files have a similar format.
    """

    infiles = ['ib6o23rsq_flt.fits',
                'ib6o23ruq_flt.fits',
                'ib6o23ryq_flt.fits',
                'ib6o23s0q_flt.fits']

    for fname in infiles:
        axetasks.gol2af(grism=fname,
                        config='G141.F140W.V4.31.conf',
                        mfwhm=3.,
                        back=False,
                        orient=True,
                        slitless_geom=True,
                        exclude=False,
                        lambda_mark=4000,
                        dmag=None,
                        out_af=None,
                        in_gol=None,
                        silent=False)


        # validate the output
        outfile = "OUTPUT/"+fname.split(".")[0] + "_2.OAF"
        assert os.path.isfile(outfile)
