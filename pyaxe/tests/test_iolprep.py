"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

CONVERSION NOTES

The original functions called iraf.wtraxy, the new ones combine those calls
with the translation of the rotation angle from the catalog and just save
out the new catalog


"""
import os
import numpy as np
import pytest

from astropy.table import Table
from numpy.testing.utils import assert_allclose

from hstaxe import axetasks


@pytest.fixture(scope='module')
def global_data():
    yield {"drizzle_image": "F140W_drz.fits",
           "input_cat": "cookbook.cat",
           "dimension_in": "183, 85, 50, 50",
           "config": "G141.F140W.V4.31.conf",
           "output_directory": "F140W/",
           "output_roots": ["ib6o23rtq", "ib6o23rwq", "ib6o23rzq", "ib6o23s2q"],
           "output_suffix": ["_flt.fits", "_flt_1.cat"],
           "reference_suffix": ".reference",
           }


def test_iolprep(global_data):
    """the iolprep task should create filter catalogs using
    the catalog made from the drizzled image of all the pointings"""

    os.chdir('F140W')  # axe expects the user to change directories
    axetasks.iolprep(drizzle_image=global_data['drizzle_image'],
                     input_cat=global_data['input_cat'],
                     dimension_in=global_data['dimension_in'])

    #  basic check for existance and non-zero size
    for rootname in global_data["output_roots"]:
        catalog = rootname + global_data["output_suffix"][1]
        assert os.path.isfile(catalog)
        stats = os.stat(catalog)
        assert stats.st_size > 0
    os.chdir('../')


def test_catalog_values(global_data):
    """check that the created catalogs match to expected.

    The catalogs get saved into the filter directory
    """
    os.chdir('F140W')
    for rootname in global_data['output_roots']:
        catalog = rootname + global_data["output_suffix"][1]
        print("Current catalog: {}".format(catalog))
        current_catalog = Table.read(catalog, format='ascii.sextractor')
        reference = rootname + global_data["output_suffix"][1]+ global_data["reference_suffix"]
        print("reference catalog: {}".format(reference))
        assert os.path.isfile(reference)
        reference_catalog = Table.read(reference, format='ascii.sextractor')
        # assert (len(reference_catalog) == len(current_catalog))

        # check that the x and y coordinates in the image are close to each other
        # the catalog number of the object is used here to limit confusion on ordering
        for number in reference_catalog["NUMBER"]:
            ref = reference_catalog[np.where(reference_catalog["NUMBER"] == number)]
            compare = current_catalog[np.where(current_catalog["NUMBER"] == number)]
            if len(compare) == 0:
                print(f"\nMissing object! {number}\nreference has: {ref}\n")
            assert_allclose(ref["THETA_IMAGE"], compare["THETA_IMAGE"], atol=1., rtol=1.)
            assert_allclose(ref["X_IMAGE"], compare["X_IMAGE"], atol=1e-3, rtol=1e-3)
            assert_allclose(ref["Y_IMAGE"], compare["Y_IMAGE"], atol=1e-3, rtol=1e-3)
            assert (ref["A_IMAGE"] == compare["A_IMAGE"])
            assert (ref["B_IMAGE"] == compare["B_IMAGE"])
            assert_allclose(ref["A_WORLD"], compare["A_WORLD"], atol=1e-3)
            assert_allclose(ref["B_WORLD"], compare["B_WORLD"], atol=1e-3)
            assert_allclose(ref["FLUX_RADIUS"], compare["FLUX_RADIUS"], atol=1e-3)
            assert_allclose(ref["X_WORLD"], compare["X_WORLD"], atol=1e-3)
            assert_allclose(ref["Y_WORLD"], compare["Y_WORLD"], atol=1e-3)

    os.chdir('../')



