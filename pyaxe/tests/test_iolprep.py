"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks
import pytest
from astropy.table import Table


@pytest.fixture(scope='module')
def global_data():
    yield {"drizzle_image": "F140W_drz.fits",
           "input_cat": "cookbook.cat",
           "dimension_in": "183, 85, 50, 50",
           "config": "G141.F140W.V4.31.conf",
           "output_directory": "F140W/",
           "output_roots": ["ib6o23rtq", "ib6o23rwq", "ib6o23rzq", "ib6o23s2q"],
           "output_suffix": ["_flt.fits", "_flt_1.cat"],
           "reference_suffix": "_reference.",
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
    for catalog, reference in global_data['output_files']:
        current_catalog = Table.read(catalog, format='ascii.sextractor')
        assert os.path.isfile(reference)
        reference_catalog = Table.read(reference, format='ascii.sextractor')
        assert all(True for item in current_catalog == reference_catalog)

    os.chdir('../')


""" CONVERSION NOTES

The original functions called iraf.wtraxy