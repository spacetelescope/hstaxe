"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
from pyaxe import axetasks
import pytest


@pytest.fixture(scope='module')
def global_data():
    yield {"drizzle_image": "F140W_drz.fits",
           "input_cat": "cookbook.cat",
           "dimension_in": "183, 85, 50, 50",
           "config": "G141.F140W.V4.31.conf",
           "background_image": "WFC3.IR.G141.sky.V1.0.fits",
           "output_files": ['ib6o23rtq_flt_1.cat',
                            'ib6o23rwq_flt_1.cat',
                            'ib6o23rzq_flt_1.cat',
                            'ib6o23s2q_flt_1.cat']
           }


def test_iolprep(global_data):
    """test the iolprep task"""

    os.chdir('F140W')
    axetasks.iolprep(drizzle_image=global_data['drizzle_image'],
                     input_cat=global_data['input_cat'],
                     dimension_in=global_data['dimension_in'])

    #  basic check for existance and non-zero size
    for f in global_data['output_files']:
        assert os.path.isfile(f)
        stats = os.stat(f)
        assert stats.st_size > 0


def test_catalog_values():
    """check that the created catalog matches to expected."""
    pass
# if (len(rsq_catalog) == len(rsq_catalog_new)):
#     diff_catalog = copy.deepcopy(rtq_catalog)
#     for i,(r1,r2) in enumerate(zip(rsq_catalog, rsq_catalog_new)):
#         for col in r1.colnames:
#             if r1['NUMBER'] == r2['NUMBER']:
#                 diff_catalog[i][col] = r1[col] - r2[col]
#                 diff_catalog[i]['NUMBER'] = r1['NUMBER']
#             else:
#                 diff_catalog[i]['NUMBER'] = 9999