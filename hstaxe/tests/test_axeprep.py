"""Licensed under a 3-clause BSD style license - see LICENSE.rst.

"""
import os
import pytest
from hstaxe import axetasks


@pytest.fixture(scope="module")
def global_data():
    yield {"axe_list": "aXe.lis",
           "config": "G141.F140W.V4.31.conf",
           "background_image": "WFC3.IR.G141.sky.V1.0.fits",
           "output_directory": "OUTPUT/",
           "output_roots": ["ib6o23rsq", "ib6o23ruq", "ib6o23ryq", "ib6o23s0q"],
           "output_suffix": ["_flt.fits", "_flt_2.cat", "_flt_2.OAF", "_flt_2.MSK.fits", "_flt_2.SGRISM.fits", "_flt_2.cat"],
           "reference_suffix": "_reference."
           }

@pytest.fixture
def axe_lis_inputs(scope="module"):
    """Write out the input list that axe is expecting."""
    axelist = """ib6o23rsq_flt.fits ib6o23rtq_flt_1.cat ib6o23rtq_flt.fits
ib6o23ruq_flt.fits ib6o23rwq_flt_1.cat ib6o23rwq_flt.fits
ib6o23ryq_flt.fits ib6o23rzq_flt_1.cat ib6o23rzq_flt.fits
ib6o23s0q_flt.fits ib6o23s2q_flt_1.cat ib6o23s2q_flt.fits

"""
    with open("aXe.lis", 'w') as fd:
            fd.write(axelist)
            fd.close()
    yield
    print("\nTearing down axe_lis_inputs")
    os.remove("aXe.lis")


def test_axeprep(axe_lis_inputs, global_data):
    """test the axeprep task"""

    axetasks.axeprep(inlist=global_data["axe_list"],
                     configs=global_data["config"],
                     backims=global_data["background_image"],
                     backgr=True,
                     norm=False,
                     mfwhm=3.0)

    # make some basic existence checks
    for image in global_data["output_roots"]:
        for suffix in global_data["output_suffix"]:
            out_name = global_data["output_directory"] + image + suffix
            assert os.path.isfile(out_name)
            stats = os.stat(out_name)
            assert stats.st_size > 0


def test_output_files(global_data):
    """Compare the catalogs that were generated."""
    for root in global_data["output_roots"]:
        cat_file = global_data["output_directory"] + root + global_data["output_suffix"][1]
        splitref = cat_file.split(".")
        reference_file = splitref[0] + global_data["reference_suffix"] + splitref[1]

        cat_fh = open(cat_file)
        ref_fh = open(reference_file)

        for l1, l2 in zip(cat_fh, ref_fh):
            assert(l1 == l2)

        cat_fh.close()
        ref_fh.close()