.. _examples:

Example Notebooks
=================
The following are example notebooks written by HSTaXe users and developers.
They are divided by their corresponding instrument.

Wide Field Camera 3 (WFC3)
--------------------------

The aXe WFC3 cookbook data can be downloaded by cloning this repository:
https://github.com/npirzkal/aXe_WFC3_Cookbook

1. git clone https://github.com/npirzkal/aXe_WFC3_Cookbook.git

2. Enter the `aXe_WFC3_Cookbook` directory and work through
   `aXe_WFC3_Cookbook.ipynb` within the conda environment created as
   part of the `hstaxe` installation.

Advanced Camera for Surveys
---------------------------

The aXe ACS cookbook and associated data can be downloaded from the
following area: https://stsci.box.com/s/eo98zjtyccnoq7z73akfrx94jog3pg7j

aXe Cookbook Examples
---------------------

There are also a couple simple test scripts in hstaxe/tests that can be
run against the example data in the aXe cookbook:

  * `run_acs_cookbook.py`: runs basic aXe against ACS data

  * `run_cookbook.py`: runs basic aXe against WFC3 data from the cookbook

  * `run_cookbook_part2.py`: runs axedrizzle against WFC3 data from the cookbook


Other Uses of aXe
-----------------
aXe has been used successfully in several large science programs, such
as GRAPES (ACS/WFC, [PIRZKAL1]_ ) and PEARS (ACS/WFC, [PIRZKAL2]_ ). The aXe software was
central in extracting ACS/G800L and, using a customized version of aXe,
NICMOS/G141 data within the corresponding Hubble Legacy Archive (HLA)
projects (see [FREUDLING]_ and [KUMMEL4]_).
