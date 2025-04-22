.. _examples:

Example Notebooks
=================
The following are example notebooks written by HSTaXe users and developers.
They are divided by their corresponding instrument.

Wide Field Camera 3 (WFC3)
--------------------------

- The HSTaXe WFC3 cookbooks are hosted in the cookbook folder on the main HSTaXe repository:
	https://github.com/spacetelescope/hstaxe/tree/main/cookbooks/WFC3 

- To run the cookbooks, you must create a virtual environment using the `cookbook_env.yml` file provided in the cookbook folder: 
	https://github.com/spacetelescope/hstaxe/blob/main/cookbooks/cookbook_env.yml

- Once the `.yml` file is in the current working directory, the (conda) environment can be built using the command below:
	`conda env create -f cookbook_env.yml`

- The data associated with a given WFC3 cookbook can be downloaded from Box:
	https://stsci.box.com/s/qtrq6frlstvb5553opbyphu4r1rnih73

- Once the environment is created and activated, and the data are downloaded from Box, each cookbook should run from start to finish without any user input, serving as a means of validating the installation of HSTaXe. More information about the cookbooks can be found in our WFC3 Instrument Science Report (ISR) 2023-07 "HSTaXe - ACS & WFC3 Cookbook Tutorials" https://ui.adsabs.harvard.edu/abs/2023wfc..rept....7K/abstract

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
