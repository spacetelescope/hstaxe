from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os


"""A Place to store global variables and constants required to
configure axe modules. Use the teal interface to ConfigObj to return a
config object cfgobj=teal.teal('axe',loadonly=True)

"""

# use these as defaults if the user has not already set them
global GLOB_VARS_SET
global AXE_IMAGE_PATH
global AXE_OUTPUT_PATH
global AXE_CONFIG_PATH
global AXE_DRIZZLE_PATH
global AXE_SIMDATA_PATH
global AXE_OUTSIM_PATH
global AXE_BINDIR

global AXE_DRZTMP_SUB
global AXE_DRZTMP_LOC

global GLOB_VARS_SET


GLOB_VARS_SET = False
AXE_IMAGE_PATH = './IMAGE'
AXE_OUTPUT_PATH = './OUTPUT'
AXE_CONFIG_PATH = './CONFIG'
AXE_DRIZZLE_PATH = './DRIZZLE'
AXE_SIMDATA_PATH = './SIMDATA'
AXE_OUTSIM_PATH = './OUTSIM'
AXE_BINDIR = './BINDIR'

AXE_DRZTMP_SUB = 'tmp'
AXE_DRZTMP_LOC = os.path.join(AXE_DRIZZLE_PATH, AXE_DRZTMP_SUB)

GLOB_VARS_SET = True
