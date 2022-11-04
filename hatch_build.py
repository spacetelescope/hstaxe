#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# install cfitsio, gsl, wcstools, astropy, drizzlepac
#
#
import sys
import os
import pathlib

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from subprocess import check_call, CalledProcessError

if (sys.version_info < (3, 7)):
    sys.stderr.write("ERROR: hstaxe requires Python 3.7 or later\n")
    sys.exit(1)


AXELIB_DIR = "cextern/"
CONF_H_NAME = os.path.join(AXELIB_DIR, "config.h")

def clean_executables():
    current_env = sys.prefix + "/bin/"
    axe_bins = ["aXe_GOL2AF",
                "aXe_AF2PET",
                "aXe_BE",
                "aXe_PET2SPC",
                "aXe_STAMPS",
                "aXe_DRZPREP",
                "aXe_PETCONT",
                "aXe_PETFF",
                "aXe_DRZ2PET",
                "aXe_GPS",
                "aXe_FILET",
                "aXe_FRIGEN",
                "aXe_FRINGECORR",
                "aXe_TFIT",
                "aXe_INTPIXCORR",
                "aXe_PETIPC",
                "aXe_NICBACK",
                "aXe_TEST",
                "aXe_DIRIMAGE",
                "aXe_SEX2GOL",
                "aXe_SCALEBCK"]

    try:
        check_call(["make", "clean"], cwd=AXELIB_DIR)
    except CalledProcessError as e:
        print(e)
        exit(1)

    for file in axe_bins:
        myfile = current_env + file
        if os.access(myfile, os.F_OK):
            os.remove(myfile)
    if os.access(CONF_H_NAME, os.F_OK):
        os.remove(CONF_H_NAME)

def BuildExtWithConfigure():
    CURRENT_ENV = sys.prefix
    try:
        check_call(["sh", "autogen.sh"],
                   cwd=AXELIB_DIR)
        check_call(["sh", "./configure",
                    "--with-cfitsio="+CURRENT_ENV,
                    "--with-wcstools="+CURRENT_ENV,
                    "--with-gsl="+CURRENT_ENV,
                    "--libdir="+CURRENT_ENV,
                    "--prefix="+CURRENT_ENV,
                    "--with-wcstools-libname=wcstools"],
                   cwd=AXELIB_DIR)
        check_call(["make", "clean"], cwd=AXELIB_DIR)
        check_call(["make", "install"], cwd=AXELIB_DIR)
    except CalledProcessError as e:
        print(e)
        exit(1)

class CustomBuildHook(BuildHookInterface):
    """Configure, build, and install the aXe C code."""
    description = "Configure, build, and install the aXe C code."

    def clean(self, version):
        # Runs if -c / --clean flag was passed to the build command
        clean_executables()

    def initialize(self, version, build_data):
        # Runs before the rest of the build
        BuildExtWithConfigure()

