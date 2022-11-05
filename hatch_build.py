#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# install cfitsio, gsl, wcstools, astropy, drizzlepac
#
#
import sys
import sysconfig
import os
import pathlib

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from subprocess import check_call, CalledProcessError

if (sys.version_info < (3, 7)):
    sys.stderr.write("ERROR: hstaxe requires Python 3.7 or later\n")
    sys.exit(1)


AXELIB_DIR = "cextern/"
CONF_H_NAME = os.path.join(AXELIB_DIR, "config.h")


def find_program(name):
    paths = os.environ.get("PATH")
    if not paths:
        raise EnvironmentError("PATH is not set")

    for path in paths.split(":"):
        if not os.path.exists(path):
            continue
        for program in os.listdir(path):
            program_path = os.path.join(path, program)
            if name == program:
                return program_path
    return ""


def check_pkgconfig(name):
    pkg_configs = []
    ext = ".pc"
    target = name
    result = {
        "usable": False,
        "search_path": "",
        "file": "",
    }

    # Append file extension to target if not present
    if not target.endswith(ext):
        target += ext

    # Determine where pkg-config lives and generate a search path
    pkg_config_prefix = os.path.dirname(os.path.dirname(find_program("pkg-config")))
    if pkg_config_prefix:
        result["usable"] = True
        pkg_configs += [os.path.join(pkg_config_prefix, "lib64", "pkgconfig")]
        pkg_configs += [os.path.join(pkg_config_prefix, "lib", "pkgconfig")]
        pkg_configs += [os.path.join(pkg_config_prefix, "share", "pkgconfig")]

    # pkg-config might not be local to Python's runtime environment
    # Generate a search path for it as well
    if sys.prefix != pkg_config_prefix:
        pkg_configs += [os.path.join(sys.prefix, "lib64", "pkgconfig")]
        pkg_configs += [os.path.join(sys.prefix, "lib", "pkgconfig")]
        pkg_configs += [os.path.join(sys.prefix, "share", "pkgconfig")]
        
    # Consume user defined pkg-config search path
    pkg_configs_user = os.environ.get("PKG_CONFIG_PATH", "")
    if pkg_configs_user:
        pkg_configs = [*pkg_configs_user.split(":"), *pkg_configs]

    # Is the target file present in any of the search paths?
    for path in pkg_configs:
        if result["file"]:
            # the previous iteration found something
            break

        if not os.path.exists(path):
            continue

        for config in os.listdir(path):
            config_path = os.path.join(path, config)
            if config == target:
                result["file"] = config_path
                # only store where we found the file. let the caller
                # handle PKG_CONFIG_PATH modifications
                result["search_path"] = os.path.dirname(config_path)
                break

    return result


def find_in_file(filename, s):
    with open(filename, "r") as fp:
        for line in fp.read().splitlines():
            if s in line:
                return line
        return ""


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
    CONFIGURE_ARGS = ["--prefix="+CURRENT_ENV,
                      "--with-cfitsio="+CURRENT_ENV,
                      "--with-gsl="+CURRENT_ENV,
                      "--libdir="+sysconfig.get_config_var("LIBDIR")]

    PKG_CONFIG_PATH = []
    if os.environ.get("PKG_CONFIG_PATH", ""):
        PKG_CONFIG_PATH = os.environ["PKG_CONFIG_PATH"].split(":")

    WCSTOOLS = check_pkgconfig("wcstools")
    WCSTOOLS_USABLE = WCSTOOLS["usable"]
    WCSTOOLS_FILE = WCSTOOLS["file"]
    WCSTOOLS_SEARCH_PATH = WCSTOOLS["search_path"]

    if not WCSTOOLS_USABLE:
        # pkg-config is not installed
        # configure options manually
        CONFIGURE_ARGS += ["--with-wcstools="+CURRENT_ENV]

        if WCSTOOLS_FILE and find_in_file(WCSTOOLS_FILE, "-lwcstools"):
            # we *do* have a .pc that wants libwcstools
            CONFIGURE_ARGS += ["--with-wcstools-libname=wcstools"]
        else:
            # otherwise this must be a vanilla wcstools from upstream
            # with "libwcs", and no .pc file
            CONFIGURE_ARGS += ["--with-wcstools-libname=wcs"]
    else:
        # pkg-config is installed
        if WCSTOOLS_SEARCH_PATH:
            # prepend the search path
            PKG_CONFIG_PATH.insert(0, WCSTOOLS_SEARCH_PATH)

    # apply the search path to the environment
    os.environ["PKG_CONFIG_PATH"] = ":".join(PKG_CONFIG_PATH)
    try:
        check_call(["sh", "autogen.sh"],
                   cwd=AXELIB_DIR)
        check_call(["sh", "./configure",
                    *CONFIGURE_ARGS],
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

