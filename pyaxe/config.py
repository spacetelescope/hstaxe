from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

"""A Place to store global variables and constants required to
configure axe modules. Use the teal interface to ConfigObj to return a
config object cfgobj=teal.teal('axe',loadonly=True)

"""
import os
import tempfile
from .axeException import DirError
from astropy import config

# defaults
_user_paths = {"AXE_IMAGE_PATH": 'IMAGE',
               "AXE_OUTPUT_PATH": 'OUTPUT',
               "AXE_CONFIG_PATH": 'CONFIG',
               "AXE_DRIZZLE_PATH": 'DRIZZLE'}

__AXE_DRZTMP_SUB = 'tmp'
__AXE_BINDIR = './cextern/aXe_c_code/'
__AXE_DRZTMP_LOC = os.path.join(_user_paths["AXE_DRIZZLE_PATH"],
                                            __AXE_DRZTMP_SUB)


def set_defaults():
    """defaults for the aXe directories"""

    # override with user environment
    user_env = os.environ
    for name in _user_paths:
        if name in user_env:
            _user_paths[name] = user_env[name]
            print("{0:s} -> {1:s}".format(name, user_env[name]))
        else:
            print("{0:s} not defined, using {1:s}".format(name,
                                                          _user_paths[name]))

            if os.access(_user_paths[name], os.R_OK):
                print("{0:s} already exists".format(name))
            else:
                try:
                    os.mkdir(_user_paths[name])
                except OSError:
                    raise OSError("Problem creating {0:s} -> {1:s}".format(name, _user_paths[name]))

    __AXE_DRZTMP_LOC = os.path.join(_user_paths["AXE_DRIZZLE_PATH"],
                                                __AXE_DRZTMP_SUB)


def check_axe_dirs():
    """Check for usability of all axe directories"""
    for location in _user_paths.values():
        if not os.access(location, os.W_OK):
            raise IOError("{0:s} not writeable".format(location))


def check_axesim_dirs():
    """Check for existence of all directories"""
    user_env = os.environ
    _user_paths["AXE_SIMDATA_PATH"] = './SIMDATA'
    _user_paths["AXE_OUTSIM_PATH"] = './OUTSIM'

    for location in _user_paths:
        if not os.access(location, os.W_OK):
            raise IOError("{0:s} not writable".format(location))


def handle_drztmp_dir():
    """Set the drizzle temporary directory"""

    # delete an already existing directory
    if os.path.isdir(_user_paths['AXE_DRZTMP_LOC']):
        print("Deleting old directory...")
        shutil.rmtree(_user_paths['AXE_DRZTMP_LOC'])

    # create a new, empty directory
    os.mkdir(_user_paths['AXE_DRZTMP_LOC'])
    print('Creating temporary directory: ', _user_paths['AXE_DRZTMP_LOC'])


def axe_setup(tmpdir=False, axesim=False):
    """Setup the aXe file and pathnames"""

    # check whether we are
    # in axesim
    if axesim:
        # set the AXE_SIMDATA_PATH;
        # keep default if not explicitly given
        try:
            AXE_SIMDATA_PATH = os.environ['AXE_SIMDATA_PATH']
        except KeyError:
            print('AXE_SIMDATA_PATH not defined, using default.')

        # set the AXE_OUTSIM_PATH;
        # keep default if not explicitly given
        try:
            AXE_OUTSIM_PATH = os.environ['AXE_OUTSIM_PATH']
        except KeyError:
            print('AXE_OUTSIM_PATH not defined, using default.')

        # make sure all directories
        # do exist
        check_axesim_dirs()
    else:

        # make sure all directories
        # do exist
        check_axe_dirs()

        if tmpdir:
            # deal with the drizzle tmp directory
            handle_drztmp_dir()


def getCONF(name=None):
    # return either AXE_CONFIG_PATH or
    # the pathname to the input file
    # in AXE_CONFIG_PATH
    if name is None:
        return _user_paths['AXE_CONFIG_PATH']
    else:
        return os.path.join(_user_paths['AXE_CONFIG_PATH'], name)


def getIMAGE(name=None):
    # return either AXE_IMAGE_PATH or
    # the pathname to the input file
    # in AXE_IMAGE_PATH
    if name is None:
        return _user_paths['AXE_IMAGE_PATH']
    else:
        return os.path.join(_user_paths['AXE_IMAGE_PATH'], name)


def getOUTPUT(name=None):
    # return either AXE_OUTPUT_PATH or
    # the pathname to the input file
    # in AXE_OUTPUT_PATH
    if name is None:
        return _user_paths['AXE_OUTPUT_PATH']
    else:
        return os.path.join(_user_paths['AXE_OUTPUT_PATH'], name)


def getOUTSIM(name=None):
    # return either AXE_OUTSIM_PATH or
    # the pathname to the input file
    # in AXE_OUTSIM_PATH
    if name is None:
        return AXE_OUTSIM_PATH
    else:
        return os.path.join(AXE_OUTSIM_PATH, name)


def getSIMDATA(name=None):
    # return either AXE_SIMDATA_PATH or
    # the pathname to the input file
    # in AXE_SIMDATA_PATH
    if name is None:
        return AXE_SIMDATA_PATH
    else:
        return os.path.join(AXE_SIMDATA_PATH, name)


def getDRIZZLE(name=None):
    # return either AXE_DRIZZLE_PATH or
    # the pathname to the input file
    # in AXE_DRIZZLE_PATH
    if name is None:
        return _user_paths['AXE_DRIZZLE_PATH']
    else:
        return os.path.join(_user_paths['AXE_DRIZZLE_PATH'], name)


def getDRZTMP(name=None):
    # return either AXE_DRZTMP_LOC or
    # the pathname to the input file
    # in AXE_DRZTMP_LOC
    if name is None:
        return __AXE_DRZTMP_LOC
    else:
        return os.path.join(__AXE_DRZTMP_LOC, name)


def get_axe_names(image, ext_info):
    """Derive the name of all aXe products for a given image"""
    # make an empty dictionary
    axe_names = {}

    # get the root of the image name
    pos = image.rfind('.fits')
    root = image[:pos]

    # FILL the dictionary with names of aXe products
    #
    # the GOL:
    axe_names['GOL'] = ("{0:s}_{1:i}.cat".format(root, ext_info['axe_ext']))

    # the OAF/BAF:
    axe_names['OAF'] = ("{0:s}_{1:i}.OAF".format(root, ext_info['axe_ext']))
    axe_names['BAF'] = ("{0:s}_{1:i}.BAF".format(root, ext_info['axe_ext']))

    # the PET:
    axe_names['PET'] = ("{0:s}_{1:i}.PET.fits"
                        .format(root, ext_info['axe_ext']))
    axe_names['BCK_PET'] = ("{0:s}_{1:i}.BCK.PET.fits"
                            .format(root, ext_info['axe_ext']))

    # the DPP:
    axe_names['DPP'] = ("{0:s}_{1:i}.DPP.fits"
                        .format(root, ext_info['axe_ext']))
    axe_names['BCK_DPP'] = ("{0:s}_{1:i}.BCK.DPP.fits"
                            .format(root, ext_info['axe_ext']))

    # the SPC:
    axe_names['SPC'] = ("{0:s}_{1:i}SPC.fits"
                        .format(root, ext_info['axe_ext']))

    # the STP:
    axe_names['STP'] = ("{0:s}_{1:i}.STP.fits"
                        .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['MSK'] = ("{0:s}_{1:i}.MSK.fits"
                        .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['NBCK'] = ("{0:s}_{1:i}.NBCK.fits"
                         .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['SGRI'] = ("{0:s}_{1:i}.SGRISM.fits"
                         .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['FLX'] = ("{0:s}_{1:i}FLX.fits"
                        .format(root, ext_info['axe_ext']))

    # return the dictionary
    return axe_names


def isstringlike(item):
    """Checks whether a term is a string or not"""
    return isinstance(item, str)


def get_random_filename(dirname, ext):
    """Deliver a random file name"""
    return tempfile.mkstemp(suffix=ext, prefix=dirname+"tmp", dir=dirname)
