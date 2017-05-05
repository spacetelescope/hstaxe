from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import sys
import shutil
from astropy.io import fits
from .axeerror import aXeError

# boolean to indicate that
# the variables have been set
GLOB_VARS_SET = False

# defaults for the
# aXe directories
AXE_IMAGE_PATH = './'
AXE_OUTPUT_PATH = './'
AXE_CONFIG_PATH = './'
AXE_DRIZZLE_PATH = './'
AXE_SIMDATA_PATH = './'
AXE_OUTSIM_PATH = './'

# name for the drizzle
# tmp-directory
AXE_DRZTMP_SUB = 'tmp'
AXE_DRZTMP_LOC = None


def safe_mkdir(s):
    # I just want the directory to exist - I don't care how it got there.
    try:
        os.mkdir(s)
    except OSError:
        pass


def check_axe_dirs():
    """Check for existence of all directories"""
    safe_mkdir(AXE_IMAGE_PATH)
    safe_mkdir(AXE_OUTPUT_PATH)
    safe_mkdir(AXE_CONFIG_PATH)
    safe_mkdir(AXE_DRIZZLE_PATH)


def check_axesim_dirs():
    """Check for existence of all directories"""
    safe_mkdir(AXE_IMAGE_PATH)
    safe_mkdir(AXE_OUTPUT_PATH)
    safe_mkdir(AXE_CONFIG_PATH)
    safe_mkdir(AXE_SIMDATA_PATH)
    safe_mkdir(AXE_OUTSIM_PATH)


def handle_drztmp_dir():
    """Set the drizzle temporary directory"""
    global AXE_DRZTMP_LOC

    # get the path name of the tmp directory
    AXE_DRZTMP_LOC = os.path.join(AXE_DRIZZLE_PATH, AXE_DRZTMP_SUB)

    # delete an already existing directory
    if os.path.isdir(AXE_DRZTMP_LOC):
        shutil.rmtree(AXE_DRZTMP_LOC)

    # create a new, empty directory
    os.mkdir(AXE_DRZTMP_LOC)
    print('Creating temporary directory: ', AXE_DRZTMP_LOC)


def axe_setup(tmpdir=False, axesim=False):
    """Setup the aXe file and pathnames"""
    global GLOB_VARS_SET
    global AXE_IMAGE_PATH
    global AXE_OUTPUT_PATH
    global AXE_CONFIG_PATH
    global AXE_DRIZZLE_PATH
    global AXE_SIMDATA_PATH
    global AXE_OUTSIM_PATH
    global AXE_BINDIR

    # set the error counter
    ret = 0

    if GLOB_VARS_SET and not tmpdir:
        return 0

    # set the AXE_IMAGE_PATH;
    # keep default if not explicitly given
    try:
        AXE_IMAGE_PATH = os.environ['AXE_IMAGE_PATH']
    except KeyError:
        print('AXE_IMAGE_PATH not defined, using default.')
        ret += 1

    # set the AXE_CONFIG_PATH;
    # keep default if not explicitly given
    try:
        AXE_CONFIG_PATH = os.environ['AXE_CONFIG_PATH']
    except KeyError:
        print('AXE_CONFIG_PATH not defined, using default.')
        ret += 1

    # set the AXE_OUTPUT_PATH;
    # keep default if not explicitly given
    try:
        AXE_OUTPUT_PATH = os.environ['AXE_OUTPUT_PATH']
    except KeyError:
        print('AXE_OUTPUT_PATH not defined, using default.')
        ret += 1

    # check whether we are
    # in axesim
    if axesim:
        # set the AXE_SIMDATA_PATH;
        # keep default if not explicitly given
        try:
            AXE_SIMDATA_PATH = os.environ['AXE_SIMDATA_PATH']
        except KeyError:
            print('AXE_SIMDATA_PATH not defined, using default.')
            ret += 1

        # set the AXE_OUTSIM_PATH;
        # keep default if not explicitly given
        try:
            AXE_OUTSIM_PATH = os.environ['AXE_OUTSIM_PATH']
        except KeyError:
            print('AXE_OUTSIM_PATH not defined, using default.')
            ret += 1

        # make sure all directories
        # do exist
        check_axesim_dirs()
    else:
        # set the AXE_DRIZZLE_PATH;
        # keep default if not explicitly given
        try:
            AXE_DRIZZLE_PATH = os.environ['AXE_DRIZZLE_PATH']
        except KeyError:
            print('AXE_DRIZZLE_PATH not defined, using default.')
            ret += 1

        # make sure all directories
        # do exist
        check_axe_dirs()

        if tmpdir:
            # deal with the drizzle tmp directory
            handle_drztmp_dir()

    # define the path to the binaries
    AXE_BINDIR = get_axebindir()

    # set the boolean
    GLOB_VARS_SET = True

    # return number of defaults
    return ret


def get_axebindir():
    """Define the path to the aXe binaries"""

    if 'axesrc' in sys.modules:
        modfile = sys.modules['axesrc'].__file__
        axebindir = os.path.abspath(os.path.join(os.path.dirname(modfile),
                                                 '../bin/'))
    return axebindir


def getBINDIR(name=None):
    # return either AXE_BINDIR or
    # the pathname to the input file
    # in AXE_BINDIR
    if name is None:
        return AXE_BINDIR
    else:
        return os.path.join(AXE_BINDIR, name)


def getCONF(name=None):
    # return either AXE_CONFIG_PATH or
    # the pathname to the input file
    # in AXE_CONFIG_PATH
    if name is None:
        return AXE_CONFIG_PATH
    else:
        return os.path.join(AXE_CONFIG_PATH, name)


def getIMAGE(name=None):
    # return either AXE_IMAGE_PATH or
    # the pathname to the input file
    # in AXE_IMAGE_PATH
    if name is None:
        return AXE_IMAGE_PATH
    else:
        return os.path.join(AXE_IMAGE_PATH, name)


def getOUTPUT(name=None):
    # return either AXE_OUTPUT_PATH or
    # the pathname to the input file
    # in AXE_OUTPUT_PATH
    if name is None:
        return AXE_OUTPUT_PATH
    else:
        return os.path.join(AXE_OUTPUT_PATH, name)


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
        return AXE_DRIZZLE_PATH
    else:
        return os.path.join(AXE_DRIZZLE_PATH, name)


def getDRZTMP(name=None):
    # return either AXE_DRZTMP_LOC or
    # the pathname to the input file
    # in AXE_DRZTMP_LOC
    if name is None:
        return AXE_DRZTMP_LOC
    else:
        return os.path.join(AXE_DRZTMP_LOC, name)


def get_ext_info(image, conf):
    """Determines the extension information on an image"""
    # initialize a dictionary
    ext_info = {}

    # set default values
    ext_info['ext_name'] = None
    ext_info['axe_ext'] = None
    ext_info['fits_ext'] = None
    ext_info['ext_version'] = None

    # open the image
    fits_image = fits.open(image, 'readonly')

    # check the keyword for string-like
    if isstringlike(conf.get_gvalue('SCIENCE_EXT')):
        # set the extension name
        ext_info['ext_name'] = conf.get_gvalue('SCIENCE_EXT')

        # check whether OPTKEY1 and OPTVAL1 are set
        if ((conf.get_gvalue('OPTKEY1') is not None) and
            (conf.get_gvalue('OPTVAL1') is not None)):

            # store the key
            optkey = conf.get_gvalue('OPTKEY1')

            # store the value, convert
            # CCDCHIP value to integer type
            if optkey == 'CCDCHIP':
                optval = int(conf.get_gvalue('OPTVAL1'))
            else:
                optval = conf.get_gvalue('OPTVAL1')

            # go over all fits extensions
            for index in range(len(fits_image)):
                # check whether the extension name fits
                if (('EXTNAME' in fits_image[index].header) and
                    (fits_image[index].header['EXTNAME'] is ext_info['ext_name'])):

                    # check whether OPTKEY1 and OPTKEYVAL1 fits
                    if ((optkey in fits_image[index].header) and
                        (fits_image[index].header[optkey] is optval)):

                        # set the extension numbers
                        ext_info['fits_ext'] = index
                        ext_info['axe_ext'] = index + 1
        else:
            for index in range(len(fits_image)):
                if (('EXTNAME' in fits_image[index].header) and
                    (fits_image[index].header['EXTNAME'] is ext_info['ext_name'])):
                    ext_info['fits_ext'] = index
                    ext_info['axe_ext'] = index + 1
    else:
        # set the axe extension number an axe extension number
        ext_info['axe_ext'] = int(conf.get_gvalue('SCIENCE_EXT'))
        ext_info['fits_ext'] = int(conf.get_gvalue('SCIENCE_EXT')) - 1

        # get extension name, if possible
        if 'EXTNAME' in fits_image[ext_info['fits_ext']].header:
            ext_info['ext_name'] = fits_image[ext_info['fits_ext']].header['EXTNAME']

    # check for success, complain and out if not
    if ext_info['axe_ext'] is None:
        err_msg = 'Unable to find the specified extensions in image %s!' % image
        raise aXeError(err_msg)

    # get extension version, if possible
    if 'EXTVER' in fits_image[ext_info['fits_ext']].header:
        ext_info['ext_version'] = fits_image[ext_info['fits_ext']].header['EXTVER']

    # close the fits image
    fits_image.close

    # return the dictionary
    return ext_info


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
    axe_names['GOL'] = '%s_%i.cat' % (root, ext_info['axe_ext'])

    # the OAF/BAF:
    axe_names['OAF'] = '%s_%i.OAF' % (root, ext_info['axe_ext'])
    axe_names['BAF'] = '%s_%i.BAF' % (root, ext_info['axe_ext'])

    # the PET:
    axe_names['PET'] = '%s_%i.PET.fits' % (root, ext_info['axe_ext'])
    axe_names['BCK_PET'] = '%s_%i.BCK.PET.fits' % (root, ext_info['axe_ext'])

    # the DPP:
    axe_names['DPP'] = '%s_%i.DPP.fits' % (root, ext_info['axe_ext'])
    axe_names['BCK_DPP'] = '%s_%i.BCK.DPP.fits' % (root, ext_info['axe_ext'])

    # the SPC:
    axe_names['SPC'] = '%s_%i.SPC.fits' % (root, ext_info['axe_ext'])

    # the STP:
    axe_names['STP'] = '%s_%i.STP.fits' % (root, ext_info['axe_ext'])

    # the background mask:
    axe_names['MSK'] = '%s_%i.MSK.fits' % (root, ext_info['axe_ext'])

    # the background mask:
    axe_names['NBCK'] = '%s_%i.NBCK.fits' % (root, ext_info['axe_ext'])

    # the background mask:
    axe_names['SGRI'] = '%s_%i.SGRISM.fits' % (root, ext_info['axe_ext'])

    # the background mask:
    axe_names['FLX'] = '%s_%i.FLX.fits' % (root, ext_info['axe_ext'])

    # return the dictionary
    return axe_names


def isstringlike(item):
    """Checks whether a term is a string or not"""
    ret = 1
    try:
        float(item)
        ret = 0
    except ValueError:
        pass
    return ret


def is_quant_contam(contam_model):
    """Get the flag for quantitative contamination"""
    # the list of quantitative models
    quant_models = ['GAUSS', 'FLUXCUBE']

    # set the default value
    isquantcont = True

    # check whether the contamination is not quantitative
    if not contam_model.upper() in quant_models:
        # re-set the flag
        isquantcont = False

    # return the flag
    return isquantcont


def get_random_filename(dirname, ext):
    """Deliver a random file name"""
    import random

    # assure a first go in the while loop
    found = 1

    # do until you find a unique name
    while found:

        # get a random int number
        str_num = str(random.randint(10000, 99999))

        # compose a random name
        fname = dirname + 'tmp' + str_num + ext

        # check whether the file exists
        if not os.path.isfile(fname):
            found = 0

    return fname
