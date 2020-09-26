"""A Place to store global variables and constants required to
configure axe modules.

"""
import os
import shutil

from astropy.io import fits
from hstaxe.axeerror import aXeError
from hstaxe.utils import set_logging

_log = set_logging(filename='axe_output.log')  # defaults to INFO

# defaults
__user_paths = {"AXE_IMAGE_PATH": 'DATA',
                "AXE_OUTPUT_PATH": 'OUTPUT',
                "AXE_CONFIG_PATH": 'CONF',
                "AXE_DRIZZLE_PATH": 'DRIZZLE',
                "AXE_SIMDATA_PATH": 'SIMDATA',
                "AXE_OUTSIM_PATH": 'OUTSIM',
                }

__AXE_DRZTMP_SUB = 'tmp'
__AXE_DRZTMP_LOC = os.path.join(__user_paths["AXE_DRIZZLE_PATH"],
                                __AXE_DRZTMP_SUB)
__user_paths['AXE_DRZTMP_LOC'] = __AXE_DRZTMP_LOC


# notification of python only axe
welcome_string="* Welcome to hstaxe!\nThis version is independent of IRAF and PyRAF. *"
print(f"\n{len(welcome_string)*'*'}")
print(welcome_string)
print(f"{len(welcome_string)*'*'}\n")


def set_defaults():
    """defaults for the aXe directories"""
    global __AXE_DRZTMP_LOC
    global __AXE_DRZTMP_SUB
    global __user_paths

    # override with user environment
    user_env = os.environ
    for name in __user_paths:
        if name in user_env:
            __user_paths[name] = user_env[name]
            print("{0:s} -> {1:s}".format(name, user_env[name]))
        else:
            if os.access(__user_paths[name], os.R_OK):
                print("{0:s} already exists, using.".format(name))
            else:
                print("{0:s} not defined, using {1:s}".format(name,
                                                              __user_paths[name]))
                try:
                    os.mkdir(__user_paths[name])
                except OSError:
                    raise OSError("Problem creating {0:s} -> {1:s}"
                                  .format(name, __user_paths[name]))
            os.environ[name] = __user_paths[name]


def check_axe_dirs():
    """Check for usability of all axe directories"""
    for location in __user_paths.values():
        if not os.access(location, os.W_OK):
            raise IOError("{0:s} not writeable".format(location))


# TODO: add the environment check here
def check_axesim_dirs():
    """Check for existence of all directories"""
    global __user_paths
    __user_paths["AXE_SIMDATA_PATH"] = './SIMDATA'
    __user_paths["AXE_OUTSIM_PATH"] = './OUTSIM'

    for location in __user_paths:
        if not os.access(location, os.W_OK):
            raise IOError("{0:s} not writable".format(location))


def handle_drztmp_dir():
    """Set the drizzle temporary directory"""

    # delete an already existing directory
    if os.path.isdir(__user_paths['AXE_DRZTMP_LOC']):
        print("Deleting old directory...")
        shutil.rmtree(__user_paths['AXE_DRZTMP_LOC'])

    # create a new, empty directory
    os.mkdir(__user_paths['AXE_DRZTMP_LOC'])
    print('Creating temporary directory: ', __user_paths['AXE_DRZTMP_LOC'])


# TODO: Update this for axesim
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


# def getCONF(name=None):
#     # return either AXE_CONFIG_PATH or
#     # the pathname to the input file
#     # in AXE_CONFIG_PATH
#     if name is None:
#         return __user_paths['AXE_CONFIG_PATH']
#     elif isinstance(name, list):
#         outlist=[]
#         for cname in name:
#             if 'CONF' not in cname:
#                 outlist.append(os.path.join(__user_paths['AXE_CONFIG_PATH'], cname))
#             else:
#                 outlist.append(cname)
#         return outlist
#     elif (isinstance(name,str) and ',' in name):
#         return [os.path.join(__user_paths['AXE_CONFIG_PATH'], cname) for cname in name.split(',')]
#     else:
#         return [os.path.join(__user_paths['AXE_CONFIG_PATH'], name)]

def getCONF(name=None):
    fullpath = __user_paths['AXE_CONFIG_PATH']
    if fullpath in name:
        return name
    if name is None:
        raise aXeError("No name passed to getCONF")
    if len(name.split(',')) > 1:
        namesplit=name.split(',')        
        newstring = [os.path.join(i,j) for i,j in zip([fullpath]*len(namesplit), namesplit) if fullpath not in j]
        return ','.join(newstring)
    else:
        return  os.path.join(__user_paths['AXE_CONFIG_PATH'], name)


def getDATA(name=None):
    # return either AXE_IMAGE_PATH or
    # the pathname to the input file
    # in AXE_IMAGE_PATH
    if name is None:
        return __user_paths['AXE_IMAGE_PATH']
    elif len(os.path.split(name)) > 1:
        return os.path.join(__user_paths['AXE_IMAGE_PATH'],
                            os.path.split(name)[-1])
    else:
        return os.path.join(__user_paths['AXE_IMAGE_PATH'], name)


def getOUTPUT(name=None):
    # return either AXE_OUTPUT_PATH or
    # the pathname to the output file
    # in AXE_OUTPUT_PATH
    if name is None:
        return __user_paths['AXE_OUTPUT_PATH']
    elif len(os.path.split(name)) > 1:
        return os.path.join(__user_paths['AXE_OUTPUT_PATH'],
                            os.path.split(name)[-1])
    else:
        print("adding output path to name")
        return os.path.join(__user_paths['AXE_OUTPUT_PATH'], name)


def getOUTSIM(name=None):
    # return either AXE_OUTSIM_PATH or
    # the pathname to the output file
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
        return __user_paths['AXE_DRIZZLE_PATH']
    else:
        return os.path.join(__user_paths['AXE_DRIZZLE_PATH'], name)


def getDRZTMP(name=None):
    # return either AXE_DRZTMP_LOC or
    # the pathname to the input file
    # in AXE_DRZTMP_LOC
    global __AXE_DRZTMP_LOC
    if name is None:
        return __AXE_DRZTMP_LOC
    else:
        return os.path.join(__AXE_DRZTMP_LOC, name)

def get_ext_info(image, conf):
    """Determines the extension information on an image.

    Parameters
    ----------
    conf: configfile.ConfigFile

    Returns
    -------
    ext_info: int
        Extension information

    """
    # initialize a dictionary
    ext_info = {}

    # set default values
    ext_info['ext_name'] = None
    ext_info['axe_ext'] = None
    ext_info['fits_ext'] = None
    ext_info['ext_version'] = None

    # open the image
    with fits.open(image, 'readonly') as fits_image:

        # check the keyword for string-like
        if isinstance(conf.get_gvalue('SCIENCE_EXT'), str):
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
                        (fits_image[index].header['EXTNAME'] == ext_info['ext_name'])):

                        # check whether OPTKEY1 and OPTKEYVAL1 fits
                        if optkey in fits_image[index].header and fits_image[index].header[optkey] == optval:

                            # set the extension numbers
                            ext_info['fits_ext'] = index
                            ext_info['axe_ext'] = index + 1
            else:
                for index in range(len(fits_image)):
                    if (('EXTNAME' in fits_image[index].header) and
                        (fits_image[index].header['EXTNAME'] == ext_info['ext_name'])):
                        ext_info['fits_ext'] = index
                        ext_info['axe_ext']  = index + 1
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
    axe_names['GOL'] = ("{0:s}_{1:d}.cat".format(root, ext_info['axe_ext']))

    # the OAF/BAF:
    axe_names['OAF'] = ("{0:s}_{1:d}.OAF".format(root, ext_info['axe_ext']))
    axe_names['BAF'] = ("{0:s}_{1:d}.BAF".format(root, ext_info['axe_ext']))

    # the PET:
    axe_names['PET'] = ("{0:s}_{1:d}.PET.fits"
                        .format(root, ext_info['axe_ext']))
    axe_names['BCK_PET'] = ("{0:s}_{1:d}.BCK.PET.fits"
                            .format(root, ext_info['axe_ext']))

    # the DPP:
    axe_names['DPP'] = ("{0:s}_{1:d}.DPP.fits"
                        .format(root, ext_info['axe_ext']))
    axe_names['BCK_DPP'] = ("{0:s}_{1:d}.BCK.DPP.fits"
                            .format(root, ext_info['axe_ext']))

    # the SPC:
    axe_names['SPC'] = ("{0:s}_{1:d}SPC.fits"
                        .format(root, ext_info['axe_ext']))

    # the STP:
    axe_names['STP'] = ("{0:s}_{1:d}.STP.fits"
                        .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['MSK'] = ("{0:s}_{1:d}.MSK.fits"
                        .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['NBCK'] = ("{0:s}_{1:d}.NBCK.fits"
                         .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['SGRI'] = ("{0:s}_{1:d}.SGRISM.fits"
                         .format(root, ext_info['axe_ext']))

    # the background mask:
    axe_names['FLX'] = ("{0:s}_{1:d}FLX.fits"
                        .format(root, ext_info['axe_ext']))

    # return the dictionary
    return axe_names


def isstringlike(item):
    """Checks whether a term is a string or not"""
    return isinstance(item, str)


def get_random_filename(dirname, ext):
    """Deliver a random file name."""
    import uuid
    return dirname + 'tmp' + uuid.uuid4().hex + ext


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
