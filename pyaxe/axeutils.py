from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import sys
import shutil
import tempfile
from astropy.io import fits
from .axeerror import aXeError


def get_random_filename(dirname, ext):
    """Deliver a random file name"""
    fname = tempfile.mkstemp(suffix=ext, prefix=dirname+"tmp", dir=dirname)
    return fname[1]

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
