import os
import logging
from astropy.io import fits

from hstaxe.axeerror import aXeError
import hstaxe.config as config_util

from . import configfile
from . import axelowlev
from . import axeinputs


# make sure there is a logger
_log = logging.getLogger(__name__)

class DPPdumps(object):
    """Class to intitially handle all DPP files"""
    def __init__(self, inima, confterm, back=False):
        self.inima = inima
        self.confterm = confterm

        # find all ddp file names
        self.dpp_list = self._get_dpp_list(inima, confterm, back)

        # make sure all DPP's exist
        self._check_files(self.dpp_list)

    def _get_dpp_list(self, inima, confterm, back):
        """Determine the name of all DPP files"""
        DPP_list = []

        # generate the input list for aXe
        axe_inputs = axeinputs.aXeInput(inima, confterm)

        # go over the list of all inputs
        for an_input in axe_inputs:
            # load the configuration file
            conf = configfile.ConfigFile(config_util.getCONF(an_input['config']))

            # get the image extensions
            ext_info = config_util.get_ext_info(config_util.getDATA(an_input['grisim']), conf)

            # get the name of all axe files
            axe_names = config_util.get_axe_names(an_input['grisim'], ext_info)

            # if requested,
            # append the background DPP file
            if back:
                DPP_list.append(axe_names['BCK_DPP'])
            else:
                # append the 'normal' DPP name to the list
                DPP_list.append(axe_names['DPP'])

        # return the DPP list
        return DPP_list

    def _check_files(self, dpp_list):
        """check for the existence of all DPP's"""

        for one_dpp in dpp_list:
            # make the full path name
            full_path = config_util.getOUTPUT(one_dpp)

            # check for existence
            if not os.path.isfile(full_path):
                # complain and out
                err_msg = 'Can not find the DPP file: %s!' % full_path
                raise aXeError(err_msg)

    def _get_contam_model(self):
        """Get the contamination model from in DPP"""

        # make a default return
        contam_model = None

        # open the fits and get the header
        fits_img = fits.open(config_util.getOUTPUT(self.dpp_list[0]), 'readonly')
        fits_head = fits_img[0].header

        # transfer value, if possible
        if 'CONTAM' in fits_head:
            contam_model = fits_head['CONTAM']

        # close the fits
        fits_img.close()

        # return the contamination model
        return contam_model

    def is_quant_contam(self):
        """Get the flag for quantitative contamination"""
        # set the default value
        isquantcont = True

        # do nothing if there are not contributors
        if len(self.dpp_list) < 1:
            return isquantcont

        # get the contamination model for the first contributor
        contam_model = self._get_contam_model()

        # check whether the model is quantitative or not
        isquantcont = config_util.is_quant_contam(contam_model)

        # return the flag
        return (contam_model, isquantcont)

    def filet_dpp(self, opt_extr=False):
        """Dump all DPP files"""
 
        drztmp = config_util.getDRZTMP()

        # go over all DPP files
        for one_dpp in self.dpp_list:

            # get the root-dir name
            root_dir = one_dpp.split('.DPP.fits')[0]
            root_dir_path = os.path.join(drztmp, root_dir)
            os.mkdir(root_dir_path)

            # create a filet objects and run the program
            filet = axelowlev.aXe_FILET(one_dpp, opt_extr=opt_extr, drztmp=root_dir_path)
            filet.run()
            del filet
