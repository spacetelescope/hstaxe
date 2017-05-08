from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import os
from stsci.tools import wcsutil

from .. import axe_asciidata
from ..axeerror import aXeError
from .. import axeutils
from . import configfile
from . import axeiol


class Sex2GolPy(object):
    """ Class for the SEX2GOL task """
    def __init__(self, grisim, config,
                 in_sex=None,
                 dirname=None,
                 out_sex=None,
                 spec_hdu=None,
                 dir_hdu=None):

        # store some parameters
        self.grisim = grisim
        self.config = config
        self.in_sex = in_sex
        self.dirname = dirname

        # determine the grism image extensions
        self.grism_extinfo = self._get_grism_ext_info(grisim, config, spec_hdu)

        # determine the grism image extensions
        self.dirname, self.dirname_extinfo = self._get_dirname_information(dirname, config, grisim, self.grism_extinfo, dir_hdu)

        # get information on the input and output lists
        self.in_sex, self.out_sex = self._resolve_list_names(self.dirname,
                                                             self.dirname_extinfo,
                                                             self.grisim,
                                                             self.grism_extinfo,
                                                             in_sex, out_sex)

        # save a name for stdout
        self.stdout = axeutils.getOUTPUT("pysex2gol.stdout")

    def __str__(self):
        """String method for the class"""
        # define the prefix
        prefix = "py_SEX2GOL: "

        # compose the feedback
        big_str = "{0:s}  Setup:\n".format(prefix)
        big_str += "{0:s}  Input g/prism image:      {0:s} \n".format(prefix, self.grisim)
        big_str += "{0:s}  Configuration file name:  {0:s} \n".format(prefix, self.config)
        big_str += "{0:s}  Direct image:             {0:s} \n".format(prefix, self.dirname)
        big_str += "{0:s}  G/Prism extension:        {0:s} \n".format(prefix, self.grism_extinfo['axe_ext'])
        big_str += "{0:s}  Direct image extension:   {0:s} \n".format(prefix, self.dirname_extinfo['axe_ext'])
        big_str += "{0:s}  Input catalog name:       {0:s} \n".format(prefix, self.in_sex)
        big_str += "{0:s}  Output catalog name:      {0:s}   ".format(prefix, self.out_sex)

        # return the string
        return big_str

    def _cleanup(self):
        """Clean up files created for stdout

        This is a usual cleaning procedure in case nothing bad happened.
        """
        # delete stdout/stderr
        if os.path.isfile(self.stdout):
            os.unlink(self.stdout)

    def _get_grism_ext_info(self, grisim, config, spec_hdu=None):
        """Determine the extension information on the grism image"""
        # check for an explicit extension
        if spec_hdu is None:
            # load the configuration file;
            # determine the extension information
            conf = configfile.ConfigFile(axeutils.getCONF(config))
            ext_info = axeutils.get_ext_info(axeutils.getIMAGE(grisim), conf)
            del conf

        else:
            # make by hand the extension information
            ext_info = {'axe_ext': spec_hdu, 'fits_ext': spec_hdu-1}

        # return the extension info
        return ext_info

    def _get_dirname_information(self,
                                 dirname,
                                 config,
                                 grisim,
                                 grism_extinfo,
                                 dir_hdu=None):
        """Determine the direct image information"""
        # check whether ANY direct image information exists
        if ((dirname is None) and (dir_hdu is None)):
            # set the grism image as direct image
            dirname = grisim
            dirname_extinfo = grism_extinfo

        elif ((dirname is not None) and (dir_hdu is None)):
            # load the configuration file;
            # determine the extension information
            conf = configfile.ConfigFile(axeutils.getCONF(config))
            dirname_extinfo = axeutils.get_ext_info(axeutils.getIMAGE(grisim),
                                                    conf)
            del conf

        elif ((dirname is not None) and (dir_hdu is not None)):
            # make by hand the extension information
            dirname_extinfo = {'axe_ext': dir_hdu, 'fits_ext': dir_hdu-1}

        else:
            # error and out
            err_msg = ("Specifying NO direct image but a direct image HDU: "
                       "{0:d} makrs NO sense!".format(dir_hdu))
            raise aXeError(err_msg)

        # return the name and the extension info
        return dirname, dirname_extinfo

    def _resolve_list_names(self,
                            dirname,
                            dirname_extinfo,
                            grisim,
                            grism_extinfo,
                            in_sex,
                            out_sex):
        """Determine the lists for input and output"""
        # compose the default name
        if in_sex is None:
            # compose the filename from the direct image name
            in_sex = dirname.replace(".fits", "_{0:d}.cat".format(dirname_extinfo['axe_ext']))

        # check whether the explicitly given filename exists
        if not os.path.isfile(in_sex):
            err_msg = "The Input Object List: {0:s}  does not exist!".format(in_sex)
            raise aXeError(err_msg)

        if out_sex is None:
            # compose the name for the output GOL
            out_sex = os.path.basename(grisim).replace(".fits", "_{0:d}.cat".format(grism_extinfo['axe_ext']))

        # return the IOL and the GOL names
        return in_sex, out_sex

    def _transfer_fromIOL_toGOL(self):
        """Copies relevant data from IOL to GOL"""

        # load the IOL
        iol = axeiol.InputObjectList(self.in_sex)

        # check for an empty table
        if iol.nrows < 1:
            # return only none's
            return None, None

        # determine the number of columns
        ncols_gol = len(iol.mand_cols) + len(iol.wav_cols)

        # create an empty GOL
        gol = axe_asciidata.create(ncols_gol, iol.nrows)
        gol.toSExtractor()

        # go over all
        # mandatory columns
        c_index = 0
        for one_col in iol.mand_cols:
            # rename the GOL column
            gol[c_index].rename(one_col['name'])

            # go over all rows
            for r_index in range(iol.nrows):

                # transfer the value
                gol[one_col['name']][r_index] = iol[one_col['name']][r_index]

            # enhance the index
            c_index += 1

        # go over all wavelength columns
        for one_col in iol.wav_cols:

            # rename the GOL column
            gol[c_index].rename(one_col['name'])

            # go over all rows
            for r_index in range(iol.nrows):

                # transfer the value
                gol[one_col['name']][r_index] = iol[one_col['name']][r_index]

            # enhance the index
            c_index += 1

        # return the GOL
        return iol, gol

    def _transfer_coos(self, iol, gol):
        """Transfer coordinates from the IOL to the GOL"""
        # compose the WCS-term for the direct and grism images
        dir_term = axeutils.getIMAGE("{0:s} [{1:d}]".format(self.dirname, self.dirname_extinfo['fits_ext']))
        gri_term = axeutils.getIMAGE("{0:s} [{1:d}]".format(self.grisim, self.grism_extinfo['fits_ext']))

        # generate the WCS objects
        dir_wcs = wcsutil.WCSObject(dir_term)
        gri_wcs = wcsutil.WCSObject(gri_term)

        # go over each row
        for index in range(iol.nrows):

            # make a position tuple
            xy_dirname = (iol['X_IMAGE'][index], iol['Y_IMAGE'][index])

            # convert to RADEC
            radec_pos = dir_wcs.xy2rd(xy_dirname)

            # convert to XY on grism
            xy_grism = gri_wcs.rd2xy(radec_pos)

            # store projected vals in the GOL
            gol['X_IMAGE'][index] = float(xy_grism[0])
            gol['Y_IMAGE'][index] = float(xy_grism[1])

    def _treat_NULL_table(self, out_name):
        """Transfer an empty table

        The header of the empty table is written to the
        GOL file. However then no check is done whether
        the column is complete.
        """
        # open the "GOL"
        out_file = open(out_name, 'w+')

        # go over the "IOL"
        for a_line in self.in_sex:
            # transfer everyhing
            out_file.write(a_line)

        # close the file
        out_file.close()

    def run(self, silent=False):
        """Make the SEX2GOL transformations"""
        # give feedback
        if not silent:
            print(self)
            print("py_SEX2GOL:  Start processing ...",)
        else:
            # open stdout/stderr
            sout = open(self.stdout, 'w+')
            sout.write(str(self)+"\n")
            sout.write("py_SEX2GOL:  Start processing ...")

        # copy the relevant data to the GOL
        iol, gol = self._transfer_fromIOL_toGOL()

        # check whether something can be done
        if ((iol is not None) and (gol is not None)):
            # transfer the coordinates
            self._transfer_coos(iol, gol)

            # store the GOL
            gol.writeto(axeutils.getOUTPUT(self.out_sex))

        else:
            # if there are no objects, just copy the empty table
            # header to the GOL
            self._treat_NULL_table(axeutils.getOUTPUT(self.out_sex))

            # give feedback
            if not silent:
                print(self)
                print("py_SEX2GOL:  Warning! Empty table copied to GOL")
            else:
                # open stdout/stderr
                sout.write("py_SEX2GOL:  Warning! Empty table copied to GOL")

        # give feedback
        if not silent:
            print("     Done")
        else:
            sout.write("     Done\n")
            sout.close()

    def runall(self, silent=False):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        @param silent: boolean for silent mode
        @type silent: boolean

        @return: the return code of the C-executable
        @rtype: int
        """
        # run the executable
        retcode = self.run(silent=silent)

        # check whether the run was good
        if silent:
            # do the cleaning
            self._cleanup()

        return retcode
