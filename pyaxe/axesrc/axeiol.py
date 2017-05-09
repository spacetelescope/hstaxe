from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .. import axe_asciidata
from ..axeerror import aXeError


class InputObjectList(axe_asciidata.AsciiData):
    """Subclass of the AsciiData class for aXe"""
    def __init__(self, filename):
        # intialize via the superclass
        super(InputObjectList, self).__init__(filename=filename)

        if self.nrows > 0:
            # check for the mandatory columns
            # in the table and save the indices
            self.mand_cols, self.wav_cols = self._find_columns()

    def _find_columns(self):
        """Identify all important columns"""
        # initialize the
        # dict for mandatory columns
        mand_cols = []

        # the list of mandatory columns
        mand_colnames = ["NUMBER", "X_IMAGE", "Y_IMAGE", "A_IMAGE",
                         "B_IMAGE", "THETA_IMAGE", "X_WORLD", "Y_WORLD",
                         "A_WORLD", "B_WORLD", "THETA_WORLD"]

        # the list of optional columns
        opt_colnames = ["MODSPEC", "MODIMAGE"]

        # go over all mandatory columns
        for one_column in mand_colnames:

            # find the column
            col_num = self.find(one_column)

            # check for existence
            if col_num < 0:
                # complain and out
                err_msg = 'Input Object List: %s does not contain column "%s" !' % (self.filename, one_column)
                raise aXeError(err_msg)
            else:
                # store the column index
                mand_cols.append({'name': one_column, 'index':col_num})

        # go over all optional columns
        for one_column in opt_colnames:

            # find the column
            col_num = self.find(one_column)

            # check for existence
            if col_num > -1:
                mand_cols.append({'name': one_column, 'index':col_num})

        # check for the MAG_AUTO-column
        mauto_col = self.find("MAG_AUTO")

        # get the columns with the
        # encoded wavelength
        wav_cols = self.search_mcols()

        # check whether there is no mag-column
        if mauto_col < 0 and len(wav_cols) < 1:
            # complain and out
            err_msg = 'Catalogue: %s does not contain any magnitude column!' % self.filename
            raise aXeError(err_msg)

        # check whether there is only MAG_AUTO
        elif mauto_col > -1 and len(wav_cols) < 1:

            # Store the column index
            wav_cols = [{'name': 'MAG_AUTO', 'index':mauto_col, 'lambda': None}]

        # check whether there is MAG_AUTO and magval columns
        elif mauto_col > -1 and len(wav_cols) > 1:
            err_msg = 'Catalogue: %s contains "MAG_AUTO" and %i other magnitude columns!' % (self.filename, len(wav_cols))
            raise aXeError(err_msg)

        # return both column lists
        return mand_cols, wav_cols

    def find_magcol(self, mag_cols, mag_wave):
        """
        Input:
            mag_cols - the list with infos on magnitude columns
            mag_wave - the target wavelength

        Return:
            min_ind - the index of the closes column within mag-cols

        Description:
            The method analyses all magnitude columns and finds the one
            which is closest to a wavelength given in the input.
            The index of the closest wavelength in the input
            list is returned.
        """
        import math

        # define a incredible large difference
        min_dist = 1.0e+30

        # define a non-result
        min_ind  = -1

        # go over al magnitude columns
        for index in range(len(mag_cols)):
            # check wehether a new minimum distance is achieved
            if math.fabs(mag_cols[index][1]-mag_wave) < min_dist:
                # transport the minimum and the index
                min_ind  = index
                min_dist = math.fabs(mag_cols[index][1]-mag_wave)

        # return the index
        return min_ind

    def search_mcols(self):
        """
        Input:
            -

        Return:
            mag_cols - a list with tuples

        Description:
            The method collects all magnitude columns
            with an encoded wavelength in the column name.
            For each such column the column index and the
            wavelength is stored in a list, and the list
            of all columns is returned.
        """

        # initialize the list with the result
        mag_cols = []

        # go over all columns
        for index in range(self.ncols):

            # get the column name
            colname = self[index].colname

            # try to decode the wavelength
            wave = self.get_wavelength(colname)

            # if a wavelength is encoded
            if wave:
                # compose and append the info
                # to the resulting list
                mag_cols.append({'name': colname, 'index':index, 'lambda': wave})

        # return the result
        return mag_cols


    def get_wavelength(self, colname):
        """
        Input:
            colname - the column name

        Return:
            wave - the wavelength encoded in the
                    column name, or 0

        Description:
            The method tries to extract the wavelength
            encoded into a column name. The encoding
            format is "MAG_<C><WAVE>*" with <C> a
            single character, <WAVE> an integer number
            and anything (*) afterwards.
            in case that there is not wavelength encoded,
            the value 0 id given back.
        """
        # set the value for 'nothing found'
        wave = 0

        # check for the start string
        if colname.find('MAG_') == 0:

            # copy the rest to a substring
            rest_name = colname.split('MAG_')[1][1:]

            # prepare to analyse the whole substring
            for index in range(len(rest_name)):

                # make a progressively longer
                # substring, starting from the beginning
                cand_name = rest_name[0:index+1]

                # try to convert the substirng to
                # and integer, set a new wavelength
                # if it is possible
                try:
                    wave = int(cand_name)
                # as soon as the substring can NOT
                # be transferred to an int
                # return the current, best wavelength
                except ValueError:
                    return wave

        # return the best result
        return wave
