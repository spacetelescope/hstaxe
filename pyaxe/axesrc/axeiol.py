import re
from astropy.table import Table
from ..axeerror import aXeError


class InputObjectList(object):
    """The input object list should be column selected SExtractor cat"""
    def __init__(self, filename):
        self.filename = filename

        # Read and validate the catalog
        # this expects an ascii file with atleast one header line
        # a full sextrctor catalog is acceptable
        self.catalog = Table.read(filename, format='ascii.sextractor')

        # check for the mandatory columns
        # self.mand_cols, self.wav_cols = self._find_columns()
        self._validate_columns()

    def _validate_columns(self):
        """Validate all required columns"""
        # the list of mandatory columns
        mand_colnames = ["NUMBER", "X_IMAGE", "Y_IMAGE", "A_IMAGE",
                         "B_IMAGE", "THETA_IMAGE", "X_WORLD", "Y_WORLD",
                         "A_WORLD", "B_WORLD", "THETA_WORLD"]

        # the list of optional columns
        # opt_colnames = ["MODSPEC", "MODIMAGE"]

        # go over all mandatory columns
        cat_columns = self.catalog.colnames
        for colname in mand_colnames:
            if colname not in cat_columns:
                err_msg = ("Input Object List: {0:s} does not contain column "
                           "{1:s}".format(self.filename, colname))
                raise aXeError(err_msg)

        # check for the MAG_AUTO-column
        mauto_col = "MAG_AUTO" in cat_columns

        # get the columns with an encoded wavelength
        wav_cols = self.search_mcols()

        # check whether there is no mag-column
        if ((not mauto_col) and (len(wav_cols) == 0)):
            # complain and out
            err_msg = ("Catalogue: {0:s} does not contain any magnitude "
                       "column!".format(self.filename))
            raise aXeError(err_msg)

        # check whether there is only MAG_AUTO
        elif ((mauto_col) and (len(wav_cols) == 0)):
            wav_cols = [{'name': 'MAG_AUTO', 'lambda': None}]

        # check whether there is MAG_AUTO and magval columns
        elif ((mauto_col) and (len(wav_cols) > 1)):
            err_msg = ("Catalogue: {0:s} contains 'MAG_AUTO' and {1:i} other "
                       "magnitude columns!"
                       .format(self.filename, len(wav_cols)))
            raise aXeError(err_msg)

    # def find_magcol(self, mag_cols, mag_wave):
    #     """
    #     Input:
    #         mag_cols - the list with infos on magnitude columns
    #         mag_wave - the target wavelength
    #
    #     Description:
    #         The method analyses all magnitude columns and finds the one
    #         which is closest to a wavelength given in the input.
    #     """
    #     # define a incredibly large difference
    #     min_dist = 1.0e+30
    #
    #     # define a non-result
    #     min_ind  = -1
    #
    #     # go over al magnitude columns
    #     for index in range(len(mag_cols)):
    #         # check wehether a new minimum distance is achieved
    #         if math.fabs(mag_cols[index][1]-mag_wave) < min_dist:
    #             # transport the minimum and the index
    #             min_ind  = index
    #             min_dist = math.fabs(mag_cols[index][1]-mag_wave)
    #
    #     # return the index
    #     return min_ind

    def search_mcols(self):
        """

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
        for colname in self.catalog.colnames:
            if (("MAG" in colname) and ("AUTO" not in colname)):
                # try to decode the wavelength
                wave = self.get_wavelength(colname)

                # if a wavelength is encoded
                if wave:
                    # compose and append the info
                    # to the resulting list
                    mag_cols.append({'name': colname, 'lambda': wave})

        # return the result
        return mag_cols

    def get_wavelength(self, colname):
        """Return the wavelength or None.

        The method tries to extract the wavelength
        encoded into a column name. The encoding
        format is "MAG_<C><WAVE>*" with <C> a
        single character, <WAVE> an integer number
        and anything (*) afterwards.

        Input
        -----
            colname - the column name

        Returns
        -------
            wave - the wavelength encoded in the
                    column name, or None

        """
        # check for the start string
        check = re.compile("^MAG_([A-Z,a-z]{1})([0-9]*)")
        found = check.split(colname)
        if len(found) > 2:
            return int(found[2])
        else:
            return None
