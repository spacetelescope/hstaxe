import re
from astropy.table import Table

from hstaxe.axeerror import aXeError


class InputObjectList(object):
    """The input object list should be column selected SExtractor cat"""
    def __init__(self, filename):
        self.filename = filename

        # Read and validate the sextrator catalog
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
        for colname in mand_colnames:
            if colname not in self.catalog.colnames:
                err_msg = ("Input Object List: {0:s} does not contain column "
                           "{1:s}".format(self.filename, colname))
                raise aXeError(err_msg)

        # check for the MAG_AUTO-column
        mauto_col = "MAG_AUTO" in self.catalog.colnames

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
        single character, <WAVE> an integer number.

        Input
        -----
        colname: str
            the column name

        Returns
        -------
        wave: int
            the wavelength encoded in the
            column name, or None

        """
        # check for the start string
        check = re.compile("^MAG_([A-Z,a-z]{1})([0-9]*)")
        found = check.search(colname)
        if found is not None:
            return int(found[2])
        else:
            return None
