import os
import sys
import math
from astropy.io import fits
from stwcs.wcsutil import HSTWCS
from hstaxe.axeerror import aXeSIMError
from hstaxe.axesrc.axeiol import InputObjectList

class MagColList(InputObjectList):
    """Subclass for lists with magnitude columns

    This class loads ASCII tables and identifies the magnitude colunmns,
    columns which have AB magnitudes as their content.
    The magnitude columns have column names "MAG_?<number>*", with
    "<number>" the associated wavelength of the AB magnitude.

    Moreover there are mechanisms to search for column names and store
    the column index in a dictionary.

    TODO: FINISH UPDATING THIS CLASS FOR AXEIOL BASE
    """
    def __init__(self, filename, mag_wavelength=None):
        """Initializer for the class

        Parameters
        ----------
        filename: str
            name of the magnitude column list
        mag_wavelength: float
            special wavelength
        """
        super(InputObjectList, self).__init__(filename=filename)

        # initialize the dictionary
        # with indices of required columns
        self.reqColIndex = {}

        # check whether rows do exist
        if self.nrows > 0:
            # check for the mandatory columns
            # in the table
            self.find_magnitude_columns(mag_wavelength)

        # do some fundamental schecks
        self._basic_checks()

    def _basic_checks(self):
        """Do some fundamental checks"""
        # this check should be removed
        # at some time.....
        # check that there is at leas one row
        if self.nrows < 1:
            err_msg = ("\nThe table: {0:s} must contain at least one row!"
                       .format(self.filename))
            raise aXeSIMError(err_msg)

        # this check should be removed
        # at some time.....
        # check that there is at leas one row
        if (len(self.header) > 0):
            err_msg = ("\nThe table: {0:s} must not contain header comments!"
                       .format(self.filename))
            raise aXeSIMError(err_msg)

    def _find_required_columns(self, columnList):
        """Search and store the index of a list of required columns

        The method searches for a list of column names given in the input.
        In case that the column does not ecist, an exception is thrown.
        The index of existing columns is written into a dictionary with the
        column name as key.

        @param columnList: list with required columns
        @type columnList: []
        """
        # go over all columns
        for aColumn in columnList:
            # find it
            colIndex = self.find(aColumn)

            # complain if not found
            if colIndex < 0:
                err_msg = ("\nThe table: {0:s} does not contain column: {1}"
                           .format(self.filename, aColumn))
                raise aXeSIMError(err_msg)
            else:
                # add the index to the dictionary
                self.reqColIndex[aColumn] = colIndex

    def _identify_magcol(self, mag_cols, mag_wave):
        """Identify the meagnitude column closest to a characteristic value

        The method analyses all magnitude columns and finds the one
        which is closest to a wavelength given in the input.
        The index of the closest wavelength in the input
        list is returned.

        Parameters
        ----------
        mag_cols: list
            list of [index, wavelength] with information on all magnitude cols
        mag_wave: float
            characteristic wavelength

        Returns
        -------
        min_ind: int
            the index of the magnitude column closest to mag_wave
        """
        # define a incredible large difference
        min_dist = 1.0e+30

        # define a non-result
        min_ind = -1

        # go over al magnitude columns
        for index in range(len(mag_cols)):
            # check wehether a new minimum distance is achieved
            if math.fabs(mag_cols[index][1]-mag_wave) < min_dist:
                # transport the minimum and the index
                min_ind = index
                min_dist = math.fabs(mag_cols[index][1]-mag_wave)

        # return the index
        return min_ind

    def _search_mcols(self):
        """Search the magnitude columns

        The method collects all magnitude columns
        with an encoded wavelength in the column name.
        For each such column the column index and the
        wavelength is stored in a list, and the list
        of all columns is returned.

        Returns
        -------
        mag_cols: list
            a list if lists of column number-wavelength pairs
        """
        # initialize the list with the result
        mag_cols = []

        # go over all columns
        for index in range(self.ncols):

            # get the column name
            colname = self[index].colname

            # try to decode the wavelength
            wave = self._get_wavelength(colname)

            # if a wavelength is encoded
            if wave:
                # compose and append the info
                # to the resulting list
                mag_cols.append([index, wave])

        # return the result
        return mag_cols

    def _get_wavelength(self, colname):
        """Determines the wvelength from a column name

        The method tries to extract the wavelength
        encoded into a column name. The encoding
        format is "MAG_<C><WAVE>*" with <C> a
        single character, <WAVE> an integer number
        and anything (*) afterwards.
        in case that there is not wavelength encoded,
        the value 0 id given back.

        Parameters
        ----------
        colname: str
            the column name

        Returns
        -------
        wave: float
            the wavelength encoded in the column name
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
                # as soon as the substirng can NOT
                # be transferred to an int
                # return the current, best wavelength
                except ValueError:
                    return wave

        # return the best result
        return wave

    def find_magnitude_columns(self, mag_wavelength=None):
        """Identify all magnitude columns

        The method identifiy all magnitude columns and can select the
        one with the wavelength closest to a given input wavelength.
        An exception is thrown in case that no magnitude column is found.

        @param mag_wavelength: characteristic wavelength
        @type mag_wavelength: float
        """
        # search for
        # magnitude columns in general
        self.mag_cols = self._search_mcols()

        # if magnitude columns exist,
        # search the one closest to the
        # desired wavelength
        if len(self.mag_cols) > 0:
            if (mag_wavelength is not None):
                mag_index = self._identify_magcol(self.mag_cols,
                                                  mag_wavelength)
            else:
                mag_index = 0

            # set the column number and the wavelength
            self.magwave = float(self.mag_cols[mag_index][1])
            self.reqColIndex['MAGNITUDE'] = self.mag_cols[mag_index][0]
        else:
            # enhance the error counter
            err_msg = ("\nModel spectrum list: {0:s} does not contain any "
                       "magnitude column!".format(self.filename))
            raise aXeSIMError(err_msg)


class ModelSpectrumList(MagColList):
    """Subclass for model spectra list"""
    def __init__(self, filename, mag_wavelength=None):
        """Initializes the class

        Parameters
        ----------
        filename: str
            name of the model spectra list
        mag_wavelength: float
            characteristic wavelength
        """
        super(ModelSpectrumList, self).__init__(filename=filename,
                                                mag_wavelength=mag_wavelength)

        # find the columsn reqested for a
        # model spectrum list
        self._find_modspeclist_columns()

        # re-set all values in the column "MODSPEC"
        self._reset_modspec()

    def _find_modspeclist_columns(self):
        """Find the requested columns

        The method finds the columns which are mandatory for a model
        spectra list.
        """
        # the required columns
        reqColumnsNames = ['NUMBER', 'Z', 'SPECTEMP']

        # search for all required columns
        self._find_required_columns(reqColumnsNames)

    def _reset_modspec(self):
        """Reset the values in the column "MODSPEC"

        The method sets all values in the column "MODSPEC" to
        the value 0.
        """
        # go over all rows
        for index in range(self.nrows):
            # set the entries to 0
            self["MODSPEC"][index] = 0

    def find_tempmax(self):
        """Find the largest entry in the spectral template column

        The method identifies and returns the largest value found in the
        column "SPECTEMP".

        Returns
        -------
        stemp_max: int
            The largest entry in "SPECTEMP"
        """
        # initialize the value
        stemp_max = 0

        # go over all rows
        for index in range(self.nrows):
            # transfer a new maximum number, if necessary
            stemp_max = max(self[self.reqColIndex['SPECTEMP']][index],
                            stemp_max)

        # return the maximum number
        return stemp_max


class ModelObjectTable(MagColList):
    """Subclass for model object table"""
    def __init__(self, filename, model_spectra=None, model_images=None):
        """
        Initializes the class

        Parameters
        ----------
        filename: str
            name of the model object table
        model_spectra: str
            the model spectra file
        model_images: str
            the model image file
        """
        super(ModelObjectTable, self).__init__(filename=filename)

        # find the columsn reqested for a
        # model spectrum list
        # the required columns
        reqColumnsNames = ['NUMBER', 'X_IMAGE', 'Y_IMAGE',
                           'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE']

        if model_spectra is not None:
            reqColumnsNames.extend(['MODSPEC'])
        if model_images is not None:
            reqColumnsNames.extend(['MODIMAGE'])

        # search for all required columns
        self._find_required_columns(reqColumnsNames)

    def fill_columns(self, WCSimage, WCSext=None):
        """Fill up column information to be ready for aXe

        The method completes the model object tables with columns
        requested by aXe.

        Parameters
        ----------
        WCSimage: str
            the reference image with WCS
        WCSext: str
            the extension to use


        ra_dec = iraf.xy2rd(infile=WCSimage, x=self['X_IMAGE'][index],\
                            y=, hms='NO', Stdout=1)
        """

        # check that the image exists
        if not os.path.isfile(WCSimage):
            err_msg = "The WCS image: {0:s} does not exist!".format(WCSimage)
            raise aXeSIMError(err_msg)

        # append the image extesnion
        if WCSext is not None:
            WCSimage = fits.getheader(WCSimage, ext=WCSext)
        else:
            WCSimage = fits.getheader(WCSimage)

        # go over all rows
        for index in range(self.nrows):

            # just copy some information
            # later it would be reasonable
            # to give more reasonable values
            self['A_WORLD'][index] = self['A_IMAGE'][index]
            self['B_WORLD'][index] = self['B_IMAGE'][index]
            self['THETA_WORLD'][index] = self['THETA_IMAGE'][index]

            # transform x,y to ra and dec
            ra, dec = WCSimage.all_world2pix(self['Y_IMAGE'][index],
                                             self['X_IMAGE'][index])

            # store ra and dec
            self['X_WORLD'][index] = float(ra)
            self['Y_WORLD'][index] = float(dec)

        # save the changes
        self.flush()
