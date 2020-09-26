import copy
import math
import logging

from astropy.io import fits
from .. import axe_asciidata
from ..axeerror import aXeError

from hstaxe.utils import set_logging

# make sure there is a logger
_log = logging.getLogger(__name__)
_log = set_logging()  # defaults to INFO

class Interpolator:
    """General class to get interpolated values"""
    def __init__(self, input_file=None, indep=None, depen=None):
        """The method initializes an interpolator objects from various input.

        Parameters
        ----------
        input_file: str
            data file name
        indep: list
            with independent data
        depen: list
            with dependent data
       """
        # check whether a file is given
        if input_file is not None:
            # store the file name
            self.input_file = input_file

            # check whether the name ends with '.fits'
            if input_file.rfind('.fits') == len(input_file)-len('.fits'):
                # load the data from the fits file
                self._indep_data, self._depen_data = self._load_interp_fromfits(input_file)
            else:
                # load the data from the file
                self._indep_data, self._depen_data = self._load_interp_fromfile(input_file)

        # check whether two arrays are given
        elif ((indep is not None) and (depen is not None)):
            # set the file name to 'None'
            self.input_file = None

            # load the data from two lists
            self._indep_data, self._depen_data = self._load_interp_fromlist(indep, depen)

        else:
            # report that there is no data
            error_message = '\nNo data given for interpolator!'
            raise aXeSIMError(error_message)

        # store minimum and maximum of
        # independent values
        if len(self) > 0:
            self.ind_min = self._indep_data[0]
            self.ind_max = self._indep_data[len(self)-1]
        else:
            self.ind_min = None
            self.ind_max = None

        # initialize an accelerator
        self.accelerator = 0

    def __getitem__(self, value):
        """The index operator for the class

        The method returns the interpolated value for the indpendent value
        given in the input. The method of interpolation is
        'linear interpolation'>

        Parameters
        ----------
        value: float
            independent value to derive dependent value for

        Returns
        -------
        float: the interpolated value at the position of the input
        """
        # check the input against the extremies
        if ((value < self.ind_min) or (value > self.ind_max)):
            error_message = ("\nThe requested value: {0:s} is outside the "
                             "covered range!".format(str(value)))
            raise aXeSIMError(error_message)

        # find the right location within the
        # independent values
        index = self._get_indep_index(value)

        # get the fraction covered on the independend interval
        factor = ((value - self._indep_data[index-1]) /
                  (self._indep_data[index] - self._indep_data[index-1]))

        # calculate and return the ionterpolated value
        return (self._depen_data[index-1] + factor *
                (self._depen_data[index] - self._depen_data[index-1]))

    def __len__(self):
        """The length operator for the class

        The method determines and returns the length of an interpolator
        object. The length is defined to be the number of interpolator
        data points.

        Returns
        -------
        int: the length of a class instance
        """
        return len(self._depen_data)

    def __str__(self):
        """The sting method for the class

        The method implements a string method for the class. The string just
        adds all independent dependent value pairs.

        Returns
        -------
        str: the string representation of the interpolator instance
       """
        # initialize the string
        bigstr = ""

        # go over all data
        for index in range(len(self)):
            # add to the string
            bigstr += "{0:e} {1:e}\n".format(self._indep_data[index],
                                             self._depen_data[index])
        # return the large string
        return bigstr

    def __deepcopy__(self, memo):
        """The deep copy method of the class instance

        The method defines a deep copy method to create and return an
        identical instance of an interpolation object.

        Parameters
        ----------
        memo: a required parameter for a deep copy method

        Returns: Interpolator
            returns the deep copy of the instance
        """
        # make deep copies of the ingredients
        indep = copy.deepcopy(self._indep_data)
        depen = copy.deepcopy(self._depen_data)

        # create and return a new object
        return Interpolator(indep=indep, depen=depen)

    def __mul__(self, in_object):
        """Defines multiplicator for the class

        The method defines the multiplication of two interpolator class
        instances. The multiplication of two interpolators is an interpolator
        object. The independent data of this interpolator unification of the
        independent data of the multiplicants. The dependent data is the
        product of interpolated values.

        Parameters
        ----------
        in_object: Interpolator
            the multiplicant

        Returns: the product of two Interpolator
        """

        # determine the min and max for the product
        ind_min = max(self.ind_min, in_object.ind_min)
        ind_max = min(self.ind_max, in_object.ind_max)

        # compose the independent data for the mult
        indep = copy.deepcopy(self._indep_data)
        indep.extend(in_object._indep_data)
        indep.sort()

        # initialize the backward index
        r_index = len(indep)-1

        # go over the new array;
        # clean the independent vector
        for index in range(1, len(indep)):
            # check the boundaries
            if indep[r_index] > ind_max:
                del indep[r_index]
            elif indep[r_index] < ind_min:
                del indep[r_index]
            # check for repetition
            elif indep[r_index] == indep[r_index-1]:
                del indep[r_index]

            # adjust the index
            r_index -= 1

        # check the last element independently
        if indep[0] > ind_max or indep[0] < ind_min:
            del indep[0]

        # build up the array with
        # the dependent values
        depen = []
        for item in indep:
            depen.append(self[item]*in_object[item])

        # return a new array
        return Interpolator(indep=indep, depen=depen)

    def __div__(self, in_object):
        """Defines division for the class

        The method defines the multiplication of two interpolator class
        instances. The multiplication of two interpolators is an interpolator
        object. The independent data of this interpolator unification of the
        independent data of the multiplicants. The dependent data is the
        product of interpolated values.

        Parameters
        ----------
        in_object: Interpolator
            the multiplicant

        Returns: the product of two Interpolators
        """
        import copy

        # determine the min and max for the product
        ind_min = max(self.ind_min, in_object.ind_min)
        ind_max = min(self.ind_max, in_object.ind_max)

        # compose the independent data for the mult
        indep = copy.deepcopy(self._indep_data)
        indep.extend(in_object._indep_data)
        indep.sort()

        # initialize the backward index
        r_index = len(indep)-1

        # go over the new array;
        # clean the independent vector
        for index in range(1, len(indep)):
            # check the boundaries
            if indep[r_index] > ind_max:
                del indep[r_index]
            elif indep[r_index] < ind_min:
                del indep[r_index]
            # check for repetition
            elif indep[r_index] == indep[r_index-1]:
                del indep[r_index]

            # adjust the index
            r_index -= 1

        # check the last element independently
        if indep[0] > ind_max or indep[0] < ind_min:
            del indep[0]

        # build up the array with
        # the dependent values
        depen = []
        for item in indep:

            # check for 0.0
            if in_object[item] == 0.0:
                err_msg = 'Value 0.0 in divisor!'
                raise aXeSIMError(err_msg)

            # make the division
            depen.append(self[item] / in_object[item])

        # return a new array
        return Interpolator(indep=indep, depen=depen)

    def _load_interp_fromfile(self, input_file):
        """Load the interpolator from an ASCII file

        The method loads independent and dependent data from an ascii file
        and creates an interpolator object.

        Parameters
        ----------
        input_file: str
            the data file name

        Returns: the interpolator with the data
        """
        # load the data file
        indata = axe_asciidata.open(input_file)

        # check the input
        self._check_input(indata, input_file)

        # return the column data
        return indata[0]._data, indata[1]._data

    def _load_interp_fromfits(self, fits_file):
        """Load the interpolator from an ASCII file

        The method loads independent and dependent data from a fits file
        and creates an interpolator from this data.

        Parameters
        ----------
        fits_file: str
            the fits file name with the data

        Returns: the interpolator with the data
        """
        # open the fits file
        infits = fits.open(fits_file, mode='readonly')

        # go to the table data
        tabdata = infits[1].data

        # create a new asciidata structure
        indata = axe_asciidata.create(2, len(tabdata.field(0)))

        # go over the new table
        for index in range(indata.nrows):
            # transfer data from the fits table to the ascii table
            try:
                indata[0][index] = float(tabdata.field(0)[index])
            except ValueError:
                msg = ("\nColumn 0 of file: {0:s} contains wrong type!"
                       .format(fits_file))
                raise aXeSIMError(msg)
            try:
                indata[1][index] = float(tabdata.field(1)[index])
            except ValueError:
                msg = ("\nColumn 1of file: {0:s} contains wrong type!"
                       .format(fits_file))
                raise aXeSIMError(msg)

        # check the input
        self._check_input(indata, fits_file)

        # return the column data
        return indata[0]._data, indata[1]._data

    def _load_interp_fromlist(self, indep_data, depen_data):
        """Creates an interpolator from data lists

        The method creates an interpolator from lists with the
        dependent and independent data values given in the input.

        @param indep_data: the list with the independent data
        @type indep_data: string
        @param depen_data: the list with dependent data
        @type depen_data: string

        @return: the interpolator with the data
        @rtype: Interpolator
        """
        # initialize new vectors
        out_indep = []
        out_depen = []

        # check the length of the lists
        if len(indep_data) != len(depen_data):
            error_message = '\nData length is inhomogeneous!'
            raise aXeSIMError(error_message)

        # go over all rows
        for index in range(len(indep_data)):

            # check for None-entries
            if ((indep_data[index] is None) or (depen_data[index] is None)):
                raise aXeSIMError("Data contains NULL entries!")

            # test conversion to float
            try:
                tmp = float(indep_data[index])
            except ValueError:
                msg = ("\nValue: {0:s} is not convertible to float!"
                       .format(str(indep_data[index])))
                raise aXeSIMError(msg)

            # test conversion to float
            try:
                tmp = float(depen_data[index])
            except ValueError:
                msg = ("\nValue: {0:s} is not convertible to float!"
                       .format(str(depen_data[index])))
                raise aXeSIMError(msg)

            # check for ascending order
            if index > 0:
                # check that independent data is rising
                if indep_data[index] <= indep_data[index-1]:
                    raise aXeSIMError("\nIndependent column data is not"
                                      "monotonically rising!")

            # append the values to the arrays
            out_indep.append(float(indep_data[index]))
            out_depen.append(float(depen_data[index]))

        # return the new arrays
        return out_indep, out_depen

    def _get_indep_index(self, value):
        """Locate the correct index

        The method locates the 'index' for the independent data which
        has "indep_data[index] <= value indep_data[index+1]".
        The approximate index postition "self.accelerator" is used as
        a starting point.

        Parameters
        ----------
        value: str
            the list with dependent data

        Returns: index position for the independent data
        """
        # check whether you have to search upwards or downwards
        if value > self._indep_data[self.accelerator]:
            # in case that you search upwards, go up
            # the independent values until you find the right interval
            self.accelerator += 1
            while value > self._indep_data[self.accelerator]:
                self.accelerator += 1
        else:
            # in case that you search downwards, go down
            # the independent values  until you find the right interval
            #    while(wavelength < resp->spec[nact-1].lambda_mean)
            while value < self._indep_data[self.accelerator-1]:
                self.accelerator -= 1

        # return the position
        return self.accelerator

    def _check_input(self, indata, input_file):
        """Check whether the input valid

        The method performs some basic checks on data which is supposed
        to form an interpolator. Checks for the correct data type and
        for rising independent data values and against NULL entries are done.

        Parameters
        ----------
        indata: axe_asciidata.AsciiData()
            the list with dependent data
        input_file: str
            the input file name
        """
        # check the type of first column
        if isinstance(float, indata[0].get_type()):

            # check whether it is int
            if isinstance(int, indata[0].get_type()):
                # go over the column
                for index in range(indata.nrows):
                    # convert to float
                    indata[0][index] = float(indata[0][index])
            else:
                msg = ("\nColumn 0 of file: {0:s}} contains wrong "
                       "type!".format(input_file))
                raise aXeSIMError(msg)

        # check the type of second column
        if isinstance(float, indata[1].get_type()):

            # check whether it is int
            if (isinstance(int, indata[1].get_type())):
                # go over the column
                for index in range(indata.nrows):
                    # convert to float
                    indata[1][index] = float(indata[1][index])
            else:
                msg = ("\nColumn 1 of file: {0:s} contains wrong type!"
                       .format(input_file))
                raise aXeSIMError(msg)

        # go over all rows
        for index in range(indata.nrows):

            # check for None-entries
            if ((indata[0][index] is None) or (indata[1][index] is None)):
                msg = ("\nData in file: {0:s} contains NULL"
                       "entries!".format(input_file))
                raise aXeSIMError(msg)

            # check for ascending order
            if index > 0:
                # check that independent data is rising
                if indata[0][index] <= indata[0][index-1]:
                    msg = ("\nIndependent column data in file: {0:s}"
                           " is not monotonically rising!".format(input_file))
                    raise aXeSIMError(msg)

    def _get_fits_name(self, fits_name):
        """Determine the proper fits name

        If not explicitly given, the method determines a proper fits name
        for an interpolator. As a base serves the file name the
        interpolator was built from.

        Parameters
        ----------
        fits_name: str
            a fits file name

        Returns
        -------
        fits_name: str
            a fits name for the interpolator data
        """
        # check whether an explicit name
        # is given
        if fits_name is not None:
            # return that explicit name
            return fits_name

        # check whether a root for the
        # fits name is there
        if self.input_file is None:
            # give an error and out
            error_message = ("\nCan not derive a proper name for the fits "
                             "file.\nPlease specify a name explicitly!")
            raise aXeSIMError(error_message)
        else:
            # find the position of the
            # last dot
            pos = self.input_file.rfind('.')

            # check whether there is a dot
            if pos > -1:
                # compose the new name
                fits_name = self.input_file[:pos] + '.fits'
            else:
                # compose the new name
                fits_name = self.input_file + '.fits'

        return fits_name

    def writetofits(self, fits_name=None, colname1=None, colname2=None):
        """Write the interplator values to a fits file

        The method writes the data of an interpolator to a binary fits table.
        The fiter name as well as the column names of independent and
        dependent data are specified explicitly.

        Parameters
        ----------
        fits_name: str
            the fits file name
        colname1: str
            column name for independent data
        colname2: str
            column name for dependent data

        Returns
        -------
        fits_name: str
            the fits file name
        """
        # get the fits name
        out_name = self._get_fits_name(fits_name)

        # create a new table
        out_tab = axe_asciidata.create(2, len(self))

        # go over all data
        for index in range(len(self)):
            # store the data in the table
            out_tab[0][index] = self._indep_data[index]
            out_tab[1][index] = self._depen_data[index]

        # rename the columns, if
        # there are names
        if colname1 is not None:
            out_tab[0].rename(colname1)
        if colname2 is not None:
            out_tab[1].rename(colname2)

        # write the table to fits
        out_tab.writetofits(out_name)

        # return the fits name
        return out_name

    def writeto(self, filename):
        """Write the interplator values to an ASCII

        The method writes the data of an interpolator to a binary fits table.
        The fiter name as well as the column names of independent and
        dependent data are specified explicitly.

        Parameters
        ----------
        filename: str
            the output file name

        Returns
        -------
        filename: str
            the file name
        """

        # create a new table
        out_tab = axe_asciidata.create(2, len(self))

        # go over all data
        for index in range(len(self)):
            # store the data in the table
            out_tab[0][index] = self._indep_data[index]
            out_tab[1][index] = self._depen_data[index]

        # write the table to fits
        out_tab.writeto(filename)

        # return the fits name
        return filename

    def integrate(self):
        """Evaluate the integral over the dependent values

        The method computes and returns the integral over the
        interpolator values. Only the fixed independent and dependent
        data values are used for computing the integral by simple
        adding of the differnetial elements.

        @return: the integral over the interpolator values
        @rtype: float
        """
        # initialize the
        # integrated value
        integral = 0.0

        # go over all data
        for index in range(1, len(self)-1):
            # add the increment
            integral += self._depen_data[index] * (self._indep_data[index+1] -
                                                   self._indep_data[index-1])

        # dont forget the start and end piece
        integral += self._depen_data[0] * (self._indep_data[1] -
                                           self._indep_data[0])
        integral += self._depen_data[len(self)-1] * (self._indep_data[len(self)-1] - self._indep_data[len(self)-2])

        # return the integral,
        # correcting the bin extension
        return integral/2.0

    def toSensitivity(self, A=None):
        """Transfer the bandpass to sensitivity

        The method performs all steps to transform an interpolator
        representing a total passband curve to a sensitivity curve.

        Parameters
        ----------
        A: float
            collecting area of telescope
        """
        h_erg = 6.6260693E-27
        c_cm = 2.99792458E+10

        # give a default for A
        if A is None:
            A = math.pi * 120.0 * 120.0

        # go over all data
        for index in range(len(self)):

            # compute the conversion factor
            factor = A / (h_erg * c_cm / (1.0E-08 * self._indep_data[index]))

            # apply the conversion factor
            self._depen_data[index] = self._depen_data[index] * factor

    def toThroughput(self, A=None):
        """Transfer the sensitivity to a passband

        The method performs all steps to transform an interpolator
        representing a sensitivity curve to a total passband for a given
        collecting area.

        Parameters
        ----------
        A: float
            collecting area of telescope
        """
        h_erg = 6.6260693E-27
        c_cm = 2.99792458E+10

        #  give a default for A
        if A is None:
            A = math.pi * 120.0 * 120.0

        # go over all data
        for index in range(len(self)):

            # compute the conversion factor
            factor = A / (h_erg * c_cm / (1.0E-08 * self._indep_data[index]))

            # apply the conversion factor
            self._depen_data[index] = self._depen_data[index] / factor

    def tonm(self):
        """Transfer the independent column to unit [nm]

        Provided the unit of the independent column is [AA], this method
        transforms the unit to [nm]
        """
        self.mult_indep(0.1)

    def tofits(self, colname1=None, colname2=None):
        """Transfer the data to a fits extension

        The method writes the data of an interpolator to a binary fits
        table extension, which is returned to the calling routine.
        The fiter name as well as the column names of independent and
        dependent data are specified explicitly.

        Parameters
        ----------
        colname1: str
            column name for independent data
        colname2: str
            column name for dependent data

        Returns
        -------
        table.hdu()
            the interpolator as fits table extension
        """
        # create a new ascii table
        new_table = axe_asciidata.create(2, len(self))

        # go over all data
        for index in range(len(self)):
            # transfer the values
            new_table[0][index] = self._indep_data[index]
            new_table[1][index] = self._depen_data[index]

        # re-name the table column
        if colname1 is not None:
            new_table[0].rename(colname1)
        if colname2 is not None:
            new_table[1].rename(colname2)

        # return the fits version
        return new_table.tofits()

    def mult_indep(self, factor):
        """Multipy the independent column

        The method multiplies all independent data values with
        a given factor.

        @param factor: multiplication factor
        @type factor: float

        @return: None
        @rtype: None
        """
        # go over all data
        for index in range(len(self)):
            # apply the conversion factor
            self._indep_data[index] = self._indep_data[index] * factor

        # also adjust the min and max
        self.ind_min = self.ind_min * factor
        self.ind_max = self.ind_max * factor

    def mult_depen(self, factor):
        """Multipy the dependent column

        The method multiplies all dependent data values with
        a given factor.

        Parameter
        ---------
        factor: float
            multiplication factor
        """
        # go over all data
        for index in range(len(self)):
            # apply the conversion factor
            self._depen_data[index] = self._depen_data[index] * factor

    def pivot(self):
        """Compute the pivot wavelength

        Provided the interpolator represents a total passband, this
        method computes the pivot wavelength of the passband.

        Returns
        -------
        pivot: float
            the pivot wavelength
        """
        import copy
        import math

        # create two copies
        nomin = copy.deepcopy(self)
        denom = copy.deepcopy(self)

        # go over all points
        for index in range(len(nomin)):
            # build the nominator and denominator interpolator
            nomin._depen_data[index] = nomin._depen_data[index] * nomin._indep_data[index]
            denom._depen_data[index] = denom._depen_data[index] / denom._indep_data[index]

        # compute the pivot wavelength
        pivot = math.sqrt(nomin.integrate()/denom.integrate())

        # return the pivot wavelength
        return pivot

    def toflambda(self):
        """Converts a mag_AB value at a wavelength to f_lambda

        Given an interpolator that represents a spectrum in
        lambda [Angstrom] - f_nu, the method converts to
        lambda [Angstrom] - f_lambda. NO unit scaling in
        f_lambda is applied!

        """
        # go over all data
        for index in range(len(self)):

            # extract the f_nu value
            fnu = self._depen_data[index]

            # extract the wavelength, converted to nm
            wlength = self._indep_data[index] / 10.0

            # compute the corresponding f_lambda value
            flambda = 2.99792458e+16*fnu/(wlength*wlength)

            # set the f_lambda value
            self._depen_data[index] = flambda
