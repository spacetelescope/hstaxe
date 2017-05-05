"""
Fits related classes

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-01-08 18:13:38 +0100 (Tue, 08 Jan 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciifits.py $
"""
__version__ = "Version 1.0 $LastChangedRevision: 329 $"

import string
import os
import pyfits

class AsciiFits(object):
    """
    Class for all fits-related things
    """
    def __init__(self, asciiData):
        """
        Initializes the class

        @param asciiData: an AsciiData object
        @type asciiData: AsciiData
        """
        # point to using numpy
        self.use_numpy = self._get_numpy_indicator()

        # transform the AsciiData object to a
        # pyfits table hdu object
        self.tabhdu    =  self._to_tabhdu(asciiData)

    def _get_numpy_indicator(self):
        """
        Determine whether using numpy or numarray

        @return: indicator for numpy
        @rtype: int
        """
        # set the default to numpy
        use_numpy = 1

        try:
            # get the 'numerix' value from pyfits
            numerix = pyfits.numerix

            # set the indicator according
            # to the 'numerix' value
            if numerix == 'numarray':
                use_numpy = 0
            else:
                use_numpy = 1

        # catch if 'numerix' does not exist
        except AttributeError:
            # make a default
            use_numpy = 0

            # get the pyfits version
            pyfits_version = pyfits.__version__.split('.')

            # determine the version number
            main_version = int(pyfits_version[0])
            if main_version < 1:
                # early versions have
                # only numarray
                use_numpy = 0
            elif main_version == 1:
                # get the next version number
                sub_version = int(pyfits_version[1])

                # 1.0 still gets numarray
                if sub_version < 1:
                    use_numpy = 0
                # higher versions get numpy
                else:
                    use_numpy = 1
            elif main_version > 1:
                # later versions use numpy
                use_numpy = 1

        # return the indicator
        return use_numpy

    def _to_tabhdu(self, asciiData):
        """
        Create a pyfits table HDU from an AsciiData object

        @param asciiData: an AsciiData object
        @type asciiData: AsciiData

        @return: the table HDU instance
        @rtype: tabHDU
        """

        # create the columns
        fits_cols = self._create_fits_cols(asciiData)

        # create the table HDU instance
        tabhdu = pyfits.new_table(fits_cols)

        # transfer the header of the AsciiData instance
        # to the fits header
        self._append_header(tabhdu, asciiData)

        # return the table HDU
        return tabhdu

    def _append_header(self, tabhdu, asciiData):
        """
        Store the AsciiData header in the fits HDU

        @param asciiData: an AsciiData object
        @type asciiData: AsciiData
        @param: tabhdu: the table HDU instance
        @type: tabhdu: tabHDU
        """

        # extract the header of the AsciiData instance
        asciiHeader = asciiData.header

        # extract the header of the table HDU
        theader = tabhdu.header

        # go through all lines in the header
        for line in asciiHeader.hdata:

            # check whether there is indeed some content
            if len(string.strip(line)) > 0:

                # add the header in a history line
                theader.add_history(line)

    def _create_fits_cols(self, asciiData):
        """
        Create pyfits columns from AsciiColumns

        @param asciiData: an AsciiData object
        @type asciiData: AsciiData

        @return: list of pyfits columns
        @rtype: [pyfits.Column]
        """

        # initialize the column list
        fits_cols = []

        # go over all AsciiColumns
        for ii in range(asciiData.ncols):

            # determine the column format
            format = self._get_fcolformat(asciiData[ii])

            # create and append a fits column to the list,
            # using either numpy of numarray
            if self.use_numpy:
                fits_cols.append(pyfits.Column(name=asciiData[ii].colname,
                                               format = format,
                                               array=asciiData[ii].tonumpy()))
            else:
                fits_cols.append(pyfits.Column(name=asciiData[ii].colname,
                                               format = format,
                                               array=asciiData[ii].tonumarray()))

        # return the column list
        return fits_cols


    def _get_fcolformat(self, asciiColumn):
        """
        Determine the fits column format

        @param asciiColumn: an AsciiColumn object
        @type asciiColumn: AsciiColumn

        @return: fits column format
        @rtype: string
        """

        # determine the column type
        ctype = asciiColumn.get_type()

        # check whether the column format is float
        if ctype == type(1.0):
            format = 'E'
        # check whether the column format is int
        elif ctype == type(1):
            format = 'J'
        # check whether the column format is string
        elif ctype == type('a'):
            maxlen = self._get_maxlen(asciiColumn)
            format = str(maxlen) + 'A'
        else:
            err_msg = 'Column type: ' + str(ctype) + 'not known in fits!'
            raise Exception(err_msg)

        # return the format
        return format

    def _get_maxlen(self, asciiColumn):
        """
        Determine the maximum string length

        @param asciiColumn: an AsciiColumn object
        @type asciiColumn: AsciiColumn

        @return: maximum string length
        @rtype: int
        """

        # initialize the maximum length
        maxlen = 0

        # go over all elements
        for ii in range(asciiColumn.get_nrows()):

            if asciiColumn[ii] == None:
                continue

            # check and perhaps replace tthe
            # maximum length
            if len(string.strip(asciiColumn[ii])) > maxlen:
                maxlen = len(string.strip(asciiColumn[ii]))

        # return the maximum lengt
        return maxlen

    def flush(self, fits_name='asciiData.fits'):
        """
        Write the object to a fits file
        """

        # check whether a file with this name
        # does just exist
        if os.path.isfile(fits_name):
            # delete the old file
            os.unlink(fits_name)

        # write the fits file
        self.tabhdu.writeto(fits_name)
