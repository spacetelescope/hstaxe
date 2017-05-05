"""
Main class of the asciidata module

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-01-08 18:17:08 +0100 (Tue, 08 Jan 2008) $
$LastChangedRevision:  $
$HeadURL: $
"""
__version__ = "Version 1.1 $LastChangedRevision: 330 $"

import string, sys, os, types,copy
from asciiheader import *
from asciicolumn import *
from asciisorter import *
from asciierror  import *
from asciiutils  import *

class NullData(object):
    """
    Null class as a parent class for the AsciiData class

    This parent classs of the AsciiData class offers to create
    a new AsciiData instance without a file to read from.
    All elements are set to None, but of course can later
    be filled by the user.
    """
    def __init__(self, ncols, nrows, null=None):
        """
                Constructor for the NullData Class

        Creates an empty AsciiData instance with columns and
        rows as specified. All entries are 'None'.

        @param ncols: the number of columns to be created
        @type ncols: integer
        @param nrows: the number of rows to be created
        @type nrows: integer
        @param null: string to be interpretet as NULL
        @type null: string
        """
        # set the default null string
        if null:
            self._null = [string.strip(null)]
        else:
            self._null  = ['Null']

        # create the colum list
        self.columns = []
        for index in range(ncols):
            # get the column name
            colname = self._def_colname(index)
            # create and append an empty column
            self.columns.append(AsciiColumn(nrows=nrows, colname=colname,
                                            null=self._null))
    def _def_colname(self, index):
        """
        Gives the default column name.

        The method composes and returns the
        default column name for a column at a
        given index.

        @param index: the index of the column
        @type index: integer
        """
        return 'column'+str(index+1)

class AsciiData(NullData):
    """
    Basic class in the AstroAsciiData project

    This class and its methods forms the complete API for the
    for the
    """
    def __init__(self, filename=None, ncols=0, nrows=0, null=None,
                 delimiter=None, comment_char=None, columnInfo=0, headerComment=1):
        """
        Constructor for the AsciiData Class

        The data is taken from  a file specified in the input.
        As addition, a NULL string, a delimiter string and a comment_char
        string can be specified. The ascii data is read in from the
        file and stored in a list of Columns

        @param filename: the filename to create the AsciiData from
        @type filename: string
        @param ncols: the number of columns to be created
        @type ncols: integer
        @param nrows: the number of rows to be created
        @type nrows: integer
        @param null: string to be interpretet as NULL
        @type null: string
        @param delimiter: string to be used as delimiter
        @type delimiter: string
        @param comment_char: string to be used as comment character
        @type comment: string
        """
        self.ncols = 0
        self.nrows = 0


        # set the default comment_char
        if comment_char:
            self._comment_char = comment_char
        else:
            self._comment_char = '#'

        # set the default null string
        if null:
            self._null = [string.strip(null)]
        else:
            self._null  = ['Null', 'NULL', 'None', '*']

        # set the delimiter
        self._delimiter = delimiter

        # set the separator
        self._separator = Separator(delimiter)

        # create the header
        self.header = Header(filename, self._comment_char)

        # check whether a filename is given
        if filename != None:
            # check whether the file exists
            if os.path.exists(filename):
                self.filename = filename
            else:
                err_msg = "Filename: " + filename + " does not exist!"
                raise Exception(err_msg)

            # set public output flags
            if self.header.SExtractorFlag:
                self.columnInfo = 1
                self.headerComment = 1
            else:
                self.columnInfo = 0
                self.headerComment = 1

           # load in all data from the files
            self.columns = self._load_columns(filename, self._null,
                                              self._comment_char, self._separator)
        else:
            # set the filename to none
            self.filename = None

            # check whether valid numbers where given
            if nrows > 0 and ncols > 0:
                # create the empty instance
                super(AsciiData, self).__init__(ncols, nrows, null)
            else:
                err_msg = "Number of columns, rows: " \
                      + str(ncols) + str(nrows)  + " are not reasonable!"
                raise Exception(err_msg)

            # set the public output flags
            # as the corresponding parameters
            self.columnInfo    = columnInfo
            self.headerComment = headerComment

        # find the number of undefined columns
        self._undef_cols = self._find_undefined_cols(self.columns)

        # find the number of columns and rows
        self.ncols = len(self.columns)
        if self.ncols:
            self.nrows = self.columns[0].get_nrows()

    def __getitem__(self, element):
        """
        Defines the list operator for indexing

        This method returns the index or indices as specified
        in the input. In the current class therefore returns
        either a column or a column slice as specified in the input.

        @param element: either column index or slice or name
        @type element: string/integer

        @return: a column
        @rtype: AsciiColumn(s)
        """
        # this part deals with slices
        if type(element) == types.SliceType:
            # FIXME this must be possible to do more elegantly
            start,stop,step = element.indices(self.ncols)
            newAD = copy.deepcopy(self)
            all = range(self.ncols)
            inclusive = [x for x in all[start:stop:step]]
            while all:
                idx = all.pop()
                if not idx in inclusive:
                    del newAD[idx]
            return newAD

        # this part deals with individual
        # columns, specified by index or name
        try:
            index = self._loc_column(element)
        except ColumnError:
            index = self.append(element)

        # return the desired column
        return self.columns[index]

    def __setitem__(self, element, column):
        """
        Defines the list operator for indexed assignement

        The method inserts a column to the class at the
        specified index. As of now, it is not possible
        to create extra columns with this method,
        only existing columns can be overwritten.

        @param element: either column index or name
        @type element: string/integer
        @param column: the column to assign to an index
        @type column: AsciiColumn
        """

        index = self._loc_column(element)

        # check whether the column does have the same number
        # of rows as the class
        # raise an error if not
        if column.get_nrows() != self.nrows:
            err_msg = 'Nrows: '+str(column.get_nrows())+' different than nrows: '\
        +str(self.nrows)+'!!'
            raise Exception(err_msg)

        # check whether the column has a name
        if not column.colname:
            # give it a default name
            column.colname = self._def_colname(index)

        # assign the new column
        self.columns[index] = column.copy()

        # transfer the null element to the new column
        self.columns[index]._null[0] = self._null[0]

    def __delitem__(self, element):
        """
        Deletes an index.

        The method deletes a column specified in the input.
        The column can be specified either by the column
        name or the index.

        @param element: either column index or name
        @type element: string/integer
        """
        # get the index from the input
        index = self._loc_column(element)

        # delete the column
        del self.columns[index]

        # adjust the number of columns
        self.ncols -= 1

    def __iter__(self):
        """
        Provide an iterator object.

        The function provides and returns an interator object
        for the AstroAsciiData class. Due to this iterator object
        sequences like:
        for column  in ascii_data_object:
            <do something with column>
        are possible.
        """
        return AsciiLenGetIter(self)

    def __len__(self):
        """
        Defines a length method for the object

        @return: the length of the object
        @rtype: integer
        """
        return self.ncols

    def str(self):
        """
        Defines a string method for the object.

        Gives a simple string method such that str(AsciiData)
        does work. The formatting is close to the formatting
        for the output to files.

        @return: the string representation of the object
        @rtype: string
        """
        bigstring = ''

        # take the object delimiter or ' '
        if not self._delimiter:
            delim = ' '
        else:
            delim = self._delimiter

        # add the header to the string
        bigstring = bigstring + str(self.header)

        # go over each row
        for ii in range(self.nrows):
            # create the string list
            strlist = self._row_tostring(ii)

            # treat the first line different
            if ii:
                # transform the listing to one string and append it
                # put a linefeed at the beginning
                bigstring = bigstring + '\n' + delim.join(strlist)
            else:
                # transform the listing to one string and append it
                bigstring = bigstring + delim.join(strlist)

        return bigstring

    def __str__(self):
        """
        Defines a string method for the object.

        Gives a simple string method such that str(AsciiData)
        does work. The formatting is close to the formatting
        for the output to files.

        @return: the string representation of the object
        @rtype: string
        """
        bigstring = ''

        # take the object delimiter or ' '
        if not self._delimiter:
            delim = ' '
        else:
            delim = self._delimiter

        # print the column information
        if self.columnInfo:
            for n, col in enumerate(self.columns):
                bigstring += str(col.collheader(n,self._comment_char))

        # print the header
        if self.headerComment:
            bigstring += str(self.header)

       # go over each row
        for ii in range(self.nrows):
            # create the string list
            strlist = self._row_tostring(ii)

            # treat the first line different
            if ii:
                # transform the listing to one string and append it
                # put a linefeed at the beginning
                bigstring = bigstring + '\n' + delim.join(strlist)
            else:
                # transform the listing to one string and append it
                bigstring = bigstring + delim.join(strlist)

       # return the string
        return bigstring


    def flush(self):
        """
        Prints the current status to the file.

        The methods gives the opportunity to replace the data in
        the AsciiData with the current version in memory.
        """
        if self.filename != None:
            # well, that an easy job
            self.writeto(self.filename)
        else:
            raise Exception('No filename given. Use "writeto()" instead.')

    def writeto(self, filename, colInfo=None, headComment=None):
        """
        Prints the AsciiData to a new file.

        The method prints the current status of the
        object to a new file. The name of the file
        is given in the input. An already existing
        file is replaced.

        @param filename: the filename to write the object to
        @type filename: string
        """
        # check whether the parameter is set
        if colInfo==None:
            # if not, take the class variable
            colInfo = self.columnInfo

        # check whether the parameter is set
        if headComment == None:
           # if not, take the class calue
           headComment = self.headerComment

        # open the file
        fstream = file(filename,'w+')

        # open a printstream
        nprinter = NicePrinter(fstream, delimiter=self._delimiter)

        # print everything to the stream
        self._print_tostream(nprinter, colInfo, headComment)

        #close the file
        fstream.close()

        # use the given name as class filename
        # if no one is yet defined
        if self.filename == None:
            self.filename = filename


    def tofits(self):
        """
        Transforms the AsciiData object to fits

        @return: pointer to the fits object
        @rtype: binary table HDU
        """
        import asciifits

        # create an AsciiFits object
        asciiFits = asciifits.AsciiFits(self)

        # return the table HDU
        return asciiFits.tabhdu


    def writetofits(self, fits_name=None):
        """
        Prints the AsciiData to a new file.

        @param fits_name: the name for the fits file
        @type fits_name: string

        @return: the name of the fits file
        @rtype: string
        """
        import asciifits

        # check whether a file name is given
        if fits_name == None:
            # check wheter the instance has a filename
            if self.filename == None:
                # no automatic filename possible; raise error
                raise Exception('Please specify a name for the fits-file!')
            else:
                # determine a filename for the fits
                fits_name = self._get_fitsname(self.filename)

        # create an AsciiFits object
        asciiFits = asciifits.AsciiFits(self)

        # write out the object onto disk
        asciiFits.flush(fits_name)

        # return the name of the fits object
        return fits_name

    def writetohtml(self, html_name=None, tr_attr=None, td_attr=None):
        """
        Prints the AsciiData object as table in a html-file

        @param filename: the filename to write the object to
        @type filename: string
        @param tr_attr: the attributes for the tr-tag
        @type tr_att: string
        @param td_attr: the attributes for the td-tag
        @type td_att: string

        @return: the name of the html-file
        @rtype: string
        """
        # check whether a file name is given
        if html_name == None:
            # check wheter the instance has a filename
            if self.filename == None:
                # no automatic filename possible; raise error
                raise Exception('Please specify a name for the html-file!')
            else:
                # determine a filename for the html-file
                html_name = self._get_htmlname(self.filename)

        # determine the line start, element delimiter and the line end
        l_start, l_delim, l_end = self._get_lineparams(tr_attr, td_attr)

        # open the file
        fstream = file(html_name,'w+')

        # open a printstream
        nprinter = NicePrinter(fstream, delimiter=l_delim,
                               linestart=l_start, linend=l_end)

        # print the data
        # go over each row
        for ii in range(self.nrows):
            # create the string list
            strlist = self._row_tostring(ii)

            # send the list to the printer
            nprinter.print_list(strlist)

        #close the file
        fstream.close()

        # return the filename
        return html_name


    def writetolatex(self, latex_name=None):
        """
        Prints the AsciiData object as table in a latex-file

        @param filename: the filename to write the object to
        @type filename: string

        @return: the name of the latex-file
        @rtype: string
        """
        # check whether a file name is given
        if latex_name == None:
            # check wheter the instance has a filename
            if self.filename == None:
                # no automatic filename possible; raise error
                raise Exception('Please specify a name for the latex-file!')
            else:
                # determine a filename for the latex-file
                latex_name = self._get_latexname(self.filename)

        # open the file
        fstream = file(latex_name,'w+')

        # open a printstream with the correct parameters
        # please note that each '\' must be protected by
        # another '\' to be interpreted as string
        nprinter = NicePrinter(fstream, delimiter='&', linend='\\\\\n')

        # print the data
        # go over each row
        for ii in range(self.nrows):
            # create the string list
            strlist = self._row_tostring(ii)

            # send the list to the printer
            nprinter.print_list(strlist)

        #close the file
        fstream.close()

        # return the filename
        return latex_name


    def info(self):
        """
        Print class info to the screen.

            The method gives some basic information on the
            class. The output is directly written onto
            the screen.

        @return: the string representing the information
        @rtype: string
        """

        # define the return string
        bigstring = ''

        # assemble the basic table information
        bigstring += 'File:       ' + str(self.filename) +'\n'
        bigstring += 'Ncols:      ' + str(self.ncols) + '\n'
        bigstring += 'Nrows:      ' + str(self.nrows) + '\n'
        bigstring += 'Delimiter:  ' + str(self._delimiter) + '\n'
        bigstring += 'Null value: ' + str(self._null) + '\n'
        bigstring += 'comment_char:    ' + str(self._comment_char) + '\n'

            # go over each column and add
            # the individual column info
        for col in self.columns:
            bigstring += col.info()

        # return the result
        return bigstring


    def append(self, colname):
        """
            Appends a new column to the object.

            This method creates and appends a new column to the
            object. The new column must be specified with a name.
            The new column doe have only Null entries.

        @param colname: the name of the column
        @type colname: string
            """
            # check whether the column name does exist
        # raise a warning if yes
        if self.find(colname) > -1:
            err_msg = 'Column with name: '+colname+' does just exist!'
            raise Exception(err_msg)

        # get the index of the new column
        index = self.ncols

        # create and append the new column
        self.columns.append(AsciiColumn(nrows=self.nrows, colname=colname,
                                        null=self._null))

        # adjust the number of columns
        self.ncols +=1

        #return the index of the column
        return index

    def find(self, colname):
        """
        Finds the column number for a name.

        The method looks through all columns of the instance
        for a matching column name. In case the column name exists,
        the column index is returned. If the column name does
        not exist, -1 is returned.

        @param colname: the name of the column
        @type colname: string

        @return: the index of the column, or -1
        @rtype: integer
        """
        for index in range(len(self.columns)):
            if self.columns[index].colname == colname:
                return index
        return -1

    def delete(self, start, end=None):
        """
        Deletes a row slice or element from all columns.

        The method deletes one or several rows from all columns.
        It uses the __delelte__ or __delitem__ operators
        in the AsciiColumn class.

        @param start: the starting row index
        @type start: integer
        @param end: the end row index
        @type end: integer
        """
        if end:
            if start < self.nrows:
                # go over each column
                for col in self.columns:
                    # delete the row
                    del col[start: end]

                # adjust the number of rows
                self.nrows -= end-start

        else:
            # go over each column
            for col in self.columns:
                        # delete the row
                del col[start]

            # adjust the number of rows
            self.nrows -= 1

        # make a lower limit to the number of rows
        if self.nrows < 0:
                   self.nrows = 0


    def newcomment_char(self, comment_char):
        """
        Define a new comment_char string

        @param comment_char: the new null string
        @type comment_char: string
        """
        # store the new null element
        self._comment_char = comment_char
        self.header.set_comment_char(comment_char)

    def newnull(self, newnull):
        """
        Define a new null string

        @param newnull: the new null string
        @type newnull: string
        """
        # store the new null element
        self._null[0] = newnull

        # store the new null in the columns
        for column in self.columns:
            column._null[0] = newnull

    def newdelimiter(self, delimiter):
        """
        Set a new delimiter string

        @param delimiter: the new delimiter string
        @type delimiter: string
        """
            # set the new delimiter
        self._delimiter = delimiter

            # set the separator
        self._separator = Separator(delimiter)

    def insert(self, nrows, start=0):
        """
        Inserts one or several rows

        The method inserts one or several rows into all
        columns of the class instance. The number of rows
        as well as the positioning of the new rows are
        specified in the input. The parameter 'start'
        gives the index which the first inserted row
        will have.
        Setting "start=-1" means appending the rows at
        the end of the columns

        @param nrows: the number of rows to add
        @type nrows: integer
        @param start: the position of the inserted rows
        @type start: integer
        """

        # go over all columns
        for col in self.columns:
            # add none elements at the end
            for ii in range(nrows):
                col.add_element(None)

        # check whether the new rows are inserted inside
        # the old rows, then the elements must be moved
        if start < self.nrows and start != -1:

            # go over all columns
            for col in self.columns:

                # repeat over rows to be inserted
                for ii in range(self.nrows-start):
                    # reorder the column elements
                    index = self.nrows - ii - 1
                    col[index+nrows] = col[index]

                # repeat over rows to be inserted
                for ii in range(nrows):
                    # insert None in the new rows
                    index = ii + start
                    col[index] = None


        # update the number of rows
        self.nrows = self.columns[0].get_nrows()


    def sort(self, colname, descending=0, ordered=0):
        """
        Sorts the entries along the values in one column

        The method sorts all columns of the AsciiData object according
        to the order in one specified column. Both, sorting in ascending
        and descending order is possible.

        @param colname: the column to use for sorting
        @type colname: string/integer
        @param descending: indicates ascending (=0) or descending (=1) sorting
        @type descending: integer
        @param ordered: indicates ordered (1) or non-ordered sorting
        @type ordered: integer
        """
        # initialize a temporary array
        sort_data = []

        # transfer the data from the sort column
        # to the temporary array
        for index in range(self.nrows):
            sort_data.append(self[colname][index])

        # create the sorting index
        sorter = ColumnIndex()

        # sort according to the data in the temporary array
        sorter.sort(sort_data, descending, ordered)

        # go over all colums
        for index in range(self.ncols):
            # reorder the data in the column according
            # to the sorting order
            self[index]._data = sorter.enindex(self[index]._data)

    def rstrip(self,x=None):
        '''
        Removes trailing rows which contain the value of x
        null is default (and the only value which really works)
        syntactic sugar for _strip(-1,x)
        @param x: Data value in rows to strip of - defaults to Null
        @type x: any legal asciidata type
        '''
        self._strip(-1,x)


    def lstrip(self,x=None):
        '''
        Removes leading rows which contain the value of x
        null is default (and the only value which really works)
        syntactic sugar for _strip(0,x)
        @param x: Data value in rows to strip of - defaults to Null
        @type x: any legal asciidata type
        '''
        self._strip(0,x)


    def strip(self,x=None):
        '''
        Removes both leading and trailing rows which contain the value of x
        null is default (and the only value which really works)
        syntactic sugar for _strip
        @param x: Data value in rows to strip of - defaults to Null
        @type x: any legal asciidata type
        '''
        self._strip(-1,x)
        self._strip(0,x)

    def toSExtractor(self):
        """
        convenience function to set the ouput to be in SEextractor style
        """
        self.headerComment = 1
        self.columnInfo = 1
        self.newcomment_char('#')
        self.newdelimiter(' ')

    def toplain(self):
         """
         convenience procedure to toggle to plain ACSII output
         delimiters are not changed
         """
         self.headerComment = 1
         self.columnInfo = 0

    def _get_fitsname(self, filename):
        """
        Determines the fitsname for a given file name

        @param filename: the input filename
        @type filename: string

        @return: the name of the fits file
        @rtype: string
        """

        # search for the extension
        dot_pos =  filename.rfind('.')

        # if an extension exists
        if dot_pos > -1:
            # replace the old extension with '.fits'
            fits_name = filename[:dot_pos] + '.fits'

        else:
            # append the extension '.fits'
            fits_name = filename + '.fits'

        # return the fits name
        return fits_name

    def _get_htmlname(self, filename):
        """
        Determines the html name for a given file name

        @param filename: the input filename
        @type filename: string

        @return: the name for the html file
        @rtype: string
        """

        # search for the extension
        dot_pos =  filename.rfind('.')

        # if an extension exists
        if dot_pos > -1:
            # replace the old extension with '.html'
            html_name = filename[:dot_pos] + '.html'

        else:
            # append the extension '.html'
            html_name = filename + '.html'

        # return the html name
        return html_name


    def _get_latexname(self, filename):
        """
        Determines the latex filename for a given file name

        @param filename: the input filename
        @type filename: string

        @return: the name for the latex file
        @rtype: string
        """

        # search for the extension
        dot_pos =  filename.rfind('.')

        # if an extension exists
        if dot_pos > -1:
            # replace the old extension with '.html'
            latex_name = filename[:dot_pos] + '.tex'

        else:
            # append the extension '.html'
            latex_name = filename + '.tex'

        # return the html name
        return latex_name


    def _get_lineparams(self, tr_attr=None, td_attr=None):
        """
        Prints the AsciiData object as table in html-file

        @param tr_attr: attributes for the tr-tag
        @type tr_attr: string
        @param td_attr: attributes for the td-tag
        @type td_attr: string

        @return: the html-table linestart, delimiter and lineend
        @rtype: string, string, string
        """
        # form the string for the tr-attributes
        if tr_attr == None:
            str_tr_add = ''
        else:
            str_tr_add = ' ' + tr_attr

        # form the string for the td-attributes
        if td_attr == None:
            str_td_add = ''
        else:
            str_td_add = ' ' + td_attr

        # compose linestart, delimiter and lineend
        lstart = '<tr'+str_tr_add+'><td'+str_td_add+'>'
        delim  = '</td><td'+str_td_add+'>'
        lend   = '</td></tr>\n'

        # return linestart, delimiter, lineend
        return lstart, delim, lend

    def _loc_column(self, element):
        """
        Localizes a column

        The method localizes the column from any possible input.
        Possible input is either the column name or column index.
        Basic checks are done whether the column exists.

        @param element: either column index or name
        @type element: string/integer

        @return: the column index
        @rtype: integer
        """
        # create an element
        elem = Element(element)

        # check the types and derive the column index
        if elem.get_type() == types.IntType:
            # check for -1, which indicates the last column
            if element == -1:
                # set the index of the last column
                index = self.ncols-1
            else:
                # set the index to the input index
                index = element
        elif elem.get_type() == types.StringType:
            index = self.find(element)

        # check whether the column index exists
        # raise an error if not
        if index > self.ncols-1:
            err_msg = 'Index: '+str(index)+' is larger than ncols: ' +str(self.ncols)+'!!'
            raise Exception(err_msg)
        elif index < 0:
            raise ColumnError('Column name: "'+element+'" does not exist!')

        # return the index
        return index

    def _load_columns(self, filename, null, comment_char, separator):
        """
            Transforms the content of a file into columns

        Opens the file, defines the columns, adds all data rows,
        and returns the columns.

        @param filename: the filename to create the AsciiData from
        @type filename: string
        @param null: string to be interpreted as NULL
        @type null: string
        @param separator: string to be used as delimiter
        @type separator: string
        @param comment_char: string to be used as comment character
        @type comment_char: string

        @return: the columns loaded
        @rtype: [AsciiColumn]
        """

        undef_cols = []
        collist    = []

        # open the file, and parse through all rows
        for line in file(filename, 'r'):

            # throw away trailing and leading whitespaces
            str_line = string.strip(line)
            if len(str_line) < 1 or str_line[0] == comment_char:
                continue

                # if collumns exist, add a row
            if collist:
                self._add_row(collist, line,  null, separator)
                # if columns do not exist, define them
            else:
                collist = self._define_cols(line,  null, separator)


        # return the column list
        return collist

    def _find_undefined_cols(self, collist):
        """
        Finds undefined columns

        The method finds undefined columns in a column list.
        An undefined column is a column with the flag "self._defined"
        not set. This means that column type and column format
        are not specified, and the column elements are Null.
        The indices of the undefined columns is returned as a list

        @param collist: the list of existing columns
        @type collist: list of AsciiColumns

        @return: a list with the indices of undefined columns
        @rtype: [integer]
        """
        undefined = []

        # go over each column
        index=0
        for col in collist:
            # check whether the column is defined
            # append the index to the list if not
            if not col.get_defined():
                undefined.append(index)
            # increment the index
            index = index+1

        # return the list
        return undefined

    def _add_row(self, collist, line,  null, separator):
        """
        Adds a line from the file to the column list.

        The method gets a line from the input file.
        The line is split up into its items.
        Then each item is added to the column
        it belongs to. Items matching the NULL
        string are added as "None". A delimiter
        is taken into account in the splitting,
        if specified.

        @param collist: the list of existing columns
        @type collist: list of AsciiColumns
        @param line: the line to be added to the columns
        @type line: string
        @param null: string to be interpretet as NULL
        @type null: string
        @param separator: string to be used as delimiter
        @type separator: string
        """
        # split the line, either according toa whitespace,
        # or according to a specified delimiter
        items = separator.separate(line)

        # check whether there is an item for each column
        if len(collist) != len(items):
            err_msg = "Number of columns does not fit to number of items in " + line
            raise Exception(err_msg)

        # go over each item
        index = 0
        for item in items:

            # check whether the item is NULL.
            # add the item to the column,
            # using 'None' for NULL items
            if null.count(string.strip(item)) > 0:
                collist[index].add_element(None)
            else:
                collist[index].add_element(item)

            # increment the index
            index += 1

    def _define_cols(self, line,  null, separator):
        """
        Defines the columns from an input line.

        The method splits an ascii line from the input file into its
        items. For each item a new column is created and added
        to a column list. The column list is finally returned.

        @param line: the line to be added to the columns
        @type line: string
        @param null: string to be interpretet as NULL
        @type null: string
        @param separator: string to be used as delimiter
        @type separator: string

        @return: the columns created
        @rtype: [AsciiColumn]
        """
        collist = []

        # split the line, either according toa whitespace,
        # or according to a specified delimiter
        items = separator.separate(line)

        # go over each item, and create a column
        # for each. NULL items are transformed to 'None'
        index = 0
        for item in items:

            # set the default column unit and comment
            colunit = ''
            colcomment = ''

            # check whether there is column
            # information from the header
            if self.header.SExtractorFlag:
                # extract the header information
                colname,colunit,colcomment = self.header.getCollInfo(index)
            else:
                # make the default column name
                colname = self._def_colname(index)

            # check whether the element is a NULL-value
            if null.count(string.strip(item)) > 0:
                # append an undefined column
                collist.append(AsciiColumn(element=[None], colname=colname,
                                           null=null))
            else:
                # append a defined column
                collist.append(AsciiColumn(element=[item], colname=colname,
                                           null=null))
            # transfer the resto of the column information
            if colunit:
                collist[-1].set_unit(colunit)
            if colcomment:
                collist[-1].set_colcomment(colcomment)

            # increment the index
            index += 1

        # return the column list
        return collist

    def _print_tostream(self, nprinter, colInfo, headComment):
        """
        Prints the AsciiData to a stream

        The method forms for each row in the AsciiData a list
        with formated strings, each list element representing
        one element. The list is sent to a printing stream
        which is responsible for the output.

        @param nprinter: the NicePrinter object with the stream
        @type nprinter: NicePrinter
        """
        # print the column information
        if colInfo:
            for n, col in enumerate(self.columns):
                nprinter.print_string(col.collheader(n,self._comment_char))

        # print the header
        if headComment:
            nprinter.print_string(str(self.header))

        # print the data
        # go over each row
        for ii in range(self.nrows):
            # create the string list
            strlist = self._row_tostring(ii)

            # send the list to the printer
            nprinter.print_list(strlist)

    def _row_tostring(self, index):
        """
        Creates the formatted string list for one row.

        The method extracts from each column the formatted
        string representation of the element in a specified
        row. The list of strings is returned.

        @param index:
        @type index: integer

        @return: the list with formatted strings
        @rtype: [string]
        """
        # initialize the list
        strlist = []

        # go over each column
        for jj in range(self.ncols):
            # append the string of the requested
            # element to the list
            strlist.append(self.columns[jj].fprint_elem(index))

        # return the list
        return strlist

    def _strip(self,rowindex, x=None):
     '''
     Removes rows which contain the value of x
     null is default (and the only value which really works)
     @param rowindex: select if it is lstrip (0) or rstrip (-1)
     @type rowindex: int
     '''
     while self.nrows>0:
         equal = True
         for col in self.columns:
             equal = equal and (col[rowindex] == x)
         if equal:
             self.delete(rowindex)
         else:
             break
