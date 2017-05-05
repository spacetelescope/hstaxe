"""
Various column classes

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-07-03 10:27:47 +0200 (Thu, 03 Jul 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciicolumn.py $
"""
__version__ = "Version 1.0 $LastChangedRevision: 503 $"

import string
from asciielement import *
from asciierror   import *
from asciiutils   import *

class NullColumn(object):
    """
    Class for an empty column with 'None'-elements
    """
    def __init__(self, nrows):
        """
            Constructor for the class

        @param nrows: the number of rows
        @type nrows: integer
            """

            # create the data array
        dummy_list = []

        # set the row number
        self._nrows   = nrows

            # append the reuqested number of None
        #       for index in range(nrows):
        #             self._data.append(None)

        # append the reuqested number of None
        # perhaps a bit faster than the above lines
        self._data = map(dummy_list.append, range(nrows))

        # set the row number
        self._nrows   = nrows


class AsciiColumn(NullColumn):
    """
    Class to store the ascii data into.
    """
    def __init__(self, element=None, colname=None, null=None, nrows=None):
        """
        Constructor for the column class.

        Instances of this column class hold the data in
        a private list. Moreover there exist few
        attributes in addition. A column does have
        a type, which is either string/integer/float.
        The column can be undefined, which means it contains
        only 'None', but the default type is string.

        @param element: list of elements to start the data with
        @type element: string/integer/float
        @param colname: the name of the column
        @type colname: string
        """
        self.colname = colname
        self.unit = ''
        self.colcomment =''
        self._data    = []
        self._defined = 0
        self._type    = types.StringType
        self._format  = ['%10s','%10s']
        self._nrows   = 0

        # set the default null string
        if null:
            self._null = [string.strip(null[0])]
        else:
            self._null  = ['*']

        if not element and nrows:
            super(AsciiColumn, self).__init__(nrows)
        else:
            # go over each element in the list
            for item in element:
                # check for 'None'
                if item != None:
                    # check whether the column is defined
                    if not self._defined:
                        # if undefined, set the type
                        # append the element, set to defined
                        elem = ForElement(item)
                        self._type   = elem.get_type()
                        self._format = elem.get_fvalue()
                        self._data.append(elem.get_tvalue())
                        self._defined = 1
                        self._nrows += 1
                    else:
                        # if defined, add the element
                        self.add_element(item)
                else:
                            # simply append the 'None'
                            self._data.append(item)

                            # increment the number of rows
                            self._nrows += 1

    def __getitem__(self, index):
        """
        Defines the list operator for indexing

        This method returns the value at a given index.
        In the current class this means the method returns
        the column value  at the requested index.

        @param index: the index of the column to be returned
        @type index: integer

        @return: the column value
        @rtype: string/integer/float
        """
        # check whether the requested index is available.
        # raise an error if not
        # [BUG] ?  self._nrows-1: ->  self._nrows:
        if index > self._nrows-1:
            err_msg = 'Index: '+str(index)+' is larger than nrows: '\
        +str(self._nrows)+'!!'
            raise Exception(err_msg)

        # return the value at the index
        return self._data[index]

    def __setitem__(self, index, value):
        """
        Defines the list operator for indexed assignement

        The method inserts a value into the column at the
        specified index. It is not possible to create
        extra rows with this method. Only existing
        elements can be overwritten.

        @param index: the index to put the colun to
        @type index: integer
        @param value: the value to assign to an index
        @type value: string/integer/float
        """
        # check whether the indexed element does exist
        # raise an error if not
        if index > self._nrows-1:
            err_msg = 'Index: '+str(index)+' is larger than nrows: '\
                 +str(self._nrows)+'!!'
            raise Exception(err_msg)

        if value != None:
                # check whether the column is defined
            if not self._defined:
                # create an element object
                val = ForElement(value)

                # if not set the type and the define flagg
                self._type    = val.get_type()
                self._format  = val.get_fvalue()
                self._defined = 1

            else:

                # create an element object
                val = ValElement(value)

                # if defined, check whether the element
                # type matches the column type
                if self._type != val.get_type():

                    # create a transformator object if the types do not match
                    type_trans = TypeTransformator(self._type, val.get_type())

                    # check whether the element is transformable
                    if type_trans.istransf:

                        # determine the transformed value
                        # the old code uses the typed value
                        # as the basis for the transformation:
                        #trans_value = type_trans.to_higher_type(val.get_tvalue())
                        # the new code uses the string value
                        # as the basis for the transformation:
                        trans_value = type_trans.to_higher_type(val.get_value())


                        val.set_tvalue(trans_value)

                    else:
                        # change the entire column type
                        self._change_column_type(type_trans, value)

            # set the column element to the transformed value
            self._data[index] = val.get_tvalue()

        else:
            self._data[index] = None


    def _change_column_type(self, t_trans, value):
        """
        Changes the type of a column

        The method changes the type of a column. It transformes
        all element into the new type and also defines the
        new type and formats.

        @param t_trans: the transformator object
        @type t_trans:
        @param value: the template value
        @type value: string/integer/float
        """
        # create an element object
        val = ForElement(value)

        # if not set the type and the define flagg
        self._format  = val.get_fvalue()

        # set the type to the one in the transformator object
        self._type    = t_trans.higher_type

        # go over all data
        for index in range(len(self._data)):
            if self._data[index] != None:
                # transform all non-Null entries
                self._data[index] = t_trans.to_higher_type(self._data[index])

    def _get_nullformat(self, newformat):
        """
        Find the null-format

        The method finds an appropriate format for the null
        elements for a given new format and the column type.
        This null-format may be smaller than needed to fully
        represent the null element.

        @param newformat: the new column format
        @type newformat: string
        @return: the format for the null elements
        @rtype: string
        """
        if self._type == types.IntType:
            length = len(str(newformat % 1))
            return '%'+str(length)+'s'
        elif self._type == types.FloatType:
            length = len(str(newformat % 1.0))
            return '%'+str(length)+'s'
        else:
            return newformat

    def __iter__(self):
        """
        Provide an iterator object.

        The function provides and returns an interator object
        for the AstroColumnData class. Due to this iterator object
        sequences like:
        for elem  in ascii_column_object:
            <do something with elem>
        are possible.
        """
#        return AsciiColumnIter(self)
        return AsciiLenGetIter(self)


    def __str__(self):
        """
        Print the column elements to the screen.

        The method prints the column name and the elements onto
        the screen. The format is column format, each
        element is written otno a new line.

        @return: the string representation of the column
        @rtype: string
        """
        # make the header
        bigstring = 'Column: '+str(self.colname)

        # go over each row
        for index in range(self._nrows):
            # append the string repr. to the output
            bigstring += '\n'+self.fprint_elem(index)

        # return the string
        return bigstring


    def __len__(self):
        """
        Defines a length method for the object

        @return: the length of the object
        @rtype: integer
        """
        return self._nrows


    def __delitem__(self, index):
        """
        Deletes an index.

        The method deletes a column row specified in the input.
        The column is specified by the index.

        @param index: row index
        @type index: integer
        """
        # delete the column
        del self._data[index]

        # adjust the number of columns
        self._nrows -= 1

    def  __delslice__(self, start, end):
        """
        Deletes an index slice.

        The method deletes a slice from the AsciiColumn
        data. Start and end index to be deleted are specfified
        in the input. This standard method redirect calls such
        as "del gaga[i:j]".

        @param start: starting row index
        @type start: integer
        @param end: ending row index
        @type end: integer
        """
        # delete the slice from the data
        del self._data[start:end]

        # determined the length of the data element
        self._nrows = len(self._data)

    def rename(self, newname):
        """
        Rename a column

        The method renames the column. The old column
        name is simply overwritten.

        @param newname: the new column name
        @type newname: string
        """
        # set the new column name
        self.colname = newname

    def reformat(self, newformat):
        """
        Gives a new column format

        The method gives a new formar to a column.
        The old column format is simply overwritten.

        @param newformat: the new column format
        @type newformat: string
        """
        # check whether the column is defined
        if self._defined:
            # get the appropriate null-format
            nullformat = self._get_nullformat(newformat)
            # set the new formats
            self._format = [newformat, nullformat]
        else:
            # first the column type must be defined
            raise Exception('The data type of this column is not yet defined!')


    def add_element(self, element):
        """
        Adds an element to the the column

            The method adds an element at the end of the data list
            of the column object. Type cheking is performed, and
            and error is thrown if the types do not match.

        @param element: string to be interpretet as NULL
        @type element: string/integer/float
        """
            # check for 'None'
        if element != None:

            # check whether the column is defined
            if not self._defined:
                # create an element object
                elem = ForElement(element)

                # if not, set the type and the define flagg
                self._type    = elem.get_type()
                self._format  = elem.get_fvalue()
                self._defined = 1

            else:
                # create an element object
                elem = ValElement(element)

            # if defined, check whether the element
            # type matches the column type
            if self._type != elem.get_type():

                # create a transformator object if the types do not match
                type_trans = TypeTransformator(self._type,elem.get_type())

                # check whether the element is transformable
                if type_trans.istransf:

                    # determine the transformed value
                    trans_value = type_trans.to_higher_type(elem.get_tvalue())
                    elem.set_tvalue(trans_value)


                else:
                    # change the entire column type
                    self._change_column_type(type_trans, element)



            # set the column element to the given value
            self._data.append(elem.get_tvalue())
#            print elem.get_tvalue()

        else:
                # append a 'None' element
            self._data.append(element)

        # increment the number of rows
        self._nrows += 1


    def fprint_elem(self, index):
        """
        Create and return a formatted string representation for an element.

        The method creates a formatted string representation
        for an element in an AsciiColumn. The element is specified
        by the row index. The string representation is returned.

        @param index: the index of the element
        @type index: integer

        @return: the string representation of the element
        @rtype: string
        """
        # check whether the row exists
        # raise an error if not
        if index > self._nrows-1:
            err_msg = 'Index: '+str(index)+' is larger than nrows: '\
            +str(self._nrows)+'!!'
            raise Exception(err_msg)
        # check for 'None'-entry
        if self._data[index] != None:
            # valid entries are formatted with the first format
            return self._format[0] % self._data[index]
        else:
            # None entries get the second format
            return self._format[1] % self._null[0]

    def tonumarray(self):
        """
        Transforms column to a numarray

        If possible, the column data is transformed to a
        numarray object and returned. Type specific numarrays
        are created to shorten the effort in the numarray module.

        @return: the numarray representation of the data
        @rtype: numarray
        """
        import numarray

        # initialize the return
        narray = None

        if None in self._data:
            raise Exception('There are "None" elements in the column. They can not be\ntransformed to numarrays!')

        # check for string column
        if self._type == types.StringType:
            # import CharArrays
            import numarray.strings
            # transform the array to CharArrays
            narray = numarray.strings.array(self._data)

        elif self._type == types.IntType:
            # transform the data to integer numarray
            narray = numarray.array(self._data, type='Int32')

        elif self._type == types.FloatType:
            # transform the data to float numarray
            narray = numarray.array(self._data, type='Float64')

        else:
            # raise an exception in case of string column
            err_msg = 'Can not transform column type: '+str(self._type)+' to numarray!'
            raise Exception(err_msg)
        # return the result
        return narray

    def tonumpy(self):
        """
        Transforms column to a numpy

        The column data is transformed to a numpy object
        and returned.

        @return: the numpy representation of the data
        @rtype: numpy/numpy masked array
        """
        import numpy
        from numpy import ma

        # initialize the return
        narray = None

        if None in self._data:

            # define a lambda function
            # to create the mask array
            make_mask = lambda x: x == None

            # create the numpy array,
            # making on the fly the mask
            narray = numpy.ma.array(self._data, mask=map(make_mask, self._data))

        else:
            # convert the list to a numpy object
            narray = numpy.array(self._data)

        # return the numpy object
        return narray


    def copy(self):
        """
        Returns a copy of the AsciiColumn

        The method creates a deep copy of the instance
        itself. The new AsciiColumn is then returned.

        @return: the copy of the current column
        @rtype: AsciiColumn
        """
        data_copy = []

        # make a copy of the data
        for ii in range(self._nrows):
            data_copy.append(self._data[ii])

        # explicitly create a column from the data copy
        self_copy = AsciiColumn(element=data_copy, colname=self.colname,
                                null=self._null)

        # explicitly transport the format
        self_copy._format = self._format

        # return the new column
        return self_copy

    def get_nrows(self):
        """
        Returns the number of rows.

        @return: the number of rows
        @rtype: integer
        """
        return self._nrows

    def get_type(self):
        """
        Returns the column type.

        @return: the column type
        @rtype: <types>-name
        """
        return self._type

    def set_type(self, type):
        """
        Sets the column type.

        @param type: the column type
        @type type: <types>-name
        """
        self._type = type

    def get_format(self):
        """
        Returns the column format

        @return: the format of the column
        @rtype: string
        """
        return self._format[0]

    def get_defined(self):
        """
        Returns the defined flagg.

        @return: the defined status of the column
        @rtype: integer
        """
        return self._defined

    def set_defined(self):
        """
        Sets the column status to defined.
        """
        self._defined = 1


    def set_unit(self,unit):
        """
        Sets the column unit

        @param unit: the column unit
        @type unit: string
        """
        self.unit = unit

    def get_unit(self):
        """
        Returns the column unit

        @return: the unit of the column
        @rtype: string
        """
        return self.unit

    def set_colcomment(self,colcomment):
        """
        Sets the colcomment

        @param colcomment: the column comment
        @type colcomment: string
        """
        self.colcomment = colcomment

    def get_colcomment(self):
        """
        Returns the column colcomment

        @return: the comment of the column
        @rtype: string
        """
        return self.colcomment


    def info(self):
        """
        Prints some column info onto the screen.

        @return: the string representing the information
        @rtype: string
        """
        # define the return string
        bigstring = ''

        # assemble the information on the column
        bigstring += 'Column name:        ' + self.colname + '\n'
        bigstring += 'Column type:        ' + str(self._type) + '\n'
        bigstring += 'Column format:      ' + str(self._format) + '\n'
        bigstring += 'Column null value : ' + str(self._null) + '\n'
        if self.unit:
            bigstring += 'Column unit : ' + self.unit + '\n'
        if self.colcomment:
            bigstring += 'Column comment : ' + self.colcomment + '\n'

        # return the result
        return bigstring

    def collheader(self,n,commentesc):
        '''
        returns string as used by column information in the header

        @param n: Column number, 0 is first column
        @type n: int
        @param commentesc: the currently used escapesequence for comments
        @type  commentesc: string
        @return: the full line of column definition to append to the header
        @rtype: string
        '''
        outstring = commentesc + ' ' + str(n+1)
        if n>8:
            outstring += ' '
        else:
            outstring += '  '
        outstring += self.colname
        if self.colcomment:
            outstring += '  '+self.colcomment
        if self.unit:
             outstring += '  ['+self.unit +']'
        outstring += '\n'
        return outstring
