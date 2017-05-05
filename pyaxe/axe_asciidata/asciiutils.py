"""
Unspecific helper classes

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-07-03 10:27:47 +0200 (Thu, 03 Jul 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciiutils.py $
"""
__version__ = "Version 1.0 $LastChangedRevision: 503 $"

import string, sys, os, types

class NicePrinter(object):
    """
    Class to print to I/O-streams

    The class is a wrapper around an I/O stream. It offers
    methods to format strings and print to a given I/O stream.
    Linend, delimiter and linestarts are attributes of the
    class and allow a nice formatting of the print.
    """
    def __init__(self, stream=None, delimiter=None, linestart=None, linend=None):
        """
        Initializes the class
    
        A simple initializer. Most of the class attributes
        are given as parameters
    
        @param stream: I/O stream to write to
        @type stream: I/O stream
        @param delimiter: optional delimiter 
        @type delimiter: string
        @param linend: optional linenend
        @type linend: string
        """
        #set the stream
        self._stream = stream

        # set a start value
        self._start  = ''
        
        # set the delimiter
        if delimiter != None:
            #       self._delimiter = ' '+delimiter+' '
            self._delimiter = delimiter
        else:
            self._delimiter = ' '

        # set the linend
        if linend != None:
            self._linend = linend
        else:
            self._linend = '\n'

        # set the linestart
        if linestart != None:
            self._start = linestart
        else:
            self._linestart = ''


    def print_string(self, hstring):
        """
        Prints a string to the stream

        This general method prints any string
        to stream.
        
        @param hstring: the header to print 
        @type hstring: string
        """
        # that's easy up to now
        self._stream.write(hstring)

    def print_list(self, strlist):
        """
        Prints a list to the stream.

            The method combines a string list from the input
        to a string which represents a line. Delimiter,
            linend and linestart are taken into account.
        The lines is directly sent to the I/O stream.

        @param strlist: list 
        @type strlist: [string]         
            """
        self._stream.write(self._start
                           + self._delimiter.join(strlist) + self._linend)


class Separator(object):
    """
    Class to separate an ascii line into items

    Instance of this class split an ascii line into
    the different items. The methods on how to split
    a line work with a delimiter, or according to
    whitespace or according to a fixed format given
    in a file (not yet implemented.
    """
    def __init__(self, delimiter=None, file=None):
        """
        The class constructor
        """
        self._delimiter = delimiter
        self._file      = file

    def separate(self, line):
        """
        Separates a line into its items

        @param line: the ascii line to be separated
        @type line: string
        
        @return: the list of items
        @rtype: [string]
        """
        # delete the trailing newline
        if line[-1] == '\n':
            line = line[:len(line)-1]

        # separate either along a delimiter
        if self._delimiter != None:
            items = self.separate_delim(line)
        # or along whitespaces
        else:
            items = self.separate_white(line)

        return items

    def separate_white(self, line):
        """
        Separates a line along the whitespace

        The method transforms a line into the list
        of its space-separated items. The first space
        is the delimiter, any further spaces are interpreted
        to belong to the item and are preserved.
        This is advantageous to keep the item length for
        string columns with leading spaces.
        
        @param line: the ascii line to be separated
        @type line: string

        @return: the list of items
        @rtype: [string]
        """
        # create the item list
        witems = []

        # split it conventionally
        items = string.split(string.strip(line))

        # go again over the line and identify
        # the exact starting position of each
        # item, preserving the leading spaces
        start=0
        for item in items:
            pos = line.find(item,start)
            if pos > -1:
                witems.append(line[start:pos+len(item)])
                start = pos+len(item)+1

        # return the list
        return witems

    def separate_delim(self, line):
        """
        Separates a line along a delimiter

        The method transforms a line into the list
        of its delimiter separated items.
        
        @param line: the ascii line to be separated
        @type line: string

        @return: the list of items
        @rtype: [string]
        """
        # split the line
        items = string.split(line, self._delimiter)

        # return the list
        return items


class AsciiLenGetIter(object):
    """
    A general purpose iteratorfor any class with len() and get[]
    """
    def __init__(self, len_get_object):
        """
        The class constructor
        """
        # store the associated AsciiData object
        self._len_get_object = len_get_object

        # set the index of the actual column
        self._index = -1

        # set the maximum column index
        self._max_index = len(self._len_get_object) - 1

    def _iter(self):
        """
        Mandatory method for an iterator class
        """
        return self

    def next(self):
        """
        Mandatory method for an iterator class
        
        The method gives the next object in the iterator sequence.
        In case that a next object does no longer exist,
        a corresponding exception is thrown to indicate
        the end of the iterator sequence.
        """
        # check whether the next iteration does exist
        if self._index >= self._max_index:
            # no next iteration, raise exception
            raise StopIteration

        # enhance the actual index
        self._index += 1

        # return the next iteration
        return self._len_get_object[self._index]


class AsciiColumnIter(object):
    """
    An iterator class for the AsciiData class
    """
    def __init__(self, ascii_column):
        """
        The class constructor
        """
        # store the associated AsciiColumn object
        self.ascii_column = ascii_column

        # set the index of the actual row
        self._row_index = -1

        # set the maximum column index
        self._max_index = ascii_column.get_nrows() - 1

    def _iter(self):
        """
        Mandatory method for an iterator class
        """
        return self

    def next(self):
        """
        Mandatory method for an iterator class
        
        The method gives the next object in the iterator sequence.
        In case that a next object does no longer exist,
        a corresponding exception is thrown to indicate
        the end of the iterator sequence.
        """
        # check whether the next iteration does exist
        if self._row_index >= self._max_index:
            # no next iteration, raise exception
            raise StopIteration

        # enhance the actual column index
        self._row_index += 1

        # return the next iteration
        return self.ascii_column[self._row_index]
