"""
Various header classes to be part of the asciidata class

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-07-03 10:27:47 +0200 (Thu, 03 Jul 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciiheader.py $
"""
__version__ = "Version 1.1 $LastChangedRevision: 503 $"

import string
import re

from asciiutils  import *

class Header(object):
    """
    The header object

    This object offers the possibility to store additional
    information such as change comments or column information.
    This additional information may just be present at the
    beginning of the data file or later be added.
    """
    def __init__(self, filename=None, comment_char=None):
        """
        Constructor for the Header class

        @param filename: the data file
        @type filename: string
        @param comment_char: the comment_char string
        @type comment_char: string
        """
        # store the comment_char
        self._comment_char = comment_char

        # Fullhdata contains the full nonparsed header - probably
        # superflupus now
        self.Fullhdata = []

        # CollInfo is a list of column names extracted from the header
        # please note that only is current at readintime and is currently
        # not updated when columns are changed
        self.CollInfo = []

        # SexVectorColls are the known sextractor output parameters which
        # come as vectors
        self.SexVectorColls = ('MAG_APER','MAGERR_APER','FLUX_RADIUS','FLUX_APER','FLUXERR_APER','VECTOR_SOMFIT','VECTOR_ASSOC','FLUX_GROWTH','VIGNET','VIGNET_SHIFT')

        # SExtarctorFlag marks whether sextractorlike header information
        # was parsed
        self.SExtractorFlag = False

        # retrieve the comment from the data file
        # hdata is the header minus the column info lines
        # in case the header column info is invalid at loading hdata defaults to Fullhdata
        if filename == None:
            self.hdata = []
        else:
            self.hdata = self._load_header(filename, comment_char)

        # set the number of elements
        self._nentry = len(self.hdata)


    def __getitem__(self, index):
        """
        Defines the list operator for indexing

        The method returns the indexed header entry,
        if it exists. An error is raised otherwise

        @param index: the index of the header entry to be returned
        @type index: integer

        @return: a header line
        @rtype: string
        """
        if index+1 > self._nentry:
            err_msg = 'Index: '+str(index)+' does not exist! The header contains '\
             + str(self._nentry) + ' items!'
            raise Exception(err_msg)

        # return the desired header entry
        return self.hdata[index]


    def __setitem__(self, index, hentry):
        """
        Defines the list operator for indexed assignement

        @param element: either column index or name
        @type element: string/integer
        @param column: the column to assign to an index
        @type column: AsciiColumn
        """

        # check whether the target index exists;
        # raise error if not
        if index+1 > self._nentry:
            err_msg = 'Index: '+str(index)+' does not exist! The header contains '\
             + str(self._nentry) + ' items!'
            raise Exception(err_msg)

        # split the string to lines
        hitems = string.split(string.strip(hentry),'\n')

        # check whether more than one line
        # wants to be added
        if len(hitems) > 1:
            raise Exception('Only one line can be set!')

        # replace the header entry,
        # add a newline if necessary
        if hentry[-1] != '\n':
            self.hdata[index] = hentry + '\n'
        else:
            self.hdata[index] = hentry


    def __delitem__(self, index):
        """
        Deletes an index.

        @param index: the index of the header item to be deleted
        @type index: integer
        """
        # check whether the target index exists;
        # raise error if not
        if index+1 > self._nentry:
            err_msg = 'Index: '+str(index)+' does not exist! The header contains '\
             + str(self._nentry) + ' items!'
            raise Exception(err_msg)

        # delete the column
        del self.hdata[index]

        # adjust the number of entries
        self._nentry -= 1


    def __str__(self):
        """
        Defines a string method for the object

        @return: the string representation
        @rtype: string
        """

        # start the string
        hstring = ''

        # add the different items
        for line in self.hdata:
            if len(line) > 0:
                hstring += self._comment_char + line
            else:
                hstring += self._comment_char + '\n'

        # return the string
        return hstring

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
        The length operator

        @param length: the length of the instance
        @type length: integer
        """
        # thats rather trivial
        length = self._nentry

        # return the length
        return length


    def append(self, hlist):
        """
        Append something to the header data

        @param hlist: the string to append
        @type hlist: string
        """
        # split the string to lines
        hitems = string.split(hlist,'\n')

        # for each line
        for item in hitems:
            # append the new content
            # to the header content
            self.hdata.append(item+'\n')
            self._nentry += 1



    def _load_header(self, filename, comment_char):
        """
        Loads the header from the data file

        @param filename: the data file
        @type filename: string
        @param comment_char: the comment_char string
        @type comment_char: string
        """

        # start the item list
        data = []
        lastcoll,currcoll =0,0
        lastname =''
        # Define patterns for some common header formats
        commentpattern = re.compile(comment_char)
        sextractor_header = re.compile('^#\s*(\d+)\s+([+*-/()\w]+)([^\[]*)(\[\w+\])?(.*)\n')
        # open the data file and go over its rows
        for line in file(filename, 'r'):
            if commentpattern.match(line):
                #append everything after the comment_char separator to Fullhdata
                line_with_comment_char_stripped_off = commentpattern.sub('',line,count=1)
                self.Fullhdata.append(line_with_comment_char_stripped_off)
                SEmatch = sextractor_header.match(line)
                if SEmatch:  #sextractor_header.match(line):
                    # seems we have a SExtractorheader
                    if not self.SExtractorFlag:
                        self.SExtractorFlag = True
                    groups = SEmatch.groups()
                    currcoll = int(groups[0])
                    name = groups[1]
                    if currcoll <= lastcoll:
                        #ignore multiple and definitions out of order
                        continue
                    if currcoll >  (lastcoll +1):
#                        print currcoll,lastcoll
                        # we jumped some lines, pad CollInfo
                        vcounter = 1
                        while (lastcoll +1) < currcoll:
                            if lastname in self.SexVectorColls:
                                self.CollInfo.append({'NAME':lastname+str(vcounter)})
                                vcounter +=1
                            else:
                                self.CollInfo.append(None)
                            lastcoll +=1
                    self.CollInfo.append({'NAME':name})
                    lastcoll =  currcoll
                    lastname = name
                    if groups[3]:
                        # a unit was extracted
                        self.CollInfo[-1]['UNIT'] = str(groups[3].strip('[]'))
                    if groups[2] or groups[4]:
                            self.CollInfo[-1]['COMMENT'] =''
                            self.CollInfo[-1]['COMMENT'] += groups[2].strip()
                            if groups[2] and groups[4]:
                                self.CollInfo[-1]['COMMENT'] += ' '
                            self.CollInfo[-1]['COMMENT'] += groups[4].strip()
                else:
                    data.append(line_with_comment_char_stripped_off)
            else:
                # leave the file at the first
                # non-comment line
                break
        return data


    def reset(self):
        """
        Reset the header
        """
        self.hdata   = []
        self._nentry = 0

    def set_comment_char(self, comment_char):
        """
        Set the comment_char string

        @param comment_char: the new comment_char string
        @type comment_char: string
        """
        self._comment_char = comment_char


    def getCollInfo(self,index):
        """
        Robustly return column info from header
        returns (columnname,unit,comment)

        @param index: The column index
        @type index: int
        """
        #default values
        name = 'column' + str(index+1)
        unit = None
        comment = None
        if index < len(self.CollInfo):
            if self.CollInfo[index]:
                if 'NAME' in self.CollInfo[index]:
                    name =  str(self.CollInfo[index]['NAME'])
                if 'UNIT' in self.CollInfo[index]:
                    unit = str(self.CollInfo[index]['UNIT'])
                if 'COMMENT' in self.CollInfo[index]:
                    comment = str(self.CollInfo[index]['COMMENT'])
        else:
            # is the very last column in the list a known vector?
            if self.CollInfo[-1]['NAME'] in self.SexVectorColls:
                name = self.CollInfo[-1]['NAME']+str(index-len(self.CollInfo)+1)

        # return name, unit, comment of the column
        return name, unit, comment
