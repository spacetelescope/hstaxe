"""
All functions for the ascidata module  

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: jhaase $ 
$LastChangedDate: 2007-07-11 08:45:32Z $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciifunction.py $
"""
__version__ = "Version 1.1 $LastChangedRevision: 234 $"

from asciidata import *

def open(filename, null=None, delimiter=None, comment_char=None):
    """
    Constructor for the AsciiData class
    
    @param filename: the filename to create the AsciiData from
    @type filename: string
    @param null: string to be interpretet as NULL
    @type null: string
    @param delimiter: string to be used as delimiter
    @type delimiter: string
    @param comment_char: string to be used as comment_char
    @type comment_char: string

    @return: the created AsciiData instance
    @rtype: AsciiData
    """
    return AsciiData(filename=filename, null=null, delimiter=delimiter,\
                     comment_char=comment_char)


def create(ncols, nrows, null=None, delimiter=None):
    """
    Constructor for the empty AsciiData class
    
    @param ncols: number of columns to be created
    @type ncols: integer
    @param nrows: number of columns to be created
    @type nrows: integer
    @param null: string to be interpretet as NULL
    @type null: string
    @param delimiter: string to be used as delimiter
    @type delimiter: string

    @return: the created AsciiData instance
    @rtype: AsciiData
    """
    return AsciiData(ncols=ncols, nrows=nrows, null=null, delimiter=delimiter,
                     columnInfo=0, headerComment=1)


def createSEx(ncols, nrows, null=None, delimiter=None):
    """
    Constructor for the empty class in the SExtractor catalogue style
    
    @param ncols: number of columns to be created
    @type ncols: integer
    @param nrows: number of columns to be created
    @type nrows: integer
    @param null: string to be interpretet as NULL
    @type null: string
    @param delimiter: string to be used as delimiter
    @type delimiter: string

    @return: the created AsciiData instance
    @rtype: AsciiData
    """
    return AsciiData(ncols=ncols, nrows=nrows, null=null, delimiter=delimiter,
                     columnInfo=1, headerComment=1)
