"""
Table element classes for input

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: mkuemmel $
$LastChangedDate: 2008-07-03 10:27:47 +0200 (Thu, 03 Jul 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciielement.py $
"""
__version__ = "Version 1.0 $LastChangedRevision: 503 $"

import string, sys, os, types
from asciierror import *

class Element(object):
    """
    Class to analyze a data element. The data element is
    given as a string, and the class methods find the type
    of the element. The type of the element is different to
    string if it can be transformed to this type
    (such as e.g. "1.0" --> 1.0).
    """
    def __init__(self, item):
        """
        Constructor for the Element class.

        @param item: the element to be analyzed
        @type item: string/integer/float
        """
        self._item   = item
        self._type   = self._find_type(item)

    def get_value(self):
        """
        Returns the elements value as its type.

        @return: the element as string
        @rtype: string/integer/float
        """
        return self._item

    def get_type(self):
        """
        Returns the element type.
        
        @return: the element type
        @rtype: <types> name
        """
        return self._type

    def _find_type(self, item):
        """
        Finds the proper type of an element.

        @param item: the item to analyze
        @type item: string/integer/float

        @return: the element type
        @rtype: <types> name
        """
        # check for int
        if self._isint(item):
            return types.IntType
        # check for float
        elif self._isfloat(item):
            return types.FloatType
        # the default is string
        else:
            return types.StringType        

    def _isint(self, item):
        """
        Checks whether an element is of type integer.

        @param item: the element to check
        @type item: any type

        @return: 1/0
        @rtype: integer
        """
        try:
            int(item)
        except:
            return 0
        return 1

    def _isfloat(self, item):
        """
        Checks whether an element is of type float.

        @param item: the element to check
        @type item: any type

        @return: 1/0
        @rtype: float
        """
        try:
            float(item)
        except:
            return 0
        return 1

class ValElement(Element):
    """
    Derived class from the Element class. In addition
    this class fills attributes with the element value
    in its proper type.
    """
    def __init__(self, item):
        """
        Constructor for the ValElement class.

        @param item: the element to be analyzed
        @type item: string/integer/float
        """
        # check whether it is a string
        if isinstance(item, type("a")):
            # if yes, initialize it in the super class
            super(ValElement, self).__init__(item)
            
            # get the typed value
            self._tvalue = self._get_tvalue(item, self.get_type())
        else:
            # if no string, determine the type.
            # store the typed value and the
            # value as string
            self._item   = str(item)
            self._type   = type(item)
            self._tvalue = item

    def get_tvalue(self):
        """
        Returns the elements value as its type.

        @return: the element value
        @rtype: string/integer/float
        """
        return self._tvalue

    def set_tvalue(self, tvalue):
        """
        Sets the typed value

        @param tvalue: the element to transform
        @type tvalue: string
        """
        self._tvalue = tvalue

    def _get_tvalue(self, item, type):
        """
        Transforms and returns the typed value.
        
        For a string element with a type different from
        string, the string is transformed into this type
        (e.g. "  1", int ----> 1).

        @param item: the element to transform
        @type item: string
        @param type: the type to transform into
        @type type: <types> name        

        @return: the typed element value
        @rtype: string/integer/float
        """
        if type == types.IntType:
            return int(item)
        elif type == types.FloatType:
            return float(item)
        else:
            return self._item

class ForElement(ValElement):
    """
    Derived class from the ValElement class. In addition
    this class fills attributes with the proper format of
    the element.
    """
    def __init__(self, item):
        """
        Constructor for the ForElement class.

        @param item: the element to be analyzed
        @type item: string/integer/float
        """
        # if yes, initialize it in the super class
        super(ForElement, self).__init__(item)

        # check whether it is a string
        if isinstance(item, type("a")):
            
            # get the format
            self._fvalue = self._get_fvalue(self.get_type())
        else:
            # get the default format
            self._fvalue = self._get_fdefaults(self.get_type())

    def get_fvalue(self):
        """
        Returns the element format.

        @return: the element value
        @rtype: string/integer/float
        """
        return self._fvalue

    def _get_fvalue(self, type):
        """
        Determines and returns the element format.

        The proper format for the element is derived from
        it string representation. This string representation
        originates directly from the input data.
        
        @param type: the type to transform into
        @type type: <types> name        

        @return: the format string
        @rtype: [string]
        """
        # check for te data type
        if type == types.IntType:

            # get the length of the stripped string version
            svalue = string.strip(self._item)
            flength = len(svalue)

            # correct for a sign in the stripped string
            if self._tvalue < 0 or svalue[0]=='+':
                flength -= 1

            # there should always be five digits
            if flength < 5:
                # return the minimum format
                return ['%5i','%5s']
            else:
                # return the format
                return ['% '+str(flength)+'i','%'+str(flength+1)+'s']

        elif type == types.FloatType:
            # store the stripped string
            svalue = string.strip(self._item)

            # check for an exponent
            epos = string.find(svalue, 'E')
            if epos < 0:
                epos = string.find(svalue, 'e')

            # get the floating point format
            if epos > -1:

                # compute the accuracy, to say the number
                # of digits after '.', taking into account
                # a possible sign
                if self._tvalue < 0.0 or svalue[0] == '+':
                    accuracy = epos-3
                else:
                    accuracy = epos-2

                # check whether there is a '.'
                if string.find(svalue, '.') < 0:
                    # correct for missing dot
                    accuracy += 1
                
                # just for security:
                if accuracy < 0:
                    accuracy = 0
                    
                # compute the total length
                tlength = accuracy+6

                # return the format
                return ['% '+str(tlength)+'.'+str(accuracy)+'e', \
                        '%'+str(tlength+1)+'s']

            # get the fixed point format
            else:

                # find the position of the '.' and the total length
                dpos = string.find(svalue, '.')
                tlength = len(svalue)
                
                # compute the accuracy, to say the number
                # of digits after '.'
                accuracy = tlength-dpos-1

                # correct the length for possible signs
                if self._tvalue < 0.0 or svalue[0] == '+':
                    tlength -=1

                # return the format
                return ['% '+str(tlength)+'.'+str(accuracy)+'f', \
                        '%'+str(tlength+1)+'s']
                       
        else:
            # default format for strings
            flength = str(len(self._item))
            return ['% '+flength+'s', '%'+flength+'s']

    def _get_fdefaults(self, type):
        """
        Determines and returns the default format
        
        @param type: the type to find the format for
        @type type: <types> name        

        @return: the list of format strings
        @rtype: [string]
        """
        if type == types.IntType:
            # default format for integers
            return ['%5i','%5s']
        elif type == types.FloatType:
            # default format for floats
            return ['% 12.6e', '%13s']
        else:
            # default format for strings
            flength = str(len(self._item))
            return ['% '+flength+'s', '%'+flength+'s']

class TypeTransformator(object):
    """
    The class contains all rules about the transformation
    of the different possible column types. It determines
    whether a transformation is possible or not. It also
    performs the transformation on elements.
    """
    def __init__(self, orig_type, new_type):
        """
        Constructor for the  TypeTransformator class.

        @param orig_type: the element to be analyzed
        @type orig_type: <types>-name
        @param new_type: the element to be analyzed
        @type new_type: <types>-name
        """
        self.istransf = self._analyze_types(orig_type, new_type)

        if self.istransf:
            self.higher_type = orig_type
        else:
            self._check_type(new_type)
            self.higher_type = new_type

    def _analyze_types(self, orig_type, new_type):
        """
        Analyzes two types for transformability

        The method analyzes whether a new type can be
        transformed into the original type. An integer
        is returned which gives the result.
        
        @param orig_type: the element to be analyzed
        @type orig_type: <types>-name
        @param new_type: the element to be analyzed
        @type new_type: <types>-name

        @return: booleans to show transformability
        @rtype: integer
        """

        # initialize the return value
        istransf = 0

        # check whether the original type is string
        if orig_type == types.StringType:
            # everything can be transformed to a string
            istransf = 1

        # check whether the original type is float
        elif orig_type == types.FloatType:
            # integers can be transformed to float
            if new_type == types.IntType:
                istransf = 1

        # check whether the original type is integer
        elif orig_type == types.IntType:
            # NOTHING can be transformed to an integer 
            istransf = 0

        else:
            raise ColTypeError('Column type: "'+str(orig_type)+'" is not a valid column type!')

        return istransf

    def _check_type(self, in_type):
        """
        Checks the validity of a type

        Only a limited number of types are admited. The method
        checks whether a certain type is valid or not.
        An exception is thrown for invalid types.
        
        @param in_type: the element to be checked
        @type in_type: <types>-name
        """

        # check whether the type is not of any valid type
        if in_type != types.StringType \
           and in_type != types.FloatType \
           and in_type != types.IntType:

            # raise an exception for invalid types
            raise ColTypeError('Column type: "'+str(in_type)+'" is not a valid column type!')

    def to_higher_type(self, tvalue):
        """
        Transforms an element to a higher type

        @param tvalue: the element to be analyzedtype to check
        @type tvalue: a value of any accepted type

        @return: the value tranformed to a higher type
        @rtype: value of any accepted type
        """
        # check if higher type is string
        if self.higher_type == types.StringType:
            # transform to string
            try:
                rvalue = str(tvalue)
            except:
                raise TypeTransError('Element: "' + str(tvalue) + '" can not be transformed to ' + str(types.StringType) + '!')
            
        # check if higher type is float
        elif self.higher_type == types.FloatType:
            # transform to float
            try:
                rvalue = float(tvalue)
            except:
                raise TypeTransError('Element: "' + str(tvalue) + '" can not be transformed to ' + str(types.FloatType) + '!')

        # check if higher type is integer
        elif self.higher_type == types.IntType:
            # transform to integer
            try:
                rvalue = int(tvalue)
            except:
                raise TypeTransError('Element: "' + str(tvalue) + '" can not be transformed to ' + str(types.IntType) + '!')

        # it should never come to this
        else:
            raise ColTypeError('Column type: "'+str(in_type)+'" is not a valid column type!')

        return rvalue

