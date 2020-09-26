from __future__ import (absolute_import, unicode_literals, division,
                        print_function)


class axeException(Exception):
    """Exception class for the axe package.

    This class is just a simple extension to the general exception class.
    All errors specific to 'axe' are thrown using this class to be able
    to distinguish them from different errors.

    A message is provided by the code that raised the exception
    """
    def __init__(self, message='Unspecified aXe exception'):
        """Initializer for the exception class.

        @param message: message associated to the exception
        @type message: string
        """
        Exception.__init__(self, message)
        self.name = "Unspecified aXe Exception"
        self.message = message

    def __str__(self):
        """
        String method for the class

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class FileError(axeException):
    """Problem using an aXe specific file"""

    def __init__(self, filename="aXe File Exception"):
        """Initializer for the class.

        @param message: message associated to the exception
        @type message: string
        """
        message = "File: {0} not found".format(filename)
        axeException.__init__(self, message)
        self.name = "aXe File Exception"

    def __str__(self):
        """String method for the class.

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class DirError(axeException):
    """Problem using an aXe specific directory"""

    def __init__(self, dirname):
        """
        Initializer for the class

        @param message: message associated to the exception
        @type message: string
        """
        message = "Dir: {0} not found".format(dirname)
        axeException.__init__(self, message)
        self.name = "aXe Directory Exception"

    def __str__(self):
        """String method for the class.

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class ParamError(axeException):
    """Problem using an aXe specific directory."""

    def __init__(self, message):
        """Initializer for the class

        @param message: message associated to the exception
        @type message: string
        """
        axeException.__init__(self, message)
        self.name = "aXe Parameter Exception"

    def __str__(self):
        """String method for the class.

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class EnvError(axeException):
    """Required axe environment variable not found"""
    def __init__(self, env_var=None):
        """
        Initializer for the class

        @param message: message associated to the exception
        @type message: string
        """
        message = "Regquired envrionment variable {0} not found".format(env_var)
        axeException.__init__(self, message)
        self.name = "aXe Environment Error"

    def __str__(self):
        """
        String method for the class

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class axesimException(axeException):
    """General Error in aXeSIM.

    This class is just a simple extension to the general exception class.
    All errors specific to 'axesim' are thrown using this class to be able
    to distinguish them from different errors.
    """
    def __init__(self, message="aXeSIM Exception"):
        """
        Initializer for the class

        @param message: message associated to the exception
        @type message: string
        """
        axeException.__init__(self, message)
        self.name = "aXeSIM Exception"

    def __str__(self):
        """
        String method for the class

        @return: the string representation of the class
        @rtype: string
        """
        return self.message
