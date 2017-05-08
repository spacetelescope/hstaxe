from __future__ import (absolute_import, unicode_literals, division,
                        print_function)


class aXeError(Exception):
    """Error class for aXe

    This class is just a simple extension to the general exception class.
    All errors specific to 'axe' are thrown using this class to be able
    to distinguish them from different errors.
    """
    def __init__(self, message, *args, **kwargs):
        """
        Parameters
        ----------
        message: str
            message associated to the exception
        """
        super(aXeError, self).__init__(self, message, *args, **kwargs)
        self.message = message

    def __str__(self):
        """String method for the class

        @return: the string representation of the class
        @rtype: string
        """
        return self.message


class aXeSIMError(Exception):
    """General Error in aXeSIM

    This class is just a simple extension to the general exception class.
    All errors specific to 'axesim' are thrown using this class to be able
    to distinguish them from different errors.
    """
    def __init__(self, message):
        """Initializer for the class

        @param message: message associated to the exception
        @type message: string
        """
        self.message = message

    def __str__(self):
        """String method for the class

        @return: the string representation of the class
        @rtype: string
        """
        return self.message
