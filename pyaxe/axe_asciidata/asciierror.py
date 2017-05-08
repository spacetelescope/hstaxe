
import string

class AsciiDataError(Exception):
    """
    A general exception class for the AsciiData object

    This class is the parent class for all specific
    Errors branched off it.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ColumnError(AsciiDataError):
    """
    Exception if a column does not exist
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ColTypeError(AsciiDataError):
    """
    Exception if a column type is not valid
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class TypeTransError(AsciiDataError):
    """
    Exception if a column type is not valid
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
