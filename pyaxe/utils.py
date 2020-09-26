import sys
import logging

# To guide any import *
__all__ = ["set_logging"]


# Set up logging ability for the user
# consider making a private logging level for data retension
def set_logging(filename=None, on=True, level=logging.INFO):
    """Turn on or off logging to file or stdout.
    Parameters
    ----------
    filename: str, optional
        name of the file for logging information
    on: bool, optional
        turn logging on or off, will close the file output
        in order to turn off stdout logging, set the level
        to logging.CRITICAL
    level: logging, optional
        the logging level to be recorded or displayed
    """

    formatter = logging.Formatter('\n%(funcName)s \n%(message)s')
    root = logging.getLogger(__name__)
    root.setLevel(level)
    stream_attached = False
    file_attached = False

    if on:
        #  Try to avoid adding duplicate file handlers
        if len(root.handlers) > 0:
            for handler in root.handlers:
                if isinstance(handler, logging.StreamHandler):
                    stream_attached = True
                    handler.setLevel(level)
                if isinstance(handler, logging.FileHandler):
                    file_attached = True
                    raise ValueError("File for logging already specified,\
                                      turn off logging first.")
        else:
            # to prevent warning in unhandled contexts and messages to stderr
            root.addHandler(logging.NullHandler())

        if isinstance(filename, str) and not file_attached:
            file_handler = logging.FileHandler(filename=filename,
                                               mode='a',
                                               delay=True)
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(formatter)
            print("Saving hstaxe logging to {0:s}".format(repr(filename)))
            root.addHandler(file_handler)

        if not stream_attached:
            # set the stdout stream handler
            stdout_handler = logging.StreamHandler(stream=sys.stdout)
            stdout_handler.setLevel(logging.INFO)
            # stdout_handler.setFormatter(formatter)
            root.addHandler(stdout_handler)

    #  turning the logging off to the file and set level on stream handler
    else:
        for handler in root.handlers:
            # close the file logger
            if isinstance(handler, logging.FileHandler):
                handler.close()
                root.removeHandler(handler)
            # set stream logging level
            if isinstance(handler, logging.StreamHandler):
                handler.setLevel(level)
    return root
