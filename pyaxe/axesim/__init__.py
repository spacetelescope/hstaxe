from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .axesimtasks import *


def straighten_string(in_string):
    if in_string is not None:
        in_string = in_string.strip()
        if len(in_string) < 1:
            in_string = None
    return in_string
