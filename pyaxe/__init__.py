"""
import the toplevel functions and check the directory setup
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from . import config
from .axesrc import axetasks

config.set_defaults()
config.check_axe_dirs()
