"""
import the toplevel functions and check the directory setup
"""
from . import config
from .axesrc import axetasks

config.set_defaults()
config.check_axe_dirs()
