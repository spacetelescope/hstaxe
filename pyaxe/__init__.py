"""
import the toplevel functions and check the directory setup
"""
from pkg_resources import get_distribution, DistributionNotFound

from . import config
from .axesrc import axetasks
from .utils import set_logging

config.set_defaults()
config.check_axe_dirs()


try:
    release = get_distribution('hstaxe').version
    __version__ = '.'.join(release.split('.')[:4])
except DistributionNotFound:
    # package is not installed
    __version__ = 'unknown'
    __githash__ = ''