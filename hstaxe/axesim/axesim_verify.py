from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
from ..axeerror import aXeSIMError

# excutables needed for axesim
_execfiles = ['aXe_AF2PET', 'aXe_GOL2AF', 'aXe_PET2SPC',
              'aXe_PETCONT', 'aXe_STAMPS',
              'aXe_DIRIMAGE']

_bin_dir = os.environ('AXESIMBIN')


# From PyDrizzle.fileutil...
def findFile(input):
    """ Search a directory for full filename with optional path. """

    _fdir, _fname = os.path.split(input)

    if _fdir == '':
        _fdir = os.curdir

    flist = os.listdir(_fdir)

    found = False
    for name in flist:
        if not name.find(_fname):
            found = True

    return found


#
# Check to see that all the executables are present
#
for binary in _execfiles:
    if not findFile(_bin_dir + binary):
        raise aXeSIMError("Missing binary: {0:s}".format(_bin_dir + binary))
