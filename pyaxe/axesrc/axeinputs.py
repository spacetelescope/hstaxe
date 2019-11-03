import os
import numpy as np

from astropy.io import fits
from astropy.table import Table, Column
from astropy.io.registry import IORegistryError

from pyaxe.axeerror import aXeError
from pyaxe import config as config_utils


class aXeInput(object):
    """Class for the administration of all aXe input."""
    def __init__(self, inlist, configterm="", fringeterm=None):
        # load the Input Image List
        self._inimlist = Table.read(inlist, format='ascii.no_header')

        # find out what data has been read and assign column names
        self._validate_columns()

        # add default columns if necessary
        self._complete_columns(configterm, fringeterm)

        # reject subarray images
        self._check_subarray()

        # Now that all the basic columns have been Added
        # expand the Input Image List if necessary
        if ',' in self._inimlist['objcat'][0]:
            self._extend_asciidata()

    def __iter__(self):
        """Allow the class to iterate over the table of inputs."""
        return iter(self._inimlist)

    def _validate_columns(self):
        """Identify columns according to the Input Image List format"""

        columns = self._inimlist.colnames

        # check for number of columns
        if len(columns) < 2:
            err_msg = ("Expected at least 2 columns in {0:s}!"
                       .format(self._inimlist))
            raise aXeError(err_msg)

        columns.reverse()
        name = columns.pop()
        if ((np.issubsctype(self._inimlist[name], np.str)) and
                ('fits' in self._inimlist[name][0])):
            self._inimlist.rename_column(name, 'grisim')
        else:
            err_msg = ("Column 1 should be the names of the grism"
                       " fits files: {0}".format(self._inimlist))
            raise aXeError(err_msg)

        # check whether second column has type string
        name = columns.pop()
        if np.issubsctype(self._inimlist[name], np.str):
            self._inimlist.rename_column(name, 'objcat')
        else:
            err_msg = ("Column 2 should be the names of the object catalogs"
                       "{0}".format(self._inimlist))
            raise aXeError(err_msg)
        # do a check on the first object catalog to make sure it's a format
        # that we can read
        try:
            __ = Table.read(config_utils.getDATA(self._inimlist[0]['objcat']), format='ascii.sextractor')
        except IORegistryError:
            raise aXeError("Catalog format not recognized , checked for: {0:s}"
                           .format(self._inimlist[name][0]))

        # go over all remaining rows to find DMAG and DIRIM
        while columns:
            name = columns.pop()

            # assume if it's a number it's DMAG
            if np.issubsctype(self._inimlist[name], np.float):
                self._inimlist.rename_column(name, 'dmag')
            # assume if it's a string its the direct image
            elif np.issubsctype(self._inimlist.columns[0], np.str):
                if '.fits' in self._inimlist.columns[0][0]:
                    self._inimlist.rename_column(name, 'dirim')
            else:
                err_msg = ("Problem identifying last columns in: {0}"
                           .format(self._inimlist))
                raise aXeError(err_msg)

    def _complete_columns(self, configterm='', fringeterm=None):
        """Adds columns to get standard format"""

        # check whether a column for the configuration
        # files exist
        if 'config' not in [x.lower() for x in self._inimlist.colnames]:
            col = Column(name='config',
                         data=[configterm]*len(self._inimlist))

            # add the configuration file column
            self._inimlist.add_column(col)

        # files exist
        if 'dmag' not in [x.lower() for x in self._inimlist.colnames]:
            col = Column(name='dmag',
                         data=[0]*len(self._inimlist))
            self._inimlist.add_column(col)

        # check whether there are fringe configs
        if fringeterm is not None:
            col = Column(name='fringe',
                         data=[fringeterm] * len(self._inimlist))
            # add a column for the fringes
            self._inimlist.add_column(col)

    def _extend_asciidata(self):
        """Extend the table to one data set per line."""

        # get a booolean for fringe correction

        for row in self._inimlist:
            clist = row['config'].split(',')
            olist = row['objlist'].split(',')
            if (len(clist) != len(olist)):
                err_msg = ("Number of object cats in {0:s} "
                           "is different from number of "
                           "configs in {1:s}"
                           .format(row['objlist'], row['config']))
                raise aXeError(err_msg)

            if ',' in row['fringe']:
                flist = row['fringe'].split(',')
                if (len(clist) != len(flist)):
                    err_msg = ("Number of fringe configs in: {0:s} "
                               "is different from number of aXe configs "
                               "in {1:s}".format(row['fringe'], row['config']))
                    raise aXeError(err_msg)

            # remove the original row and replace it with a row for each image
            self._inimlist.delete_row(row.index)

            while (clist and olist):
                conf = clist.pop()
                objcat = olist.pop()
                new_row = [row['grisim'], row['dirim'],
                           row['dmag'], objcat, conf]
                if flist:
                    new_row.append(flist.pop())
                self._inimlist.add_row(new_row)

    def _tolist(self, inimlist):
        """Transforms the Input Image List to a list of dictionaries."""
        imagedict = []

        # check whether a column with direct images exists
        dirim = False
        if 'dirim' in [x.lower() for x in self._inimlist.colnames]:
            dirim = True

        # check whether a column with fringe info exists
        fringe = False
        if 'fringe' in [x.lower() for x in self._inimlist.colnames]:
            fringe = True

        # check whether a column with dmag info exists
        dmag = False
        if 'dmag' in [x.lower() for x in self._inimlist.colnames]:
            dmag = True

        # go over all rows in the list
        for row in self._inimlist:
            ndict = {}
            ndict['GRISIM'] = row['GRISIM']
            ndict['OBJCAT'] = row['OBJCAT']
            ndict['CONFIG'] = row['CONFIG']

            if dirim:
                ndict['dirim'] = row['dirim']
            else:
                ndict['dirim'] = None

            if fringe:
                ndict['fringe'] = row['fringe']
            else:
                ndict['fringe'] = None

            if dmag:
                ndict['dmag'] = row['dmag']
            else:
                ndict['dmag'] = 0.

            imagedict.append(ndict)

        return imagedict

    def _check_grisms(self):
        """Check the existence of the grism files."""
        # go over all rows in the list
        for row in self._inimlist:
            # check the existence of the Input Image List
            if not os.path.isfile(row['grisim']):
                err_msg = ("Grism image: {0:s} does not exist!"
                           .format(row['grisim']))
                raise aXeError(err_msg)

    def _check_subarray(self):
        """ check for and reject subarray images """
        # go over all rows in the list
        for row in self._inimlist:

            # check the existence of the Input Image List
            image = config_utils.getDATA(row['grisim'])
            subarray = fits.getval(image, "SUBARRAY", ext=0)
            if subarray:
                err_msg = ("Grism image: {0:s} is a subarray"
                           " which is not supported".format(image))
                raise aXeError(err_msg)

            if 'dirim' in self._inimlist.colnames:
                image = config_utils.getDATA(row['dirim'])
                subarray = fits.getval(image, "SUBARRAY", ext=0)
                if subarray:
                    err_msg = ("Direct image: {0:s} is a subarray"
                               " which is not supported".format(image))
                    raise aXeError(err_msg)

    def writeto(self, filename):
        """Write list to a file"""

        # check and delete an existing file
        # of that name
        if os.path.isfile(filename):
            os.unlink(filename)

        # open/create the file
        ofile = open(filename, 'w')

        # write the content to this file
        ofile.write(str(self))

        # close the file
        ofile.close()
