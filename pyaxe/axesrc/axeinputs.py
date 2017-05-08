from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .. import axe_asciidata
from ..axeerror import aXeError


class aXeInputList(object):
    """Class for the administration of all aXe input"""
    def __init__(self, inlist, configterm, fringeterm=None):
        # load the Input Image List and identify
        # all columns therein
        self._inimlist = self._identify_columns(inlist)

        # add columns if necessary
        self._complete_columns(self._inimlist, configterm, fringeterm)

        # check the existence of the grism files
        # self._check_grisms(self._inimlist)

        #check for subarray images
        self._check_subarray(self._inimlist)

        # get the multiplicity, which means the number
        # of Input Object Lists per grism image
        multiplicity = len(self._resolve_term(self._inimlist['OBJCAT'][0]))

        # expand the Input Image List if necessary
        if multiplicity > 1:
            self._inimlist = self._extend_asciidata(self._inimlist, multiplicity)

        # transform the Input Image List to a
        # dictionary list
        self.axeinputs = self._tolist(self._inimlist)


    def __str__(self):
        """String method for the class"""
        # intitialize the string
        rstring = ''

        # go over all items
        for index in range(len(self)):

            # compose the string for one item
            single = self[index]['GRISIM'] + ' ' + self[index]['OBJCAT'] + ' '

            # add the direct image
            if self[index]['DIRIM']:
                single += self[index]['DIRIM'] + ' '

            # add the direct image
            if self[index]['FRINGE']:
                single += self[index]['FRINGE'] + ' '

            # add dmag and config file
            single += str(self[index]['DMAG']) + ' ' + self[index]['CONFIG'] + '\n'

            # add the item string
            rstring += single

        # return the final string
        return rstring


    def __len__(self):
        """Defines and returns the length of a class instance"""
        return len(self.axeinputs)


    def __getitem__(self,index):
        """The index operator for the class"""
        return self.axeinputs[index]

    def __iter__(self):
        """Provide an iterator object.

        The function provides and returns an interator object
        for the AstroColumnData class. Due to this iterator object
        sequences like:
        for elem  in ascii_column_object:
            <do something with elem>
        are possible.
        """
        return LenGetIter(self)

    def _identify_columns(self, inlist):
        """Identify columns according to the Input Image List format"""
        # load the file as an AsciiData object,
        # complain in case of error
        try:
            inimlist = axe_asciidata.open(inlist)
        except:
            err_msg = 'Problems loading file: ' + inlist + ' into axe_asciidata!'
            raise aXeError(err_msg)

        # check for number of columns
        if inimlist.ncols < 2:
            err_msg = 'Not enough columns in Input Image List: ' + inlist + '!'
            raise aXeError(err_msg)

        # check whether first column has type string
        if inimlist[0].get_type() is type('a'):
            # rename the column
            inimlist[0].rename('GRISIM')
        else:
            err_msg = 'Column 1 in Input Image List: ' + inlist + ' must contain string!'
            raise aXeError(err_msg)

        # check whether second column has type string
        if inimlist[1].get_type() == type('a'):
            # rename the column
            inimlist[1].rename('OBJCAT')
        else:
            err_msg = 'Column 2 in Input Image List: ' + inlist + ' must contain string!'
            raise aXeError(err_msg)

        # go over all remaining rows
        for cindex in range(2, inimlist.ncols):

            # store the column type
            ctype = inimlist[cindex].get_type()

            # check whether column has number type
            if ctype == type(1.0) or ctype == type(1):
                # rename it as the dmag-column
                inimlist[cindex].rename('DMAG')

            # check whether column has string type
            elif ctype ==  type('a'):
                # check whether it contains the name of fits files
                if inimlist[cindex][0].find('.fits') > -1:

                    # rename the column as direct-image column
                    inimlist[cindex].rename('DIRIM')

                else:
                    # rename the column as direct-image column
                    inimlist[cindex].rename('CONFIG')


        return inimlist


    def _complete_columns(self, inimlist, configterm, fringeterm=None):
        """Adds columns to get standard format"""

        # check whether a column for the configuration
        # files exist
        if inimlist.find('CONFIG') < 0:

            # add the configuration file column
            conf_col = inimlist.append('CONFIG')

            # go over all rows
            for index in range(inimlist.nrows):

                # insert the term for configuration files
                inimlist[conf_col][index] = configterm

        # check whether a column for the configuration
        # files exist
        if inimlist.find('DMAG') < 0:

            # add the configuration file column
            dmag_col = inimlist.append('DMAG')

            # go over all rows
            for index in range(inimlist.nrows):

                # insert the term for configuration files
                inimlist[dmag_col][index] = 0.0

        # check whether there are fringe configs
        if fringeterm:

            # add a column for the fringes
            fringe_col = inimlist.append('FRINGE')

            # go over all rows
            for index in range(inimlist.nrows):

                # insert the term for configuration files
                inimlist[fringe_col][index] = fringeterm


    def _resolve_term(self, interm):
        """
        Resolves a comma separated term
        """
        # initialize the return list
        items = []

        # split the term at the commas
        terms = interm.split(',')

        # go over the list
        for term in terms:
            # add the stripped terms to the return ;ist
            items.append(term.strip())

        # return the list
        return items


    def _extend_asciidata(self, inimlist, multiplicity):
        """
        Extend the Input Image List such that there is one data set per line
        """
        # create a new, larger asciidata object
        new_inimlist = axe_asciidata.create(inimlist.ncols,
                                            multiplicity*inimlist.nrows)

        # rename the column of the new object
        for index in range(inimlist.ncols):
            new_inimlist[index].rename(inimlist[index].colname)


        # get a booolean for fringe correction
        if inimlist.find('FRINGE') > -1:
            fringe = 1
        else:
            fringe = 0

        # go over all wows of the old list
        new_index = 0
        for index in range(inimlist.nrows):

            # resolve the terms for configuration file
            # and Input Object List
            cterm = self._resolve_term(inimlist['CONFIG'][index])
            oterm = self._resolve_term(inimlist['OBJCAT'][index])

            # check the number of items
            if len(cterm) != len(oterm):
                err_msg = 'Number of object cats in: '+ inimlist['OBJCAT'][index]\
                          +' is different from number of configs in: '\
                          +inimlist['CONFIG'][index]+' !'
                raise aXeError(err_msg)

            # check for the right length
            if len(cterm) != multiplicity:
                err_msg = 'Numer of configs in: '+ inimlist['CONFIG'][index]\
                          +' must be: '+multiplicity+'!'
                raise aXeError(err_msg)

            if fringe:
                fterm = self._resolve_term(inimlist['FRINGE'][index])

                # check the number of items
                if len(cterm) != len(fterm):
                    err_msg = 'Number of fringe configs in: '+ inimlist['FRINGE'][index]\
                              +' is different from number of aXe configs in: '\
                              +inimlist['CONFIG'][index]+' !'
                    raise aXeError(err_msg)

            # go over every term
            for ii in range(multiplicity):
                new_index += ii

                # transfer the 'constant' data
                new_inimlist['GRISIM'][new_index] = inimlist['GRISIM'][index]
                new_inimlist['DIRIM'][new_index]  = inimlist['DIRIM'][index]
                new_inimlist['DMAG'][new_index]   = inimlist['DMAG'][index]

                # transfer the 'variable' data
                new_inimlist['OBJCAT'][new_index] = oterm[ii]
                new_inimlist['CONFIG'][new_index] = cterm[ii]

                if fringe:
                    new_inimlist['FRINGE'][new_index] = fterm[ii]

            # enhance the index for the new list
            new_index += 1


        # return the new list
        return new_inimlist


    def _tolist(self, inimlist):
        """
        Transforms the Input Image List to a list of dictionaries
        """

        # initialize the list
        axeinputs = []

        # check whether a column with
        # direct images exists
        if inimlist.find('DIRIM') > -1:
            dirim = 1
        else:
            dirim = 0

        # check whether a column with
        # direct images exists
        if inimlist.find('FRINGE') > -1:
            fringe = 1
        else:
            fringe = 0

        # go over all rows in the list
        for index in range(inimlist.nrows):

            # create a new, empty dictionary
            axeitem= {}

            # fill the dictionary
            axeitem['GRISIM'] = inimlist['GRISIM'][index]
            axeitem['OBJCAT'] = inimlist['OBJCAT'][index]
            axeitem['DMAG']   = inimlist['DMAG'][index]
            axeitem['CONFIG'] = inimlist['CONFIG'][index]

            # check whether there are direct images
            # add the image name or 'None'
            if dirim:
                axeitem['DIRIM'] = inimlist['DIRIM'][index]
            else:
                axeitem['DIRIM'] = None

            # check whether there are direct images
            # add the image name or 'None'
            if fringe:
                axeitem['FRINGE'] = inimlist['FRINGE'][index]
            else:
                axeitem['FRINGE'] = None

            # append the dictionary to the list
            axeinputs.append(axeitem)

        # return the dictionary list
        return axeinputs


    def _check_grisms(self, inimlist):
        """
        Check the existence of the grism files
        """
        from . import axeutils

        # go over all rows in the list
        for index in range(inimlist.nrows):

            # check the existence of the Input Image List
            if  not os.path.isfile(axeutils.getIMAGE(inimlist['GRISIM'][index])):
                err_msg = 'Grism image: ' + axeutils.getIMAGE(inimlist['GRISIM'][index])\
                          + ' does not exist!'
                raise aXeError(err_msg)

    def _check_direct(self, inimlist):
        """
        Check the existence of the direct image files
        Just look in the first row since all rows should
        have the same format
        """
        from . import axeutils
        # check the existence of the Input Image List
        if inimlist['DIRIM'][0]:
            return True
        else:
            return False



    def _check_subarray(self,inimlist):
        """ check for and reject subarray images """
        from . import axeutils
        from astropy.io import fits as pyfits
        # go over all rows in the list
        for index in range(inimlist.nrows):

            # check the existence of the Input Image List
            image=axeutils.getIMAGE(inimlist['GRISIM'][index])
            subarray=pyfits.getval(image,"SUBARRAY",ext=0)
            if subarray:
                err_msg = 'Grism image: ' + image\
                          + ' is a subarray - which is not supported'
                raise aXeError(err_msg)
            if self._check_direct(inimlist):
                image=axeutils.getIMAGE(inimlist['DIRIM'][index])
                subarray=pyfits.getval(image,"SUBARRAY",ext=0)
                if subarray:
                    err_msg = 'Direct image: ' + image\
                              + ' is a subarray - which is not supported'
                    raise aXeError(err_msg)

    def writeto(self, filename):
        """
        Write list to a file
        """

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


class LenGetIter(object):
    """
    A general purpose iterator for any class with len() and get[]
    """
    def __init__(self, len_get_object):
        """
        The class constructor
        """
        # store the associated AsciiData object
        self._len_get_object = len_get_object

        # set the index of the actual column
        self._index = -1

        # set the maximum column index
        self._max_index = len(self._len_get_object) - 1

    def _iter(self):
        """
        Mandatory method for an iterator class
        """
        return self

    def __next__(self):
        """
        Mandatory method for an iterator class

        The method gives the next object in the iterator sequence.
        In case that a next object does no longer exist,
        a corresponding exception is thrown to indicate
        the end of the iterator sequence.
        """
        # check whether the next iteration does exist
        if self._index >= self._max_index:
            # no next iteration, raise exception
            raise StopIteration

        # enhance the actual index
        self._index += 1

        # return the next iteration
        return self._len_get_object[self._index]

    def next(self):
        """
        For Python 2 Compatibility
        """
        return self.__next__()
