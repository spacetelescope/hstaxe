import os
import math
import logging

from hstaxe import config as config_util
from hstaxe.axeerror import aXeError

# make sure there is a logger
_log = logging.getLogger(__name__)

class ConfigList:
    """Configuration File Object"""
    def __init__(self, keylist, header=None):
        """
        Initializes the ConfigList object by tranfsforming
        a list of keywords into a structured list including
        beams descriptions

        keylist: list
            List of configuration keys
        header: str
            the header string
        """
        # beam indices which might be found the file
        idents = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                  'L', 'M', 'N', 'O', 'P', 'Q']

        # create the (visible) dictionary
        self.beams = {}

        # create the hidden beam list
        self._beams = []

        # store the header
        self.header = header

        # load the general required keywords
        self.gkeys = self._find_gkeys(keylist)
        # try to load beams as long as there
        # are keywords and as long as there
        # are candidate beam numbers
        iindex = 0
        while (len(keylist) > 0 and iindex < len(idents)):
            try:
                # try to load a beam
                self._beams.append(ConfigBeam(idents[iindex], keylist))
                self.beams[idents[iindex]] = self._beams[iindex]

            except BeamNotFound:
                # no information on this beam is in the file
                pass

            # enhance the counter
            iindex += 1

        # inform about the useless keywords
        if len(keylist) > 0:
            _log.info('\nDispensable Keywords: ')
            for key in keylist:
                _log.info(key)

    def __str__(self):
        """String method for the class

        The method transforms the configuration
        file object into its string representation.

        Returns
        -------
        a string representation of the object
        """
        # take the string of the header
        rstring = str(self.header) + '\n'

        # add the strings for the global keys
        for key in self.gkeys:
            rstring += str(key)

        for beam in self._beams:
            rstring += str(beam)

        # return the total string
        return rstring

    def __delitem__(self, item):

        # find the index of the requested item
        index = self._find_gkey(item)

        # check whether the item was found
        if index > -1:
            del self.gkeys[index]

    def __getitem__(self, item):

        # find the index of the requested item
        index = self._find_gkey(item)

        # check whether the item was found
        if index > -1:
            # return the identified item
            return self.gkeys[index].keyvalue
        else:
            if item in self.beams.keys():
                return self.beams[item]
            else:
                # return NULL
                return None

    def _find_gkey(self, item):

        # set the default return value
        found = -1
        # go over all items
        for index in range(len(self.gkeys)):
            # check whether it is the right item
            if self.gkeys[index].keyword == item:
                # set the return value to the index
                found = index

        # return the result
        return found

    def _load_file(self):
        """Configuration file --> keyword list

        The method to load a configuration file and
        extract all valid keyword-keyvalue-comment information
        from it. The keyword-keyvalue pairs are
        organized and returned as a list of
        configuration key objects.

        Returns
        -------
        keylist: list
        list of ConfKey's

        Notes
        -----
        Double check reading two files here
        """
        # initialize the list
        keylist = []
        # open the file and parse through it
        with open(self.filename, 'r') as fopen:
            for line in fopen:
                # strip the line
                str_line = line.strip()

                # check whether the line contains a keyword
                if len(str_line) and str_line[0] != '#':
                    # create and append the keyword
                    keylist.append(self._key_from_line(str_line))

        # return the list
        return keylist

    def _get_gkey_index(self, keyword):
        """Retrieve the index of a global keyword

        The method searches for the index of
        a requested keyword in the list of global
        keywords. If the keyword does not exists,
        the index -1 is returned

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        index: int
            the index of the keyword
        """
        # initialize the return value
        kindex = -1

        # go over all keys
        for index in range(len(self.gkeys)):
            # check whether the current key matches
            if self.gkeys[index].keyword == keyword:
                # return it if it matches
                return index

        # return the default
        return kindex

    def _key_from_line(self, line):
        """Creates a keyword from a line

        The method extracts the konfiguration keyword,
        the associated value and, if present,
        a comment from a line in the configuration file.
        A configuration key object representing the extracted
        keyword is created and returned.

        Parameters
        ----------
        line: list
            line to analyze

        Returns
        -------
        configuration key  object
        """
        # split the line into items
        items = line.split()

        # for more than one item the
        # first item is the keyword
        if len(items) > 1:
            keyword = items[0].strip()

            # check for a comment
            cpos = line.rfind(';')
            if cpos < 0:
                # evaluate the keyvalue
                keyvalue = line[line.find(keyword)+len(keyword):].strip()
                comment = None
            else:
                # evalute keyvalue and comment
                tmp_val = line[line.find(keyword)+len(keyword):].strip()
                keyvalue = tmp_val.split(';')[0].strip()
                comment = tmp_val.split(';')[1].strip()
        else:
            # something's wrong here
            err_msg = 'Only one item in: ' + line + ' !'
            raise aXeError(err_msg)

        # create and return the keyword
        return ConfKey(keyword, keyvalue, comment)

    def _find_gkeys(self, keylist):
        """Finds and extracts the global keywords

        The method finds the all predefined global keywords in
        a keyword list. The list of global keywords is
        returned. Their counterparts in the input keyword list
        are deleted.

        Parameters
        ----------
        keylist: list
            list of keywords

        Returns
        -------
        keys: list
            global keywords
        """
        gkeywords = ['INSTRUMENT', 'CAMERA', 'TELAREA',
                     'SCIENCE_EXT', 'ERRORS_EXT',
                     'DQ_EXT', 'OPTKEY1', 'OPTVAL1', 'FFNAME', 'DQMASK',
                     'DRZRESOLA', 'DRZSCALE', 'DRZLAMB0', 'DRZXINI',
                     'DRZROOT', 'EXPTIME', 'WEIGHT_EXT', 'DRZPFRAC',
                     'DRZPSCALE', 'DRZKERNEL', 'MODEL_EXT', 'VARIANCE_EXT',
                     'RDNOISE', 'PSFCOEFFS', 'PSFRANGE', 'IPIXFUNCTION',
                     'POBJSIZE', 'SMFACTOR']

        # initialize the global keylist
        # and the list with indices to be deleted
        gkeys = []
        dindex = []

        # go over the keylist read in,
        # keeping and index variable
        iindex = 0
        for key in keylist:

            # identify the current keyword in the
            # list of possible ones
            if key.keyword in gkeywords:

                # store the index
                dindex.append(iindex)

                # create and append the new keyword
                gkeys.append(ConfKey(key.keyword, key.keyvalue, key.comment))

            iindex += 1

        # delete the input keywords which
        # have been 'used'
        dindex.sort()
        dindex.reverse()
        for index in dindex:
            del keylist[index]

        # return the list of global keys
        return gkeys

    def _check_gfiles(self):
        """Checks whether all files exist

        The method checks whether the files whose names
        are within the class data do exist or not.
        An error is reported in case that the files
        do not exist.
        """
        # list of the root of all
        # global keys indicating a file
        fkeys = ['FFNAME']

        # go over all file keywords
        for key in fkeys:

            # identify the keyword in the list
            index = self._get_gkey_index(key)

            # check for existence
            if index > -1:

                # extract the keyvalue
                kvalue = self.gkeys[index].keyvalue

                # if the keyvalue is NOT None but the file does not exist
                if ((kvalue.upper() is not 'NONE') and
                    (not os.path.isfile(config_util.getCONF(kvalue)))):
                    # report an error
                    err_msg = ("The file: {0:s} does not exist!"
                               .format(config_util.getCONF(kvalue)))
                    raise aXeError(err_msg)

    def get_gkey(self, keyword):
        """Retrieve a requested global keyword

        The method searches the list of global keywords
        for a fitting keyword. In case that the requested
        keyword exists, it is returned.
        If not 'None' is returned

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        key: str or None
            the requested keyword or 'None'
        """
        # initialize the return value
        rkey = None

        # search for the index in the keyword list
        index = self._get_gkey_index(keyword)

        # check whether the keyword exists
        if index > -1:
            # return the keyword
            return self.gkeys[index]
        else:
            # return the default
            return rkey

    def add_gkey(self, keyword, keyvalue, comment=None):
        """Add global keyword

        The method adds a keyword to the list of global
        keywords. In case that the keyword just exists,
        it is overwritten, otherwise it is appended
        to the global keyword list.

        Parameters
        ----------
        keyword: str
            name of the requested keyword
        keyvalue: any
            value of the requested keyword
        comment: str
            comment for the keyword
        """

        # search for the index in the keyword list
        index = self._get_gkey_index(keyword)

        if index > -1:
            # if it matches, copy the data
            self.gkeys[index].keyvalue = keyvalue
            self.gkeys[index].comment = comment
        else:
            # the keyword does not yet exist, just create and add it
            self.gkeys.append(ConfKey(keyword, keyvalue, comment))

    # def drizzle_check(self):
    #     """Check for drizzle keywords

    #     The method assures that all necessary drizzle keywords
    #     are present. Nonexisting keywords are added with default
    #     values. Finally the value for the drizzle kernel is checked
    #     against all valid values.

    #     Returns
    #     -------
    #     bool:  True if the drizzle kernel is valid
    #     """
    #     # list with all valid kernels
    #     kernels = ['square', 'point', 'turbo', 'gaussian', 'tophat',
    #                'lanczos2', 'lanczos3']

    #     # make sure that some important drizzle keywords are there

    #     pself = self.setdefault('DRZPSCALE', 1.0)
    #     pfrac = self.setdefault('DRZPFRAC', 1.0)
    #     dkernel = self.setdefault('DRZKERNEL', 'square')
    #     droot = self.setdefault('DRZROOT', 'aXedrizzle')

    #     # check for valid drizzle kernel
    #     if dkernel not in kernels:
    #         return False
    #     return True

    # def setdefault(self, keyword, keyvalue, comment=None):
    #     """Add global keyword

    #     The method mimics the setdefault method for dictionary
    #     objects. A keyword is added with the given value and
    #     comment, but only in case that it does not yet exist.
    #     If it exists, nothing is done

    #     Parameters
    #     ----------
    #     keyword: str
    #         name of the requested keyword
    #     keyvalue: any
    #         value of the requested keyword
    #     comment: str
    #         comment for the keyword

    #     Returns
    #     -------
    #     The keyword value
    #     """

    #     # search for the index in the keyword list
    #     index = self._get_gkey_index(keyword)

    #     if index < 0:
    #         # the keyword does not yet exist, just create and add it
    #         self.gkeys.append(ConfKey(keyword, keyvalue, comment))
    #         # extract the keyvalue
    #         value = self.gkeys[-1].keyvalue
    #     else:
    #         # extract the keyvalue
    #         value = self.gkeys[index].keyvalue

    #     # return the keyvalue
    #     return value

    def get_gvalue(self, keyword):
        """Retrieve a requested global keyword value

        The method returns the value of the keyword
        which matches the requested value.
        If there is no matching keyword, 'None'
        is returned.

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        The keyword value
        """

        # set the default return value
        rvalue = None

        # search for the keyword
        key = self.get_gkey(keyword)

        # check whether it is non-NULL
        if key:
            # extract the value
            rvalue = key.keyvalue

        # return the value
        return rvalue

    def writeto(self, filename):
        """Save the object to a file

        The method saves the object to a file
        with name specified in the input.

        Parameters
        ----------
        filename: str
            name of the file
        """
        # destroy the old file
        if os.path.isfile(filename):
            os.unlink(filename)

        # open the new file
        ofile = open(filename, 'w')

        # write the string to the file
        ofile.write(str(self))

        # close the file
        ofile.close()

    def flush(self):
        """Save the object back to file

        The method saves the object back to a file
        with the identical filename it was read from.
        """
        # just use the more general method
        self.writeto(self.filename)

    def check_files(self, check_glob=True):
        """Checks whether all files exist

        The method checks whether the files whose names
        are within the class data do exist or not.
        An error is reported in case that the files
        do not exist.
        """
        n_sens = 0

        # check global files if desired
        if check_glob:
            self._check_gfiles()

        # create the (visible) dictionary
        for bkey in self.beams.keys():
            n_sens += self.beams[bkey].check_files()

        # return the number
        # of existing sensitivity files
        return n_sens


class ConfigFile(ConfigList):
    """Configuration File Object"""
    def __init__(self, filename=None):
        """
        Initializes the ConfigFile object either
        by reading in a configuration file
        or by creating a default configuration file

        Parameters
        ----------
        filename: str
            name of the configuration file
        """
        _log.info(f"Initializing configfile with {filename}")
        # check if a filename is given
        if filename is None:
            # load the default
            _log.info('No file given, can do nothing!!')
        else:
            # save the file name
            self.filename = filename  # list(filename.split(','))
            # create a keyword list
            keylist = self._load_file()

            # load the header
            header = ConfHeader(self.filename)

            super(ConfigFile, self).__init__(keylist, header)

    def _get_simul_name(self):
        """Get the filename used in aXeSIM"""
        # just add '.simul' and return the result
        return self.filename + '.simul'

    def confirm_extrkeys(self):
        """Confirm that all keywords for the extraction exist"""
        # default is true!
        extr_ready = 1

        # check existence of 'POBJSIZE'
        if self['POBJSIZE'] is None:
            extr_ready = 0
        # check for reasonable value
        elif float(self['POBJSIZE']) < 0.0:
            extr_ready = 0

        # check existence of 'SMFACTOR'
        if self['SMFACTOR'] is None:
            extr_ready = 0
        # check for reasonable value
        elif float(self['SMFACTOR']) < 0.0:
            extr_ready = 0

        # return the value
        return extr_ready

    def confirm_lambda_psf(self):
        """Check whether a 'lambda_psf' value is needed, provide one"""
        # check whether 'lambda_psf' is needed
        if ((self['PSFCOEFFS'] is not None) and
            (self['PSFRANGE'] is not None)):
            # split the term
            psf_range = self['PSFRANGE'].split()

            # extract the defined range as float
            lambda_min = float(psf_range[0])
            lambda_max = float(psf_range[1])

            # make 'lambda_psf' to the mean value
            lambda_psf = 0.5 * (lambda_max + lambda_min)
        else:
            # leave it at None
            lambda_psf = None

        # return the value
        return lambda_psf

    def axesim_prep(self):
        """Removes or modifies some keywords"""
        # derive the new configuration file name
        new_name = self._get_simul_name()

        # check whether the science extension has other
        # than the allowed values
        if self['SCIENCE_EXT'] != 'SCI' and self['SCIENCE_EXT'] not in ['1','2']:

            # find the index of the sciecne extension
            index = self._find_gkey('SCIENCE_EXT')

            # check whether the item was found
            if index > -1:
                # set it to the allowed value
                self.gkeys[index].keyvalue = 'SCI'
            else:
                raise aXeError("No SCI index found")

        # check whether the telesocpe are is known
        if self['TELAREA'] is None:
            # set the telescope are to the
            # Hubble default
            self.add_gkey('TELAREA', 45238.93)

        index = 1
        while self['OPTKEY'+str(index)] is not None:
            del self['OPTKEY'+str(index)]
            del self['OPTVAL'+str(index)]
            index += 1

        # just make sure that
        # the error=- and dq-
        # extensions are set
        self.add_gkey('ERRORS_EXT', 'ERR')
        self.add_gkey('DQ_EXT', 'DQ')

        # write the file back
        self.writeto(new_name)

        # return the baseic filename of the
        # simulation configuration file
        return os.path.basename(new_name)

class ConfigBeam:
    """Configuration Beam object"""
    def __init__(self, ident=None, keylist=None):
        """
        A configuration beam object is intialized. This is done
        by either extracting the relevant keywords for a certain
        beam from a keyword list or creating a default beam.

        Parameters
        ----------
        ident: char
            beam identification
        keylist: list
            list of keywords
        """

        # check if a filename is given
        if ident is None or keylist is None:
            # load the default
            _log.info('No ID or no keywords given, can do nothing!!')
        else:
            # try to load the beam keywords
            try:
                # store the ident
                self.ident = ident

                # load the general beam keywords
                self.beamkeys = self._find_beamkeys(ident, keylist)

                # load the trace keywords
                self.trace = ConfigTrace(ident, keylist)

                # load the dispersion keywords
                self.disp = ConfigDisp(ident, keylist)

            # catch a pure CKeyNotFound exception
            # which is raised if a beam is competely
            # absent in the keyword list
            except CKeyNotFound:
                raise BeamNotFound(ident)

    def __str__(self):
        """String method for the class

        The method transforms theconfiguration
        beam object into its string representation.
        """
        # initialize the return string
        rstring = ("\n#-----------\n#\n# Beam {0:s}:\n#\n#-----------\n"
                   .format(str(self.ident)))

        # add the strings for the global keys
        for key in self.beamkeys:
            rstring += str(key)

        # add the string for the trace
        rstring += str(self.trace)

        # add the string for the dispersion
        # solution
        rstring += str(self.disp)

        # return the total string
        return rstring

    def __getitem__(self, item):

        full_item = item + self.ident

        rvalue = self.get_bvalue(full_item)

        return rvalue

    def __setitem__(self, item, value):

        full_item = item + self.ident

        index = self._get_bkey_index(full_item)

        if index > -1:
            self.beamkeys[index].keyvalue = value

    def _find_beamkeys(self, ident, keylist):
        """Load the global beam keywords

        The method extracts all global beam keywords
        from a keyword list. The extracted keywords are returned
        as a list. They are removed from the input list.

        Parameters
        ----------
        ident: char
            beam identification
        keylist: list
            list of keywords
        """

        # list of the root of all globale
        # beamword keys
        bkeys = ['BEAM', 'MMAG_EXTRACT_', 'MMAG_MARK_', 'XOFF_',
                 'YOFF_', 'SENSITIVITY_']

        # list of optional keywords
        okeys = ['PSF_OFFSET_']

        # appen the beam identifier to the
        # keyword roots to get a list of keywords
        # to search for
        id_keys = []
        for key in bkeys:
            id_keys.append(key + ident)

        # initiate and fill
        # collect a list of optional keywords
        opt_keys = []
        for key in okeys:
            opt_keys.append(key + ident)

        # here is some kind of extra
        # keyword
        # ekey = 'DLD1P_' + ident + '_PRANGE'
        opt_keys.append('DLD1P_' + ident + '_PRANGE')

        # initialize the global keylist
        # and the list with indices to be deleted
        bkeys = []
        dindex = []

        # go over the keylist read in,
        # keeping and index variable
        iindex = 0
        nfound = 0
        for key in keylist:

            # identify the current keyword in the
            # list of possible ones
            if key.keyword in id_keys:

                # store the index
                dindex.append(iindex)

                # create and append the new keyword
                bkeys.append(ConfKey(key.keyword,
                                     key.keyvalue, key.comment))

                # enhance the nuber of keywords found
                nfound += 1

            elif key.keyword in opt_keys:
                # store the index
                dindex.append(iindex)

                # create and append the new keyword
                bkeys.append(ConfKey(key.keyword,
                                     key.keyvalue, key.comment))

            # enhance the index
            iindex += 1

        # check whether all keywords were found
        if nfound < len(id_keys):
            # raise an exeption if not
            raise CKeyNotFound('general')

        # delete the input keywords which
        # have been 'used'
        dindex.sort()
        dindex.reverse()
        for iindex in dindex:
            del keylist[iindex]

        # return the list of global keys
        return bkeys

    def _get_bkey_index(self, keyword):
        """Retrieve the index of a beam keyword

        The method searches for the index of
        a requested keyword in the list of beam
        keywords. If the keyword does not exists,
        the index -1 is returned

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        index: int
            the index of the keyword
        """
        # initialize the return value
        bindex = -1

        # go over all keys
        for index in range(len(self.beamkeys)):
            # check whether the current key matches
            if self.beamkeys[index].keyword == keyword:
                # return it if it matches
                return index

        # return the default
        return bindex

    def get_bkey(self, keyword):
        """Retrieve a requested beam keyword

        The method searches the list of beam keywords
        for a fitting keyword. In case that the requested
        keyword exists, it is returned.
        If not 'None' is returned

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        key: str or None
            the requested keyword or 'None'
        """
        # initialize the return value
        rkey = None

        # search for the index in the keyword list
        index = self._get_bkey_index(keyword)

        # ckeck whehter the keyword exists
        if index > -1:
            # return the keyword
            return self.beamkeys[index]
        else:
            # return the default
            return rkey

    def get_bvalue(self, keyword):
        """Retrieve a requested beam-keyword value

        The method returns the value of the keyword
        which matches the requested value.
        If there is no matching keyword, 'None'
        is returned.

        Parameters
        ----------
        keyword: str
            name of the requested keyword

        Returns
        -------
        key: str or None
            the requested keyword or 'None'
        """
        # set the default return value
        rvalue = None

        # search for the keyword
        key = self.get_bkey(keyword)

        # check whether it is non-NULL
        if key:
            # extract the value
            rvalue = key.keyvalue

        # return the value
        return rvalue

    def check_files(self):
        """Checks whether all files exist

        The method checks whether the files whose names
        are within the class data do exist or not.
        An error is reported in case that the files
        do not exist.

        """
        n_sens = 0

        # list of the root of all
        # beamword keys indicating a file
        fkeys = ['SENSITIVITY_']

        # append the beam identifier to the
        # keyword roots to get the full keyname
        for key in fkeys:
            full_keyword = key + self.ident

            # go over all beam keys
            for bkey in self.beamkeys:

                # check whether the current keyword is right
                # and whether the keyvalue is not 'None'
                if ((bkey.keyword is full_keyword) and
                     (bkey.keyvalue.upper() is not 'NONE')):
                    # check for the file
                    if not os.path.isfile(config_util.getCONF(bkey.keyvalue)):
                        # report an error
                        err_msg = ("The file: {0:s} does not exist!"
                                   .format(config_util.getCONF(bkey.keyvalue)))
                        raise aXeError(err_msg)
                    else:
                        n_sens += 1

        return n_sens


class TwoDimPolyN:
    """Object for a polynomial with 2D variance"""
    def __str__(self):
        """The method transforms the 2D polynomial object into its str
        representation.

        Returns
        -------
        object: str
            string representation of the object
        """
        # initialize the return string
        rstring = str(self.norder)

        for key in self.twodkeys:
            rstring += str(key)

        # return the total string
        return rstring

    def __getitem__(self, index):
        """Getindex method for the class

        The operator method which is called
        when an index is requested on a
        class instace
        test = kkk[0]

        Parameters
        ----------
        index: int
            the index to address
        Returns
        -------
        key : ConfListKey
            the indexed object
        """
        # check whether the index exists
        if index > len(self.twodkeys)-1:
            # raise an exception
            err_msg = "Index: {0:s} does not exist!".format(str(index))
            raise aXeError(err_msg)

        # return the indexed object
        return self.twodkeys[index]

    def __setitem__(self, index, obj):
        """Setindex method for the class

        The operator method which is called
        when the index of a class instance is
        set to a value.
        kkk[0] = test

        Parameters
        ----------
        index: int
            the index to address
        obj: ConfListKey
            description of the object content
        """
        # check whether the index exists
        if (index > (len(self.twodkeys))-1):
            # raise an exception
            err_msg = 'Index ' + str(index) + ' does not exist!'
            raise aXeError(err_msg)
        # check whether the input type is correct
        elif (not isinstance(type(self[0]), obj)):
            # raise an exception
            err_msg = ("Object: {0:s} has wrong type: {1:s}!"
                       .format(str(obj), str(type(obj))))
            raise aXeError(err_msg)

        # set the index to the input object
        self.twodkeys[index] = obj

    def _find_order(self, prefix, ident, keylist):
        """Find the keyword with the polynomial order

        The method finds and extracts the keyword
        indicating the polynomial degree from
        a keyword list. The keyword is returned.

        Parameters
        ----------
        prefix: str
            keyword prefix
        ident: char
            beam identification
        keylist: list
            list of keywords

        Returns
        -------
        keyword: str
            keyword with number of orders
        """
        # create the name of the keyword with the
        # polynomial order
        order_key = prefix + 'ORDER_' + ident

        # extract and return the keyword from the
        # keyword list
        return self._find_key(order_key, keylist)

    def _find_twodkeys(self, prefix, ident, keylist):
        """Find the all 2D polynomial keywords

        Given a prefix and a beam identifier the method
        extracts all orders of the 2D polynomial which
        describes the trace or dispersion. The number
        of orders expected is taken from the object data.

        Parameters
        ----------
        prefix: str
            keyword prefix
        ident: char
            beam identification
        keylist: list
            list of keywords

        Returns
        -------
        keys: list
            list of keywords
        """

        # initialize an empty list
        twodkeys = []

        # for each expected keyword
        for ii in range(int(self.norder.keyvalue)+1):
            # form the keyword name
            twodkey = prefix + ident + '_' + str(ii)

            # extract the new keyword
            newkey = self._find_key(twodkey, keylist, 1)

            if self._check_twodkey(newkey):
                # extract the keyword and append it to the list
                twodkeys.append(newkey)
            else:
                raise CKeyLengthWrong(ident, twodkey)

        # return the list
        return twodkeys

    def _find_key(self, keyword, keylist, lkey=0):
        """Extract a certain keyword from the list

        The methods searches for a particular keyword
        in a keyword list. If found, the keyword is
        copied and destroied in the input list.
        If not found, an exception is fired.

        Parameters
        ----------
        keyword: str
            the keyword name
        keylist: list
            list of keywords

        Returns
        -------
        keyword: str
            the extracted keyword
        """

        # initialize the index
        iindex = 0

        # set indicator to "not found"
        found = -1

        # go over all keys in the list
        for key in keylist:

            # checke whether the keyword is the desired one
            if key.keyword == keyword:
                # create a list keyword if desired
                if lkey:
                    nkey = ConfListKey(key.keyword, key.keyvalue, key.comment)
                else:
                    nkey = ConfKey(key.keyword, key.keyvalue, key.comment)

                # store the index
                found = iindex

            # enhance the index
            iindex += 1

        # fire an exception if nothing was found
        if found < 0:
            raise CKeyNotFound(keyword)
        # delete the keyword from the inlist
        else:
            del keylist[found]

        # return the keyword
        return nkey

    def _check_twodkey(self, inkey):
        """Check the length of the a field dependent keyword

        Field dependent keywords such as the polynimial
        coefficients in the trace description and dispersion
        solution must have a certain number of values,
        which is:
        n = m^2/2 + m/2
        The method checks whether the number of values
        is in agreement with this.

        @param inkey: the keyword name
        @type inkey: ConfListKey

        @return: 1/0
        @rtype: int
        """

        # determine the length of the list
        n = float(len(inkey.kvallist))

        # compute the 'order' of the xy-dependence
        m = (-1.0 + math.sqrt(1.0+8.0*n))/2.0

        # chech whether the 'order' is integer
        if math.fabs(m-int(m)) > 1.0e-16:
            # no integer -> key length wrong
            return 0

        # integer -> key length correct
        return 1

    def str_header(self, description):
        """Create a header string

        The method offers to the subclasses the possibility
        to have a meaningful string header before the
        actual data string.

        Parameters
        ----------
        @param description: description of the object content
        @type description: string
        @return: the header string
        @rtype: string
        """
        # pre-decoration
        rstring = '\n#\n# '

        # add description
        rstring += description

        # add post-decoration
        rstring += ':\n#\n'

        # return the result
        return rstring


class ConfigTrace(TwoDimPolyN):
    """Configuration Beam object"""
    def __init__(self, ident=None, keylist=None):
        """The method initializes a configuration beam

        object for a given beam identifier.
        All necessary keywords are extracted from
        an input keyword list.
        In case of missing keywords an exception
        is fired.

        Parameters
        ----------
        ident: char
            beam identification
        keylist: list
            list of keywords
        """

        # try to read in the keywords
        try:
            self.ident = ident
            self.norder = self._find_order('DYDX_', ident, keylist)
            self.twodkeys = self._find_twodkeys('DYDX_', ident, keylist)

        # raise an exception if keywords are missing
        except CKeyNotFound as e:
            raise TraceNotFound(ident, e.keyword)
        except CKeyLengthWrong as e:
            _log.info('Field dependent keyword: ' + e.keyword)

    def __str__(self):
        """Returns string representation of the object"""
        # create the label or description
        description = 'Trace description for Beam ' + str(self.ident)

        # get the string header
        rstring = super(ConfigTrace, self).str_header(description)

        # get the data string
        rstring += super(ConfigTrace, self).__str__()

        # return the result
        return rstring


class ConfigDisp(TwoDimPolyN):
    """Configuration Beam object"""
    def __init__(self, ident=None, keylist=None):
        """The method initializes a configuration dispersion

        object for a given beam identifier.
        All necessary keywords are extracted from
        an input keyword list.
        In case of missing keywords an exception
        is fired.

        Parameters
        ----------
        ident: char
            beam identification
        keylist: list
            list of keywords
        """
        # try to read in the keywords
        try:
            self.ident = ident
            self.norder = self._find_order('DISP_', ident, keylist)
            self.twodkeys = self._find_twodkeys('DLDP_', ident, keylist)
        # raise an exception if keywords are missing
        except CKeyNotFound as e:
            try:
                self.twodkeys = self._find_twodkeys('DLD1P_', ident, keylist)
            # raise an exception if keywords are missing
            except CKeyNotFound as e:
                raise DispNotFound(ident, e.keyword)
            except CKeyLengthWrong as e:
                _log.info('\nField dependent keyword: {0:s} has wrong length!'
                      .format(e.keyword))
                raise DispNotFound(ident, e.keyword)
        except CKeyLengthWrong as e:
            _log.info('\nField dependent keyword: {0:s} has wrong length!'
                  .format(e.keyword))
            raise DispNotFound(ident, e.keyword)

    def __str__(self):
        """return string representation of the object"""
        # create the label or description
        description = 'Dispersion solution for Beam ' + str(self.ident)

        # get the string header
        rstring = super(ConfigDisp, self).str_header(description)

        # get the data string
        rstring += super(ConfigDisp, self).__str__()

        # return the result
        return rstring


class DefConfHeader:
    """Default header for a configuration file"""
    def __init__(self):
        self.header = []
        self.header.append("#-----------------------------------------------"
                           "------------\n# Default configuration file for aXe"
                           "\n#\n#-------------------------------------------"
                           "---------------")

    def __str__(self):
        """returns string representation of the object"""
        rstring = ''
        for line in self.header:
            rstring += line

        return rstring


class ConfHeader(DefConfHeader):
    """Header class for the configuration file"""
    def __init__(self, filename=None):
        """Initializes the configuration header class

        The method extracts the header from a configuration
        file. If no filename is provided, a default
        header is created.

        Parameters
        ----------
        filename: str
            name of the configuration file
        """
        # no filename -> default header
        if filename is None:
            super(ConfHeader, self).__init__()
        else:
            # initialize the data list
            self.header = []

            # intialize the start pointer
            start = 1

            # open and parse through the file
            with open(filename, 'r') as fopen:
                for line in fopen:
                    # check whether the start pointer is still set
                    if start:
                        # strip the line
                        str_line = line.strip()

                        # check whether the first character
                        # is a comment, which qualifies
                        # the line as part of the header
                        if ((len(str_line) > 0) and (str_line[0] is '#')):
                            # append the line to the header data
                            self.header.append(line.strip()+'\n')
                        else:
                            # set the starter pointer to 0,
                            # thus indicating the end of the header
                            start = 0


class ConfKey:
    """Class for a keyword in a configuration file

    This keyword class is a light, but yet versatile
    and important class to strore a keyword entry in a
    configuration file. All important values are
    directly read from the object attributes.
    """
    def __init__(self, keyword, keyvalue, comment=None):
        """Constructor for the keyword class

        The keyword instance is created using
        all input values.

        Parameter
        ---------
        keyword: str
            the keword name
        keyvalue: str
            the keyword value
        comment: str
            the keyword comment
        """
        self.keyword = keyword
        self.keyvalue = keyvalue
        self.comment = comment

    def __str__(self):
        """String method for the class

        The method creats and returns
        the string representation of the
        keyword.

        Returns
        -------
        obj: str
            string representation of the object
        """
        rstring = self.keyword + ' ' + str(self.keyvalue)
        if self.comment is not None:
            rstring = rstring + ' ; ' + self.comment

        rstring += '\n'

        return rstring


class ConfListKey(ConfKey):
    """Class for a keyword list

    The keyword list class is a subclass derived from the
    keyword class. In the keyword list class has as an
    additional attribute the keyvalues transformed to a list
    of floats.
    """
    def __init__(self, keyword, keyvalue, comment=None):
        """Constructor for the keyword list class

        Initializer for the keyword list class.
        The keyword instance is created using
        all input values.

        Parameters
        ----------
        keyword: str
            the keword name
        keyvalue: str
            the keyword values
        comment: str
            the keyword comment
        """
        # initialize the keyvalue list
        self.kvallist = []

        # create a traditional keyword instance
        super(ConfListKey, self).__init__(keyword, keyvalue, comment)

        # split the string keyvalue
        vlist = self.keyvalue.split()
        for value in vlist:
            # append the floats to the list
            self.kvallist.append(float(value))

    def __getitem__(self, index):
        """Getindex method for the class

        The operator method which is called
        when an index is requested on a
        class instace
        test = kkk[0]

        Parameters
        ----------
        index: int
            the index to address

        Returns
        -------
        obj: float
            the indexed object
        """
        # check whether the index exists
        if index > len(self.kvallist)-1:
            # raise an exception
            err_msg = 'Index: ' + str(index) + ' does not exist!'
            raise aXeError(err_msg)

        # return the indexed object
        return self.kvallist[index]

    def __setitem__(self, index, obj):
        """Setindex method for the class

        The operator method which is called
        when the index of a class instance is
        set to a value.
        kkk[0] = test

        Parameters
        ----------
        index: int
            the index to address
        obj: list
            description of the object content
        """
        # check whether the index exists
        if index > len(self.kvallist)-1:
            # raise an exception
            err_msg = 'Index ' + str(index) + ' does not exist!'
            raise aXeError(err_msg)
        # check whether the input type is correct
        elif not (isinstance(type(self[0]), obj)):
            # raise an exception
            err_msg = ("Object: {0:s} has wrong type: {1:s}!"
                       .format(str(obj), str(type(obj))))
            raise aXeError(err_msg)

        # set the index to the input object
        self.kvallist[index] = obj

    def __str__(self):
        """returns the string representation of the keyword."""
        # first comes the keyword
        rstring = self.keyword

        # append the keyvalues using a default format
        for value in self.kvallist:
            rstring = rstring + ' %12.6g' % value

        # append the comment
        if self.comment is not None:
            rstring = rstring + ' ; ' + self.comment

        # append a linefeed
        rstring += '\n'

        # return the complete string
        return rstring


class ConfError(Exception):
    """Base class for exceptions in this module"""
    pass


class CKeyNotFound(ConfError):
    """Error for missing keyword"""
    def __init__(self, keyword):
        self.keyword = keyword


class BeamNotFound(ConfError):
    """Error for unknown beam """
    def __init__(self, ident):
        self.ident = ident


class TraceNotFound(ConfError):
    """Error for unknown trace"""
    def __init__(self, ident, keyword=None):
        self.ident = ident
        self.keyword = keyword


class DispNotFound(ConfError):
    """Error for unknown dispersion"""
    def __init__(self, ident, keyword=None):
        self.ident = ident
        self.keyword = keyword


class CKeyLengthWrong(ConfError):
    """Error for wrong lengt in KeywordList"""
    def __init__(self, ident, keyword=None):
        self.ident = ident
        self.keyword = keyword
