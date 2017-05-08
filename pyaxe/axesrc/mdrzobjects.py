from __future__ import (absolute_import, unicode_literals, division,
                        print_function)


import sys
import math
from .. import axeutils
from . import drizzleobjects


class MulDrzObjList(drizzleobjects.DrizzleObjectList):
    """List class for all objects to be drizzled"""
    def __init__(self,
                 drizzle_params,
                 mult_drizzle_par,
                 cont_info=None,
                 opt_extr=False,
                 back=False,
                 drztmp_dir=None,
                 drizzle_dir=None):

        # load the drizzle parameters
        self.drizzle_params = drizzle_params.copy()

        # save the multidrizzle parameters
        self.mult_drizzle_par = mult_drizzle_par.copy()

        # store the quantitative contamination flag
        self.cont_info = cont_info

        # store the optimal extraction flag
        self.opt_extr = opt_extr

        # store the background flag
        self.back = back

        # save the drizzle tmp-directory;
        # use the default if not explicitly given
        if drztmp_dir is not None:
            self.drztmp_dir = drztmp_dir
        else:
            self.drztmp_dir = axeutils.getDRZTMP()

        # save the drizzle directory;
        # use the default if not explicitly given
        if drizzle_dir is not None:
            self.drizzle_dir = drizzle_dir
        else:
            self.drizzle_dir = axeutils.getDRIZZLE()

        # get the identifier for drizzle objects
        self.regexp = self._get_regexp(back=False)

        # get all drizzle objects
        objectlist = self._find_drizzle_objects(self.drztmp_dir,
                                                self.regexp,
                                                back=False)

        # convert the objects list to a list of objects
        self.drizzle_objects = self._objlist_to_drzobjects(objectlist,
                                                           self.drizzle_params,
                                                           self.mult_drizzle_par,
                                                           self.cont_info,
                                                           self.opt_extr,
                                                           self.back,
                                                           self.drztmp_dir,
                                                           self.drizzle_dir)


    def _objlist_to_drzobjects(self,
                               objectlist,
                               drizzle_params,
                               mult_drizzle_par,
                               cont_info,
                               opt_extr,
                               back,
                               drztmp_dir,
                               drizzle_dir):
        """
        Converts the object list into drizzle objects
        """
        # create an empty list
        drzobjects = []

        # split the dictionary
        # into a key-value list
        item_list = objectlist.items()

        # go over the list
        for an_item in item_list:
            # create a drizzle object and append it to the list
            drzobjects.append(MultDrzObj(an_item[0],
                                         an_item[1],
                                         drizzle_params,
                                         mult_drizzle_par,
                                         cont_info,
                                         opt_extr,
                                         back,
                                         drztmp_dir,
                                         drizzle_dir))

        # return the list
        return drzobjects

    def store_reject_info(self):
        """Retrieve and store the info on the rejection process"""
        # go over all drizzle object
        for drizzleObject in self.drizzle_objects:

            # determine the information on rejected pixels
            reject_info = drizzleObject.get_reject_info()

            # store the information on the images
            drizzleObject.update_reject_info(reject_info)

    def multidrizzle(self):
        """Drizzle all objects"""
        # go over all drizzle object
        for drizzleObject in self.drizzle_objects:

            # do the individual drizzles
            drizzleObject.singdrizzle()

            # combine the individual drizzles
            drizzleObject.mediancombine()

            # blot the combined image
            drizzleObject.blot()

            # identify the cosmic rays
            drizzleObject.drzrej()

        # do the final drizzle
        self.drizzle()

        # get information in the
        # rejection process
        self.store_reject_info()


class MultDrzObj(drizzleobjects.DrizzleObject):
    """List class for all objects to be drizzled"""
    def __init__(self,
                 objID,
                 file_list,
                 drizzle_params,
                 mult_drizzle_par,
                 cont_info,
                 opt_extr,
                 back,
                 drztmp_dir,
                 drizzle_dir):

        # save the object number
        self.ID = int(objID[2:])

        # define and save the object ID
        self.objID = objID

        # save the drizzle parameters
        self.drizzle_params = drizzle_params.copy()

        # save the multidrizzle params
        self.mult_drizzle_par = mult_drizzle_par.copy()

        # store the quantitative contamination flag
        self.cont_info = cont_info

        # save the optimal extraction flag
        self.opt_extr = opt_extr

        # save the background flag
        self.back = back

        # save the drizzle directory
        # and the drizzle tmp-directory
        self.drztmp_dir = drztmp_dir
        self.drizzle_dir = drizzle_dir

        # define the name of the object directory
        self.objID_dir = self._get_objID_dirname(self.objID,
                                                 self.drztmp_dir,
                                                 back=False)

        # generate the list of contributors
        self.contrib_list = self._make_contrib_list(self.objID,
                                                    file_list,
                                                    opt_extr,
                                                    drztmp_dir=drztmp_dir)

        # determine all relevant names
        self.ext_names = self._get_ext_names(self.drizzle_params['ROOT'],
                                             self.objID, back=False,
                                             drizzle_dir=self.drizzle_dir)

        # get the number of contributors
        self.ncontrib = self._get_ncontrib()

    def __str__(self):
        """Defines a string representation"""
        return "{0:s}: {1:d} image contributions.\n".format(self.objID,
                                                            len(self))

    def _make_contrib_list(self, objID, file_list, opt_extr, drztmp_dir):
        """Generates a list of contributing objects"""
        # make an empty list
        contrib_list = []

        # go over all contributing files
        for a_file in file_list:
            # generate an object and append it to the list
            contrib_list.append(drizzleobjects.DrizzleObjectContrib(a_file,
                                                                    objID,
                                                                    opt_extr,
                                                                    back=False,
                                                                    drztmp_dir=drztmp_dir))

        # return the llist of contributors
        return contrib_list

    def _migrate_crr(self, wht_image, crr_image):
        from astropy.io import fits as pyfits

        # open the two images
        wht_img = pyfits.open(wht_image, 'update')
        crr_img = pyfits.open(crr_image, 'readonly')

        # migrsate the crr's
        wht_img[0].data = wht_img[0].data * crr_img[0].data

        # close everything
        wht_img.flush()
        wht_img.close()
        crr_img.close()

    def singdrizzle(self):
        """MultiDrizzle all contributors together"""

        msg = "MultiDrizzling object : {0:10s} ... ".format(self.objID)
        print(msg)
        sys.stdout.flush()

        # create a drizzle object
        drizzleObject = dither.Drizzle()

        # go over all contributing objects
        for one_contrib in self.contrib_list:

            img_nx = int(one_contrib.info['LENGTH'])
            img_ny = 2*int(math.ceil(one_contrib.info['OWIDTH'])) + 10

            # run drizzle for the object data
            drizzleObject.run(one_contrib.ext_names['FLT'],
                              one_contrib.ext_names['WHT'],
                              one_contrib.ext_names['SING_SCI'],
                              one_contrib.ext_names['SING_WHT'],
                              one_contrib.ext_names['CFF'],
                              one_contrib.info['EXPTIME'],
                              self.drizzle_params,
                              img_nx,
                              img_ny)

        # give feedback
        print('Done!')
        sys.stdout.flush()

    def mediancombine(self):
        """
        Combine the individual drizzles
        """
        from . import dither

        mcomb = dither.MedianCombine(self.contrib_list,
                                     self.drizzle_params,
                                     self.mult_drizzle_par,
                                     self.ext_names)
        mcomb.run()

    def blot(self):
        """Blot the median image back"""
        from . import dither

        blotObject = dither.Blot()

        # go over all contributing objects
        for one_contrib in self.contrib_list:

            # blot the median image back
            blotObject.run(self.ext_names['MED'],
                           one_contrib.ext_names['BLT'],
                           one_contrib.ext_names['CFF'],
                           one_contrib.info['NAXIS1'],
                           one_contrib.info['NAXIS2'],
                           self.drizzle_params,
                           self.mult_drizzle_par)

    def drzrej(self):
        """Do the CR-rejection"""

        # make a deriv object
        derivObject = dither.Deriv()

        # go over all contributing objects
        for one_contrib in self.contrib_list:

            # generate a derivate of the blotted back image
            derivObject.run(one_contrib.ext_names['BLT'],
                            one_contrib.ext_names['DER'])

        # make an identification object
        cridentObject = dither.CRIdent(self.drizzle_params,
                                       self.mult_drizzle_par)

        # go over all contributing objects
        for one_contrib in self.contrib_list:

            # identify the CR's on each contributor
            cridentObject.run(one_contrib.ext_names['FLT'],
                              one_contrib.ext_names['BLT'],
                              one_contrib.ext_names['DER'],
                              one_contrib.info['EXPTIME'],
                              one_contrib.info['SKY_CPS'],
                              one_contrib.ext_names['CRR'])

            # apply the crr information onto the wht-image
            self._migrate_crr(one_contrib.ext_names['WHT'],
                              one_contrib.ext_names['CRR'])
