import os
from astropy.io import fits

from . import interpolator
from . import modspeclist

from .. import axe_asciidata
from ..axeerror import aXeSIMError
from ..axeutils import getDATA, getSIMDATA, getOUTSIM


class TemplateSpectrumList:
    """Class for the template spectra"""

    def __init__(self, inlist, incat, tpass_flux):
        """Initializer for the class

        @param inlist: name of list with ascii spectra
        @type inlist: string
        @param incat: model object list
        @type incat: string
        @param tpass_flux: file or expression for the total flux passband
        @type tpass_flux: string
        """
        # initialize the list of spectra
        self._speclist = []

        # load the spectrum template list
        ilist = axe_asciidata.open(inlist)

        # load the total passband
        self.tpass_flux = tpass_flux
        tpass = self._load_tpass_flux(tpass_flux)
        self.piv_wav = tpass.pivot()

        # load the model spectrum list
        self.modspeclist = modspeclist.ModelSpectrumList(getDATA(incat),
                                                         self.piv_wav)
        maxnum = self.modspeclist.find_tempmax()

        # check that the spectral template index
        # and the number of template spectra match
        if maxnum > ilist.nrows:
            err_msg = ("There are not enough templates in the spectrum "
                       "template list {0:s} !".format(inlist))
            raise aXeSIMError(err_msg)

        # initialize the
        # extension counter
        spec_index = 1
        # go through the list
        for index in range(self.modspeclist.nrows):
            # check whether a hig-res model
            # is desired, skip oterwise
            if (self.modspeclist[self.modspeclist.reqColIndex['SPECTEMP']][index] < 1):
                print("Nothing to be done for object: {0:s} in row {1:s}"
                      .format(str(self.modspeclist[self.modspeclist.reqColIndex['NUMBER']][index]), str(index+1)))

            # extract the template name
            spectemp = ilist[0][self.modspeclist[self.modspeclist.reqColIndex['SPECTEMP']][index]-1]

            # extract and compose the catalogue magnitude
            maginfo = [self.modspeclist[self.modspeclist.reqColIndex['MAGNITUDE']][index], self.modspeclist.magwave]

            # extract the redshift
            zvalue = self.modspeclist[self.modspeclist.reqColIndex['Z']][index]

            # compose and add a new shift model
            self._speclist.append(TemplateSpectrum(getSIMDATA(spectemp),
                                                   tpass, maginfo, zvalue))

            # insert the column value for the
            # aXe simulations
            self.modspeclist["MODSPEC"][index] = spec_index

            # enhance the
            # extension counter
            spec_index += 1

    def _load_tpass_flux(self, tpass_flux):
        """
        Create/load the total passband

        The method loads an interpolator from the expression given in the
        input parameter. Te input either specifies an existing file in fits
        or in ASCII format, or an expression "lambda_min, lambda_max"
        with lambda_min and lambda_max min and max of the independent data,
        respectively. The corresponding dependent values are taken as
        1.0 everywhere.


        @param tpass_flux: expression for the total flux passband
        @type tpass_flux: string

        @return: the interpolator for the total passband
        @rtype: Interpolator()
        """
        # check whether a file exists
        if os.path.isfile(getSIMDATA(tpass_flux)):
            # load the passband from file
            ip = interpolator.Interpolator(getSIMDATA(tpass_flux))

            # the file is in Angstrom,
            # convert to [nm]
            ip.tonm()

        else:
            # resove the expression
            # to get the independent data
            indep = self._get_tpass_fromexpr(tpass_flux)

            # the dependent data is
            # always constant
            depen = [0.0, 1.0, 1.0, 0.0]

            # load the interpolator from the list
            ip = interpolator.Interpolator(indep=indep, depen=depen)

        # return the interpolator
        return ip

    def _get_tpass_fromexpr(self, tpass_flux):
        """
        Resolve an input expression in a passband

        The method analyzes the expression given in the input. The expression
        is supposed to be a comma separated pair 'wav_min, wax_max'.
        On both sides a point with the distance 0.01*(max-min) is added
        for a smooth transgression (though perhaps not really needed)

        @param tpass_flux: file or expression for the total flux passband
        @type tpass_flux: string

        @return: independent values derived from the expression
        @rtype: []
        """
        wavs = []
        indep = []

        # check whether there is a comma
        if tpass_flux.find(',') < 0:
            err_msg = ("Can not resolve the expression: {0:s} "
                       "into a total passband lambda_min,lambda_max:\n"
                       .format(str(tpass_flux)))
            if (AXE_SIMDATA_PATH == './'):
                err_msg += ("\nOr the file: {0:s} does not exist!"
                            .format(getSIMDATA(tpass_flux)))
                raise aXeSIMError(err_msg)

        # split the expression
        # according to the commas
        tp_list = tpass_flux.split(',')

        # go over the comma separated list
        for item in tp_list:
            # convert the items
            # to float
            try:
                wavs.append(float(item))
            except ValueError:
                err_msg = ("\nCan not resolve item: {0:s} in expression {1}"
                           "to a wavelength value!".
                           format(str(item), tpass_flux))
                raise aXeSIMError(err_msg)

        # sorted the array
        wavs.sort()

        # build up the array with extensions
        indep.append(wavs[0] - 0.01*(wavs[1] - wavs[0]))
        indep.append(wavs[0])
        indep.append(wavs[1])
        indep.append(wavs[1] + 0.01*(wavs[1] - wavs[0]))

        # return the array
        return indep

    def beamto(self):
        """Transform every list member to the right redshift and flux level

        The method passes over the whole list of template spectra
        and invokes that they are shifted in redshift and flux space.
        """
        # go over all spectra
        for spec in self._speclist:
            # beam them to the desired magnitude
            # and redshift
            spec.beamto()

    def tofits(self, fitsname, indata_copy=0):
        """Converts and stores the spectra in a fits file

        The method stores the redshift- and flux- shifted template spectra
        in a model spetra file, which is a multi-extension fits table
        with one model spectrum on each extension.
        A flagg indicates whether, besides the normal output to
        AXE_OUTSIM_PATH, a copy to AXE_IMAGE_PATH is desired.

        Parameters
        ----------
        fitsname: str
            name for the MEX fits file
        indata_copy: int
            flag to save also a copy
        """
        # create a HDU list
        hdulist = fits.HDUList()

        # create an empty primary HDU
        phdu = fits.PrimaryHDU()

        # put the primary to the list
        hdulist.append(phdu)

        # give a linefeed
        print()

        # initialize an index
        index = 0

        # go over all spectra
        for spec in self._speclist:
            # enhance the counter
            index += 1

            # convert the spectrum to a fits extension
            hdu = spec.tofits()

            # fill in some header keywords
            hdu.header['SPECNAME'] = (os.path.basename(spec.specname),
                                      "Name of spectrum")
            hdu.header['EXTNAME'] = (os.path.basename(spec.specname),
                                     "Name of spectrum")
            hdu.header['PASSBAND'] = (self.tpass_flux,
                                      "Total passband for flux scaling")
            hdu.header['PB_PIVOT'] = ("{0:.2f}".format(self.piv_wav),
                                      "Pivot wavelength of total passband")
            hdu.header['MAG_AB'] = (spec.mag_info[0],
                                    "AB-magnitude for flux scaling")
            hdu.header['WAV_AB'] = (spec.mag_info[1],
                                    "Wavelength for AB magnitude [nm]")
            hdu.header['Z'] = (spec.z, "Redshift of spectrum")

            # append the fits created from
            # the spectrum
            print("Adding: {0:s}, z={1:s}, mag={2:s} to {3:s}, ext: {4:s}"
                  .format(os.path.basename(spec.filename), str(spec.z),
                          str(spec.mag_info[0]), fitsname,  str(index)))
            hdulist.append(hdu)

        # delete older versions
        # of the fits name
        if os.path.isfile(getOUTSIM(fitsname)):
            os.unlink(getOUTSIM(fitsname))

        print("\nWriting model spectra to file: {0:s} ...."
              .format(getOUTSIM(fitsname)))
        # write it to fits
        hdulist.writeto(getOUTSIM(fitsname))
        print("Done")

        # check whether a copy
        # to AXE_IMAGE_PATH is desired
        if indata_copy:
            # delete older versions
            # of the fits name
            if os.path.isfile(getDATA(fitsname)):
                os.unlink(getDATA(fitsname))

            print("Writing model spectra to file: {0:s} ..."
                  .format(getDATA(fitsname)))
            # write it to fits
            hdulist.writeto(getDATA(fitsname))
            print('Done')

        # add an extra linefeed
        print('')

    def save_modspec_list(self, mspec_name_new=None):
        """
        Save the modified model spectrum list

        The class model spectrum list is modified to contain the references
        to the correct extensions in the model spectra file created.
        The method saves the model spectrum list to the disk with a name
        given in the input

        @param mspec_name_new: name for the modified model spectra file
        @type mspec_name_new: string
        """
        # look whether a new name is given
        if (mspec_name_new is not None):
            # compose the new full name
            new_mspec = getDATA(mspec_name_new)
        else:
            # compose the new full name
            new_mspec = getDATA(self.modspeclist.filename)

        # write the model spectrum list out
        print("Writing model spectrum list to file: {0:s} ..."
              .format(new_mspec))
        self.modspeclist.writeto(new_mspec)
        print("Done\n")

    def check(self):
        """Invokes a check on all spectra

        The method passes over all spectra in the list., A basic
        check is invoked on each individual spectrum.
        """
        # go over all spectra
        for spec in self._speclist:
            # check the current spectrum
            spec.check()


class TemplateSpectrum:
    """General class to shift spectra in flux-redshift range
    """
    def __init__(self, filename, tpass,  mag_info, z=0.0):
        """Initializes the class

        filename: str
            file name of the template spectrum
        tpass: Interpolator
            interpolator to represent the total passband
        mag_info: list of floats
            information on target AB magnitude and corresponding weavelength
        z: float
            target redshift for the template spectrum
        """
        # check for the spectru file
        if not os.path.isfile(filename):
            err_msg = "File: {0:s} does not exist!".format(filename)
            raise aXeSIMError(err_msg)

        # store the template file
        self.filename = filename

        # store the passband
        self.tpass = tpass

        # store the AB magnitude to shift to
        self.mag_info = mag_info

        # store the redshift
        self.z = z

        # load in the spectrum
        self.data = interpolator.Interpolator(filename)
        self.data.tonm()

        # report on what you did
        print("File: {0:s} loaded for z={1:s}, mag={2:s}"
              .format(filename, str(z), str(mag_info[0])))

        # define the specturm name,
        # which is the file name until
        # the last '.', e.g. 'sbc.dat' --> 'sbc'
        self.specname = filename[:filename.rfind('.')]

    def check(self):
        """Check whether the shift in redshift is possible

        The method does a basic check and compares the boundaries of the
        flux normalization passband and the boundaries of the redshifted
        template spectrum. In case that the redshifted template spectrums
        is not totally included in the boundaries of the passband, an
        exception is issued.
        """
        # compute the interval where the redshifted spectra is defined
        ind_min_redshift = (1.0 + self.z) * self.data.ind_min
        ind_max_redshift = (1.0 + self.z) * self.data.ind_max

        # check whether the total passband for normalization
        # is completely inside this interval
        if ((self.tpass.ind_min < ind_min_redshift) or
           (self.tpass.ind_max > ind_max_redshift)):
            err_msg = ("The redhifted template spectrum: {0:s} "
                       "is defined on the interval:".format(self.filename))
            err_msg += ("[{0:.2f}, {1:.2f}]AA.\nThe total passpand is defined"
                        "on [{2:.2f}, {3:.2f}]AA.\n"
                        .format(ind_min_redshift*10.0,
                                ind_max_redshift*10.0,
                                self.tpass.ind_min*10.0,
                                self.tpass.ind_max*10.0))
            err_msg += ("The redshifted template spectrum does not entirely"
                        "cover the total passband, hence an accurate shift"
                        "in flux space is not possible!\n"
                        "Please extend the template spectrum or choose a"
                        "different passband!")
            raise aXeSIMError(err_msg)

    def beamto(self):
        """Beam the spectrum to a certain redshift and AB-magnitude

        The method perfomrs the shift in redshift and the normalization
        of the flux in the given passband and to the given
        AB-magnitude value.
        """
        # shift to the redshift
        self.data.mult_indep(1.0+self.z)

        # get the flux the spectrum shall be scaled to
        flux_AB = self.get_flambda_from_magab(self.mag_info[0],
                                              self.mag_info[1])

        # determine the average flux value over
        # the desired wavelength range
        flux_tpass = self.tpass * self.data

        # check for spectral overlap
        if (len(flux_tpass) == 0):
            err_msg = ("\nThere is no spectral overlap between the total "
                       "passband and the redhsifted template!")
            raise aXeSIMError(err_msg)

        # integrate over the passband;
        # check for zero integral
        tpass_int = self.tpass.integrate()
        if tpass_int <= 0.0:
            err_msg = "\nThe integral over the total passband is <= 0.0!"
            raise aXeSIMError(err_msg)

        # compute the average flux of the spectrum in the passband
        # check for zero average
        flux_ave = flux_tpass.integrate() / tpass_int
        if flux_ave <= 0.0:
            err_msg = ("\nThe average flux of the spectrum in the "
                       "passband is <= 0.0!")
            raise aXeSIMError(err_msg)

        # compute the scaling factor
        flux_scale = flux_AB / flux_ave

        # apply the scaling factor
        self.data.mult_depen(flux_scale)

    def tofits(self):
        """Transforms the template spectrum to fits

        The method transforms the model spectrum to
        a fits extension.
        """
        return self.data.tofits('WAV_NM', 'FLUX')

    def get_flambda_from_magab(self, mag, wlength):
        """Converts a mag_AB value at a wavelength to f_lambda

        @param mag: mag_AB value
        @type mag: float
        @param wlength: wavelength value
        @type wlength: float

        @return: the f_lambda value
        @rtype: float
        """
        """
        /**
        *
        * Function: get_flambda_from_magab
        * The subroutine calculates the flambda value for a
        * mag_AB value given with its wvavelength as input
        * parameters.
        *
        * Parameters:
        * param  mag     - the mag_AB value
        * param  lambda  - the wavelength for mag_AB [nm]
        *
        * Returns:
        * return flambda - the calculated flambda value
        */
        double
        get_flambda_from_magab(double mag, double lambda)
        {
        double flambda=0.0;
        double fnu=0.0;

        fnu     = pow(10.0, -0.4*(mag+48.6));
        flambda = 1.0e+16*LIGHTVEL*fnu/(lambda*lambda);

        return flambda;
        };
        """
        import math

        fnu = math.pow(10.0, -0.4*(mag+48.6))
        if (wlength != 0):
            flambda = 2.99792458e+16 * fnu / (wlength * wlength)
        else:
            flambda = 0.
        return flambda
