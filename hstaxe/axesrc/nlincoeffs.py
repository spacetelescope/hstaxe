import os
import numpy as np

from astropy.io import fits

from stwcs.wcsutil import HSTWCS
from hstaxe.axeerror import aXeError


class NonLinCoeffs:
    def __init__(self, image, ext_info):
        """Initializes the class"""
        # dummy example
        # self.default_coeffs_data ="""# Dummy polynomial distortion coefficients
        # # file, created by aXe
        # cubic
        # 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

        # 0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        # """
        self.detector=1
        self.image = image
        self.ext_info = ext_info
        self._set_coefficients()
        self._set_coeff_filename()
        self._make_header()

    def _set_coefficients(self):
        """Extracts the drizzle coefficients."""
        try:
            self.imwcs = HSTWCS(self.image, ext=self.detector)
            #self.xcoeffs, self.ycoeffs = coeff_converter.sip2idc(self.imwcs)
        except ValueError:
            raise aXeError("Could not determine distorsion coefficients from header.")

        # dictionary with the fixed names for the orders
        fixed_orders = {1: 'constant', 2: 'linear', 3: 'quadratic', 4: 'cubic', 5: 'quintic'}
        self.xcoeffs = np.array([[0.0],[1.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]])
        self.ycoeffs = np.array([[0.0],[0.0],[1.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]])

        # return the order
        # for HST the a and b orders should always match
        # sip.ap_order is the inverse
        # if (self.imwcs.sip.a_order != self.imwcs.sip.b_order):
        #     raise aXeError("SIP a and b orders don't match!")
        # self.order = fixed_orders[self.imwcs.sip.a_order]
        try:
            self.order = fixed_orders[len(self.xcoeffs)/2]
        except KeyError:
            raise aXeError(f'Coefficient order not available: {self.order} in {fixed_orders}')

    def _set_coeff_filename(self):
        """Get the missing information about the file name"""
        # compose the name of the coefficients file
        self.coeff_filename = self.image.replace('.fits', '_coeffs'+str(self.imwcs.chip)+'.dat')
        


    # def _create_dummy_coeffs(self):
    #     """Make default coefficients"""

    #     # get the instrumental setup
    #     # instr_setup = self._get_spect_setup()

    #     # get the default coefficients
    #     # self._get_default_coeffs(instr_setup, detector)

    #     # get the default coefficient
    #     # name and destroy any existing
    #     # file with this name
    #     if os.path.isfile(self.coeff_filename):
    #         os.unlink(self.coeff_filename)

    #     # create, fill and close the
    #     # default coefficients file
    #     cfile = open(self.coeff_filename, 'w+')
    #     cfile.write(self.default_coeffs_data)
    #     cfile.close()

    def _make_header(self):
        """Generate a header"""
        # make an empty header
        self.coeff_fileheader = []

        # make some specific header phrases
        self.coeff_fileheader.append('-----------------------------------------')
        self.coeff_fileheader.append('coefficients file generated for aXedrizzle')
        self.coeff_fileheader.append(f'from keywords in image: {os.path.basename(self.image)}')

    def write_file(self):
        """Generate the non-linear coefficients file"""
        print("\nWriting non-linear coefficients to file\n")

        # check whether the file exists
        with open(self.coeff_filename, 'w+') as cfile:
            cfile.write("# created by axe\n")
            cfile.write(f"# Coefficients generated from SIP in {self.image} and {self.imwcs.idctab}\n")
            cfile.write("refpix " + " ".join(map(str,self.imwcs.wcs.crpix)) + f"\n{self.order}\n")
            for i in self.xcoeffs:
                cfile.write(("\t".join(str(j) for j in i) + " "))
            cfile.write("\n\n")
            for i in self.ycoeffs:
                cfile.write(("\t".join(str(j) for j in i) + " "))       
        

    def store_coeffs(self):
        """Store coeff information in the image header"""
        print("\nStoring non-linear coefficients: ")
        print(f"{self.imwcs.instrument}: {self.imwcs.detector}\n")

        # open the fits and go to the first header
        flt_imag = fits.open(self.image, mode='update', memmap=0)
        flt_head = flt_imag[0].header

        
        # translate instrument and detector to
        # a linear pixel scale;
        # store this pixel scale
        flt_head['DRZSCALE'] = (self.imwcs.idcscale, 'Scale for drizzling')
        
        print(f"DRZSCALE set to: {flt_head['DRZSCALE']}")

        # store the number of coefficients
        flt_head['DRZCNUM'] = (len(self.xcoeffs.flatten()/len(self.xcoeffs)),
                               'Number of coefficients per xy coordinate')
        # store the x-coefficients
        if (len(self.xcoeffs) == 0):
            raise ValueError("No x-coefficients recorded")
        else:
            for idx, val in enumerate(self.xcoeffs.flatten()):
                keyname = 'DRZ{0:1d}X{1:02d}'.format(int(self.ext_info['axe_ext']), idx+1)
                keycomm = 'Drizzle coefficient {0:02d} in X'.format(idx+1)
                flt_head[keyname] = (val, keycomm)

        # store the y-coefficients
        if (len(self.ycoeffs) == 0):
            raise ValueError("No y-coefficients recorded")
        else:
            for idx, val in enumerate(self.ycoeffs.flatten()):
                keyname = 'DRZ{0:1d}Y{1:02d}'.format(int(self.ext_info['axe_ext']), idx+1)
                keycomm = 'Drizzle coefficient {0:02d} in Y'.format(idx+1)
                flt_head[keyname] = (val, keycomm)

        # close the fits
        flt_imag.close()
