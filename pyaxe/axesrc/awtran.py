import sys
import math
from drizzlepac import astrodrizzle
from stsci.tools import fileutil
"""

WTRAN - a wrapper for wtraxy/wtranback to transform a position (x,y) in pixels
       on an input, distorted image to/from an output image.

There are two methods: wtran.f for forward transforms and wtran.b for
the reverse.

The syntax is:

 wtran.f(original_image,drizzled_image,x,y)
  --or--
 wtran.f(original_image,drizzled_image,List='list')

and

 wtran.b(drizzled_image,original_image,x,y)
  --or--
 wtran.b(drizzled_image,original_image,List='list')

In the 'list' case the list is a normal text file with two, free-
format columns of X,Y pixel positions.

All the information is extracted from the header by searching for
the drizzle header records. The coefficients file must be present and
have the same name. This is in "drizzle" format, as produced by PyDrizzle,
and not the IDCTAB.

Note - the 'original_image' name must match the string written to
the header by drizzle, exactly.

It is assumed that this script is invoked from Pyraf and that the
wtraxy and wtranback IRAF tasks are available. They are in the dither
package of STSDAS.

Example:

--> import wtran

--forwards---

--> wtran.f('j8c0c1011_crj.fits[sci,1]','f606w_z61c1_drz.fits[1]',136.109,371.455)
-Reading drizzle keywords from header...
-Image  j8c0c1011_crj.fits[sci,1]  was # 3
-to be drizzled onto output  f606w_z61c1_drz.fits[1]
-Transforming position...
 Xin,Yin:    136.109   371.455 Xout,Yout:    123.000   432.000

--backwards---

--> wtran.b('f606w_z61c1_drz.fits[1]','j8c0c1011_crj.fits[sci,1]',123,432)
-Reading drizzle keywords from header...
-Image  j8c0c1011_crj.fits[sci,1]  was # 3
-to be drizzled onto output  f606w_z61c1_drz.fits[1]
-Transforming position...
 Xin,Yin:    136.109   371.455 Xout,Yout:    123.000   432.000

Richard Hook, ST-ECF/STScI, April 2003
Added "List" feature and other small improvements, November 2003

Added trap for error messages in "list" form, May 2004.

Version 0.12 for STSDAS 3.3 release, October 2004
Added more robust handling of wavelengths and DGEO image support.

Version 0.20
PyDrizzle will be automatically run to generate coeffs files if not already
present for input image.

Version 0.21
Syntax for calling PyDrizzle updated for new 'bits' syntax.

Created wtran to use wtraxy and wtranback, and hence be MUCH faster in list mode.
Richard Hook, ST-ECF, April 2008

Comments: rhook@eso.org
"""
"""
Howard Bushouse, STScI, 08-Mar-2001, version 1.3
Updated 'raise' statements to use Exception class, to conform to latest python
syntax rules.
"""

# Some convenient definitions
yes=iraf.yes
no=iraf.no
MaxImages=999
TRUE=1
FALSE=0

__version__ = '0.1 (Apr2008)'

# A class for drizzle geometrical parameters
class DrizGeoPars:

    # Constructor, set to drizzle default values
    def __init__(self,image=None,inimage=None):

        if image == None:
            self.scale=1.0
            self.coeffs=None
            self.lam=555.0
            self.xsh=0.0
            self.ysh=0.0
            self.rot=0.0
            self.shft_un="input"
            self.shft_fr="input"
            self.align="center"
            self.xgeoim=""
            self.ygeoim=""
            self.d2xscale=0.0
            self.d2yscale=0.0
            self.d2xsh=0.0
            self.d2ysh=0.0
            self.d2rot=0.0
            self.d2shft_fr="output"
        else:

            # Read geometric parameters from a header using an image name as
            # the key

            found=FALSE

            # First search for the entry for this image
            i=1
            while i < MaxImages:
                datkey = 'D%3iDATA' % i
                datkey=datkey.replace(' ','0')

                iraf.keypar(image,datkey,silent='yes')

                # If we can't read this no point considering
                if iraf.keypar.value == '':
                    break

                # If we have a match set flag and leave
                if iraf.keypar.value == inimage:
                    found=TRUE
                    break

                i += 1

            if not found:
                raise Exception("Failed to get keyword information from header.")

            # Now we know that the selected image is present we can
            # get all the other parameters - we don't check whether this
            # succeeds, if it doesn't let it crash
            stem=datkey[:4]

            iraf.keypar(image,stem+"SCAL",silent='yes')
            self.scale=float(iraf.keypar.value)

            iraf.keypar(image,stem+"COEF",silent='yes')
            self.coeffs=iraf.keypar.value
            # Check for existence
            if fileutil.findFile(self.coeffs) == FALSE:
               try:
                  print('\n-Coeffs file not found.  Trying to reproduce them using PyDrizzle...')
                  # Try to generate the coeffs file automatically
                  indx = inimage.find('[')
                  p = astrodrizzle.adrizzle.drizzle(inimage[:indx], bits_single=None, bits_final=None, updatewcs=FALSE)
                  del p
               except:
                  print("! Cannot access coefficients file. (",self.coeffs,")")
                  raise Exception("File missing or inaccessible.")

            iraf.keypar(image,stem+"LAM",silent='yes')
            if iraf.keypar.value != '':
               self.lam=float(iraf.keypar.value)
            else:
               self.lam=555.0

            iraf.keypar(image,stem+"XSH",silent='yes')
            self.xsh=float(iraf.keypar.value)

            iraf.keypar(image,stem+"YSH",silent='yes')
            self.ysh=float(iraf.keypar.value)

            iraf.keypar(image,stem+"ROT",silent='yes')
            self.rot=float(iraf.keypar.value)

            iraf.keypar(image,stem+"SFTU",silent='yes')
            self.shft_un=iraf.keypar.value

            iraf.keypar(image,stem+"SFTF",silent='yes')
            self.shft_fr=iraf.keypar.value

            iraf.keypar(image,stem+"XGIM",silent='yes')
            self.xgeoim=iraf.keypar.value
            indx = self.xgeoim.find('[')
            # Check for existence
            if fileutil.findFile(self.xgeoim[:indx]) == FALSE and self.xgeoim != '':
               print("! Warning, cannot access X distortion correction image")
               print(" continuing without it. (",self.xgeoim,")")
               self.xgeoim=''

            iraf.keypar(image,stem+"YGIM",silent='yes')
            self.ygeoim=iraf.keypar.value
            indx = self.ygeoim.find('[')
            # Check for existence
            if fileutil.findFile(self.ygeoim[:indx]) == FALSE and self.ygeoim != '':
               print("! Warning, cannot access Y distortion correction image")
               print(" continuing without it. (",self.ygeoim,")")
               self.ygeoim=''

            # The case of the "align" parameter is more tricky, we
            # have to deduce it from INXC keyword
            iraf.keypar(image,stem+"INXC",silent='yes')
            inxc=float(iraf.keypar.value)

            # Need the X and Y dimensions as well - both input and
            # output
            iraf.keypar(inimage,'i_naxis1',silent='yes')
            xdim=int(iraf.keypar.value)
            iraf.keypar(inimage,'i_naxis2',silent='yes')
            ydim=int(iraf.keypar.value)

            self.nxin=xdim
            self.nyin=ydim

            iraf.keypar(image,'i_naxis1',silent='yes')
            xdim=int(iraf.keypar.value)
            iraf.keypar(image,'i_naxis2',silent='yes')
            ydim=int(iraf.keypar.value)

            self.nxout=xdim
            self.nyout=ydim

            if abs(inxc-float(xdim/2)-0.5) < 1e-4:
                self.align='corner'
            else:
                self.align='center'

            # Check for the presence of secondary parameters
            iraf.keypar(image,stem+"SECP",silent='yes')
            if iraf.keypar.value == "yes":
                raise Exception("Sorry, this version does NOT support secondary parameters.")
            else:
                self.secp=FALSE

# Main TRAN methods - f for forward and b for back
#
# inimage - the input image which is to have its WCS updated
# drizimage - the reference image, assumed to contain the drizzle parameters
#         in its header
#
# x,y - a single position for transformation
#
# List - a text file name containing x y pairs
#
def f(origimage,drizimage,x=None,y=None,List=None):
    # define the output
    all_out = []


    # Get the parameters from the header
    GeoPar=DrizGeoPars(drizimage,origimage)

    iraf.unlearn('wtraxy')

    # Use wtraxy, along with all the parameters specified above, to
    # transform to the output image
    iraf.wtraxy.nxin=GeoPar.nxin
    iraf.wtraxy.nyin=GeoPar.nyin
    iraf.wtraxy.nxout=GeoPar.nxout
    iraf.wtraxy.nyout=GeoPar.nyout
    iraf.wtraxy.scale=GeoPar.scale
    iraf.wtraxy.xsh=GeoPar.xsh
    iraf.wtraxy.ysh=GeoPar.ysh
    iraf.wtraxy.rot=GeoPar.rot
    iraf.wtraxy.coeffs=GeoPar.coeffs
    iraf.wtraxy.shft_un=GeoPar.shft_un
    iraf.wtraxy.shft_fr=GeoPar.shft_fr
    iraf.wtraxy.align=GeoPar.align
    iraf.wtraxy.lam=GeoPar.lam
    iraf.wtraxy.xgeoim=GeoPar.xgeoim
    iraf.wtraxy.ygeoim=GeoPar.ygeoim
    iraf.wtraxy.geomod='user'

    if List != None:

            x=1.0
            y=1.0

            iraf.wtraxy.xylist=List

            str=iraf.wtraxy(x,y,mode='h',Stdout=1)

            all_out.append("   Xin         Yin         Xout        Yout")

            # Just show the lines of interest
            for line in str:
                if line[0:1] == '!':
                    print(line)
                    sys.exit()

                if line[0:3] == ' Xi':
                    xin = float(line.split()[1])
                    yin = float(line.split()[2])
                    xout = float(line.split()[4])
                    yout = float(line.split()[5])
                    all_out.append("%10.3f %10.3f %10.3f %10.3f" % (xin,yin,xout,yout))

    else:

        iraf.wtraxy.xylist=''

        # Transform and display the result
        print("-Transforming position...")
        str=iraf.wtraxy(x,y,mode='h',Stdout=1)

        # Just show the lines of interest
        for line in str:

            if line[0:1] == '!':
                all_out.appen(line)

            if line[0:3] == ' Xi':
                all_out.appen(line)

    # return the output
    return all_out

def b(drizimage,origimage,x=None,y=None,List=None):
    # define the output
    all_out = []

    # Get the parameters from the header
    GeoPar=DrizGeoPars(drizimage,origimage)

    iraf.unlearn('wtranback')

    # Use wtranback, along with all the parameters specified above, to
    # transform to the output image
    iraf.wtranback.nxin=GeoPar.nxin
    iraf.wtranback.nyin=GeoPar.nyin
    iraf.wtranback.nxout=GeoPar.nxout
    iraf.wtranback.nyout=GeoPar.nyout
    iraf.wtranback.scale=GeoPar.scale
    iraf.wtranback.xsh=GeoPar.xsh
    iraf.wtranback.ysh=GeoPar.ysh
    iraf.wtranback.rot=GeoPar.rot
    iraf.wtranback.coeffs=GeoPar.coeffs
    iraf.wtranback.shft_un=GeoPar.shft_un
    iraf.wtranback.shft_fr=GeoPar.shft_fr
    iraf.wtranback.align=GeoPar.align
    iraf.wtranback.lam=GeoPar.lam
    iraf.wtranback.xgeoim=GeoPar.xgeoim
    iraf.wtranback.ygeoim=GeoPar.ygeoim
    iraf.wtranback.geomode='user'

    if List != None:
            iraf.wtranback.xylist=List

            x=1.0
            y=1.0

            str=iraf.wtranback(x,y,mode='h',Stdout=1)

            all_out.append("   Xin         Yin         Xout        Yout")

            # Just show the lines of interest
            for line in str:
                if line[0:1] == "!":
                    print(line)
                    sys.exit()

                if line[0:3] == ' Xi':
                    xin = float(line.split()[1])
                    yin = float(line.split()[2])
                    xout = float(line.split()[4])
                    yout = float(line.split()[5])
                    all_out.append("%10.3f %10.3f %10.3f %10.3f" % (xin,yin,xout,yout))

    else:

        iraf.wtranback.xylist=''

        # Transform and display the result
        str=iraf.wtranback(x,y,mode='h',Stdout=1)

        # Just show the lines of interest
        for line in str:
            if line[0:1] == '!':
                all_out.appen(line)
                sys.exit()

            if line[0:3] == ' Xi':
                all_out.append(line)

    # return the output
    return all_out
