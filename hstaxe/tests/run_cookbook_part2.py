import os
import shutil
from drizzlepac import astrodrizzle
from hstaxe import axetasks

cwd = os.getcwd()
print("We are in %s" % (cwd))

dimension_info = "183,85,50,50"


os.chdir(cwd)

if os.path.isdir("F125W"):
    shutil.rmtree("F125W")

os.mkdir("F125W")

os.system("cp cookbook_data/F125W/*flt.fits F125W/")
os.system("cp cookbook_data/F125W/F125W.lis F125W/")
os.chdir("F125W")

ref = "../G141/G141_drz.fits[1]"
astrodrizzle.AstroDrizzle("@F125W.lis",output="F125W",in_memory=False,skysub="yes",build=True,driz_cr_corr=True,driz_cr=True,final_wcs=True,driz_separate=True,driz_sep_wcs=True,driz_sep_refimage=ref,final_refimage=ref)

os.chdir(cwd)

if os.path.isdir("FLX"):
    shutil.rmtree("FLX")

os.mkdir("FLX")

os.chdir("FLX")

os.system("cp ../F125W/F125W_drz.fits ./")
os.system("cp ../F140W/F140W_drz.fits ./")
os.system("cp ../G141/G141_drz.fits ./")

os.system("cp ../DATA/*flt.fits .")

os.system("cp ../cookbook_data/catalog/seg.fits .")

import glob

dir_images = []
for dir_image in glob.glob("F*drz.fits"):
    print (dir_image)
    fname = dir_image.split("_")[0]
    dir_images.append(dir_image)

from astropy.io import fits
import numpy as np
s = []
for dir_image in dir_images:
    print (dir_image)
    PHOTPLAM = fits.open(dir_image)[0].header["PHOTPLAM"] # Wavelength of filter in A
    PHOTFLAM = fits.open(dir_image)[0].header["PHOTFLAM"] # Wavelength of filter in A
    ABZP = -48.60 - 2.5*np.log10(PHOTFLAM * PHOTPLAM**2/3e8/1e10 )
    ss = "%s %f %f\n" % (dir_image, PHOTPLAM/10., ABZP)
    s.append(ss)
open("cubelist.lis","w").writelines(s)


axetasks.fcubeprep(grism_image = os.path.join("G141_drz.fits"),
				segm_image = os.path.join("seg.fits"),
				filter_info = "cubelist.lis",
				AB_zero = "yes",
				dim_info = dimension_info)

os.system("cp ib6o23*FLX.fits ../DATA/")

os.chdir(cwd)
os.system("cp G141/*flt.fits DATA/")

os.system("cat aXe.lis")
axetasks.axeprep(inlist="aXe.lis",
                     configs="G141.F140W.V4.31.conf",
                     backgr=True,
                     backims="WFC3.IR.G141.sky.V1.0.fits",
                     norm=False,
                     mfwhm=3.0)

print("Finished axeprep\n")
print( "sky: ",fits.open("DATA/ib6o23rsq_flt.fits")[1].header["SKY_CPS"],"e/s")
print( "sky: ",fits.open("DATA/ib6o23ruq_flt.fits")[1].header["SKY_CPS"],"e/s")
print( "sky: ",fits.open("DATA/ib6o23ryq_flt.fits")[1].header["SKY_CPS"],"e/s")
print( "sky: ",fits.open("DATA/ib6o23s0q_flt.fits")[1].header["SKY_CPS"],"e/s")

print("Running axecore...\n")
axetasks.axecore('aXe.lis',
                 "G141.F140W.V4.31.conf",
                 extrfwhm=4.,
                 drzfwhm=3.,
                 backfwhm=0.,
                 orient=False,
                 weights=True,
                 slitless_geom=False,
                 cont_model='fluxcube',
                 sampling='drizzle',
                 exclude=True)

print("Running drzprep....\n")
axetasks.drzprep(inlist = "aXe.lis",
				configs =  "G141.F140W.V4.31.conf",
				back = True)

print("Running axecrr....\n")
axetasks.axecrr(inlist = "aXe.lis",
                configs = "G141.F140W.V4.31.conf",
                infwhm = 4.0,
                outfwhm = 3.0,
                back = False,
                driz_separate = 'yes'
                )



