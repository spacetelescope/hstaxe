"""
Run the commands from nor's axe cookbook example
but setup for hstaxe
"""
import os
import shutil
from drizzlepac import astrodrizzle
from hstaxe import axetasks

cwd = os.getcwd()
print("We are in %s" % (cwd))
os.chdir(cwd)
if os.path.isdir("G141"):
    shutil.rmtree("G141")
os.mkdir("G141")

os.system("cp cookbook_data/G141/*flt.fits G141/")
os.system("cp cookbook_data/G141/G141.lis G141/")
os.chdir(cwd)
os.chdir("G141")
#!cat G141.lis
astrodrizzle.AstroDrizzle("@G141.lis", output="G141", build=True)
os.chdir(cwd)

if os.path.isdir("F140W"):
    shutil.rmtree("F140W")

os.mkdir("F140W")

os.system("cp cookbook_data/F140W/*flt.fits F140W/")
os.system("cp cookbook_data/F140W/F140W.lis F140W/")

os.chdir(cwd)
os.chdir("F140W")
#!cat F140W.lis
ref = "../G141/G141_drz.fits[1]"

astrodrizzle.AstroDrizzle("@F140W.lis",output="F140W",in_memory=False,skysub="yes",
                          build=True,driz_cr_corr=True,driz_cr=True,final_wcs=True,driz_separate=True,
                          driz_sep_wcs=True,driz_sep_refimage=ref,final_refimage=ref)

os.system("cp ../cookbook_data/cookbook.cat .")
#!cat cookbook.cat

os.chdir(cwd)

if os.path.isdir("CONF"):
    shutil.rmtree("CONF")
os.mkdir("CONF")

os.system("cp cookbook_data/CONF/* CONF/")

dimension_info = "183,85,50,50"
os.system("cp G141/*flt.fits DATA/")
os.system("cp F140W/*flt.fits DATA/")

os.chdir(cwd)
os.chdir("F140W")

axetasks.iolprep(drizzle_image='F140W_drz.fits',
                     input_cat='cookbook.cat',
                     dimension_in=dimension_info)
os.system("cp ib6o23*_1.cat ../DATA/")
os.chdir(cwd)
os.system("cp cookbook_data/aXe.lis .")
#!cat aXe.lis

axetasks.axeprep(inlist="aXe.lis",
                     configs="G141.F140W.V4.31.conf",
                     backgr=True,
                     backims="WFC3.IR.G141.sky.V1.0.fits",
                     norm=False,
                     mfwhm=3.0)

axetasks.axecore('aXe.lis',
                 "G141.F140W.V4.31.conf",
                 extrfwhm=4.,
                 drzfwhm=3.,
                 backfwhm=0.,
                 orient=False,
                 weights=True,
                 slitless_geom=False,
                 cont_model='gauss',
                 sampling='drizzle',
                 exclude=True)
#!ls -altr OUTPUT

