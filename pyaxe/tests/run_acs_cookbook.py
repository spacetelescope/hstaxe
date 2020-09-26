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


os.system("cp ./cookbook_data/F775W/F775W.cat .")
os.system("rm ./DATA/jdql01jxq_flc.fits")
os.system("cp ./cookbook_data/G800L/jdql01jxq_flc.fits ./DATA/")

#!cat cookbook.cat

# os.chdir(cwd)


# if os.path.isdir("G800L"):
#     shutil.rmtree("G800L")
# os.mkdir("G800L")

# os.system("cp cookbook_data/G800L/*flc.fits G800L/")
# os.system("cp aXe.lis G800L/")

# if os.path.isdir("F775W"):
#     shutil.rmtree("F775W")
#os.mkdir("F775W")

# os.system("cp cookbook_data/F775W/F775W_drc.fits F775W/")
# os.system("cp cookbook_data/F775W/F775W.cat F775W/")
# os.system("cp cookbook_data/F775W/j*flc.fits F775W/")

#os.system("cp ../cookbook_data/F775W.cat .")

#os.chdir(cwd)

# if os.path.isdir("CONF"):
#     shutil.rmtree("CONF")
# os.mkdir("CONF")


# os.system("cp cookbook_data/CONF/* CONF/")

# os.chdir(cwd)

# if os.path.isdir("DATA"):
#     shutil.rmtree("DATA")
# os.mkdir("DATA")

dimension_info = "183,85,50,50"
# os.system("cp cookbook_data/G800L/*flc.fits DATA/")
# os.system("cp cookbook_data/F775W/*flc.fits DATA/")

dimension_info = "0,0,0,0"
#os.chdir(cwd)
os.chdir("F775W")

axetasks.iolprep(drizzle_image='F775W_drc.fits',
                     input_cat='F775W.cat',
                     dimension_in=dimension_info)
os.system("cp j*_[12].cat ../DATA/")
os.chdir(cwd)

# os.system("cp cookbook_data/aXe1.lis .")
axetasks.axeprep(inlist="aXe.lis",
                     configs="ACS.WFC.CHIP1.Cycle13.5.conf,ACS.WFC.CHIP2.Cycle13.5.conf",
                     backgr=True,
                     backims="ACS.WFC.CHIP1.msky.1.smooth.fits,ACS.WFC.CHIP2.msky.1.smooth.fits",
                     norm=True,
                     mfwhm=3.0)



# os.system("cp cookbook_data/aXe2.lis .")

# axetasks.axeprep(inlist="aXe2.lis",
#                      configs="ACS.WFC.CHIP2.Cycle13.5.conf",
#                      backgr=True,
#                      backims=""ACS.WFC.CHIP1.msky.1.smooth.fits,ACS.WFC.CHIP2.msky.1.smooth.fits",
#                      norm=True,
#                      mfwhm=3.0)

# axetasks.axecore('aXe2.lis',
#                  "ACS.WFC.CHIP2.Cycle13.5.conf",
#                  extrfwhm=4.,
#                  drzfwhm=3.,
#                  backfwhm=0.,
#                  orient=False,
#                  weights=True,
#                  slitless_geom=False,
#                  cont_model='gauss',
#                  sampling='drizzle',
#                  exclude=True)

print(f"starting axecore\n")
axetasks.axecore('aXe.lis',
                 "ACS.WFC.CHIP1.Cycle13.5.conf,ACS.WFC.CHIP2.Cycle13.5.conf",
                 fconfterm="ACS.WFC.CHIP1.msky.1.smooth.fits,ACS.WFC.CHIP2.msky.1.smooth.fits",
                 extrfwhm=4.,
                 drzfwhm=3.,
                 backfwhm=0.,
                 orient=False,
                 weights=True,
                 slitless_geom=False,
                 cont_model='gauss',
                 sampling='drizzle',
                 exclude=True)


# opt_extr=False

# axetasks.drzprep(inlist = "aXe.lis",
#                  configs="ACS.WFC.CHIP1.Cycle13.5.conf,ACS.WFC.CHIP2.Cycle13.5.conf",
#                  back = False,
#                  opt_extr=opt_extr)

opt_extr=True

axetasks.drzprep(inlist = "aXe.lis",
                 configs="ACS.WFC.CHIP1.Cycle13.5.conf,ACS.WFC.CHIP2.Cycle13.5.conf",
                 back = True,
                 opt_extr=opt_extr)

axetasks.axecrr(inlist="aXe1.lis",
    configs="ACS.WFC.CHIP1.Cycle13.5.conf,ACS.WFC.CHIP2.Cycle13.5.conf",
    infwhm = 4.0,
    outfwhm = 3.0,
    back = False,
    driz_separate = 'yes',
    opt_extr=opt_extr,
    )


opt_extr=True

axetasks.axecrr(inlist="aXe2.lis",
    configs="ACS.WFC.CHIP2.Cycle13.5.conf",
    infwhm = 4.0,
    outfwhm = 3.0,
    back = False,
    driz_separate = 'yes',
    opt_extr=opt_extr
    )




