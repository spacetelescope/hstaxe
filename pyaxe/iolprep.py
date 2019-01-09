
import os

from pyaxe.axesrc import iolmaking
from stsci.tools import parseinput
from stsci.tools import teal

__taskname__ = "iolprep"


def iolprep(drizzle_image='',
            input_cat='',
            dim_info='0,0,0,0'):
    """Function for the aXe task IOLPREP

    drizzle_image: string
        the astrodrizzled image
    input_cat: string
        the master catalog from source extractor
    dim_info: string
        the extra padding dimensions to add to the output images
    """

    # run the main run object;
    # execute the run;
    # delete the object
    iol_maker = iolmaking.IOLMaker(drizzle_image,
                                    input_cat,
                                    dim_info)
    iol_maker.run()
    del iol_maker


def help(file=None):
    helpstr = getHelpAsString(docstring=True)
    if file is None:
        print(helpstr)
    else:
        if os.path.exists(file):
            os.remove(file)
        f = open(file, mode='w')
        f.write(helpstr)
        f.close()


def getHelpAsString(docstring=False):
    """Returns documentation on the 'iolprep' function. Required by TEAL."""

    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir, 'htmlhelp', __taskname__ + '.html')
    helpfile = os.path.join(install_dir, __taskname__ + '.help')
    helpstring = "NA"
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        if os.path.exists(helpfile):
            helpString = teal.getHelpFileAsString(__taskname__, __file__)
        else:
            helpString = 'file://' + htmlfile

    return helpString


iolprep.__doc__ = getHelpAsString(docstring=True)


def run(configobj=None):
    """TEAL interface for the 'iolprep' function"""
    iolprep(configobj['input'])
