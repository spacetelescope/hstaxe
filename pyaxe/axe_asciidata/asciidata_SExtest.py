"""
Unittest classes for the asciidata module

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: jhaase $
$LastChangedDate: 2007-07-11 08:45:32Z $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciidata_SExtest.py $
"""
__version__ = "Version 1.1 $LastChangedRevision: 234 $"

import unittest
import asciidata, asciifunction
import os, string


class Test_SExtractCat(unittest.TestCase):
    """
    A test class for SExtractor catalogue format
    """
    def setUp(self):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # define the data
        data = """#   1 NUMBER          Running object number
#   2 FLUXERR_ISO     RMS error for isophotal flux                    [count]
#   3 FLUX_APER       Flux vector within fixed circular aperture(s)   [count]
#   6 FLUXERR_APER    RMS error vector for aperture flux(es)          [count]
#   9 MAG_APER        Fixed aperture magnitude vector                 [mag]
#  12 MAGERR_APER     RMS error vector for fixed aperture mag.        [mag]
#  15 FLUX_AUTO       Flux within a Kron-like elliptical aperture     [count]
#  16 FLUXERR_AUTO    RMS error for AUTO flux                         [count]
#  17 X_IMAGE         Object position along x                         [pixel]
#  18 Y_IMAGE         Object position along y                         [pixel]
#  19 FLAGS           Extraction flags
         1      385.918      700.488      9628.46      39347.2      74.2671      149.333      298.248  -7.1135  -9.9589 -11.4873   0.1151   0.0168   0.0082      52533.1      580.708    379.715     72.461   3
         2      732.112       1369.2      6577.91      21262.8      74.5715      149.143      298.172  -7.8412  -9.5452 -10.8191   0.0591   0.0246   0.0152       171543      2014.45    341.365    320.621  19
         3      1020.77      1827.83      7057.36      25838.5      74.7232      149.105      298.153  -8.1548  -9.6216 -11.0307   0.0444   0.0229   0.0125       267764      1844.97    379.148    196.397   3
         4      658.615      6044.48      26384.8      60815.5      74.6474      148.763      298.001  -9.4534 -11.0534 -11.9600   0.0134   0.0061   0.0053       178541      1290.41    367.213    123.803   3
         5      777.319      1330.99      4771.96      16097.7      74.6474      149.105      298.343  -7.8104  -9.1967 -10.5169   0.0609   0.0339   0.0201       131343      1648.34    305.545    307.027   3
         6      1545.78      73556.4       193961       448466      74.1908      149.257      298.267 -12.1666 -13.2193 -14.1293   0.0011   0.0008   0.0007  2.22738e+06      1938.01    258.692    260.341   3
         7      576.549      3410.04      13602.5      34175.4      74.6474      149.181      298.248  -8.8319 -10.3340 -11.3343   0.0238   0.0119   0.0095       111597      1525.78    336.462     97.060   3
         8      565.642      3699.42        13845      48987.7      74.7232      149.219      298.229  -8.9203 -10.3532 -11.7252   0.0219   0.0117   0.0066       129934      917.641    177.377    199.843   3
         9      478.308      1034.29      4120.77      13578.7      74.6474      149.029      298.001  -7.5366  -9.0374 -10.3321   0.0784   0.0393   0.0238      72761.7      1603.15     94.196    131.380   3
        10      200.515      523.002      2356.02      6720.24      74.7232      149.143      298.096  -6.7963  -8.4304  -9.5685   0.1552   0.0687   0.0482        14072      895.465    265.404     46.241   3"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testBasics(self):
        """
        Basic tests on the SExtractor table

        Check the dimension of the table,
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 19)
        self.assertEqual(self.tdata.nrows, 10)

    
    def testSExColumns(self):
        """
        Basic tests on the SExtractor table

        Determine some column numbers by their name,
        then check whether the number is all right.
        """

        # find the column numbers
        numb = self.tdata.find('NUMBER')
        mape = self.tdata.find('MAG_APER')
        feis = self.tdata.find('FLUXERR_ISO')
        xcol = self.tdata.find('X_IMAGE')
        ycol = self.tdata.find('Y_IMAGE')

        # check the column numbers
        self.assertEqual(numb, 0)
        self.assertEqual(feis, 1)
        self.assertEqual(mape, 8)
        self.assertEqual(xcol, 16)
        self.assertEqual(ycol, 17)


    def testSExNewCol(self):
        """
        Basic tests on the SExtractor table

        """
        import tempfile

        # create a new column and fill in values
        for index in range(self.tdata.nrows):
            self.tdata['NEW_COLUMN'][index] = float(index)*float(index)

        # create an new temp file
        new_tmpfile = tempfile.NamedTemporaryFile()

        # write the modified table to the tmp file
        self.tdata.writeto(new_tmpfile.name)

        # read in the tem file
        new_ascii = asciifunction.open(new_tmpfile.name)

        # find the column numbers
        newc = self.tdata.find('NEW_COLUMN')

        # check the column numbers
        self.assertEqual(newc, 19)


class Test_SExtractCatII(unittest.TestCase):
    """
    A test class for SExtractor catalogue format
    """
    def setUp(self):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # define the data
        data = """#   1 NUMBER          Running object number
#   2 FLUXERR_ISO     RMS error for isophotal flux
#   3 FLUX_APER       Flux vector within fixed circular aperture(s)
#   6 FLUXERR_APER    RMS error vector for aperture flux(es)
#   9 MAG_APER        Fixed aperture magnitude vector
#  12 MAGERR_APER     RMS error vector for fixed aperture mag.
#  15 FLUX_AUTO       Flux within a Kron-like elliptical aperture
#  16 FLUXERR_AUTO    RMS error for AUTO flux
#  17 X_IMAGE         Object position along x
#  18 Y_IMAGE         Object position along y
#  19 FLAGS           Extraction flags
         1      385.918      700.488      9628.46      39347.2      74.2671      149.333      298.248  -7.1135  -9.9589 -11.4873   0.1151   0.0168   0.0082      52533.1      580.708    379.715     72.461   3
         2      732.112       1369.2      6577.91      21262.8      74.5715      149.143      298.172  -7.8412  -9.5452 -10.8191   0.0591   0.0246   0.0152       171543      2014.45    341.365    320.621  19
         3      1020.77      1827.83      7057.36      25838.5      74.7232      149.105      298.153  -8.1548  -9.6216 -11.0307   0.0444   0.0229   0.0125       267764      1844.97    379.148    196.397   3
         4      658.615      6044.48      26384.8      60815.5      74.6474      148.763      298.001  -9.4534 -11.0534 -11.9600   0.0134   0.0061   0.0053       178541      1290.41    367.213    123.803   3
         5      777.319      1330.99      4771.96      16097.7      74.6474      149.105      298.343  -7.8104  -9.1967 -10.5169   0.0609   0.0339   0.0201       131343      1648.34    305.545    307.027   3
         6      1545.78      73556.4       193961       448466      74.1908      149.257      298.267 -12.1666 -13.2193 -14.1293   0.0011   0.0008   0.0007  2.22738e+06      1938.01    258.692    260.341   3
         7      576.549      3410.04      13602.5      34175.4      74.6474      149.181      298.248  -8.8319 -10.3340 -11.3343   0.0238   0.0119   0.0095       111597      1525.78    336.462     97.060   3
         8      565.642      3699.42        13845      48987.7      74.7232      149.219      298.229  -8.9203 -10.3532 -11.7252   0.0219   0.0117   0.0066       129934      917.641    177.377    199.843   3
         9      478.308      1034.29      4120.77      13578.7      74.6474      149.029      298.001  -7.5366  -9.0374 -10.3321   0.0784   0.0393   0.0238      72761.7      1603.15     94.196    131.380   3
        10      200.515      523.002      2356.02      6720.24      74.7232      149.143      298.096  -6.7963  -8.4304  -9.5685   0.1552   0.0687   0.0482        14072      895.465    265.404     46.241   3"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testBasics(self):
        """
        Basic tests on the SExtractor table

        Check the dimension of the table,
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 19)
        self.assertEqual(self.tdata.nrows, 10)

    
    def testSExColumns(self):
        """
        Basic tests on the SExtractor table

        Determine some column numbers by their name,
        then check whether the number is all right.
        """

        # find the column numbers
        numb = self.tdata.find('NUMBER')
        mape = self.tdata.find('MAG_APER')
        feis = self.tdata.find('FLUXERR_ISO')
        xcol = self.tdata.find('X_IMAGE')
        ycol = self.tdata.find('Y_IMAGE')

        # check the column numbers
        self.assertEqual(numb, 0)
        self.assertEqual(feis, 1)
        self.assertEqual(mape, 8)
        self.assertEqual(xcol, 16)
        self.assertEqual(ycol, 17)


    def testSExNewCol(self):
        """
        Basic tests on the SExtractor table

        """
        import tempfile

        # create a new column and fill in values
        for index in range(self.tdata.nrows):
            self.tdata['NEW_COLUMN'][index] = float(index)*float(index)

        # create an new temp file
        new_tmpfile = tempfile.NamedTemporaryFile()

        # write the modified table to the tmp file
        self.tdata.writeto(new_tmpfile.name)

        # read in the tem file
        new_ascii = asciifunction.open(new_tmpfile.name)

        # find the column numbers
        newc = self.tdata.find('NEW_COLUMN')

        # check the column numbers
        self.assertEqual(newc, 19)


class Test_SExtractCatIII(unittest.TestCase):
    """
    A test class for SExtractor catalogue format
    """
    def setUp(self):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # define the data
        data = """#   1 NUMBER
#   2 FLUXERR_ISO
#   3 FLUX_APER
#   6 FLUXERR_APER
#   9 MAG_APER
#  12 MAGERR_APER
#  15 FLUX_AUTO
#  16 FLUXERR_AUTO
#  17 X_IMAGE
#  18 Y_IMAGE
#  19 FLAGS
         1      385.918      700.488      9628.46      39347.2      74.2671      149.333      298.248  -7.1135  -9.9589 -11.4873   0.1151   0.0168   0.0082      52533.1      580.708    379.715     72.461   3
         2      732.112       1369.2      6577.91      21262.8      74.5715      149.143      298.172  -7.8412  -9.5452 -10.8191   0.0591   0.0246   0.0152       171543      2014.45    341.365    320.621  19
         3      1020.77      1827.83      7057.36      25838.5      74.7232      149.105      298.153  -8.1548  -9.6216 -11.0307   0.0444   0.0229   0.0125       267764      1844.97    379.148    196.397   3
         4      658.615      6044.48      26384.8      60815.5      74.6474      148.763      298.001  -9.4534 -11.0534 -11.9600   0.0134   0.0061   0.0053       178541      1290.41    367.213    123.803   3
         5      777.319      1330.99      4771.96      16097.7      74.6474      149.105      298.343  -7.8104  -9.1967 -10.5169   0.0609   0.0339   0.0201       131343      1648.34    305.545    307.027   3
         6      1545.78      73556.4       193961       448466      74.1908      149.257      298.267 -12.1666 -13.2193 -14.1293   0.0011   0.0008   0.0007  2.22738e+06      1938.01    258.692    260.341   3
         7      576.549      3410.04      13602.5      34175.4      74.6474      149.181      298.248  -8.8319 -10.3340 -11.3343   0.0238   0.0119   0.0095       111597      1525.78    336.462     97.060   3
         8      565.642      3699.42        13845      48987.7      74.7232      149.219      298.229  -8.9203 -10.3532 -11.7252   0.0219   0.0117   0.0066       129934      917.641    177.377    199.843   3
         9      478.308      1034.29      4120.77      13578.7      74.6474      149.029      298.001  -7.5366  -9.0374 -10.3321   0.0784   0.0393   0.0238      72761.7      1603.15     94.196    131.380   3
        10      200.515      523.002      2356.02      6720.24      74.7232      149.143      298.096  -6.7963  -8.4304  -9.5685   0.1552   0.0687   0.0482        14072      895.465    265.404     46.241   3"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testBasics(self):
        """
        Basic tests on the SExtractor table

        Check the dimension of the table,
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 19)
        self.assertEqual(self.tdata.nrows, 10)

    
    def testSExColumns(self):
        """
        Basic tests on the SExtractor table

        Determine some column numbers by their name,
        then check whether the number is all right.
        """

        # find the column numbers
        numb = self.tdata.find('NUMBER')
        mape = self.tdata.find('MAG_APER')
        feis = self.tdata.find('FLUXERR_ISO')
        xcol = self.tdata.find('X_IMAGE')
        ycol = self.tdata.find('Y_IMAGE')

        # check the column numbers
        self.assertEqual(numb, 0)
        self.assertEqual(feis, 1)
        self.assertEqual(mape, 8)
        self.assertEqual(xcol, 16)
        self.assertEqual(ycol, 17)


    def testSExNewCol(self):
        """
        Basic tests on the SExtractor table

        """
        import tempfile

        # create a new column and fill in values
        for index in range(self.tdata.nrows):
            self.tdata['NEW_COLUMN'][index] = float(index)*float(index)

        # create an new temp file
        new_tmpfile = tempfile.NamedTemporaryFile()

        # write the modified table to the tmp file
        self.tdata.writeto(new_tmpfile.name)

        # read in the tem file
        new_ascii = asciifunction.open(new_tmpfile.name)

        # find the column numbers
        newc = self.tdata.find('NEW_COLUMN')

        # check the column numbers
        self.assertEqual(newc, 19)

class Test_SExtractCatIV(unittest.TestCase):
    """
    A test class for SExtractor catalogue format
    """
    def setUp(self):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # define the data
        data = """#   1 NUMBER          Running object number
#   2 XWIN_IMAGE      Windowed position estimate along x              [pixel]
#   3 ERRX2WIN_IMAGE  Variance of windowed pos along x                [pixel**2]
#   4 YWIN_IMAGE      Windowed position estimate along y              [pixel]
#   5 ERRY2WIN_IMAGE  Variance of windowed pos along y                [pixel**2]
#   6 VIGNET          Pixel data around detection                     [count]
#  10 AWIN_IMAGE      Windowed profile RMS along major axis           [pixel]
#  11 ERRAWIN_IMAGE   RMS windowed pos error along major axis         [pixel]
#  12 BWIN_IMAGE      Windowed profile RMS along minor axis           [pixel]
#  13 ERRBWIN_IMAGE   RMS windowed pos error along minor axis         [pixel]
#  14 XWIN_WORLD      Windowed position along world x axis            [deg]
#  15 YWIN_WORLD      Windowed position along world y axis            [deg]
#  16 AWIN_WORLD      Windowed profile RMS along major axis (world un [deg]
#  17 BWIN_WORLD      Windowed profile RMS along minor axis (world un [deg]
#  18 THETAWIN_IMAGE  Windowed position angle (CCW/x)                 [deg]
#  19 ERRTHETAWIN_IMA Windowed error ellipse pos angle (CCW/x)        [deg]
#  20 THETAWIN_WORLD  Windowed position angle (CCW/world-x)           [deg]
#  21 ERRTHETAWIN_WOR Windowed error ellipse pos. angle (CCW/world-x) [deg]
#  22 MAG_AUTO        Kron-like elliptical aperture magnitude         [mag]
#  23 MAGERR_AUTO     RMS error for AUTO magnitude                    [mag]
#  24 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR                 [mag]
#  25 MAGERR_BEST     RMS error for MAG_BEST                          [mag]
#  26 MAG_APER        Fixed aperture magnitude vector                 [mag]
#  29 MAGERR_APER     RMS error vector for fixed aperture mag.        [mag]
#  32 CLASS_STAR      S/G classifier output
#  33 FLAGS           Extraction flags
#  34 FLUX_APER       Aperture flux
#
# Believe it or not, a comment!!
#
         1    100.523 4.7898588934e-03     11.911 4.7395137260e-03     2.624896     3.205168     3.886365     4.170891     2.783   0.0693     2.078   0.0688 1.2917689651e+02 9.0495749600e-01 0.0001565747 0.0001174314 -87.9 -21.3  31.2 -50.4  -5.3246   0.0416  -5.0007   0.0685  -3.0908  -4.0255  -4.4398   0.0422   0.0353   0.0360  0.00  19 100.0 110.0
         2    100.660 1.5909142490e-02      4.872 9.7755908555e-03     1.135931     1.280251     2.080871     1.316482     7.005   0.1261     3.742   0.0989 1.2917655690e+02 9.0475187327e-01 0.0003907711  0.000218085   6.1   1.1 -61.6 -58.5  -6.4538   0.0214  -5.9503   0.0387  -2.1055  -3.7328  -4.4567   0.1052   0.0476   0.0335  0.00  27 100.0 110.0
         3    131.046 4.5551252175e-03     10.382 4.4816703457e-03    0.7003819      1.52853     2.039873     3.597899     1.965   0.0681     1.714   0.0663 1.2917590593e+02 9.0636513151e-01 0.0001080815 9.955581e-05  31.6 -36.3 -89.1 -33.9  -4.6836   0.0524  -4.0672   0.0752  -2.7938  -3.6714  -4.0605   0.0403   0.0391   0.0429  0.00  17 100.0 110.0
         4    338.959 2.8985337809e-02      4.966 2.1065494058e-02     2.566542     1.813654     2.990188     3.295038    11.439   0.1704     4.337   0.1450 1.2916939683e+02 9.1610428539e-01 0.0006269945 0.0002865449  11.6  -3.9 -65.5 -55.7  -7.1747   0.0173  -6.6982   0.0218  -2.4454  -3.9800  -4.8584   0.0759   0.0372   0.0241  0.00  25 100.0 110.0
         5    166.280 6.5862141396e-03      3.956 5.9103673705e-03    0.6696419    0.8948278    0.9609098     1.225031     1.801   0.0812     1.665   0.0769 1.2917454038e+02 9.0784961355e-01 9.998853e-05 9.544059e-05  36.1   0.4  82.9 -58.7  -4.0865   0.0621  -3.6963   0.1088  -1.6781  -3.1326  -3.5869   0.1229   0.0610   0.0597  0.00  25 100.0 110.0
         6    161.446 9.1287947340e-03      3.515 6.9659061424e-03     0.925561     1.324872     1.902285    0.8756383     2.470   0.0956     1.820   0.0835 1.2917466471e+02 9.0760606082e-01 0.0001384499  0.000103526   5.0   1.3 -61.0 -58.8  -4.2483   0.0577  -3.4544   0.1239  -1.8848  -3.0356  -3.5832   0.0954   0.0671   0.0589  0.00  25 100.0 110.0
         7    116.389 1.1765660120e-02      4.199 8.9645960633e-03    0.4616966     1.656009     1.185955     1.233691     3.396   0.1085     2.358   0.0947 1.2917605196e+02 9.0548069727e-01 0.0001916109 0.0001324951  -1.8   1.5 -56.8 -58.9  -5.0263   0.0443  -4.2159   0.0943  -1.9416  -3.3972  -4.0421   0.1191   0.0635   0.0512  0.00  25 100.0 110.0"""
        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)


    def testBasics(self):
        """
        Basic tests on the SExtractor table

        Check the dimension of the table,
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 35)
        self.assertEqual(self.tdata.nrows, 7)


    def testColumns(self):
        """
        Test for column names

        Check the expanded column names
        """
        # check for the 'normal' column name
        cindex = self.tdata.find('MAG_APER')
        self.assertEqual(cindex, 25)

        # check for the expanded column name
        cindex = self.tdata.find('MAG_APER2')
        self.assertEqual(cindex, 27)

        # check for the expanded column name
        cindex = self.tdata.find('VIGNET3')
        self.assertEqual(cindex, 8)

        # check for the expanded column name
        cindex = self.tdata.find('FLUX_APER1')
        self.assertEqual(cindex, 34)

    def testToPlain(self):
        """
        Test the change to plain format
        """
        self.tdata.toplain()
        
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 35)
        self.assertEqual(self.tdata.nrows, 7)

        # delete the header comments
        self.tdata.header.reset()

        # convert the object to a string
        bstring = str(self.tdata)

        # find a comment in this string
        commpos = bstring.find('#')
        
        # there should be none, since
        # a header comment did not exist
        # and the column info is not printed
        self.assertEqual(commpos, -1)

    def testwriteto(self):
        """
        Test the writeto method
        
        This module tests the writeto method with all
        it's different options.
        """
        import tempfile

        # create an open test file
        temp_file = tempfile.NamedTemporaryFile()

        #----------------------
        # Part I
        #
        
       # write the object to the tmp-file,
        # using the default settings
        self.tdata.writeto(temp_file.name)

        # create an object from the tmp-file
        adata = asciifunction.open(temp_file.name)
        
        # check for the expanded column name
        # in the re-loaded object
        cindex = adata.find('VIGNET3')
        self.assertEqual(cindex, 8)

        # check for the comment 
        # in the re-loaded object
        self.assertEqual(adata.header[1],  " Believe it or not, a comment!!\n")

        #----------------------
        # Part II
        #
        
        # write the object to the tmp-file,
        # using non-default
        self.tdata.writeto(temp_file.name, colInfo=0)

        # create an object from the tmp-file
        adata = asciifunction.open(temp_file.name)
        
        # check for the expanded column name
        # in the re-loaded object
        cindex = adata.find('VIGNET3')
        # now it should not find the column, since
        # the column names were not saved
        self.assertEqual(cindex, -1)

        # the header must still be there
        # and have the old length
        self.assertEqual(len(adata.header), 3)

        #----------------------
        # Part III
        #

        # write the object to the tmp-file,
        # using non-default
        self.tdata.writeto(temp_file.name, colInfo=0, headComment=0)

        # create an object from the tmp-file
        adata = asciifunction.open(temp_file.name)
        
        # check for the expanded column name
        # in the re-loaded object
        cindex = adata.find('VIGNET3')
        # now it should not find the column, since
        # the column names were not saved
        self.assertEqual(cindex, -1)

        # the header has length=0
        # since it did not survive
        self.assertEqual(len(adata.header), 0)
      
class Test_SExtractCatV(unittest.TestCase):
    """
    A test class for SExtractor catalogue format
    """
    def setUp(self):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # define the data
        data = """# 1 Sp
# 2 Mv
# 3 B-V
# 4 U/B
# 5 V-R
# 6 R+I
# 7 Teff
# 8 1.05*bla
# please ignore the silly column names - its just a test after all
O5     -5.7    -0.33    -1.19    -0.15    -0.32    42000    -4.40
O9    -4.5    -0.31    -1.12    -0.15    -0.32    34000    -3.33
B0    -4.0    -0.30    -1.08    -0.13    -0.29    30000    -3.16
B2    -2.45    -0.24    -0.84    -0.10    -0.22    20900    -2.35
B5    -1.2    -0.17    -0.58    -0.06    -0.16    15200    -1.46
B8    -0.25    -0.11    -0.34    -0.02    -0.10    11400    -0.80
A0    0.65    -0.02    -0.02    0.02    -0.02    9790    -0.30
A2    1.3    0.05    0.05    0.08    0.01    9000    -0.20
A5    1.95    0.15    0.10    0.16    0.06    8180    -0.15
F0    2.7    0.30    0.03    0.30    0.17    7300    -0.09
F2    3.6    0.35    0.00    0.35    0.20    7000    -0.11
F5    3.5    0.44    -0.02    0.40    0.24    6650    -0.14
F8    4.0    0.52    0.02    0.47    0.29    6250    -0.16
G0    4.4    0.58    0.06    0.50    0.31    5940    -0.18
G2    4.7    0.63    0.12    0.53    0.33    5790    -0.20
G5    5.1    0.68    0.20    0.54    0.35    5560    -0.21
G8    5.5    0.74    0.30    0.58    0.38    5310    -0.31
K0    5.9    0.81    0.45    0.64    0.42    5150    -0.31
K2    6.4    0.91    0.64    0.74    0.48    4830    -0.42
K5    7.35    1.15    1.08    0.99    0.63    4410    -0.72
M0    8.8    1.40    1.22    1.28    0.91    3840    -1.38
M2    9.9    1.49    1.18    1.50    1.19    3520    -1.89
M5    12.3    1.64    1.24    1.80    1.67    3170    -2.73

"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)


    def testBasics(self):
        """
        Basic tests on the SExtractor table

        Check the dimension of the table,
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 8)
        self.assertEqual(self.tdata.nrows, 23)


    def testColumns(self):
        """
        Test for column names with arithmetic operator in the name

        """
        #
        cindex = self.tdata.find('Teff')
        self.assertEqual(cindex, 6)
         
        cindex = self.tdata.find('B-V')
        self.assertEqual(cindex, 2)

        cindex = self.tdata.find('R+I')
        self.assertEqual(cindex, 5)
        
        cindex = self.tdata.find('1.05*bla')
        self.assertEqual(cindex, 7)
              
        cindex = self.tdata.find('U/B')
        self.assertEqual(cindex, 3)
        
if __name__ == '__main__':

    suite = unittest.makeSuite(Test_SExtractCat)
    unittest.TextTestRunner(verbosity=2).run(suite)

#    unittest.main()
