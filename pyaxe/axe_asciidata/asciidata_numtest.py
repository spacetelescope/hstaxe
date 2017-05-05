"""
Unittest classes for the asciidata module

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: $
$LastChangedDate: $
$HeadURL: $
"""
__version__ = "Version 1.0 $LastChangedRevision: 113 $"

import unittest
import asciidata, asciifunction
import os, string

class Test_AsciiNumpy(unittest.TestCase):
    """
    A test class for the conversion to numpy
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
        data = """1  20.0 15.0 aaa
 13  10.2   7.0 bb
  1  30.33  7.0 cc
 26  10.44  1.0 dddd"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testBasicsNumpy(self):
        """
        Test the working of the numpy conversion
        """
        # create a numpy object from the first column
        numpy_col = self.tdata[0].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col[0], 1)
        self.assertEqual(numpy_col[1], 13)
        self.assertEqual(numpy_col[2], 1)
        self.assertEqual(numpy_col[3], 26)

        # create a numpy object from the second row
        numpy_col = self.tdata[1].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col[0], 20.0)
        self.assertEqual(numpy_col[1], 10.2)
        self.assertEqual(numpy_col[2], 30.33)
        self.assertEqual(numpy_col[3], 10.44)

        # create a numpy object from the second row
        numpy_col = self.tdata[3].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col[0], 'aaa')
        self.assertEqual(numpy_col[1], 'bb')
        self.assertEqual(numpy_col[2], 'cc')
        self.assertEqual(numpy_col[3], 'dddd')


class Test_AsciiNumpyNone(unittest.TestCase):
    """
    A test class for the conversion to numpy
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
        data = """1  Null 15.0 aaa
 Null  10.2   7.0 bb
  1  Null  7.0 cc
 26  10.44  1.0 Null"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testReturnType(self):
        """
        Test the basic equivalience of the elements
        """
        import numpy
        from numpy import ma

        # create a numpy object from the first column
        numpy_col = self.tdata[0].tonumpy()

        # check against the template
        self.assert_(isinstance(numpy_col, numpy.ma.MaskedArray))

        # create a numpy object from the second column
        numpy_col = self.tdata[1].tonumpy()

         # check against the template
        self.assert_(isinstance(numpy_col, numpy.ma.MaskedArray))

        # create a numpy object from the second column
        numpy_col = self.tdata[3].tonumpy()

        # check against the template
        self.assert_(isinstance(numpy_col, numpy.ma.MaskedArray))


    def testBasicsNumpy(self):
        """
        Test the basic equivalience of the elements
        """
        # create a numpy object from the first column
        numpy_col = self.tdata[0].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col[0], 1)
        self.assertEqual(numpy_col._mask[1], True)
        self.assertEqual(numpy_col[2], 1)
        self.assertEqual(numpy_col[3], 26)

        # create a numpy object from the second row
        numpy_col = self.tdata[1].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col._mask[0], True)
        self.assertEqual(numpy_col[1], 10.2)
        self.assertEqual(numpy_col._mask[2], True)
        self.assertEqual(numpy_col[3], 10.44)

        # create a numpy object from the second row
        numpy_col = self.tdata[3].tonumpy()
        # check the number of columns and rows
        self.assertEqual(numpy_col[0], 'aaa')
        self.assertEqual(numpy_col[1], 'bb')
        self.assertEqual(numpy_col[2], 'cc')
        self.assertEqual(numpy_col._mask[3], True)


class Test_AsciiNumarray(unittest.TestCase):
    """
    A test class for the conversion to numarray
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
        data = """1  20.0 15.0 aaa
 13  10.2   7.0 bb
  1  30.33  7.0 cc
 26  10.44  1.0 dddd"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testBasicsNumarray(self):
        """

        """
        import numarray

        # create a numpy object from the first column
        numarray_col = self.tdata[0].tonumarray()

        # check the type
        self.assertEqual(type(numarray_col), type(numarray.array([1,2], type='Int64')))

        # check the number of columns and rows
        self.assertEqual(numarray_col[0], 1)
        self.assertEqual(numarray_col[1], 13)
        self.assertEqual(numarray_col[2], 1)
        self.assertEqual(numarray_col[3], 26)


        # create a numpy object from the second row
        numarray_col = self.tdata[1].tonumarray()

        # check the type
        self.assertEqual(type(numarray_col), type(numarray.array([1.0,2.0], type='Float64')))

        # check the number of columns and rows
        self.assertEqual(numarray_col[0], 20.0)
        self.assertEqual(numarray_col[1], 10.2)
        self.assertEqual(numarray_col[2], 30.33)
        self.assertEqual(numarray_col[3], 10.44)


        # create a numpy object from the second row
        numarray_col = self.tdata[3].tonumarray()

        # check the type
        self.assertEqual(type(numarray_col), type(numarray.strings.array(['a','b'])))

        # check the number of columns and rows
        self.assertEqual(numarray_col[0], 'aaa')
        self.assertEqual(numarray_col[1], 'bb')
        self.assertEqual(numarray_col[2], 'cc')
        self.assertEqual(numarray_col[3], 'dddd')

if __name__ == '__main__':

    suite = unittest.makeSuite(Test_AsciiNumpy)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiNumpyNone)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiNumarray)
    unittest.TextTestRunner(verbosity=2).run(suite)

