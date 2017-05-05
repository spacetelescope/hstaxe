"""
Unittest classes for the asciidata module

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: jhaase $
$LastChangedDate: 2008-02-27 16:45:36 +0100 (Wed, 27 Feb 2008) $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciidata_test.py $
"""
__version__ = "Version 1.1 $LastChangedRevision: 389 $"

import unittest
import asciidata, asciifunction
import os, string

class Test_AsciiData(unittest.TestCase):
    """
    A test class for the asciidata module
    """
    def setUp(self):
        """
        Automatic set up for the class

        set up data used in the tests.
        setUp is called before each test function execution.
        """
            # define the data
        self.data = """ 105.2  323.4   star    20
 102.4  529.0 galaxy    21
 834.1  343.7 galaxy    23"""

        # define a test file
            # delete in case it just exists
        self.testfile = 'test_file.tmp'
        if os.path.isfile(self.testfile):
               os.unlink(self.testfile)

        # open the test file
        tfile = open(self.testfile, 'w')

        # fill data into the test file
        tfile.write(self.data)

        #close the test file
        tfile.close()

        # create the test instance
        self.tdata = asciifunction.open(self.testfile)


    def tearDown(self):
        """
        Automatic destruction after test
        """
        # explicitly destroy important class data
        del self.data
        del self.tdata

        # remove the file
        if os.path.isfile(self.testfile):
                   os.unlink(self.testfile)


    def testSimplePrint(self):
        """
        Simple check against input

        The test transforms the asciidata instance back into
        an ascii string. This string is then compared against
        the input string, which is the basis for the
        asciidata instance.
        """

        # transform the instance back into
        # a string
        out_data = str(self.tdata)

        # compare against the original string
        self.assertEqual(out_data, self.data)


    def testDimension(self):
        """
        Check the number of columns and rows

        This test checks whether the instance reports
        the correct number of rows and columns.
        """

        # check the number of rows
        self.assertEqual(self.tdata.nrows, 3)

        # check the number of columns
        self.assertEqual(self.tdata.ncols, 4)


    def testElements(self):
        """
        Check the individual elements in the table
        """
        # check the first column
        self.assertEqual(self.tdata[0][0], 105.2)
        self.assertEqual(self.tdata[0][1], 102.4)
        self.assertEqual(self.tdata[0][2], 834.1)

        # check the second column
        self.assertEqual(self.tdata[1][0], 323.4)
        self.assertEqual(self.tdata[1][1], 529.0)
        self.assertEqual(self.tdata[1][2], 343.7)

        # check the third column
        self.assertEqual(self.tdata[2][0], '  star')
        self.assertEqual(self.tdata[2][1], 'galaxy')
        self.assertEqual(self.tdata[2][2], 'galaxy')

        # check the fourth column
        self.assertEqual(self.tdata[3][0], 20)
        self.assertEqual(self.tdata[3][1], 21)
        self.assertEqual(self.tdata[3][2], 23)


    def testColumnType(self):
        """
        Check the column types

        The method checks the different column types
        against the types which should be there.
        """

        # just check all column type against
        # their proper types
        self.assertEqual(self.tdata[0].get_type(), type(1.0))
        self.assertEqual(self.tdata[1].get_type(), type(1.0))
        self.assertEqual(self.tdata[2].get_type(), type('test'))
        self.assertEqual(self.tdata[3].get_type(), type(1))


    def testValueInput(self):
        """
        Test the correct change of table values

        Go a representative number of table elements. Change the
        element content and check whether the correct number is
        read out.
        """

        # go over the first two columns
        for c_index in range(2):

            # go over each row
            for r_index in range(self.tdata.nrows):

                # compute an arbitrary float
                number =  c_index * r_index * 1.23456

                # insert the float
                self.tdata[c_index][r_index] = number

                # read and check the number
                self.assertEqual(self.tdata[c_index][r_index], number)

        # go over each row
        for r_index in range(self.tdata.nrows):

            # form a string
            string =  str(c_index * r_index +1)

            # insert the string
            self.tdata[2][r_index] = string

            # read and check the string
            self.assertEqual(self.tdata[2][r_index], string)

        # go over each row
        for r_index in range(self.tdata.nrows):

            # form an integer
            string =  c_index * r_index

            # insert the string
            self.tdata[3][r_index] = string

            # read and check the string
            self.assertEqual(self.tdata[3][r_index], string)


    def testNoneInput(self):
        """
        Test the insertion of None as table input

        The method iterates over each table element. The
        value 'None' is written in each element. Then the
        table element is read and compared to 'None'.
        """
        # go over each column
        for c_index in range(self.tdata.ncols):

            # go over each row
            for r_index in range(self.tdata.nrows):

                # insert None into the element
                self.tdata[c_index][r_index] = None

                # read and check the value against None
                self.assertEqual(self.tdata[c_index][r_index], None)


    def testTypeConversion(self):
        """
        Test the automatic type conversion of table columns

        Change the column type by insering lower type elements
        in to a table column. Check if the column type
        is really adapted.
        """

        # test float --> string
        self.assertEqual(self.tdata[0].get_type(), type(1.0))
        self.tdata[0][0] = 'tostring'
        self.assertEqual(self.tdata[0].get_type(), type('a'))

        # test int --> float
        self.assertEqual(self.tdata[3].get_type(), type(1))
        self.tdata[3][1] = 1.5
        self.assertEqual(self.tdata[3].get_type(), type(1.0))

        # test float --> string
        self.assertEqual(self.tdata[3].get_type(), type(1.0))
        self.tdata[3][1] = 'change again'
        self.assertEqual(self.tdata[3].get_type(), type('a'))


    def testValueConversion(self):
        """
        Test the automatic type conversion of table entries

        Insert a table value into a lower type column.
        Then the type of the value has to adapt to the
        column type.
        """

        # integer element into float column --> float element
        self.tdata[0][0] = 1
        self.assertEqual(self.tdata[0][0], 1.0)

        # integer element into string column --> string element
        self.tdata[2][0] = 1
        self.assertEqual(string.strip(self.tdata[2][0]), '1')

        # float element into string column --> string element
        self.tdata[2][1] = 1.0
        self.assertEqual(string.strip(self.tdata[2][1]), '1.0')


    def testColumnCreation(self):
        """
        Test the creation of a new column

        Create new columns of different type. Inser values.
        Check the inserted values as well as the type of the
        inserted values.
        """

        # store the initial number of columns
        ncols = self.tdata.ncols

        # create a new column with floats
        self.tdata['new_float'][0] = 0.0
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['new_float'][index] = float(index)

        # check the column type
        self.assertEqual(self.tdata['new_float'].get_type(), type(1.0))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['new_float'][index], float(index))
            self.assertEqual(type(self.tdata['new_float'][index]), type(1.0))


        # create a new column with integers
        self.tdata['new_int'][0] = 0
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['new_int'][index] = index

        # check the column type
        self.assertEqual(self.tdata['new_int'].get_type(), type(1))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['new_int'][index], index)
            self.assertEqual(type(self.tdata['new_int'][index]), type(1))

        # create a new column with integers
        self.tdata['new_string'][0] = 'a'
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['new_string'][index] = str(index)

        # check the column type
        self.assertEqual(self.tdata['new_string'].get_type(), type('a'))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['new_string'][index], str(index))
            self.assertEqual(type(self.tdata['new_string'][index]), type('a'))

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols+3)


    def testAppendColumn(self):
        """
        Test to append new columns

        Append columns via the 'append()' method of the class.
        Test the value insertion and the correct type definition
        """
        # store the initial number of columns
        ncols = self.tdata.ncols

        # create a new column with integers
        self.tdata.append('testInt')

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols+1)
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['testInt'][index] = int(index)

        # check the column type
        self.assertEqual(self.tdata['testInt'].get_type(), type(1))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['testInt'][index], int(index))
            self.assertEqual(type(self.tdata['testInt'][index]), type(1))


        # create a new column with floats
        self.tdata.append('testFloat')

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols+2)
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['testFloat'][index] = float(index)

        # check the column type
        self.assertEqual(self.tdata['testFloat'].get_type(), type(1.0))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['testFloat'][index], int(index))
            self.assertEqual(type(self.tdata['testFloat'][index]), type(1.0))


        # create a new column with strings
        self.tdata.append('testString')

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols+3)
        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['testString'][index] = 'a'

        # check the column type
        self.assertEqual(self.tdata['testString'].get_type(), type('a'))

        # check the values and types
        for index in range(self.tdata.nrows):
            self.assertEqual(self.tdata['testString'][index], 'a')
            self.assertEqual(type(self.tdata['testString'][index]), type('a'))


    def testDeleteColumn(self):
        """
        Test to delete columns

        Create a new column, and then check whether
        it is possible to destroy it again.
        """
        # store the initial number of columns
        ncols = self.tdata.ncols

        # create a new column
        self.tdata.append('testColumn')

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols+1)

        # fill values into the new column
        for index in range(self.tdata.nrows):
            self.tdata['testColumn'][index] = float(index)

        # delete the new column
        del self.tdata['testColumn']

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols)

        # delete one of the old columns
        del self.tdata[0]

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols-1)

        # delete one of the old columns
        del self.tdata[0]

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols-2)

        # check the changed column numbers
        self.assertEqual(self.tdata[0][1], 'galaxy')

        # delete one of the old columns
        del self.tdata[-1]

        # check the new number of columns
        self.assertEqual(self.tdata.ncols, ncols-3)


    def testInsertRows(self):
        """
        Test to insert rows

        Several rows are inserted at the end or in the
        middle in between 'old' columns. The changing
        total number of rows is monitored. Also
        a table element is checked.
        """
        # store the initial number of rows
        nrows = self.tdata.nrows

        # insert two rows at the end
        self.tdata.insert(2,-1)

        # check the new number of rows
        self.assertEqual(self.tdata.nrows, nrows+2)

        # insert two rows at the beginning
        self.tdata.insert(2,0)

        # check the new number of rows
        self.assertEqual(self.tdata.nrows, nrows+4)

        # check a table element
        self.assertEqual(self.tdata[2][2], '  star')


    def testDeleteRows(self):
        """
        Test to delete rows

        The method deletes some rows from all the columns.
        The number of rows is tested, also one element is
        veryfied.
        """
        # store the initial number of rows
        nrows = self.tdata.nrows

        # delete two rows at the beginning
        self.tdata.delete(0,2)

        # check the new number of rows
        self.assertEqual(self.tdata.nrows, nrows-2)

        # check a table element
        self.assertEqual(self.tdata[2][0], 'galaxy')

    def testStrip(self):
        """
        test of strip,lstrip,rstrip

        The method inserts emty rows in start and end of table and strips them off again
        """
        # store the initial number of rows
        nrows = self.tdata.nrows
        # insert two rows at the end
        self.tdata.insert(2,-1)
        # insert three rows at the beginning
        self.tdata.insert(3,0)
        # check the new number of rows
        self.assertEqual(self.tdata.nrows, nrows+5)
        # strip the first three off again
        self.tdata.lstrip()
        # check they are gone
        self.assertEqual(self.tdata.nrows, nrows+2)
        # strip the trailing
        self.tdata.rstrip()
        # back to normal
        self.assertEqual(self.tdata.nrows, nrows)
        # once more with strip
        self.tdata.insert(3,-1)
        self.tdata.insert(4,0)
        self.tdata.strip()
        self.assertEqual(self.tdata.nrows, nrows)

    def testFindColumn(self):
        """
            Test to find the column number

        Derive the column number from the column name.
        Do this for default column names as well as for
        fresh columns with arbitrary column names.
        """
        # store the initial number of columns
        ncols = self.tdata.ncols

        # find and test a default column name
        self.assertEqual(self.tdata.find('column1'), 0)

        # find and test a default column name
        self.assertEqual(self.tdata.find('column3'), 2)

        # append a new column, find and test the column number
        self.tdata.append('new_column')
        self.assertEqual(self.tdata.find('new_column'), ncols)

        # create a new column, find and test the column number
        self.tdata['new_column2'][0] = 1.0
        self.assertEqual(self.tdata.find('new_column2'), ncols+1)


    def testWriteTo(self):
        """
        Test the writing to a different file

        Test to write an instance to a file. Read it in
        again and then compare it to the original string
        """

        # define a second test file
        # delete in case it just exists
        tfile = 'test_file2.tmp'
        if os.path.isfile(tfile):
            os.unlink(tfile)

        # write the data to the test file
        self.tdata.writeto(tfile)

        # read the table from the test file
        ttable = asciifunction.open(tfile)

        # delete the test file
        if os.path.isfile(tfile):
            os.unlink(tfile)

        # compare against the original string
        self.assertEqual(str(ttable), self.data)


    def testFlush(self):
        """
        Test writing a modified instance to a file

        Modify the instance, write it back to its
        original file, and test whether the content
        of this file is correct.
        """
        # delete the first row
        self.tdata.delete(0)

        # write the instance back to the original file
        self.tdata.flush()

        # read in the modified table from the old file
        ttable = asciifunction.open(self.testfile)

        # compare against the original string
        self.assertEqual(str(ttable), """ 102.4  529.0 galaxy    21
 834.1  343.7 galaxy    23""")


    def testNewNull(self):
        """
        Test providing a new Null string

        Give a new Null string for the instance.
        Set a few table elements to zero.
        Check the correct writing of the None-entries
        with the new Null string.
        """
        # fill the second row with None entries
        for index in range(self.tdata.ncols):
            self.tdata[index][1] = None

        # change the null string
        self.tdata.newnull('<*>')

        # write the instance to a string, and copmpare
        # to the exspected stirng with the new Null string
        self.assertEqual(str(self.tdata), """ 105.2  323.4   star    20
   <*>    <*>    <*>   <*>
 834.1  343.7 galaxy    23""")

    def testInfo(self):
        """
        Test the info on an instance

        The method tests the info method of the
        asciidata class. The info-string is compared
        against the hardcoded, correct string.
        """
        ref_string = """File:       test_file.tmp
Ncols:      4
Nrows:      3
Delimiter:  None
Null value: ['Null', 'NULL', 'None', '*']
Comment:    #
Column name:        column1
Column type:        <type 'float'>
Column format:      ['% 5.1f', '%6s']
Column null value : ['Null']
Column name:        column2
Column type:        <type 'float'>
Column format:      ['% 5.1f', '%6s']
Column null value : ['Null']
Column name:        column3
Column type:        <type 'str'>
Column format:      ['% 6s', '%6s']
Column null value : ['Null']
Column name:        column4
Column type:        <type 'int'>
Column format:      ['%5i', '%5s']
Column null value : ['Null']

        """

        # check against the correct string
        self.assertEqual(self.tdata.info(), self.tdata.info())

    def testColumnCopy(self):
        """
        Test to copy table columns

        Test the ability to create a deep copy
        of a table column. Check further the
        the number of rows and the preservance
        of the content.
        """

        # make a deep copy of the column
        floatCol = self.tdata[0].copy()

        # change the original values
        for index in range(self.tdata.nrows):
            self.tdata[0][index] = 0.0

        # check the number of rows
        self.assertEqual(floatCol.get_nrows(), 3)

        # check the elements in the copy
        self.assertEqual(floatCol[0], 105.2)
        self.assertEqual(floatCol[1], 102.4)
        self.assertEqual(floatCol[2], 834.1)


        # make a deep copy of the column
        intCol = self.tdata[3].copy()

        # change the original values
        for index in range(self.tdata.nrows):
            self.tdata[3][index] = 0

        # check the number of rows
        self.assertEqual(intCol.get_nrows(), 3)

        # check the elements in the copy
        self.assertEqual(intCol[0], 20)
        self.assertEqual(intCol[1], 21)
        self.assertEqual(intCol[2], 23)

        # make a deep copy of the column
        strCol = self.tdata[2].copy()

        # change the original values
        for index in range(self.tdata.nrows):
            self.tdata[2][index] = 'NN'

        # check the number of rows
        self.assertEqual(strCol.get_nrows(), 3)

        # check the elements in the copy
        self.assertEqual(strCol[0], '  star')
        self.assertEqual(strCol[1], 'galaxy')
        self.assertEqual(strCol[2], 'galaxy')


    def testReformat(self):
        """
        Test changing column formats

        Change the format of some columns and check whether
        they are reformatted correctly.
        """

        # the correct string after reformatting
        reformat_data = """105.20  323.4     star  20
102.40  529.0   galaxy  21
834.10  343.7   galaxy  23"""

        # change the column formats
        self.tdata[0].reformat('%6.2f')
        self.tdata[2].reformat('%8s')
        self.tdata[3].reformat('%3i')

        # check the string against the reference string
        self.assertEqual(str(self.tdata), reformat_data)


    def testRenameCol(self):
        """
        Test the renaming of columns

        Give columns a new name and check whether
        they are stored correctly.
        """
        # rename a column
        self.tdata[2].rename('newname')

        # check whether the new name is found
        self.assertEqual(self.tdata.find('newname'), 2)

        # rename the column again
        self.tdata['newname'].rename('newnewname')

        # check whether it is found again
        self.assertEqual(self.tdata.find('newnewname'), 2)


    def testGetFormat(self):
        """
        Test the retrieval of column formats

        Get the different column formats and check
        whether they are correct
        """

        # go over each column and check the format
        self.assertEqual(self.tdata[0].get_format(), '% 5.1f')
        self.assertEqual(self.tdata[1].get_format(), '% 5.1f')
        self.assertEqual(self.tdata[2].get_format(), '% 6s')
        self.assertEqual(self.tdata[3].get_format(), '%5i')


    def testNoneHeader(self):
        """
        Test the retrieval of column formats

        Get the different column formats and check
        whether they are correct
        """
        self.assertEqual(str(self.tdata.header), '')


    def testResetHeader(self):
        """
        Reset the header
        """

        # add something to the header
        self.tdata.header.append(' A new header entry!')

        # check whether it is still in the header
        self.assertEqual(str(self.tdata.header), '# A new header entry!\n')

        # make a new string representing the whole
        # class and test against it
#        hdata = '# A new header entry!\n' + self.data
        self.tdata.header.reset()
        self.assertEqual(str(self.tdata),self.data)

    def testAppendHeader(self):
        """
        Append something to the header
        """

        # add something to the header
        self.tdata.header.append(' A new header entry!')

        # check whether it is still in the header
        self.assertEqual(str(self.tdata.header), '# A new header entry!\n')

        # make a new string representing the whole
        # class and test against it
        hdata = '# A new header entry!\n' + self.data
        self.assertEqual(str(self.tdata),hdata)


    def testToHTML(self):
        """
        Write the instance as table to an HTML-file
        """
        # check whether you can write to html
        self.assert_(self.tdata.writetohtml())

        # do it again, just to get the filename
        html_file = self.tdata.writetohtml()

        # remove the file
        if os.path.isfile(html_file):
            os.unlink(html_file)

        # give a filename as input
        html_file = 'my_html.html'

        # check whether you can write to a dedicated file
        self.assert_(self.tdata.writetohtml(html_file))

        # remove the file
        if os.path.isfile(html_file):
            os.unlink(html_file)

        # check whether you can give attributes
        self.assert_(self.tdata.writetohtml(html_file,
                                            tr_attr='id="my_tr"',
                                            td_attr='id="my_td"'))

        # remove the file
        if os.path.isfile(html_file):
            os.unlink(html_file)

    def testToLatex(self):
        """
        Write the instance to as table a latex-file
        """
        # check whether you can write to a latex file
        self.assert_(self.tdata.writetolatex())

        # do it again, just to get the filename
        latex_file = self.tdata.writetolatex()

        # remove the file
        if os.path.isfile(latex_file):
            os.unlink(latex_file)

        # give a filename as input
        html_file = 'my_latex.tex'

        # check whether you can write to a dedicated file
        self.assert_(self.tdata.writetohtml(latex_file))

        # remove the file
        if os.path.isfile(latex_file):
            os.unlink(latex_file)


class Test_AsciiDataII(unittest.TestCase):
    """
    A second test class for the asciidata module
    """
    def tmpFileFromString(self, data):
        """
        Store a string into a temporary file

        The method creates a named temporary file and writes
        a string given on input into it.
        The file reference to the temporary file is returned
        for further use of it.
        """
        import tempfile

        # create an open test file
        tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        tfile.write(data)
        tfile.flush()

        # return the file reference
        return tfile


    def testNullDefault(self):
        """
        Test the default for 'None' in the input

        As default there exist the strings '*', 'None',
        'Null' and 'NULL' as markers for entries with
        a None-value.
        """

        # define the data
        data = """ *  323.4   star    20
 102.4  Null  galaxy    21
 834.1  343.7 NULL    None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # test the various representations of the
        # None element such as '*', 'None', ...
        self.assertEqual(tdata[0][0],None)
        self.assertEqual(tdata[1][1],None)
        self.assertEqual(tdata[2][2],None)
        self.assertEqual(tdata[3][2],None)

        # check one of the column types
        self.assertEqual(tdata[0].get_type(), type(1.0))


    def testNullInput(self):
        """
        Test a non-default None value

        The string can be given which should be interpreted
        as a None entry in the table. This is tested
        By loading in a table with a non-default
        None value.
        """

        # define the data
        data = """ !!  323.4   star    20
 102.4  123.5  galaxy    !!
 834.1  343.7 comet    25"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance, specifying
        # the null value
        tdata = asciifunction.open(testfile, null='!!')

        # test the correct loading of the None element
        self.assertEqual(tdata[0][0],None)
        self.assertEqual(tdata[3][1],None)

        # convert the element back to a string
        out_string = str(tdata)

        # find the first Null string
        first = out_string.find('!!')

        # check that it is found
        self.assert_(first > -1)

        # find the second null string
        second = out_string.find('!!', first+1)

        # check that it is found
        self.assert_(second > -1)

        # find the second null string
        third = out_string.find('!!',second+1)

        # now it should not be in any more
        self.failIf(third > -1)


    def testCommentDefault(self):
        """
        Test the default comment string

        Test whether lines marked with the default
        comment string are recognized correctly.
        Also check whether comments at the beginning
        of the file are inserted into the header,
        and ignored within the file
        """

        # define the data
        data = """# This is a comment
#
     *  323.4   star    20
# 102.4  Null  galaxy    21
 834.1  343.7 NULL    None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # check the number of columns
        # and the number of rows
        self.assertEqual(tdata.nrows, 2)
        self.assertEqual(tdata.ncols, 4)

        # test the various representations of the
        # None element such as '*', 'None', ...
        self.assertEqual(tdata[0][0],None)
        self.assertEqual(tdata[2][1],None)
        self.assertEqual(tdata[3][1],None)

        # check one of the column types
        self.assertEqual(tdata[0].get_type(), type(1.0))

        # check whether it is still in the header
        self.assertEqual(str(tdata.header), '# This is a comment\n#\n')

        # convert the instance back to a string
        out_string = str(tdata)

        # find the header in the string
        first = out_string.find('# This is a comment\n#\n')
        self.assertEqual(first, 0)

        # add something to the header
        tdata.header.append(' A new header entry!')

        # convert the instance back to a string
        out_string = str(tdata)

        # check the new, modified header in the string
        first = out_string.find('# This is a comment\n#\n# A new header entry!\n')
        self.assertEqual(first, 0)

        # check for another string character
        second = out_string.find('#', 26)

        # now it should not be in any more
        self.failIf(second > -1)


    def testNewComment(self):
        """
        Test the changing of the comment string

        The comment string is changed to a different
        value. Then the table is converted to a string.
        The correct representation of the new comment
        in the string is checked.
        """

        # define the data
        data = """# This is a comment
#
     *  323.4   star    20
# 102.4  Null  galaxy    21
 834.1  343.7 NULL    None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # check the number of columns
        # and the number of rows
        self.assertEqual(tdata.nrows, 2)

        # change the comment strng to '?'
        tdata.newcomment_char('?')

        # check whether it is still in the header
        self.assertEqual(str(tdata.header), '? This is a comment\n?\n')

        # convert the instance back to a string
        out_string = str(tdata)

        # find the header in the string
        first = out_string.find('? This is a comment\n?\n')
        self.assertEqual(first, 0)

        # add something to the header
        tdata.header.append('A new header entry!')

        # convert the instance back to a string
        out_string = str(tdata)

        # check the new, modified header in the string
        first = out_string.find('? This is a comment\n?\n?A new header entry!\n')
        self.assertEqual(first, 0)

        # check for another string character
        second = out_string.find('?', 26)

        # now it should not be in any more
        self.failIf(second > -1)


    def testCommentInput(self):
        """
        Test an input comment string

        Test whether lines marked with a non-default
        comment string are recognized correctly.
        Also check whether comments at the beginning
        of the file are inserted into the header,
        and ignored within the file
        """

        # define the data
        data = """@ This is a comment
@
     *  323.4   star    20
@     102.4  Null  galaxy    21
 834.1  343.7 NULL    None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile, comment_char='@')

        # check the number of columns
        # and the number of rows
        self.assertEqual(tdata.nrows, 2)
        self.assertEqual(tdata.ncols, 4)

        # test the various representations of the
        # None element such as '*', 'None', ...
        self.assertEqual(tdata[0][0],None)
        self.assertEqual(tdata[2][1],None)
        self.assertEqual(tdata[3][1],None)

        # check one of the column types
        self.assertEqual(tdata[0].get_type(), type(1.0))

        # check whether it is still in the header
        self.assertEqual(str(tdata.header), '@ This is a comment\n@\n')

        # convert the instance back to a string
        out_string = str(tdata)

        # find the header in the string
        first = out_string.find('@ This is a comment\n@\n')
        self.assertEqual(first, 0)

        # add something to the header
        tdata.header.append('A new header entry!')

        # convert the instance back to a string
        out_string = str(tdata)

        # check the new, modified header in the string
        first = out_string.find('@ This is a comment\n@\n@A new header entry!\n')
        self.assertEqual(first, 0)

        # check for another string character
        second = out_string.find('@', 26)

        # now it should not be in any more
        self.failIf(second > -1)


    def testDelimiterInput(self):
        """
        Test a non-default delimiter

        Test whether data given with a non-default delimiter
        is loaded correctly. Also assure that the non-default
        delimiter is written out OK.
        """

        # define the data
        data = """# This is a comment
#
     * | 323.4 | star |  20
# 102.4 | Null  |galaxy |  21
 834.1  |343.7 |NULL    | None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile, delimiter='|')

        # check the number of columns
        # and the number of rows
        self.assertEqual(tdata.nrows, 2)
        self.assertEqual(tdata.ncols, 4)

        # check one of the column types
        self.assertEqual(tdata[0].get_type(), type(1.0))

        # check whether it is still in the header
        self.assertEqual(str(tdata.header), '# This is a comment\n#\n')

        # convert the instance back to a string
        out_string = str(tdata)

        # find the header in the string
        first = out_string.find('# This is a comment\n#\n')
        self.assertEqual(first, 0)

        # add something to the header
        tdata.header.append('A new header entry!')

        # convert the instance back to a string
        out_string = str(tdata)

        # check the new, modified header in the string
        first = out_string.find('|')
        self.assert_(first > -1)


    def testNewDelimiter(self):
        """
        Test the change of the delimiter

        Read in a table, change the delimiter.
        Then check whether the delimiter appears
        as often as necessary.
        """

        # define the data
        data = """# This is a comment
#
     *  323.4   star    20
 102.4  Null  galaxy    21
 834.1  343.7 NULL    None"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # assign a new delimiter
        tdata.newdelimiter('?')

        # convert the instance back to a string
        out_string = str(tdata)

        # count how often the new delimiter
        # is in the string
        tot_num = out_string.count('?')

        # check the number
        self.assertEqual(tot_num, 9)

    def testIterator(self):
        """
        Test the iterator ofer the columns and elements
        """

        # define the data
        data = """# This is a comment
#
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # initialize the sum
        sum = 0.0

        # iterate over columns
        for columns in tdata:
            # iterate over column elements
            for elem in columns:
                sum += elem

        # check the number
        self.assertEqual(sum, 12.0)


    def testHeader(self):
        """
        Test the header methods
        """

        # define the data
        data = """# This is a comment
#
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

        # check the length of the header
        self.assertEqual(len(tdata.header), 2)

        # check the length of the header
        self.assertEqual(tdata.header[0], ' This is a comment\n')

        tdata.header.reset()
        # check the length of the header
        self.assertEqual(len(tdata.header), 0)

        # put a new, empty header item
        tdata.header.append('')

        # replace the empty item with a 'full' item
        tdata.header[0] = ' Now with some content!!'

        str_table = str(tdata)

        first = str_table.find(' Now with some content!!')
        # check the length of the header
        self.assertEqual(first, 1)

        # delete the only header entry
        del tdata.header[0]
         # check the length of the header
        self.assertEqual(len(tdata.header), 0)

    def testHeaderIterator(self):
        """
        Test the header methods
        """

        # define the data
        data = """# This is the first comment
# This is the second comment
# This is the third comment
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0
1.0  1.0   1.0    1.0"""

        # derive the test file object
        tfile = self.tmpFileFromString(data)

        # derive the name of the test file object
        testfile = tfile.name

        # create the test instance
        tdata = asciifunction.open(testfile)

         # check the length of the header
        self.assertEqual(len(tdata.header), 3)

        # iterate over the header
        for aline in tdata.header:
            # check for non-zero entry
            self.assert_(len(aline) > 0)

class Test_NullData(unittest.TestCase):
    """
    A test class for the NullData class
    """
    def setUp(self):
        """
        Automatic set up for the class

        set up data used in the tests.
        setUp is called before each test function execution.
        """
        # create a table with only null entries
        self.tdata = asciifunction.create(4, 5)

        # fill in some values
        for ii in range(self.tdata.nrows):
            # integers
            self.tdata[0][ii]  = ii
            # floats
            self.tdata[1][ii]  = float(ii)
            # a string
            self.tdata[2][ii]  = str(ii)+'a'
            # another float
            self.tdata[3][ii]  = float(ii)*float(ii)

    def testBasics(self):
        """
        Basic tests on the nulldata-table

        Check the dimension of the table,
        fill in some data and check the
        column format.
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 4)
        self.assertEqual(self.tdata.nrows, 5)

        # check the column types
        self.assertEqual(self.tdata[0].get_type(), type(1))
        self.assertEqual(self.tdata[1].get_type(), type(1.0))
        self.assertEqual(self.tdata[2].get_type(), type('test'))
        self.assertEqual(self.tdata[3].get_type(), type(1.0))

    def testToSEx(self):
        """
        Check the transformation to SExtractor format
        """
        import tempfile

        # create an open test file
        tfile = tempfile.NamedTemporaryFile()

        # change some column names
        self.tdata[0].rename('Seq')
        self.tdata[1].rename('Float')
        self.tdata[2].rename('String')
        self.tdata[3].rename('Float2')

        # change the format
        self.tdata.toSExtractor()

        # write the object to a file
        self.tdata.writeto(tfile.name)

        # read in the file to a new object
        adata = asciifunction.open(tfile.name)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Seq')
        self.assertEqual(cindex, 0)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float')
        self.assertEqual(cindex, 1)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('String')
        self.assertEqual(cindex, 2)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float2')
        self.assertEqual(cindex, 3)

    def testWriteTo(self):
        """
        Check the options for 'writeto()'
        """
        import tempfile

        # create an open test file
        tfile = tempfile.NamedTemporaryFile()

        # change some column names
        self.tdata[0].rename('Seq')
        self.tdata[1].rename('Float')
        self.tdata[2].rename('String')
        self.tdata[3].rename('Float2')

        # write the object to a file
        self.tdata.writeto(tfile.name, colInfo=1)

        # read in the file to a new object
        adata = asciifunction.open(tfile.name)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Seq')
        self.assertEqual(cindex, 0)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float')
        self.assertEqual(cindex, 1)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('String')
        self.assertEqual(cindex, 2)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float2')
        self.assertEqual(cindex, 3)

class Test_StrangeInput(unittest.TestCase):
    """
    A test class for strange input
    """
    def setUp(self):
        """
        Automatic set up for the class

        set up data used in the tests.
        setUp is called before each test function execution.
        """
        # create a table with only null entries
        self.tdata = asciifunction.create(1, 1)

    def testBasics(self):
        """
        Basic tests on the nulldata-table

        Check the dimension of the table,
        fill in some data and check the
        column format.
        """
        # fille a strange data in
        self.tdata[0][0] = '3E6'

        # check that the content is correct
        self.assertEqual(self.tdata[0].get_format(), '% 6.0e')


class Test_SExData(unittest.TestCase):
    """
    A test class for the NullData class
    """
    def setUp(self):
        """
        Automatic set up for the class

        set up data used in the tests.
        setUp is called before each test function execution.
        """
        # create a table with only null entries
        self.tdata = asciifunction.createSEx(4, 5)

        # fill in some values
        for ii in range(self.tdata.nrows):
            # integers
            self.tdata[0][ii]  = ii
            # floats
            self.tdata[1][ii]  = float(ii)
            # a string
            self.tdata[2][ii]  = str(ii)+'a'
            # another float
            self.tdata[3][ii]  = float(ii)*float(ii)

    def testBasics(self):
        """
        Basic tests on the nulldata-table

        Check the dimension of the table,
        fill in some data and check the
        column format.
        """
        # check the number of columns and rows
        self.assertEqual(self.tdata.ncols, 4)
        self.assertEqual(self.tdata.nrows, 5)

        # check the column types
        self.assertEqual(self.tdata[0].get_type(), type(1))
        self.assertEqual(self.tdata[1].get_type(), type(1.0))
        self.assertEqual(self.tdata[2].get_type(), type('test'))
        self.assertEqual(self.tdata[3].get_type(), type(1.0))

    def testOutFormat(self):
        """
        Test the write out format
        """
        import tempfile

        # create an open test file
        tfile = tempfile.NamedTemporaryFile()

        # change some column names
        self.tdata[0].rename('Seq')
        self.tdata[1].rename('Float')
        self.tdata[2].rename('String')
        self.tdata[3].rename('Float2')

        # write the object to a file
        self.tdata.writeto(tfile.name)

        # read in the file to a new object
        adata = asciifunction.open(tfile.name)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Seq')
        self.assertEqual(cindex, 0)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float')
        self.assertEqual(cindex, 1)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('String')
        self.assertEqual(cindex, 2)

        # check whether the new name
        # is loaded, too
        cindex = adata.find('Float2')
        self.assertEqual(cindex, 3)

class Test_AsciiFits(unittest.TestCase):
    """
    A test class for all fits related methods
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
        data = """ 105.2  323.4   star    20
 102.4  529.0 galaxy    21
 834.1  343.7 galaxy    23"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)


    def testBasics(self):
        """
        Basic tests on the nulldata-table

        Check the dimension of the table,
        fill in some data and check the
        column format.
        """
        import string
        import tempfile
        import pyfits

        # write the AsciiData instance to a fits file
        fits_name = self.tdata.writetofits()

        # check whether the fits file exists
        self.assert_(os.path.isfile(fits_name))

        # open the fits file
        fits = pyfits.open(fits_name)

        # check the number of fits HDU's
        self.assertEqual(len(fits), 2)

        # extract the data
        tdata = fits[1].data

        # check the number of columns's
        self.assertEqual(len(tdata.names), 4)

        # check the number of rows
        self.assertEqual(len(tdata.field(0)), 3)

        self.assertAlmostEqual(tdata.field(0)[0], 105.2, 4)
        self.assertAlmostEqual(tdata.field(0)[2], 834.1, 4)
        self.assertAlmostEqual(tdata.field(1)[1], 529.0, 4)
        self.assertEqual(string.strip(tdata.field(2)[0]), 'star')
        self.assertEqual(tdata.field(2)[1], 'galaxy')
        self.assertEqual(tdata.field(3)[0], 20)
        self.assertEqual(tdata.field(3)[1], 21)

        fits.close()


class Test_AsciiSort(unittest.TestCase):
    """
    A test class for the sorting
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
        data = """1.0  20.0 15.0 aa
 13.0  10.0  7.0 bb
  1.0  30.0  7.0 cc
 26.0  10.0  1.0 dd"""

         # define the number of rows
        self.NROWS_1 = 10000
        self.NROWS_2 = 1000

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testThroughoutAscending(self):
        """
        Test sorting a long table with radnoms
        """
        import random

        # create an empty table
        self.tdata = asciifunction.create(2, self.NROWS_1)

        # fill the table with randoms
        for index in range(self.NROWS_1):
            self.tdata[0][index] = random.random()

        # execute the sorting command
        self.tdata.sort(0, descending=0, ordered=0)

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata[0][index] >= self.tdata[0][index-1])


    def testThroughoutDescending(self):
        """
        Test sorting a long table with radnoms
        """
        import random

        # create an empty table
        self.tdata = asciifunction.create(2, self.NROWS_1)

        # fill the table with randoms
        for index in range(self.NROWS_1):
            self.tdata[0][index] = random.random()

        # execute the sorting command
        self.tdata.sort(0, 1, ordered=0)

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata[0][index] <= self.tdata[0][index-1])

    def testAscendingSort(self):
        """
        Test for sorting in ascending order
        """
        # execute the sorting command
        self.tdata.sort(0)

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata[0][index] >= self.tdata[0][index-1])

    def testDescendingSort(self):
        """
        Test for sorting in descending order
        """
         # execute the sorting command
        self.tdata.sort(0, 1)

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata[0][index] <= self.tdata[0][index-1])

    def testColnameSort(self):
        """
        Test for sorting with a column name
        """
        # execute the sorting command
        self.tdata.sort('column2')

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata['column2'][index] >= self.tdata['column2'][index-1])

    def testCharAscSort(self):
        """
        Test for ascending sort on a string column
        """
        # execute the sorting command
        self.tdata.sort('column3')

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata['column3'][index] >= self.tdata['column3'][index-1])

    def testCharDesSort(self):
        """
        Test for descending sort on a string column
        """
        # execute the sorting command
        self.tdata.sort('column3', 1)

        # go along the column, check the sorting
        for index in range(1,self.tdata.nrows):
            self.assert_(self.tdata['column3'][index] <= self.tdata['column3'][index-1])

    def testCorrAscSort(self):
        """
        Test for ascending sort on two columns
        """
        # execute the sorting command on the secondary column
        self.tdata.sort(1, 0, 1)

        # execute the sorting command on the primary column
        self.tdata.sort(0, 0, 1)

        # go along the column
        for index in range(1,self.tdata.nrows):
            value1 = self.tdata[0][index]
            value2 = self.tdata[0][index-1]
            # check the sorting on the primary column
            self.assert_(value1 >= value2)

            # in case of equal values in the primary column
            if self.tdata[0][index] == self.tdata[0][index-1]:
                value1 = self.tdata[1][index]
                value2 = self.tdata[1][index-1]
                # check the sorting on the primary column
                self.assert_(value1 >= value2)


    def testThroughCorrAscSort(self):
        """
        Test for ascending sort on two columns
        """
        import random

        # create an empty table
        self.tdata = asciifunction.create(2, self.NROWS_2)

        # fill the table with randoms
        for index in range(self.NROWS_2):
            # the second column is filled with
            # all random variable
            self.tdata[1][index] = random.random()

            # in the first column, two rows have
            # the identical random variable
            if not index % 2:
                old = random.random()
            self.tdata[0][index] = old

        # execute the sorting command on the secondary column
        self.tdata.sort(1, 0, 1)

        # execute the sorting command on the primary column
        self.tdata.sort(0, 0, 1)

        # go along the column
        for index in range(1,self.tdata.nrows):
            value1 = self.tdata[0][index]
            value2 = self.tdata[0][index-1]
            # check the sorting on the primary column
            self.assert_(value1 >= value2)

            # in case of equal values in the primary column
            if self.tdata[0][index] == self.tdata[0][index-1]:
                value1 = self.tdata[1][index]
                value2 = self.tdata[1][index-1]
                # check the sorting on the primary column
                self.assert_(value1 >= value2)

    def testCorrDesSort(self):
        """
        Test for ascending sort on two columns
        """
        # execute the sorting command on the secondary column
        self.tdata.sort(2, 1, 1)

        # execute the sorting command on the primary column
        self.tdata.sort(1, 1, 1)

        # go along the column
        for index in range(1,self.tdata.nrows):
            value1 = self.tdata[1][index]
            value2 = self.tdata[1][index-1]
            # check the sorting on the primary column
            self.assert_(value1 <= value2)

            # in case of equal values in the primary column
            if self.tdata[1][index] == self.tdata[1][index-1]:
                value1 = self.tdata[2][index]
                value2 = self.tdata[2][index-1]
                # check the sorting on the primary column
                self.assert_(value1 <= value2)

    def testThroughCorrDesSort(self):
        """
        Test for ascending sort on two columns
        """
        import random

        # create an empty table
        self.tdata = asciifunction.create(2, self.NROWS_2)

        # fill the table with randoms
        for index in range(self.NROWS_2):
            # the second column is filled with
            # all random variable
            self.tdata[1][index] = random.random()

            # in the first column, two rows have
            # the identical random variable
            if not index % 2:
                old = random.random()
            self.tdata[0][index] = old

        # execute the sorting command on the secondary column
        self.tdata.sort(1, 1, 1)

        # execute the sorting command on the primary column
        self.tdata.sort(0, 1, 1)

        # go along the column
        for index in range(1,self.tdata.nrows):
            value1 = self.tdata[0][index]
            value2 = self.tdata[0][index-1]
            # check the sorting on the primary column
            self.assert_(value1 <= value2)

            # in case of equal values in the primary column
            if self.tdata[0][index] == self.tdata[0][index-1]:
                value1 = self.tdata[1][index]
                value2 = self.tdata[1][index-1]
                # check the sorting on the primary column
                self.assert_(value1 <= value2)

class Test_AsciiStrip(unittest.TestCase):
    """
    A test class for the stripping functions
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
        data = """2.0 2.0 2.0 2.0
Null Null Null Null
 13.0  10.0  7.0 2.3
  1.0  30.0  7.0 2.3
 1.0  1.0  1.0 1.0"""

        # create an open test file
        self.tfile = tempfile.NamedTemporaryFile()

        # fill data into the test file and flush
        self.tfile.write(data)
        self.tfile.flush()

        # create the test instance
        self.tdata = asciifunction.open(self.tfile.name)

    def testRstrip(self):
       orig_nrows = self.tdata.nrows
       # strips off the last row
       self.tdata.rstrip(1.0)
       self.assertEqual(self.tdata.nrows,orig_nrows-1)
       # table should remain the same
       self.tdata.rstrip(30.0)
       self.assertEqual(self.tdata.nrows,orig_nrows-1)
       # adding two emty rows
       self.tdata.insert(2,self.tdata.nrows)
       self.assertEqual(self.tdata.nrows,orig_nrows+1)
       # and stripping them off
       self.tdata.rstrip()
       self.assertEqual(self.tdata.nrows,orig_nrows-1)
       
    def testLstrip(self):
       orig_nrows = self.tdata.nrows
       # strip off the first row
       self.tdata.lstrip(2.0)
       self.assertEqual(self.tdata.nrows,orig_nrows-1)
       self.tdata.lstrip()
       self.assertEqual(self.tdata.nrows,orig_nrows-2)
       
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_AsciiData)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiDataII)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_NullData)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiFits)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiSort)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Test_AsciiStrip)
    unittest.TextTestRunner(verbosity=2).run(suite)
#    unittest.main()
