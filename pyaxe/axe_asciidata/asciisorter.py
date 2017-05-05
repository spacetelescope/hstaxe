"""

@author: Martin Kuemmel, Jonas Haase
@organization: Space Telescope - European Coordinating Facility (ST-ECF)
@license: Gnu Public Licence
@contact: mkuemmel@eso.org
@since: 2005/09/13

$LastChangedBy: jhaase $
$LastChangedDate: 2006-10-13 18:38:13Z $
$HeadURL: http://astropy.scipy.org/svn/astrolib/trunk/asciidata/Lib/asciisorter.py $
"""
__version__ = "Version 1.0 $LastChangedRevision: 113 $"

class ColumnIndex(object):
    """
    External column index to allow variations in the index

    This class offers an external index to a column. The index
    stores the access sequence in the original column data.
    This allows to modify the column index without moving the
    data in the column itself.
    This will mainly be used in sorting the column.
    """
    def __init__(self, nrows=None):
        """
        Initializes the class
 
        @param nrows: the number of rows in the index
        @type nrows: int
        """
        # set the sort flagg
        self.sorted = 0
        
        if nrows:
             # create the initial index
            self.index_col = range(nrows)
        else:
            # define the index to none
            self.index_col = None


    def __str__(self):
        """
        Small string routine

        Converts the object into a string

        @return: the string representation of te object
        @rtype: string
        """      
        # intialize the string
        bigstring = ''

        # go over all elements
        for index in range(len(self.index_col)):

            # add the string representation
            # for one element
            bigstring += 'Row '+str(index)+': '+str(self.index_col[index])+ '\n'

        # return the string
        return bigstring


    def __getitem__(self, index):
        """
        Get method for the class

        @param index: the index to return
        @type index: int
        """
        return self.index_col[index]


    def create_index(self, nrows):
        """
        Create the index
        
        The method creates of redefines the index. The former
        is important if the class was not intialized with the
        number of columns.
        
        @param nrows: the number of rows in the index
        @type nrows: int
        """
        #create the initial index
        self.index_col = range(nrows)


    def sort(self, sort_col=None, descending=0, ordered=1):
        """
        Implementation of a sort algorithm
        
        The method is a frontend to a sorting algorithm.
        Under this method it is possible to implement
        various sorting algortihms, which are called
        according the specifications of according
        to the user request.
 
        @param sort_col: the first column to sort for
        @type sort_col: []
        @param descending: boolean to fix ascending (=0) or descending (1) sort order
        @type descending: int
        """
        # check whether the index is defined
        if not self.index_col:
            # create the initial index
            self.index_col = range(len(sort_col))

        # check whether ordered sort is requested
        if ordered:
            # use the insertion sort for ordered search,
            # which is rather slow
            if descending:
                self._insertion_sort_desc(self.index_col, sort_col)
            else:
                self._insertion_sort_asc(self.index_col, sort_col)
        else:
            # quick sort is supposed to be rather fast
            self._rand_quick_sort(self.index_col, sort_col, 0, len(sort_col)-1)
            if descending:
                self.index_col.reverse()
            

        # set the sort flagg
        self.sorted = 1
 

    def _quick_sort(self, index_col, sort_col, first, last):
        """
        Implementation of the quick sort algorithm
        
        An implementation of the classical quick sort algorithm.
       
        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        @param first: index to start sorting
        @type first: int
        @param last: index to end sorting
        @type last: int
        """
        # check whether something has to be ordered
        if first < last:

            # rearrange the pivot element
            div_index = self._partition(index_col, sort_col, first, last)

            # sort from start to the pivot elemetn
            self._quick_sort(index_col, sort_col, first  , div_index-1)

            # sort from the pivot to the end 
            self._quick_sort(index_col, sort_col, div_index+1, last  )


    def _partition(self, index_col, sort_col, first, last):
        """
        Rearrange the array using the last element as pivot
        
        The method selects the last element from the sort column
        as pivot. Then it rearranges the sort column and the index
        column such that all indices higher than the pivot have
        higher values and the indices lower than the pivot have lower
        values, respectively.

        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        @param first: index to start sorting
        @type first: int
        @param last: index to end sorting
        @type last: int

        @return: index of the pivot
        @rtype: int
        """
        # get the pivot value
        pivot_value = sort_col[last]

        # intitialize the order index
        i = first - 1

        # go along the whole array
        for index in range(first, last):

            # check whether the current element
            # is smaller than the pivot
            if sort_col[index] <= pivot_value:

                # increase the order index
                i += 1

                # exchange elements in the sort
                # and the index column
                self._exchange_elements(sort_col, index, i)
                self._exchange_elements(index_col, index, i)

        # exchange the pivot in the sort
        # and index column
        self._exchange_elements(sort_col, last, i+1)
        self._exchange_elements(index_col, last, i+1)
            
        # return the order index
        return i+1


    def _exchange_elements(self, array, index1, index2):
        """
        Exchange two elements in a vector

        @param array: the array
        @type array: []
        @param index1: first index to exchange
        @type index1: int
        @param index2: second index to exchange
        @type index2: int
        """
        # exchange the pivot in the index column
        tmp = array[index1]
        array[index1] = array[index2]
        array[index2] = tmp
       

    def _rand_quick_sort(self, index_col, sort_col, first, last):
        """
        Implementation of a randomized quick sort algorithm
        
        An implementation of the classical quick sort algorithm,
        however using a random element and NOT the last element
        as a pivot.
        
        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        @param first: index to start sorting
        @type first: int
        @param last: index to end sorting
        @type last: int
        """
        # check whether something has to be ordered
        if first < last:

            # rearrange the pivot element
            div_index = self._rand_partition(index_col, sort_col, first, last)

            # sort from start to the pivot elemetn
            self._rand_quick_sort(index_col, sort_col, first  , div_index-1)

            # sort from the pivot to the end 
            self._rand_quick_sort(index_col, sort_col, div_index+1, last  )    


    def _rand_partition(self, index_col, sort_col, first, last):
        """
        Rearrange the array using the a random element as pivot

        The method selects a random element from the sort column
        as pivot. Then it rearranges the sort column and the index
        column such that all indices higher than the pivot have
        higher values and the indices lower than the pivot have lower
        values, respectively.

        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        @param first: index to start sorting
        @type first: int
        @param last: index to end sorting
        @type last: int

        @return: index of the pivot
        @rtype: int
        """
        import random
        
        # get a random index
        rand_index = random.randint(first, last)

        # exchange the last element with the
        # one at the random index
        self._exchange_elements(sort_col, rand_index, last)
        self._exchange_elements(index_col, rand_index, last)

        # perform the usual rearranging
        # and return the index of the pivot
        return self._partition(index_col, sort_col, first, last)


    def _insertion_sort_asc(self, index_col, sort_col):
        """
        The ascending insertion sort algorithm

        An implementation of the insertion sort algortihm. This algorithm
        sort in ASCENDING order.
        It is suited for order sequences. The result of previous sortings
        is NOT unnecessarily disrupted.

        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        """
        # go along the array
        for jj in range(1, len(index_col)):
            
            # choose the current element
            sort_key  = sort_col[jj]
            index_key = index_col[jj]

            # select the position before
            # and place the current element
            # into that sub-array
            ii = jj - 1
            while ii > -1 and sort_col[ii] > sort_key:

                # move the element to the left
                sort_col[ii+1]  = sort_col[ii]
                index_col[ii+1] = index_col[ii]
                
                # decrease the index
                ii -= 1

            # place the current element
            sort_col[ii+1]  = sort_key
            index_col[ii+1] = index_key


    def _insertion_sort_desc(self, index_col, sort_col):
        """
        The descending insertion sort algorithm
    
        An implementation of the insertion sort algortihm. This algorithm
        sort in DESCENDING order.
        It is suited for order sequences. The result of previous sortings
        is NOT unnecessarily disrupted.

        @param index_col: the index column
        @type index_col: []
        @param sort_col: the sort column
        @type sort_col: []
        """
        # go along the array
        for jj in range(1, len(index_col)):
            
            # choose the current element
            sort_key  = sort_col[jj]
            index_key = index_col[jj]

            # select the position before
            # and place the current element
            # into that sub-array
            ii = jj - 1
            while ii > -1 and sort_col[ii] < sort_key:

                # move the element to the left
                sort_col[ii+1]  = sort_col[ii]
                index_col[ii+1] = index_col[ii]
                
                # decrease the index
                ii -= 1

            # place the current element
            sort_col[ii+1]  = sort_key
            index_col[ii+1] = index_key


    def deindex(self, array):
        """
        Reorders an array into the index order
        
        The method creates and returns an array in the index
        order, such that:
        new_array[index_array[index]] = old_array[index]
        This method is important to insert an external column
        into an AsciiData object.

        @param array: the input array to de-index
        @type array: []

        @return: the array in index order
        @rtype: []
        """
        # create a new array
        new_array = len(array) * [None]
        
        # go over the array
        for index  in range(len(array)):
            
            # put one element into the right place
            new_array[self.index_col[index]] = array[index]

        # return the array
        return new_array


    def enindex(self, array):
        """
        Reorders an indexed array into the normal order
 
        The method creates and returns an array in the natural
        order, such that:
        new_array[index] = old_array[index_array[index]]
        This method is important give a column of an AsciiData
        column to the outside world.

        @param array: the input array to en-index
        @type array: []

        @return: the array in normal order
        @rtype: []
        """
        # create a new array
        new_array = len(array) * [None]
        
        # go over the array
        for index  in range(len(array)):
            
            # put one element into the right place
            new_array[index] = array[self.index_col[index]]

        # return the array
        return new_array
