/**
 */
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "spc_cfg.h"
#include "aXe_grism.h"
#include "aXe_errors.h"
#include "inima_utils.h"

/**
 * Function: get_axe_inputs
 * The function reads the data in the input image list given as the first
 * parameter and extracts all axe input data from it. If the configuration
 * files are not given in the input image list, the aXe configuration files
 * given in the command line is used. All data is stored
 * and returned in an 'axe_inputs' structure.
 *
 * Parameters:
 * @param inima    - the name of the input image list
 * @param confterm - the term for the aXe configuration file(s)
 *
 * Returns:
 * @return in_list - the resolved list of aXe input
 */
axe_inputs *
get_axe_inputs(char inima[], char confterm[])
{
  char_array *conf_array;

  int nrows, ncols;
  int dirim_col, conf_col, dmag_col, nmult;

  axe_inputs *in_list;

  FILE *Filelist;

  // check whether the Input Image List exists
  Filelist = fopen (inima, "r");
  if (NULL == Filelist)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s\n",
		   inima);
  fclose(Filelist);

  // convert the configuration term
  // into a character array
  conf_array = get_items(confterm, ",");

  // get the number of rows and the number
  // of items per row in the input image list
  get_inlist_basics(inima, &nrows, &ncols);

  // there mus be at least two columns,
  // one with the grism file and one with the IOL
  if (ncols < 2)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "get_axe_inputs: less than 2 items inlist: %s\n", inima);

  // analyze the content of the input image list.
  // find the multiplicity, which is the number of
  // extensions in the grism image, and find the columns
  // for dmag, direct image and configuration files
  get_inima_cols(inima, &dirim_col, &conf_col, &dmag_col, &nmult);


  // check whether the multiplicity from the file
  // is consistent with the multiplicity from
  // the configterm-input
  if (nmult != conf_array->nitems)
    {
      // report the error and go out
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_axe_inputs: Number of configs from input (%i) conflicts with multiplicity (%i) !\n", conf_array->nitems, nmult);
      //      fprintf(stderr, "ERROR: Number of configs from input (%i) conflicts with multiplicity (%i) !\n", conf_array->nitems, nmult);
      //      exit(1);
    }

  // convert the input to an axe_inputs structure
  in_list = extract_axe_inputs(inima, confterm, nrows, nmult,
			       dirim_col, conf_col, dmag_col);

  // free the memory
  free_char_array(conf_array);

  // return the resulting structure
  return in_list;
}


/**
 * Function: copy_basic_input
 * The function fills one entry in an 'axe_items'-structure with
 * data from one line of the Input Image List. The data is just
 * transferred as an item list, column indices are already known.
 *
 * Parameters:
 * @param line_array - the item list of the line
 * @param iol_array  - the item list of the Input Object Lists
 * @param conf_array - the item list of the configuration files (from input)
 * @param arr_index  - the multiplicity index
 * @param dirim_col  - direct column index
 * @param dmag_col   - dmag column index
 * @param axe_item   - the aXe input to fill
 *
 * Returns:
 * @return -
 */
void
copy_basic_input(char_array *line_array, char_array *iol_array,
		 char_array *conf_array, int arr_index, int dirim_col,
		 int dmag_col, basic_input *axe_item)
{
  // define an entry for NULL's
  char null_entry[1] = "";

  // copy the basic input,
  // which is grism file, iol, config file
  strcpy(axe_item->grism_file,  line_array->char_items[0]);
  strcpy(axe_item->iol_cat,     iol_array->char_items[arr_index]);
  strcpy(axe_item->config_file, conf_array->char_items[arr_index]);

  // check whether a direct image exists,
  // copy it over if possible
  if (dirim_col > -1)
    strcpy(axe_item->dirima_file, line_array->char_items[dirim_col]);
  else
    // store the NULL value
    strcpy(axe_item->dirima_file, null_entry);

  // check whether a dmag value exists
  // copy it over if possible
  if (dmag_col > -1)
    axe_item->dmag = atof(line_array->char_items[dmag_col]);
  else
    // store the NULL value
    axe_item->dmag = 0.0;
}

/**
 * Function: extract_axe_inputs
 * The function converts the data given in an Input Image List and in
 * the line input (for aXe configuration files) to a list of aXe inputs.
 * the columns to find optional inputs such as direct images are
 * given as parameters.
 * The filled structure with the aXe inputs is returned.
 *
 * Parameters:
 * @param inima     - the name of the input image list
 * @param confterm  - the configuration files given in the input
 * @param nrows     - the number of rows
 * @param nmult     - the multiplicity
 * @param dirim_col - the direct column index
 * @param conf_col  - the configuration column index
 * @param dmag_col  - the dmag column index
 *
 * Returns:
 * @return in_list - the resolved list of aXe inputs
 */
axe_inputs *
extract_axe_inputs(char inima[], char confterm[], int nrows, int nmult,
		   int dirim_col, int conf_col, int dmag_col)
{
  FILE *flist;

  char Buffer[LINE_LEN_MAX];
  char input[LINE_LEN_MAX];

  axe_inputs *in_list;

  //char null_entry[1] = "";

  char_array *line_array;
  char_array *conf_array;
  char_array *iol_array;

  int index=0;
  int arr_index=0;

  in_list = alloc_axe_inputs(nrows * nmult);
  in_list->nitems = nrows * nmult;

  // check whether the Input Image List exists
  flist = fopen (inima, "r");

  // reset the counter
  index=0;

  // load each line of the Input Image List
  while (fgets (Buffer, LINE_LEN_MAX, flist) != NULL)
    {
      // leave the loop if there is a valid line
      if (!is_valid_inima_line(Buffer))
	continue;

      // trim the input line
      trim(Buffer, input);

      // tokinze the input line
      line_array = get_items(input, " ");

      // tokenize the configuration term
      iol_array = get_items(line_array->char_items[1], ",");

      // check the number of IOL's against
      // the multiplicity
      if (iol_array->nitems != nmult)
	// report the error and go out
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "extract_axe_inputs: Number of IOL's in \"%s\" conflicts with multiplicity (%i) !\n", line_array->char_items[1], nmult);

      // get the aXe configuration files
      // either from the line of from the command input
      if (conf_col > -1)
	{
	  conf_array = get_items(line_array->char_items[conf_col], ",");
	  // check the number of IOL's against
	  // the multiplicity
	  if (conf_array->nitems != nmult)
	    // report the error and go out
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "extract_axe_inputs: Number of Config files in \"%s\" conflicts with multiplicity (%i) !\n", line_array->char_items[conf_col], nmult);

	}
      else
	{
	  conf_array = get_items(confterm, ",");
	}

      // loop over the multiplicity
      for (arr_index=0; arr_index < nmult; arr_index++)
	{
	  // make a line in the inlist structure
	  copy_basic_input(line_array, iol_array, conf_array, arr_index, dirim_col,
			   dmag_col, &in_list->axe_items[index]);

	  // enhance the counter
	  index++;
	}

      free_char_array(line_array);
      free_char_array(iol_array);
      free_char_array(conf_array);
    }

  // close the file
  fclose(flist);

  return in_list;
}


/**
 * Function: get_inima_cols
 *
 * The function determines the column number with the dmag-information, the
 * direct image names and the aXe configuration files. It determines the
 * multiplicity, which is the number of input image lists and aXe configuration
 * files per grism image.
 *
 * Parameters:
 * @param inima     - the name of the input image list
 * @param dirim_col - the direct column index
 * @param conf_col  - the configuration column index
 * @param dmag_col  - the dmag column index
 * @param nmult     - the multiplicity
 *
 * Returns:
 * @return -
 */
void
get_inima_cols(char inima[], int *dirim_col,
	       int *conf_col, int *dmag_col, int *nmult)
{
  char_array *line_array;
  char_array *item_array;

  FILE *flist;
  char Buffer[LINE_LEN_MAX];
  char input[LINE_LEN_MAX];

  int index=0;

  char *t_err;

  *dirim_col = -1;
  *conf_col  = -1;
  *dmag_col  = -1;
  *nmult     = -1;

  float dmag;

  // check whether the Input Image List exists
  flist = fopen (inima, "r");

  // load each line of the Input Image List
  while (fgets (Buffer, LINE_LEN_MAX, flist) != NULL)
    {
      // leave the loop if there is a valid line
      if (is_valid_inima_line(Buffer))
	break;
    }
  // close the file
  fclose(flist);

  // get rid of whitespaces
  trim(Buffer, input);

  // tokenize the line
  line_array = get_items(input, " ");

  // tokenize the term with the input object lists
  item_array = get_items(line_array->char_items[1], ",");

  // store the number of items
  // derived from the IOL term
  *nmult = item_array->nitems;
  free_char_array(item_array);

  // go from item three to the end
  for (index=2; index < line_array->nitems; index++)
    {
      // try to convert the item to a float
      dmag = strtod(line_array->char_items[index], &t_err);

      // check whether the conversion to float was successful
      if (dmag || strcmp(line_array->char_items[index], t_err))
	// set the column for the dmag-values
	*dmag_col = index;
      // check whether there is a ".fits" in the item
      else if (strstr(line_array->char_items[index], ".fits"))
	*dirim_col = index;
      // if nothing fits, it must be the config file column
      else
	*conf_col = index;
    }

  // if there is a column with aXe configuration files
  if (*conf_col > -1)
    {
      // extract the items in an char_array
      item_array = get_items(line_array->char_items[*conf_col], ",");

      // check the multiplicity
      if (item_array->nitems != *nmult)
	// report the error and go out
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_inima_cols: Multiplicity from IOL: %i, from Configs: %i\n", *nmult, item_array->nitems);

      // release the memory
      free_char_array(item_array);
    }

  // release the memory
  free_char_array(line_array);
}


/**
 * Function:get_inlist_basics
 *
 * The function determines the basic format of the input image given as first
 * parameter, which is the number of nonzero, uncommented rows and the number
 * of items in every row. If the number of items changes from row to row,
 * an error is given.
 *
 * Parameters:
 * @param inima - the name of the input image list
 * @param nrows - the number of rows in the input image list
 * @param ncols - the multiplicity
 *
 * Returns:
 * @return -
 */
void
get_inlist_basics(char inima[], int *nrows, int *ncols)
{
  FILE *flist;
  char Buffer[LINE_LEN_MAX];

  int nitems;

  //initialize the counters
  *nrows   =  0;
  *ncols = -1;

  // check whether the Input Image List exists
  flist = fopen (inima, "r");

  // load each line of the Input Image List
  while (fgets (Buffer, LINE_LEN_MAX, flist) != NULL)
    {

      // check whether the line is valid
      nitems = is_valid_inima_line(Buffer);

      // go to the nex line
      // if the line is NOT valid
      if (!nitems)
	continue;

      // increment the row counter
      *nrows = *nrows + 1;

      // check whether the token counter
      // is not yet set
      if (*ncols < 0)
	{
	  // set the token counter
	  *ncols = nitems;
	}
      else
	{
	  // check whether the number
	  // of tokens disagrees with
	  // the previous number
	  if (*ncols != nitems)
	    {
	      // report the error and go out
	      fclose(flist);
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "get_inlist_basics: %i != %i\n", *ncols, nitems);
	    }
	}
    }

  // close the file
  fclose(flist);
}

/**
 * Function: get_nitems
 * The function determines and returns the number of items in
 * the input string when using the second parameter as token.
 *
 * Parameters:
 * @param input_line - the character string to analyze
 * @param tokens     - the tokens to split the string
 *
 * Returns:
 * @return nitems - number of items in the string
 */
int
get_nitems(char input_line[], char tokens[])
{
  char input[LINE_LEN_MAX];
  char *tok;

  int ntoken=0;

  // copy the input in a local string
  strcpy(input, input_line);

  // search for the first token
  tok = strtok(input, tokens);

  // iterate while a token exists
  while(tok)
    {
      // increase the counter
      ntoken++;

      // search for the next token
      tok = strtok(NULL, tokens);
    }

  //
  return ntoken;
}


/**
 *Funktion: get_items
 *
 * The function analyzes the input line. It splits the line using
 * the tokens given as second paramter and stores the items
 * in a 'char_array' structure. this structure is finally returned.
 *
 * Parameters:
 * @param input_line - the character string to analyze
 * @param tokens     - the tokens to split the string
 *
 * Returns:
 * @return an_array - the list of items in the input
 */
char_array *
get_items(char input_line[], char tokens[])
{
  char input[LINE_LEN_MAX];

  char_array *an_array;

  char *tok;

  int nitems = 0;
  int index    ;

  // determine the number of items
  nitems = get_nitems(input_line, tokens);

  // allocate an appropriate structure
  an_array = alloc_char_array(nitems);

  // store the number of items
  an_array->nitems = nitems;

  if (nitems > 0)
    {
      // copy the input in a local string
      strcpy(input, input_line);

      // search for the first token
      tok = strtok(input, tokens);

      // iterate while a token exists
      index = 0;
      while(tok)
	{
	  // copy the token to the array structure
	  strcpy(an_array->char_items[index], tok);
	  index ++;

	  // search for the next token
	  tok = strtok(NULL, tokens);
	}
    }

  // return the strucutre
  return an_array;
}

/**
 * Function: is_valid_inima_line
 * Checks whether a character string is a valid line in an input image
 * list, which means it has zero length and does not start with a '#'.
 * The function returns the number of items for valid lines and 0 otherwise.
 *
 * Parameters:
 * @param input_line - the character string to analyze
 *
 * Returns:
 * @return nitems - the number of items using blank as separator
 */
int
is_valid_inima_line(char input_line[])
{
  char input[LINE_LEN_MAX];
  char trim_input[LINE_LEN_MAX];

  char comment[2] = "#";

  int nitems = 0;

  // copy the input in a local string
  strcpy(input, input_line);

  // check whether there is something at all
  if (strlen(input) == 0)
    return 0;

  // trim the input line
  trim(input, trim_input);

  // check whether the line is commented
  if (!strncmp(trim_input, comment, 1))
    return 0;

  // get the number of items in the line
  nitems = get_nitems(trim_input, " ");

  // return the number of items
  return nitems;
}



/**
 * Function: alloc_char_array
 * The function allocates memory for a character array
 * structure, which is given back to the calling routine.
 *
 * Parameters:
 * @param nitems - the number of items in the charcter array
 *
 * Returns:
 * @return an_array - the allocated character array
 */
char_array *
alloc_char_array(int nitems)
{
  char_array *an_array;

  int i;

  // allocate memory for the return structure
  an_array = (char_array *)malloc(sizeof(char_array));

  // allocate the memory for the pointers to the char data
  an_array->char_items = (char **)malloc(nitems * sizeof(char *));

  // allocate memory for every cha data item
  for (i=0; i < nitems; i++)
    // allocate the item
    an_array->char_items[i] = malloc(MAXCHAR * sizeof(char));

  return an_array;
}


/**
 * Function: free_char_array
 * The function releases the memory of a character array.
 *
 * Parameters:
 * @param an_array - the character array
 *
 * Returns:
 * @return -
 */
void
free_char_array(char_array *an_array)
{
  int i;

  // go over all char data items
  for (i=0; i < an_array->nitems; i++)
    {
      // free the allocated memory
      free(an_array->char_items[i]);
    }

  // free the pointer to the data items
  free(an_array->char_items);

  // fre the whole structure
  free(an_array);

  // set the whole structure to NULL
  an_array = NULL;
}


/**
 * Function: print_char_array
 * The function prints the content of a character
 * array onto the screen.
 *
 * Parameters:
 * @param an_array - the character array
 *
 * Returns:
 * @return -
 */
void
print_char_array(char_array *an_array)
{
  int i;

  for (i=0; i < an_array->nitems; i++)
    {
      fprintf(stdout, "item %i: %s", i+1, an_array->char_items[i]);
    }
  fprintf(stdout, "\n");
}

/**
 *Function: alloc_axe_inputs
 *
 * The function allocates memory for an 'axe_inputs' structure with
 * a defined number of entries of type basic_input.
 *
 * Parameters:
 * @param nitems - the number of items in the aXe input list
 *
 * Returns:
 * @return axe_inputs - the allocated aXe input structure
 */
axe_inputs *
alloc_axe_inputs(int nitems)
{
  axe_inputs *in_list;

  //int i;

  // allocate memory for the return structure
  in_list = (axe_inputs *)malloc(sizeof(axe_inputs));

  // allocate the memory for the pointers to the axe items
  in_list->axe_items = (basic_input *)malloc(nitems * sizeof(basic_input));

  // return the allocated structure
  return in_list;
}

/**
 *Function: free_axe_inputs
 *
 * The function releases memory in an 'axe_inputs' structure.
 *
 * Parameters:
 * @param in_list - the aXe input structure
 *
 * Returns:
 * @return -
 */
void
free_axe_inputs(axe_inputs *in_list)
{
  // free the item list
  free(in_list->axe_items);

  // free the base structure
  free(in_list);

  // set the structure to NULL
  in_list = NULL;
}

/**
 *Function: print_axe_inputs
 *
 * Print the content of the structure onto the screen.
 *
 * Parameters:
 * @param in_list - the aXe input structure
 *
 * Returns:
 * @return -
 */
void
print_axe_inputs(axe_inputs *in_list)
{
  int index;

  // go over all entries
  for (index=0; index < in_list->nitems; index++)
    {
      // print the basic input onto the screen
      fprintf(stdout, "(%i) Grism: %s, IOL: %s, Conf: %s", index + 1,
	      in_list->axe_items[index].grism_file,
	      in_list->axe_items[index].iol_cat,
	      in_list->axe_items[index].config_file);

      // print the direct image name, if there is any
      if (strlen(in_list->axe_items[index].dirima_file) > 0)
	fprintf(stdout, ", Direct: %s", in_list->axe_items[index].dirima_file);

      // print the dmag value
      fprintf(stdout, ", Dmag: %f", in_list->axe_items[index].dmag);

      // print a new line
      fprintf(stdout, "\n");
    }
}

/**
 * Function: trim
 * Trim leading and trailing blank charachters from a string.
 * blank charachters are the one recognised as such by the isspace lib function
 * (Space, horizontal and vertical tab, form feed, carriage return).
 *
 * Parameters:
 * @param inString  - the input string
 * @param outString - the output string
 *
 * Returns:
 * @return -
 */
char * trim(char inString[], char outString[])
{
  int headIndex = 0;
  int tailIndex = 0;
  int outputLength = 0;

  if (inString == NULL) { return NULL; }

  tailIndex = strlen(inString) - 1; // Last char index in the input (excluding the '\0')

  // Get rid of the leading blanks:
  while (isspace(inString[headIndex]) && (headIndex <= tailIndex)) {headIndex++; }

  // Get rid of the trailing blanks:
  while (isspace(inString[tailIndex]) && (tailIndex > headIndex)) {tailIndex--; }

  // compute the total length
  // add one for the null character
  outputLength = tailIndex - headIndex + 1;

  // copy over from input to output
  strncpy(outString, &inString[headIndex], outputLength);

  // add the null character
  outString[outputLength] = '\0';

  // just for fun, return also
  // the result here
  return outString;
}
