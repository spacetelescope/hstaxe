/**
 */
#ifndef _INIMA_UTILS_H
#define _INIMA_UTILS_H

#define MAXCHAR 255
//#define LINE_LEN_MAX    1280 

/*
 * Struct: basic_input
 */
typedef struct
{
  char grism_file[MAXCHAR];
  char config_file[MAXCHAR];
  char dirima_file[MAXCHAR];
  char iol_cat[MAXCHAR];
  double dmag;
}
  basic_input;

/*
 * Struct: axe_inputs
 */
typedef struct
{
  int nitems;
  basic_input *axe_items;
}
  axe_inputs;


/*
 * Struct: char_array
 */
typedef struct
{
  int nitems;
  char **char_items;
}
  char_array;

extern axe_inputs *
get_axe_inputs(char inima[], char confterm[]);

extern void
copy_basic_input(char_array *line_array, char_array *iol_array,
		 char_array *conf_array,  int arr_index, int dirim_col,
		 int dmag_col, basic_input *axe_item);

extern axe_inputs *
extract_axe_inputs(char inima[], char confterm[], int nrows, int nmult, int dirim_col,
		   int conf_col, int dmag_col);

extern void
get_inima_cols(char inima[], int *dirim_col,
	       int *conf_col, int *dmag_col, int *nmult);

extern void
get_inlist_basics(char inima[], int *nrows, int *ncols);

extern char_array *
get_items(char input_line[], char tokens[]);

extern int
is_valid_inima_line(char input_line[]);

extern int
get_ntokens(char input_line[], char tokens[]);

extern char_array *
alloc_char_array(int nitems);

extern void
free_char_array(char_array *an_array);

extern void
print_char_array(char_array *an_array);

extern axe_inputs *
alloc_axe_inputs(int nitems);

extern void
free_axe_inputs(axe_inputs *in_list);

extern char *
trim(char inString[], char outString[]);

extern void
print_axe_inputs(axe_inputs *in_list);

#endif
