/**
 * File: spce_PET.c
 * Input/Output routines of the PET FITS binary tables containing a list
 * of ap_pixel structures
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"
#include <unistd.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_PET.h"


#define NPETCOL 19


/**
    Allocate and return a new ap_pixel structure with enough
    room for N elements

    @param N number of element
    @return a pointer to a newly allocated ap_pixel structure
    @see ap_pixel

*/
ap_pixel *
alloc_aperture_table (long N)
{
  ap_pixel *table;
  table = (ap_pixel *) malloc ((N + 1) * sizeof (ap_pixel));
  if (table == NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "alloc_aperture_table:" " Could not allocate"
		   " memory for aperture table of size %ld", N);
    }
  return table;
}

/** A function that creates a FITS file containing an empty primary header

	@param filename The name of the file to open
	@param overwrite If set to 1, then any exisiting file is deleted. Nothing is done otherwise.
*/
void
create_PET (char filename[], int overwrite)
{
  fitsfile *output;
  int f_status = 0;

  // Try to open the file
  {
    FILE *in_file;
    in_file = fopen (filename, "r");
    if ((overwrite == 1) && (in_file != NULL))
      {
	//aXe_message (aXe_M_WARN3, __FILE__, __LINE__,
	//	      "create_MEPET: File %s "
	//	      "exits. Overwriting it.", filename);
	fclose (in_file);
	unlink (filename);
      }
    if ((overwrite != 1) && (in_file != NULL))
      {
	fclose (in_file);
	return;
      }
  }
  // Open the file for creating/appending
  fits_create_file (&output, filename, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_MEPET: Could not open" " file: %s",
		   filename);
    }

  // Create empty HDU
  {
    int naxis = 0;
    long naxes[2];
    ffiimg (output, 16, naxis, naxes, &f_status);
  }
  fits_close_file (output, &f_status);
}

/** A function that creates a FITS file containing an empty primary header
    and returns a fitsfile pointer to it
	@param filename The name of the file to open
	@param overwrite If set to 1, then any exisiting file is deleted. Nothing is done otherwise.
*/
fitsfile *create_PET_opened (char filename[], int overwrite)
{
  fitsfile *output;
  int f_status = 0;

  // Try to open the file
  {
    FILE *in_file;
    in_file = fopen (filename, "r");
    if ((overwrite == 1) && (in_file != NULL))
      {
	//aXe_message (aXe_M_WARN3, __FILE__, __LINE__,
	//	      "create_MEPET: File %s "
	//	      "exits. Overwriting it.", filename);
	fclose (in_file);
	unlink (filename);
      }
    if ((overwrite != 1) && (in_file != NULL))
      {
	fclose (in_file);
	return NULL;
      }
  }
  // Open the file for creating/appending
  fits_create_file (&output, filename, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_MEPET: Could not open" " file: %s",
		   filename);
    }

  // Create empty HDU
  {
    int naxis = 0;
    long naxes[2];
    ffiimg (output, 16, naxis, naxes, &f_status);
  }

  return output;
}



/**
	A helper function that returns the number of elements in an ap_p structure

	@param ap_p A pointer to an existing ap_pixel structure
	@return number of elements
*/
long
PET_count_elements (ap_pixel * ap_p)
{
     long N;

     //  Count number of pixels
     N = 0;
     while (ap_p->p_x != -1)
       {
	    ap_p++;
	    N++;
       }


     return N;
}

/**
	A helper function which returns the collumn number of the column matching a
	given name. Function fails if collumn is not found.

	@param input A pointer to an opened FITS binary table
	@param colname The name of the desired column
	@return Number of the column.
*/
int
get_PET_colnum (fitsfile * input, char colname[])
{
  int colnum, f_status = 0;

  fits_get_colnum (input, CASEINSEN, colname, &colnum, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_PET_colnum: Could not find collumn %s ",
		   colname);
    }

  return colnum;
}


/**
 * Function: add_ALL_to_PET
 * This function populate a BINARY table with the ALL the content of an
 * ap_pixel structure. This function is meant  to be use when one wishes
 * to keep adding new extension to a multi-extension file. It need to be passed
 * a fitsfile pointer pointing to an opened FITS file. A new FITS binary table
 * extension is appended, populated with the content od the ap_p table and the
 * fitsfile pointer pointing to this new extenstion is returned.
 *
 * Parameters:
 * @param ap_p   - An existing ap_pixel structure
 * @param ID     - the name to assign to this row
 * @param input  - a pointer to an opened FITS file
 * @param update - if set, then an existing table is updated, a new one is appended otheerwise
 *
 */
void
add_ALL_to_PET (ap_pixel * ap_p, char ID[], fitsfile *input, int update)
{
  int f_status = 0;
  long N = 0;
  int colnum;
  char colname[FLEN_KEYWORD];
  int index = -1;
  int hdunum;

  struct Col_Descr FITSData[] = {
    {"ID", "60A", NULL},
    {"N", "J1", NULL},
    {"P_X", "XXXJ1", "PIXEL"},
    {"P_Y", "XXXJ1", "PIXEL"},
    {"X", "XXXE1", "PIXEL"},
    {"Y", "XXXE1", "PIXEL"},
    {"DIST", "XXXE1", "PIXEL"},
    {"XS", "XXXE1", "PIXEL"},
    {"YS", "XXXE1", "PIXEL"},
    {"DXS", "XXXE1", "PIXEL"},
    {"XI", "XXXE1", "PIXEL"},
    {"LAMBDA", "XXXE1", "ANGSTOM"},
    {"DLAMBDA", "XXXE1", "ANGSTOM"},
    {"COUNT", "XXXE1", NULL},
    {"ERROR", "XXXE1", NULL},
    {"WEIGHT", "XXXE1", NULL},
    {"CONTAM", "XXXE1", NULL},
    {"MODEL", "XXXE1", NULL},
    {"DQ", "XXXI1", NULL}
  };


  N = PET_count_elements (ap_p);

  { /* Begin column set up */
    char *ttype[NPETCOL], *tform[NPETCOL], *tunit[NPETCOL];
    int i;

    /* Prepare column description */
    for (i = 0; i < NPETCOL; i++)
      {
	ttype[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (ttype[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_ALL_to_PET: Memory allocation failed,");
	  }
	tform[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (tform[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_ALL_to_PET: Memory allocation failed,");
	  }
	tunit[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (tunit[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_ALL_to_PET: Memory allocation failed,");
	  }
	if (FITSData[i].ttype != NULL)
	  sprintf (ttype[i], "%s", FITSData[i].ttype);
	if (FITSData[i].tform != NULL)
	  {
	    if (!strncmp ("XXX", FITSData[i].tform, 3))
	      sprintf (tform[i], "%ld%s", N,FITSData[i].tform + 3);
                else
		  sprintf (tform[i], "%s", FITSData[i].tform);
	  }
	if (FITSData[i].tunit != NULL)
	  sprintf (tunit[i], "%s", FITSData[i].tunit);
	else
	  sprintf (tunit[i], "%s", " ");
      }
    if (!update)
      {
	fits_create_tbl (input, BINARY_TBL, 0, NPETCOL, ttype, tform, tunit, ID, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_ALL_to_PET: Could not create new binary table HDU in PET");
	  }
      }
    /* Clean up */
    for (i = 0; i < NPETCOL; i++)
      {
	free (tunit[i]);
	tunit[i] = NULL;
	free (ttype[i]);
	ttype[i] = NULL;
	free (tform[i]);
	tform[i] = NULL;
      }
    /* Get current HDU number */
    fits_get_hdu_num (input, &hdunum);
  } /* End column set up */


  { /* Begin write ID field */
    char **array;
    int colnum;
    array = malloc (sizeof (char *));
    if (array == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array[0] = malloc (60 * sizeof (char));
    if (array[0] == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    sprintf (array[0],"%s", ID);
    colnum = get_PET_colnum (input, "ID");
    fits_write_col (input, TSTRING, colnum, 1, 1, 1, array, &f_status);
    if (f_status) {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ALL_to_PET: Error writing ID field in PET");
    }
    free (array[0]);
    array[0] = NULL;
    free (array);
    array = NULL;
  } /* End write ID field */

  { /* Begin writing number of elements/row in table */
    colnum = get_PET_colnum (input, "N");

    fits_write_col (input, TLONG, colnum, 1, 1, 1, &N, &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_N_to_PET: Could not write N to PET");
      }

  }
  /* If there are no elements in the table to write then exit now */
  if (N==0) return;


  { /* Begin writing P_X and P_Y */
    int *array_px, *array_py;
    int i;
    ap_pixel *pt;

    array_px = malloc (N * sizeof (int));
    if (array_px == NULL)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ALL_to_PET: Out of memory");

    array_py = malloc (N * sizeof (int));
    if (array_py == NULL)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ALL_to_PET: Out of memory");

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_px[i] = pt->p_x;
	array_py[i] = pt->p_y;

	pt++;
      }
    colnum = get_PET_colnum (input, "P_X");
    index = 1;

    fits_write_col (input, TINT, colnum, index, 1, N, array_px, &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Could not write P_X to row %d, column %s (%d) in PET"
		     , index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "P_Y");

    fits_write_col (input, TINT, colnum, index, 1, N, array_py,&f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Could not write P_Y to row %d, collumn %s (%d) in PET"
		     ,index, colname, colnum);
      }

    free (array_px);
    array_px = NULL;
    free (array_py);
    array_py = NULL;

  }

  { /* Begin writing WEIGHT */
    double *array_weight;
    int i;
    ap_pixel *pt;

    array_weight = malloc (N * sizeof (double));
    if (array_weight == NULL)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ALL_to_PET: Out of memory");


    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_weight[i] = pt->weight;
	pt++;
      }
    colnum = get_PET_colnum (input, "WEIGHT");
    index = 1;

    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_weight, &f_status);
    if (f_status)
      {
	for (i = 0; i < N; i++)
	  {
	    fprintf (stdout, "your weight is: %f ",array_weight[i]);
	  }
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Could not write WEIGHT to row %d, column %s (%d) in PET"
		     , index, colname, colnum);
      }

    free (array_weight);
    array_weight = NULL;
  }

  /* Writing X, Y */
  {
    double *array_x, *array_y;
    int i;
    ap_pixel *pt;
    array_x = malloc (N * sizeof (double));
    if (array_x == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_y = malloc (N * sizeof (double));
    if (array_y == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_x[i] = pt->x;
	array_y[i] = pt->y;
	pt++;
      }
    colnum = get_PET_colnum (input, "X");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_x,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write X to row %d, collumn %s (%d) in PET"
		     , index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "Y");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_y,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write Y to row %d, collumn %s (%d) in PET"
		     , index, colname, colnum);
      }

    free (array_x);
    array_x = NULL;
    free (array_y);
    array_y = NULL;
  }

  /* Writing ,DIST, DXS, XS, YS, XI */
  {
    double *array_xs, *array_ys, *array_dxs, *array_xi, *array_dist;
    int i;
    ap_pixel *pt;
    array_xs = malloc (N * sizeof (double));
    if (array_xs == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_ys = malloc (N * sizeof (double));
    if (array_ys == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_dxs = malloc (N * sizeof (double));
    if (array_dxs == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_xi = malloc (N * sizeof (double));
    if (array_xi == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_dist = malloc (N * sizeof (double));
    if (array_dist == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_xs[i] = pt->xs;
	array_ys[i] = pt->ys;
	array_dxs[i] = pt->dxs;
	array_xi[i] = pt->xi;
	array_dist[i] = pt->dist;
	pt++;
      }
    colnum = get_PET_colnum (input, "XS");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_xs,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write XS to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "YS");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_ys,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write YS to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "DXS");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_dxs,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write DXS to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "XI");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_xi,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write XI to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "DIST");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_dist,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write DIST to row %d, collumn %s (%d) in PET ",
		     index, colname, colnum);
      }

    free (array_xs);
    array_xs = NULL;
    free (array_ys);
    array_ys = NULL;
    free (array_dxs);
    array_dxs = NULL;
    free (array_xi);
    array_xi = NULL;
    free (array_dist);
    array_dist = NULL;
  }

  /* Writing LAMBDA */
  {
    double *array_lambda;
    int i;
    ap_pixel *pt;
    array_lambda = malloc (N * sizeof (double));
    if (array_lambda == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_lambda[i] = pt->lambda;
	pt++;
      }
    colnum = get_PET_colnum (input, "LAMBDA");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_lambda,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     " Could not write LAMBDA to row %d, collumn %s (%d) in PET ",
		     index, colname, colnum);
      }
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_lambda[i] = pt->dlambda;
	pt++;
      }
    colnum = get_PET_colnum (input, "DLAMBDA");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_lambda,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     " Could not write DLAMBDA to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    free (array_lambda);
    array_lambda = NULL;
  }

  /* Writing CONTAM, MODEL */
  {
    double *array_contam, *array_model;
    int i;
    ap_pixel *pt;

    // alloc space for CONTAM
    array_contam = malloc (N * sizeof (double));
    if (array_contam == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }

    // alloc space for MODEL
    array_model = malloc (N * sizeof (double));
    if (array_contam == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }

    // transfer the data from the PET to the vector
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_contam[i] = pt->contam;
	array_model[i] = pt->model;
	pt++;
      }

    // store the CONTAM column
    colnum = get_PET_colnum (input, "CONTAM");
    //    fprintf (stdout, "contam column no: %i\n", colnum);
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_contam,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     " Could not write CONTAM to row %d, collumn %s (%d) in PET ",
		     index, colname, colnum);
      }
    free (array_contam);
    array_contam = NULL;

    // store the MODEL column
    colnum = get_PET_colnum (input, "MODEL");
    //    fprintf (stdout, "model column no: %i\n", colnum);
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_model,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     " Could not write MODEL to row %d, collumn %s (%d) in PET ",
		     index, colname, colnum);
      }
    free (array_model);
    array_model = NULL;
  }

  /* Writing DQ */
  {
    int *array_dq;
    int i;
    ap_pixel *pt;
    array_dq = malloc (N * sizeof (int));
    if (array_dq == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_dq[i] = pt->dq;
	pt++;
      }
    colnum = get_PET_colnum (input, "DQ");
    //    fprintf (stdout, "DQ column no: %i\n", colnum);
    fits_write_col (input, TINT, colnum, index, 1, N, array_dq,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     " Could not write DQ to row %d, collumn %s (%d) in PET ",
		     index, colname, colnum);
      }
    free (array_dq);
    array_dq = NULL;
  }


  /* Writing COUNT, ERROR */
  {
    double *array_error, *array_count;
    int i;
    ap_pixel *pt;
    array_error = malloc (N * sizeof (double));
    if (array_error == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    array_count = malloc (N * sizeof (double));
    if (array_count == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: Out of memory");
      }
    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	array_count[i] = pt->count;
	array_error[i] = pt->error;
	pt++;
      }

    colnum = get_PET_colnum (input, "COUNT");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_count,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write COUNT to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }
    colnum = get_PET_colnum (input, "ERROR");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_error,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_ALL_to_PET: "
		     "Could not write ERROR to row %d, collumn %s (%d) in PET",
		     index, colname, colnum);
      }

    free (array_count);
    array_count = NULL;
    free (array_error);
    array_error = NULL;
  }
}

/*
 * Function: get_ALL_from_next_in_PET
 */
ap_pixel *get_ALL_from_next_in_PET(fitsfile *input, int *aperID, int *beamID)
{
  int f_status = 0, hdutype;
  long N;
  int colnum;
  int anynull;
  void *nullval = NULL;
  ap_pixel *ap_p;

  fits_movrel_hdu (input, 1, &hdutype, &f_status);
  if (f_status)
    {
      *aperID = -1;
      *beamID = -1;
        return NULL;
    }

    /* reading the aperture ID header keyword - OBJECTID */
  {
    long tmp;
    char comment[FLEN_COMMENT];
    fits_read_key_lng (input, "OBJECTID", &tmp, comment, &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_ALL_from_next_in_PET: Error getting index keyword OBJECTID");
      }
    *aperID = (int)tmp;
    fits_read_key_lng (input, "BEAMID", &tmp, comment, &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_ALL_from_next_in_PET: Error getting index keyword OBJECTID");
      }
    *beamID = (int)tmp;

  }


  colnum = get_PET_colnum (input, "N");
  fits_read_col (input, TLONG, colnum, 1, 1, 1, &nullval, &N, &anynull,
		 &f_status);
  if (N==0) return NULL;
  ap_p = alloc_aperture_table (N);


  /* Reading P_X and P_Y */
  {
    int *array_px, *array_py;
    int i;
    ap_pixel *pt;

    array_px = (void *) malloc (N * sizeof (int));
    if (array_px == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_PX_PY_from_MPET:" " Out of memory");
      }
    array_py = (void *) malloc (N * sizeof (int));
    if (array_py == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_PX_PY_from_MPET: " " Out of memory");
      }
    colnum = get_PET_colnum (input, "P_X");
    fits_read_col (input, TINT, colnum, 1, 1, N, NULL, array_px,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "P_Y");
    fits_read_col (input, TINT, colnum, 1, 1, N, NULL, array_py,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->p_x = array_px[i];
	pt->p_y = array_py[i];
	pt++;
      }
    pt->p_x = -1;
    pt->p_y = -1;
    free (array_px);
    array_px = NULL;
    free (array_py);
    array_py = NULL;
  }


  /* Reading WEIGHT */
  {
    float *array_weight;
    int i;
    ap_pixel *pt;

    array_weight = (void *) malloc (N * sizeof (float));
    if (array_weight == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_PX_PY_from_MPET:" " Out of memory");
      }

    colnum = get_PET_colnum (input, "WEIGHT");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_weight,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->weight = array_weight[i];
	pt++;
      }

    free (array_weight);
    array_weight = NULL;

  }

  /* Reading X and Y */
  {
    float *array_x, *array_y;
    int i;
    ap_pixel *pt;

    array_x = (void *) malloc (N * sizeof (float));
    if (array_x == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_X_Y_from_MEPET: " " Out of memory");
      }
    array_y = (void *) malloc (N * sizeof (float));
    if (array_y == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_X_Y_from_MEPET: " " Out of memory");
      }
    colnum = get_PET_colnum (input, "X");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_x,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "Y");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_y,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->x = array_x[i];
	pt->y = array_y[i];
	pt++;
      }

    free (array_x);
    array_x = NULL;
    free (array_y);
    array_y = NULL;
  }

  /* Reading XS, YS, DXS, XI and DIST */
  {
    float *array_xs, *array_ys, *array_dxs, *array_xi, *array_dist;
    int i;
    ap_pixel *pt;

    array_xi = (void *) malloc (N * sizeof (float));
    if (array_xi == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_XI_XS_DIST_from_MEPET:" "Out of memory");
      }
    array_xs = (void *) malloc (N * sizeof (float));
    if (array_xs == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_XI_XS_DIST_from_MEPET:" " Out of memory");
      }
    array_ys = (void *) malloc (N * sizeof (float));
    if (array_ys == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_XI_XS_DIST_from_MEPET:" " Out of memory");
      }
    array_dxs = (void *) malloc (N * sizeof (float));
    if (array_dxs == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_XI_XS_DIST_from_MEPET:" " Out of memory");
      }
    array_dist = (void *) malloc (N * sizeof (float));
    if (array_dist == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_XI_XS_DIST_from_MEPET:" " Out of memory");
      }
    colnum = get_PET_colnum (input, "XI");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_xi,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "XS");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_xs,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "YS");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_ys,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "DXS");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_dxs,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "DIST");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_dist,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->xi = array_xi[i];
	pt->xs = array_xs[i];
	pt->ys = array_ys[i];
	pt->dxs = array_dxs[i];
	pt->dist = array_dist[i];
	pt++;
      }

    free (array_xi);
    array_xi = NULL;
    free (array_xs);
    array_xs = NULL;
    free (array_ys);
    array_ys = NULL;
    free (array_dxs);
    array_dxs = NULL;
    free (array_dist);
    array_dist = NULL;
  }

  /* Reading LAMBDA */
  {
    float *array_lambda;
    int i;
    ap_pixel *pt;

    array_lambda = (void *) malloc (N * sizeof (float));
    if (array_lambda == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_LAMBDA_from_MEPET:" " Out of memory");
      }
    colnum = get_PET_colnum (input, "LAMBDA");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_lambda,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->lambda = array_lambda[i];
	pt++;
      }
    colnum = get_PET_colnum (input, "DLAMBDA");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_lambda,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->dlambda = array_lambda[i];
	pt++;
      }
    free (array_lambda);
    array_lambda = NULL;
  }

  /* Reading CONTAM  and MODEL */
  {
    double *array_contam, *array_model;
    int i;
    ap_pixel *pt;

    array_contam = (void *) malloc (N * sizeof (double));
    if (array_contam == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_LAMBDA_from_MEPET:" " Out of memory");
      }
    colnum = get_PET_colnum (input, "CONTAM");
    fits_read_col (input, TDOUBLE, colnum, 1, 1, N, NULL, array_contam,
		   &anynull, &f_status);

    array_model = (void *) malloc (N * sizeof (double));
    if (array_model == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_LAMBDA_from_MEPET:" " Out of memory");
      }
    colnum = get_PET_colnum (input, "MODEL");
    fits_read_col (input, TDOUBLE, colnum, 1, 1, N, NULL, array_model,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->contam = array_contam[i];
	pt->model  = array_model[i];
	pt++;
      }

    free (array_contam);
    array_contam = NULL;
    free (array_model);
    array_model = NULL;
  }

  /* Reading DQ */
  {
    int *array_dq;
    int i;
    ap_pixel *pt;

    array_dq = (void *) malloc (N * sizeof (int));
    if (array_dq == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_LAMBDA_from_MEPET:" " Out of memory");
      }
    colnum = get_PET_colnum (input, "DQ");
    fits_read_col (input, TINT, colnum, 1, 1, N, NULL, array_dq,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->dq = array_dq[i];
	pt++;
      }

    free (array_dq);
    array_dq = NULL;
  }


  /* Reading COUNT and ERROR */
  {
    float *array_count, *array_error;
    int i;
    ap_pixel *pt;

    array_count = (void *) malloc (N * sizeof (float));
    if (array_count == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "Get_COUNT_ERROR_To_Pet:" " Out of memory");
      }
    array_error = (void *) malloc (N * sizeof (float));
    if (array_error == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "Get_COUNT_ERROR_To_Pet:" " Out of memory");
      }

    colnum = get_PET_colnum (input, "COUNT");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_count,
		   &anynull, &f_status);
    colnum = get_PET_colnum (input, "ERROR");
    fits_read_col (input, TFLOAT, colnum, 1, 1, N, NULL, array_error,
		   &anynull, &f_status);

    pt = ap_p;
    for (i = 0; i < N; i++)
      {
	pt->count = array_count[i];
	pt->error = array_error[i];
	pt++;
      }

    free (array_count);
    array_count = NULL;
    free (array_error);
    array_error = NULL;
  }
  return ap_p;
}

/**
    Funtion to display the content of an ap_pixel structure

    @param output a pointer to a stream
    @param ap an ap_pixel

*/
void
fprintf_ap_pixel (FILE * output, ap_pixel ap)
{
  fprintf (output, "p_x,p_y: %d %d\n", ap.p_x, ap.p_y);
  fprintf (output, "x,y: %f %f\n", ap.x, ap.y);
  fprintf (output, "dist: %f xs: %f ys: %f dxs: %f xi: %f \n", ap.dist, ap.xs, ap.ys, ap.dxs, ap.xi);
  fprintf (output, "lambda: %f dlambda: %f\n", ap.lambda, ap.dlambda);
  fprintf (output, "count: %f error: %f weight: %f\n", ap.count, ap.error, ap.weight);
  fprintf (output, "contam: %f\n",ap.contam);
  fprintf (output, "DQ: %ld\n",ap.dq);
}

/**
   Function to display the content of an ap_pixel list

   @param output a pointer to a stream
   @param ap a pointer to a list of ap_pixels

*/
void
fprintf_ap_pixel_list (FILE * output, ap_pixel * ap)
{
  int i = 0;

  while (ap[i].p_x != -1)
    {
      fprintf_ap_pixel (output, ap[i]);
      i++;
    }
}
