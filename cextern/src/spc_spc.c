/**
 * Set of functions to handle binned spectra
 */
#include "spc_spc.h"
#include "aXe_utils.h"

/**
 * Function: allocate_spectrum
 * allocates and initializes a spectrum and its subordinate data structures
 * (this is usually called by binning routines, so I've declared it static)
 *
 * Parameters:
 * @param numbin - number of spectrum data points
 *
 * Returns:
 * @return sp  - the allocated spectrum structure
 */
spectrum *
allocate_spectrum (const int numbin)
{
  spectrum *sp;
  int i;
  
  sp = malloc (sizeof (spectrum));
  if (!sp)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  sp->spec_len = numbin;
  sp->warning = 0;
  sp->spec = malloc (numbin * sizeof (spc_entry));
  if (!sp->spec)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  for (i = 0; i < sp->spec_len; i++)
    {
      sp->spec[i].count = 0.;
      sp->spec[i].lambda_max = -1e38;
      sp->spec[i].lambda_min = +1e38;
      sp->spec[i].lambda_mean = GSL_NAN;
      sp->spec[i].error = 0.;
      sp->spec[i].flux = 0.;
      sp->spec[i].ferror = 0.;
      sp->spec[i].contam = -1.0;
      sp->spec[i].dq = 0;
      
    }

  
  return sp;
}


/**
 * Function: free_spectrum
 * Frees a spectrum and its subordinate data structures
 *
 * Parameters:
 * @param spec - a pointer to the spectrum structure.
 */
void
free_spectrum (spectrum * spec)
{
  if (spec == NULL) return;
  free (spec->spec);
  spec->spec = NULL;
  free (spec);
  spec = NULL;
}

/**
 * Function: fprintf_spectrum
 * Function to display the content of a spectrum
 *
 * Parameters:
 * @param output - a pointer to a stream 
 * @param sp     - a pointer to a spectrum
 */
void
fprintf_spectrum (FILE * output, const spectrum * const sp)
{
  int i;
  double sum = 0;
  
  fprintf(output,"# "
	  "lambda_mean count error weight i lambda_max "
	  "lambda_min delta_lambda flux ferror\n");
  for (i = 0; i < sp->spec_len; i++)
    {
      if (isnan (sp->spec[i].lambda_mean))
	{
	  //  continue;
	}
      fprintf (output, "%g %g %g %g %i %g %g %g %g %g\n",
	       sp->spec[i].lambda_mean, sp->spec[i].count,
	       sp->spec[i].error, sp->spec[i].weight, i,
	       sp->spec[i].lambda_max, sp->spec[i].lambda_min,
	       sp->spec[i].lambda_max - sp->spec[i].lambda_min,
	       sp->spec[i].flux, sp->spec[i].ferror);
      
      sum += sp->spec[i].count;
    }
  fprintf (output, "#Spectrum sum: %f\n", sum);
}

/**
 * Function: subtract_spectra
 * Subtract one spectrum form sa second one 
 *
 * Parameters:
 * @param a - the original spectrum
 * @param b - the spectrum to subtract
 *
 * Returns:
 * @param res - the subtracted spectrum
 */
spectrum *
subtract_spectra (spectrum * a, spectrum * b)
{
  int i;
  spectrum *res;
  double alignement;
  
  if ( (a==NULL) || (b==NULL) )
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		   "\nsubtract_spectra: spectra empty.\n");
      return NULL;
    }
  
  if (a->spec_len != b->spec_len)
    {
      fprintf (stderr, "%d %d\n", a->spec_len, b->spec_len);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "\nsubtract_spectra: spectra not aligned.\n");
    }
  
  /* Allocate new spectrum */
  res = allocate_spectrum (a->spec_len);
  
  for (i = 0; i < a->spec_len; i++)
    {
      if (isnan (a->spec[i].lambda_mean))
	{
	  continue;
	}
      if (isnan (b->spec[i].lambda_mean))
	{
	  continue;
	}

      // compute the difference between
      // the wavelength in the object and
      // the background spectral element
      alignement= (a->spec[i].lambda_mean - b->spec[i].lambda_mean)
	/ a->spec[i].lambda_mean;

      // check whether the alignement is OK
      if (alignement > SPCTOL)
	{
	  { 
	    int j;
	    for (j=0; j<a->spec_len; j++)
	      fprintf(stderr,"%f %f\n",a->spec[j].lambda_mean,b->spec[j].lambda_mean);
	  }
	  fprintf (stderr, "lambda: %f %f; misalignement: %f\n", a->spec[i].lambda_mean,
		   b->spec[i].lambda_mean, alignement);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "subtract_spectra: spectra not aligned.");
	}
      res->spec[i].lambda_mean = a->spec[i].lambda_mean;
      res->spec[i].dlambda = a->spec[i].dlambda;
      res->spec[i].count = a->spec[i].count - b->spec[i].count;
      res->spec[i].error =
	sqrt (pow (a->spec[i].error, 2) + pow (b->spec[i].error, 2));
      
      res->spec[i].flux = a->spec[i].flux - b->spec[i].flux;
      res->spec[i].ferror =
	sqrt (pow (a->spec[i].ferror, 2) +
	      pow (b->spec[i].ferror, 2));
      
      res->spec[i].weight = a->spec[i].weight;
      res->spec[i].lambda_max = a->spec[i].lambda_max;
      res->spec[i].lambda_min = a->spec[i].lambda_min;
      res->spec[i].contam = a->spec[i].contam;
      res->spec[i].dq = a->spec[i].dq;
    }
  return res;
}

/**
 * Function: empty_counts_spectrum_copy
 * Allocate and create an empty spectrum on the basis of a
 * template spectrum
 *
 * Parameters:
 * @param a - the input spectrum
 *
 * Returns:
 * @param res - the subtracted spectrum
 */
spectrum *
empty_counts_spectrum_copy (spectrum * a)
{
  int i;
  spectrum *res;
  
  if (a==NULL) return NULL;
  
  /* Allocate new spectrum */
  res = allocate_spectrum (a->spec_len);
  
  /* All the counts and flux are set to zero, rest is */
  /* identical to the original spectrum */
  for (i = 0; i < a->spec_len; i++)
    {
      res->spec[i].lambda_mean = a->spec[i].lambda_mean;
      res->spec[i].dlambda = a->spec[i].dlambda;
      res->spec[i].count = 0.0;
      res->spec[i].error = 0.0;
      res->spec[i].flux = 0.0;
      res->spec[i].ferror = 0.0;
      res->spec[i].weight = a->spec[i].weight;
      res->spec[i].lambda_max = a->spec[i].lambda_max;
      res->spec[i].lambda_min = a->spec[i].lambda_min;
    }
  
  return res;
}


/**
 * Function: create_SPC
 * A function that creates a FITS file containing an empty primary header
 *
 * Parameters:
 * @param filename  - the name of the file to open
 * @param overwrite - if set to 1, then any exisiting file is deleted. Nothing is done otherwise.
 */
void
create_SPC (char filename[], int overwrite)
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
		//	      "create_SPC: File %s "
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
     //fprintf(stderr,"Creating file %s\n",filename);
     fits_create_file (&output, filename, &f_status);
     if (f_status)
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "create_SPC: Could not open" " file: %s", filename);
       }

     // Create empty HDU
     {
	  int naxis = 0;
	  long naxes[2];
	  ffiimg (output, 16, naxis, naxes, &f_status);
	  if (f_status)
	    {
		 aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			      "create_SPC: Error creating "
			      " empty first HDU in: %s", filename);
	    }
     }

     fits_close_file (output, &f_status);
     if (f_status)
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "create_SPC: Could not" " close file: %s", filename);
       }
}


/**
 * Function: create_SPC_opened
 * A function that creates a FITS file containing an empty primary header
 *
 * @param filename  - The name of the file to open
 * @param overwrite - If set to 1, then any exisiting file is deleted. Nothing is done otherwise.
 *
 * Returns:
 * @return output - pointer to the opened fits-file
 */
fitsfile *create_SPC_opened (char filename[], int overwrite)
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
	//	      "create_SPC: File %s "
	//	      "exits. Overwriting it.", filename);
	fclose (in_file);
	unlink (filename);
      }
    if ((overwrite != 1) && (in_file != NULL))
      {
	fits_open_file (&output, filename, READONLY, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "create_SPC_opened: Could not open SPC file: %s",
			 filename);
	  }
	return output;
      }
  }
  // Open the file for creating/appending
  //fprintf(stderr,"Creating file %s\n",filename);
  fits_create_file (&output, filename, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_SPC: Could not open" " file: %s", filename);
    }
  
  // Create empty HDU
  {
    int naxis = 0;
    long naxes[2];
    ffiimg (output, 16, naxis, naxes, &f_status);
    if (f_status)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "create_SPC: Error creating "
		     " empty first HDU in: %s", filename);
      }
  }

  return output;
}


/**
 * Function: find_ID_in_SPC
 * A  function to return the number of the HDU which has a givenextname.
 *
 * Parameters:
 * @param filename - the name of the file to open
 * @param ID       - a string containing the ID to match
 *
 * Returns:
 * @return         - the number of the HDU with matching ID. Returns -1 if none is found. A warning is issued.
 */
int
find_ID_in_SPC (char filename[], char ID[])
{
  fitsfile *output;
  int f_status = 0;
  int hdunum = 0;
  char IDindex[60];

  /** NEED TO LOOK IN HDU 1 FOR KEY ID%D%D **/
  sprintf (IDindex, "ID%s", ID);
  hdunum = get_ID_index_to_SPC (filename, IDindex);
  /** IF FOUND EXIT NOW */
  if (hdunum != -1)
    {
      return hdunum;
    }
  
  
  //  Open the file for creating/appending
  fits_open_file (&output, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "find_ID_in_SPC: Could not open" " file: %s",
		   filename);
    }
  // Move to the HDU which has the desired extname
  fits_movnam_hdu (output, BINARY_TBL, ID, 0, &f_status);
  if (f_status)
    {
      fits_close_file (output, &f_status);
      return -1;
    }
  
  /* Get current HDU number */
  fits_get_hdu_num (output, &hdunum);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "find_ID_in_SPC: Could not get"
		   " current HDU number from file: %s", filename);
    }
  fits_close_file (output, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "find_ID_in_SPO: Could not" " close file: %s",
		   filename);
    }
  return hdunum;
}


/**
 * Function: add_ID_to_SPC_opened
 * A function that add a binary table extension to an existing opened SPC FITS file.
 * The table is appended to the FITS file. No check is performed to ensure that the
 * named extension does not exist already. This function is meant to initialize
 * a FITS binary table which will be used to store the content of an spectrum structure.
 * The fitsfile pointer is returned pointing to the newly created extension.
 *
 * Parameters:
 * @param a       - fitsfile pointer to an opened FITS file
 * @param N       - The size of the vectors in this table (i.e. number of pixels
 *                  in the ap_pixel structure).
 * @param extname - The name to give to this extension
 */
#define SPCCOL 15
void
add_ID_to_SPC_opened (fitsfile *output, int N, char ID[])
{
  /*   struct Col_Descr SPCData[SPCCOL] = {
       {"ID", "60A", NULL},
       {"N", "I1", NULL},
       {"LAMBDA", "XXXE1", "ANGSTROM"},
       {"TCOUNT", "XXXE1", "COUNT"},
       {"BCOUNT", "XXXE1", "COUNT"},
       {"COUNT", "XXXE1", "COUNT"},
       {"TERROR", "XXXE1", "COUNT"},
       {"BERROR", "XXXE1", "COUNT"},
       {"ERROR", "XXXE1", "COUNT"},
       {"FLUX", "XXXE1", "PHYSICAL UNITS"},
       {"FERROR", "XXXE1", "PHYSICAL UNITS"},
       {"WEIGHT", "XXXE1", "PIXEL"},
       {"CONTAM", "XXXI1", "FLAG"}
       };
*/
  struct Col_Descr SPCData[SPCCOL] = {
    {"ID", "60A", NULL},
    {"N", "I1", NULL},
    {"LAMBDA", "E1", "ANGSTROM"},
    {"TCOUNT", "E1", "COUNT"},
    {"BCOUNT", "E1", "COUNT"},
    {"COUNT", "E1", "COUNT"},
    {"TERROR", "E1", "COUNT"},
    {"BERROR", "E1", "COUNT"},
    {"ERROR", "E1", "COUNT"},
    {"FLUX", "E1", "PHYSICAL UNITS"},
    {"FERROR", "E1", "PHYSICAL UNITS"},
    {"WEIGHT", "E1", "PIXEL"},
    {"CONTAM", "E1", "FLAG"},
    // new for dlambda column
    {"DLAMBDA", "E1", "ANGSTROM"},
    {"DQ","I1","DQ"}};
  
     
  int f_status = 0;
  {
    char *ttype[SPCCOL], *tform[SPCCOL], *tunit[SPCCOL];
    
    int i;
    
    /* Prepare column description */
    for (i = 0; i < SPCCOL; i++)
      {
	ttype[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (ttype[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "intadd_ID_to_SPC_opened: "
			 "Memory allocation failed,");
	  }
	tform[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (tform[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "intadd_ID_to_SPC_opened: "
			 "Memory allocation failed,");
	  }
	tunit[i] = (char *) malloc (FLEN_KEYWORD * sizeof (char));
	if (tunit[i] == NULL)
	  {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "intadd_ID_to_SPC_opened: "
			 "Memory allocation failed,");
	  }
	if (SPCData[i].ttype != NULL)
	  sprintf (ttype[i], "%s", SPCData[i].ttype);
	if (SPCData[i].tform != NULL)
	  {
	    if (!strncmp ("XXX", SPCData[i].tform, 3))
	      {
		sprintf (tform[i], "%d%s", N,
			 SPCData[i].tform + 3);
	      }
	    else
	      {
		sprintf (tform[i], "%s", SPCData[i].tform);
	      }
	  }
	if (SPCData[i].tunit != NULL)
	  sprintf (tunit[i], "%s", SPCData[i].tunit);
	else
	  sprintf (tunit[i], "%s", " ");
      }
    
    fits_create_tbl (output, BINARY_TBL, 0, SPCCOL, ttype, tform, tunit,
		     ID, &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "intadd_ID_to_SPC_opened: Could not create"
		     " new binary table HDU in SPC file");
      }
    
    /* Clean up */
    for (i = 0; i < SPCCOL; i++)
      {
      	free (tunit[i]);
      	tunit[i] = NULL;
      	free (ttype[i]);
      	ttype[i] = NULL;
      	free (tform[i]);
      	tform[i] = NULL;
      }
    
    {
      char **array;
      int colnum;
      array = malloc (sizeof (char *));
      if (array == NULL)
	{
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "add_new_to_SPC: Out of memory");
	}
      array[0] = malloc (60 * sizeof (char));
      if (array[0] == NULL)
	{
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "add_net_to_SPC: Out of memory");
	}
      sprintf (array[0], "%s", ID);
      colnum = get_SPC_colnum (output, "ID");
      fits_write_col (output, TSTRING, colnum, 1, 1, 1, array,
		      &f_status);
      free (array[0]);
      array[0] = NULL;
      free (array);
      array = NULL;
      
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "add_new_to_SPC: "
		       "Could not write ID to row %d ,collumn %s (%d) in "
		       "SPC file ", 1, "ID", colnum);
	} 
    }
  } 
}


/**
 * Function: add_ID_to_SPC
 * A function that add a binary table extension to an existing SPC FITS file.
 * The table is appended to the FITS file. No check is performed to ensure that the
 * named extension does not exist already. This function is meant to initialize
 * a FITS binary table which will be used to store the content of an spectrum structure.
 *
 * Parameters:
 * @param filename - The name of the file to open
 * @param N        - The size of the vectors in this table (i.e. number of pixels in the ap_pixel structure).
 * @param extname  - The name to give to this extension
 *
 * Returns:
 * @return hdunum  - the number of the HDU that was just created
 */
int add_ID_to_SPC (char filename[], int N, char ID[])
{
 
  fitsfile *output;
  int f_status = 0;
  int hdunum = 0, hdutype;

  //  Open the file for creating/appending
  fits_open_file (&output, filename, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_new_to_SPC: Could not open" " file: %s",
		   filename);
    }
  fits_get_num_hdus (output, &hdunum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_new_to_SPC: Could not get"
		   " number of HDU from:", filename);
    }
  /* Move to last HDU */
  fits_movabs_hdu (output, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_new_to_SPC: Could not mov"
		   " to HDU number %d in file: %s", hdunum, filename);
    }
  /* Get current HDU number */
  fits_get_hdu_num (output, &hdunum);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_new_to_SPC: Could not get"
		   " current HDU number from file: %s", filename);
    }
  
  add_ID_to_SPC_opened (output, N, ID);
  
  
  fits_close_file (output, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_new_to_SPC: Could not" " close file: %s",
		   filename);
    }
  return hdunum;
}

/**
 * Function: add_ID_index_to_SPC
 * Add an HDU index keyword for beam ID%d%c into the first HDU of the FITS file.
 *
 * Parameters:
 * @param filename - a pointer to an array containing the name of an existing FITS file.
 * @param ID       - a pointer to a char array containing the beam ID in the form "ID%d%c".
 * @param hdunum   - the index of the HDU containing this beam.
 */
void
add_ID_index_to_SPC (char filename[], char ID[], int hdunum)
{
  fitsfile *input;
  int f_status = 0, hdutype;
  char comment[FLEN_COMMENT];

  //  Open the file for creating/appending
  fits_open_file (&input, filename, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ID_index_to_SPC: Could not open" " file:",
		   filename);
    }
  /* Move to first hdu */
  fits_movabs_hdu (input, 1, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ID_index_to_SPC: "
		   "Could not read extention %d from file: %s", hdunum,
		   filename);
    }
  sprintf (comment, "HDU number for beam");

  fits_write_key_lng (input, ID, (long) hdunum, comment, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ID_index_to_SPC: Error adding index keyword %s ",
		   ID);
    }

  fits_close_file (input, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_ID_index_to_SPC: Error closing file: %s ",
		   filename);
    }
}


/**
 * Function: get_ID_index_to_SPC
 * Get the value of an HDU index keyword of the for "ID%d%c" from the first HDU
 * of an existing FITS file.
 *  
 * Parameters:
 * @param filename - a pointer to an array containing the name of an existing FITS file.
 * @param ID       - a pointer to a char array containing the beam ID in the form "ID%d%c".
 *
 * Returns:
 * @return hdunum - the number of the HDU containing this beam.
 */
int
get_ID_index_to_SPC (char filename[], char ID[])
{
  fitsfile *input;
  int f_status = 0, hdutype;
  char comment[FLEN_COMMENT];
  long hdunum;
  
  //  Open the file for reading 
  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_ID_index_to_SPC: Could not open" " file:",
		   filename);
    }
  /* Move to first hdu */
  fits_movabs_hdu (input, 1, &hdutype, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_ID_index_to_SPC: "
		   "Could not read extention %d from file: %s", hdunum,
		   filename);
    }
  
  fits_read_key_lng (input, ID, &hdunum, comment, &f_status);
  if (f_status)
    {
      if (f_status == 202)
	{
	  hdunum = -1;
	  f_status = 0;
	}
      else
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "get_ID_index_to_SPC: Error getting index keyword %s ",
		       ID);
	}
    }

  fits_close_file (input, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_ID_index_to_SPC: Error closing file: %s ",
		   filename);
    }
  
  return (int) hdunum;
}


/** 
 * Function: get_SPC_colnum
 * A helper function which returns the collumn number of the column matching a 
 * given name. Function fails if collumn is not found.
 *
 * Parameters:
 * @param input   - a pointer to an opened FITS binary table
 * @param colname - the name of the desired column

 * Returns:
 * @return colnum - the number of the column. 
 */
int
get_SPC_colnum (fitsfile * input, char colname[])
{
  int colnum, f_status = 0;
  
  fits_get_colnum (input, CASEINSEN, colname, &colnum, &f_status);
#ifdef SPCDEBUG
  fprintf (stderr, "get_SPC_colnum: colum of %s is %d\n", colname, colnum);
#endif
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_SPC_colnum: Could not find collumn %s ",
		   colname);
    }

  return colnum;
}


/**
 * Function: add_spectra_to_SPC
 * Function to add a new entry into a SPC file.
 * This function takes three spectra, one containing the object, one containing
 * the background, and the last one containing the object subtracted background and
 * add the information to the appropriate (i.e. create a new one or update an existing one)
 * extension binary table.
 *
 * Parameters:
 * @param filename  - a pointer to a string containing the name of the PET file to write to/
 * @param obj_spec  - a pointer to a spectrum structure containing the object+background 1D binned spectrum
 * @param bck_spec  - a pointer to a spectrum structure containing the background 1D binned spectrum
 * @param sobj_spec - a pointer to a spectrum structure containing the object  1D binned spectrum
 * @param aperID    - the numeric aperture ID of this table.
 * @param beamID    - the numeric beam ID of this table.
 */
void
add_spectra_to_SPC (char filename[], spectrum * obj_spec, spectrum * bck_spec,
		    spectrum * sobj_spec, int aperID, int beamID)
{
     int hdunum;
     char ID[60], IDindex[60];
     int spec_len;

     create_SPC (filename, 0);

     sprintf (ID, "BEAM_%d%c", aperID, BEAM (beamID));
     sprintf (IDindex, "%s", ID);

     if (bck_spec != NULL)
       {
	    /* Check that it has the proper length */
	    if (obj_spec->spec_len != bck_spec->spec_len)
	      {
		   aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
				"add_spectra_to_SPC: Total and Background spectra"
				" do not have the same lengths!: %d vs. %d\n",
				obj_spec->spec_len, bck_spec->spec_len);
	      }
       }

     if (sobj_spec != NULL)
       {
	    /* Check that it has the proper length */
	    if (obj_spec->spec_len != sobj_spec->spec_len)
	      {
		   aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
				"add_spectra_to_SPC: Total and Object spectra"
				" do not have the same lengths!: %d vs. %d\n",
				obj_spec->spec_len, sobj_spec->spec_len);
	      }
       }

    if (obj_spec!=NULL) {
        spec_len = obj_spec->spec_len;
    }
    else {
        spec_len = 0;
    }
    
#ifdef DEBUGSPC
     fprintf (stderr, "Looking for ID:%s in SPC...", ID);
#endif


    /** NEED TO LOOK IN HDU 1 FOR KEY ID%D%D **/
     hdunum = get_ID_index_to_SPC (filename, IDindex);
    /** IF FAIL THEN DO FOLLOWING LINE */
     if (hdunum == -1)
       {
	    hdunum = find_ID_in_SPC (filename, ID);
       }

#ifdef DEBUGSPC
     fprintf (stderr, "Found at HDU: %d. ", hdunum);
#endif


     if (hdunum == -1)
       {
	    hdunum = add_ID_to_SPC (filename, spec_len, ID);
      /** NEED TO ADD KEYWORK ID%D%D IN HDU 0 WITH HDUNUM AS VALUE **/
	    add_ID_index_to_SPC (filename, IDindex, hdunum);
       }


#ifdef DEBUGSPC
     fprintf (stderr, "Done.\n");
#endif
     // Write data into HDU HERE!
     add_data_to_SPC (bck_spec, "BCOUNT", "BERROR", "WEIGHT", ID, filename,
		      hdunum, spec_len);

     add_data_to_SPC (obj_spec, "TCOUNT", "TERROR", "WEIGHT", ID, filename,
		      hdunum, spec_len);

     /* warning, if sobj_spec does not contain a valiv lambda and
	weight values, then these still get written... */
     add_data_to_SPC (sobj_spec, "COUNT", "ERROR", "WEIGHT", ID, filename,
		      hdunum, spec_len);
}


/**
 * Function: add_data_to_SPC
 * Function to add a specific spectrum structure into specifically 
 * named columns of an existing SPC FITS file. If a column name
 * is set to NULL then that particular data column is not written to the file.
 *
 * Parameters:
 * @param spec          - a pointer to an existing spectrum structure
 * @param coutncolname  - the name of the column to write spec.counts into
 * @param errorcolname  - the name of the column to write spec.errors into
 * @param weightcolname - the name of the column to write spec.weights into
 * @param ID            - the aperture/beam ID of this spectrum
 * @param filename      - the name of an existing SPC FITS table
 * @param hdunum        - the number of the extension in the SPC FITS file to write the data to
 * @param N             - the number of elements of the spectrum to write to the table
 */
void
add_data_to_SPC (spectrum * spec, char countcolname[], char errorcolname[],
		 char weightcolname[], char ID[], char filename[], int hdunum,
		 long N)
{
  fitsfile *input;
  int f_status = 0, hdutype;
  int colnum;
  char colname[FLEN_KEYWORD];
  int index = -1;

  //#define DEBUGSPC
  //  Open the file for creating/appending
  fits_open_file (&input, filename, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_data_to_SPC: " "Could not open file:",
		   filename);
    }
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_data_to_SPC: "
		   "Could not read extention %d from file: %s", hdunum,
		   filename);
    }
  if (hdutype != 2)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_data_to_SPC: "
		   "Extension %i of %s is not a BINARY table", hdunum,
		   filename);
    }

  /* Row number containing the desired ID */
  index = 1;

  /* If there are no elements in the table to write then exit now */
  if (N==0) {
    fits_close_file (input, &f_status);
    if (f_status)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "Error closing file: %s ", filename);
      }
    return;	 
  }

  /* Writing count, error, weight, lambda, contam, dq */
  {
    double *array_error, *array_count, *array_weight, *array_lambda;
    double *array_flux, *array_ferror, *array_contam;
    int *array_dq;
    int i;
    spc_entry *pt = NULL;
    
    array_error = malloc (N * sizeof (double));
    if (array_error == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_count = malloc (N * sizeof (double));
    if (array_count == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_weight = malloc (N * sizeof (double));
    if (array_weight == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_lambda = malloc (N * sizeof (double));
    if (array_lambda == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_flux = malloc (N * sizeof (double));
    if (array_flux == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_ferror = malloc (N * sizeof (double));
    if (array_ferror == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_contam = malloc (N * sizeof (double));
    if (array_contam == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_dq = malloc (N * sizeof (int));
    if (array_dq == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    if (spec != NULL)
      pt = spec->spec;
    for (i = 0; i < N; i++)
      {
	if (spec != NULL)
	  {
	    array_lambda[i] = pt->lambda_mean;
	    array_count[i] = pt->count;
	    array_error[i] = pt->error;
	    array_weight[i] = pt->weight;
	    array_flux[i] = pt->flux;
	    array_ferror[i] = pt->ferror;
            array_contam[i] = pt->contam;
            array_dq[i] = pt->dq;
	  }
	else
	  {
	    array_lambda[i] = GSL_NAN;
	    array_count[i] = GSL_NAN;
	    array_error[i] = GSL_NAN;
	    array_weight[i] = GSL_NAN;
	    array_flux[i] = GSL_NAN;
	    array_ferror[i] = GSL_NAN;
            array_contam[i] = -1.0;
            array_dq[i] = 0;
	  }
	pt++;
      }

    colnum = get_SPC_colnum (input, "LAMBDA");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_lambda,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in file %s ",
		     index, "LAMBDA", colnum, filename);
      }
    
    colnum = get_SPC_colnum (input, "FLUX");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_flux,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in file %s ",
		     index, "FLUX", colnum, filename);
      }
    
    colnum = get_SPC_colnum (input, "FERROR");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_ferror,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in file %s ",
		     index, "FERROR", colnum, filename);
      }
    
    colnum = get_SPC_colnum (input, "CONTAM");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_contam,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in file %s ",
		     index, "CONTAM", colnum, filename);
      }
    
    colnum = get_SPC_colnum (input, "DQ");
    fits_write_col (input, TINT, colnum, index, 1, N, array_dq,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in file %s ",
		     index, "DQ", colnum, filename);
      }

    /* Write the spec.counts data if the name of the column is defined */
    if (countcolname != NULL)
      {
	colnum = get_SPC_colnum (input, countcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_count, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write to row %d, collumn %s (%d) in file"
			 " %s ", index, countcolname, colnum,
			 filename);
	  }
      }

    /* Write the spec.errors data if the name of the column is defined */
    if (errorcolname != NULL)
      {
	colnum = get_SPC_colnum (input, errorcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_error, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write %s to row %d, collumn %s (%d) in file"
			 " %s ", errorcolname, index, colname,
			 colnum, filename);
	  }
      }
    
    /* Write the spec.weights data if the name of the column is defined */
    if (weightcolname != NULL)
      {
	colnum = get_SPC_colnum (input, weightcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_weight, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write %s to row %d, collumn %s (%d) in file"
			 " %s ", weightcolname, index, colname,
			 colnum, filename);
	  }
      }


    free (array_count);
    array_count = NULL;
    free (array_error);
    array_error = NULL;
    free (array_weight);
    array_weight = NULL;
    free (array_lambda);
    array_lambda = NULL;
    free (array_flux);
    array_flux = NULL;
    free (array_ferror);
    array_ferror = NULL;
    free (array_contam);
    array_contam = NULL;
    free (array_dq);
    array_dq = NULL;
     }

  fits_close_file (input, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "add_data_to_SPC: " "Error closing file: %s ",
		   filename);
    }
}



/**
 *
 * Function: add_spectra_to_SPC_opened
 * Function to add a new entry into a SPC file.
 * This function takes three spectra, one containing the object, one containing
 * the background, and the last one containing the object subtracted background and
 * add the information to the appropriate (i.e. create a new one or update an existing one)
 * extension binary table.
 *
 * @param filename a pointer to a string containing the name of the PET file to write to/
 * @param obj_spec a pointer to a spectrum structure containing the object+background 1D binned spectrum
 * @param bck_spec a pointer to a spectrum structure containing the background 1D binned spectrum
 * @param sobj_spec a pointer to a spectrum structure containing the object  1D binned spectrum
 * @param aperID the numeric aperture ID of this table.
 * @param beamID the numeric beam ID of this table.
 */
void
add_spectra_to_SPC_opened (fitsfile *input, spectrum * obj_spec, spectrum * bck_spec,
		    spectrum * sobj_spec, int aperID, int beamID)
{
  char ID[60], IDindex[60];
  int spec_len;
  
  sprintf (ID, "BEAM_%d%c", aperID, BEAM (beamID));
  sprintf (IDindex, "%s", ID);
  
  if (bck_spec != NULL)
    {
      /* Check that it has the proper length */
      if (obj_spec->spec_len != bck_spec->spec_len)
	{
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "add_spectra_to_SPC: Total and Background spectra"
		       " do not have the same lengths!: %d vs. %d\n",
		       obj_spec->spec_len, bck_spec->spec_len);
	}
    }

  if (sobj_spec != NULL)
    {
      /* Check that it has the proper length */
      if (obj_spec->spec_len != sobj_spec->spec_len)
	{
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "add_spectra_to_SPC: Total and Object spectra"
		       " do not have the same lengths!: %d vs. %d\n",
		       obj_spec->spec_len, sobj_spec->spec_len);
	}
    }
  
  if (obj_spec!=NULL) 
    spec_len = obj_spec->spec_len;
  else 
    spec_len = 0;
   
  add_ID_to_SPC_opened (input, spec_len, ID);
  /** NEED TO ADD KEYWORK ID%D%D IN HDU 0 WITH HDUNUM AS VALUE **/
  //  add_ID_index_to_SPC (filename, IDindex, hdunum);
  
  // Write data into HDU HERE!
  add_data_to_SPC_opened (bck_spec, "BCOUNT", "BERROR", "WEIGHT", ID, input, spec_len);
  
  add_data_to_SPC_opened (obj_spec, "TCOUNT", "TERROR", "WEIGHT", ID, input, spec_len);
  
  /* warning, if sobj_spec does not contain a valiv lambda and
     weight values, then these still get written... */
  add_data_to_SPC_opened (sobj_spec, "COUNT", "ERROR", "WEIGHT", ID, input, spec_len);
}

/**
 * Function: add_data_to_SPC_opened
 * Function to add a specific spectrum structure into specifically 
 * named columns of an existing opened SPC FITS file. If a column name
 * is set to NULL then that particular data column is not written to the file.
 *
 * @param spec          - a pointer to an existing spectrum structure
 * @param coutncolname  - the name of the column to write spec.counts into
 * @param errorcolname  - the name of the column to write spec.errors into
 * @param weightcolname - the name of the column to write spec.weights into
 * @param ID            - the aperture/beam ID of this spectrum
 * @param input         - a pointer to an opened FITS file and HDU
 * @param N             - the number of elements of the spectrum to write to the table
 */
void
add_data_to_SPC_opened (spectrum * spec, char countcolname[], char errorcolname[],
			char weightcolname[], char ID[], fitsfile *input,
			long N)
{
  int f_status = 0;
  int colnum;
  char colname[FLEN_KEYWORD];
  int index = -1;
  

  /* Row number containing the desired ID */
  index = 1;
  
  /* If there are no elements in the table to write then exit now */
  if (N==0)  return;	 
  
  /* Writing count, error, weight, lambda */
  {
    double *array_error, *array_count, *array_weight, *array_lambda;
    double *array_flux, *array_ferror, *array_contam, *array_dlambda;
    int *array_dq;
    int i;
    spc_entry *pt = NULL;
    
    array_error = malloc (N * sizeof (double));
    if (array_error == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_count = malloc (N * sizeof (double));
    if (array_count == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_weight = malloc (N * sizeof (double));
    if (array_weight == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_lambda = malloc (N * sizeof (double));
    if (array_lambda == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_flux = malloc (N * sizeof (double));
    if (array_flux == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_ferror = malloc (N * sizeof (double));
    if (array_ferror == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_contam = malloc (N * sizeof (double));
    if (array_contam == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_dlambda = malloc (N * sizeof (double));
    if (array_dlambda == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    array_dq = malloc (N * sizeof (int));
    if (array_dq == NULL)
      {
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: Out of memory");
      }
    
    if (spec != NULL)
      pt = spec->spec;
    for (i = 0; i < N; i++)
      {
	if (spec != NULL)
	  {
	    array_lambda[i] = pt->lambda_mean;
	    array_count[i] = pt->count;
	    array_error[i] = pt->error;
	    array_weight[i] = pt->weight;
	    array_flux[i] = pt->flux;
	    array_ferror[i] = pt->ferror;
	    array_contam[i] = pt->contam;
	    array_dlambda[i] = pt->dlambda;
	    array_dq[i] = (int) pt->dq;
	  }
	else
	  {
	    array_lambda[i] = GSL_NAN;
	    array_count[i] = GSL_NAN;
	    array_error[i] = GSL_NAN;
	    array_weight[i] = GSL_NAN;
	    array_flux[i] = GSL_NAN;
	    array_ferror[i] = GSL_NAN;
	    array_contam[i] = -1.0;
	    array_dlambda[i] = GSL_NAN;
	    array_dq[i] = (int) 0;
	  }
	
	pt++;
      }
    
    colnum = get_SPC_colnum (input, "LAMBDA");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_lambda,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "LAMBDA", colnum);
      }
       
    colnum = get_SPC_colnum (input, "FLUX");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_flux,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "FLUX", colnum);
      }
       
    colnum = get_SPC_colnum (input, "FERROR");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_ferror,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "FERROR", colnum);
      }
    
    colnum = get_SPC_colnum (input, "CONTAM");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_contam,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "CONTAM", colnum);
      }
       

    // NEW for dlambda column 
    colnum = get_SPC_colnum (input, "DLAMBDA");
    fits_write_col (input, TDOUBLE, colnum, index, 1, N, array_dlambda,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "DLAMBDA", colnum);
      }
       
    colnum = get_SPC_colnum (input, "DQ");
    fits_write_col (input, TINT, colnum, index, 1, N, array_dq,
		    &f_status);
    if (f_status)
      {
	ffrprt (stderr, f_status);
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "add_data_to_SPC: "
		     "Could not write to row %d, collumn %s (%d) in SPC file ",
		     index, "DQ", colnum);
      }
       
    /* Write the spec.counts data if the name of the column is defined */
    if (countcolname != NULL)
      {
	colnum = get_SPC_colnum (input, countcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_count, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write to row %d, collumn %s (%d) in SPC file"
			 " %s ", index, countcolname, colnum);				     
	  }
      }
    
    /* Write the spec.errors data if the name of the column is defined */
    if (errorcolname != NULL)
      {
	colnum = get_SPC_colnum (input, errorcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_error, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write %s to row %d, collumn %s (%d) in  SPC file"
			 , errorcolname, index, colname,
			 colnum);
	  }
      }
       
    /* Write the spec.weights data if the name of the column is defined */
    if (weightcolname != NULL)
      {
	colnum = get_SPC_colnum (input, weightcolname);
	fits_write_col (input, TDOUBLE, colnum, index, 1, N,
			array_weight, &f_status);
	if (f_status)
	  {
	    ffrprt (stderr, f_status);
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "add_data_to_SPC: "
			 "Could not write %s to row %d, collumn %s (%d) in SPC file"
			 , weightcolname, index, colname,
			 colnum);
	  }
      }


    free (array_count);
    array_count = NULL;
    free (array_error);
    array_error = NULL;
    free (array_weight);
    array_weight = NULL;
    free (array_lambda);
    array_lambda = NULL;
    free (array_flux);
    array_flux = NULL;
    free (array_ferror);
    array_ferror = NULL;
    free (array_contam);
    array_contam = NULL;
    free (array_dlambda);
    array_dlambda = NULL;    
    free (array_dq);
    array_dq = NULL;
  }
}

/**
 * Function: trim_spectrum
 * A function which trims the beginning and ending INDEF values 
 * out of a spectrum structure (left by the bin_naive routine.
 * The spectrum is replace in place.
 *
 * Parameters:
 * @param spc - a pointer to an existing spectrum structure
 *
 * Returns:
 * new_spc - pointer to a trimmed spectrum structure
 */
spectrum *
trim_spectrum (spectrum * spc)
{
  int i, is = 0, ie = 0;
  spectrum *new_spc;

  //fprintf(stderr,"here\n");
  for (i = 0; i < spc->spec_len; i++)
    {
      if (!(isnan (spc->spec[i].lambda_mean)))
	{
	  is = i;	/* Beginning of data found */
	  break;
	}
    }
  //fprintf(stderr,"here2\n");
  
  for (i = spc->spec_len - 1; i >= 0; i--)
    {
      if (!(isnan (spc->spec[i].lambda_mean)))
	{
	  ie = i;	/* End of data found */
	  break;
	}
    }
  //fprintf(stderr,"old length: %d\n",spc->spec_len);
  //fprintf(stderr,"Begin at: %d, end at %d\n",is,ie);

  /* Allocate new spectrum structure and transfer content from
     the old one into it. */
  new_spc = allocate_spectrum (ie - is + 1);
  for (i = is; i <= ie; i++)
    {
      new_spc->spec[i - is].lambda_mean = spc->spec[i].lambda_mean;
      new_spc->spec[i - is].dlambda = spc->spec[i].dlambda;
      new_spc->spec[i - is].lambda_min = spc->spec[i].lambda_min;
      new_spc->spec[i - is].lambda_max = spc->spec[i].lambda_max;
      new_spc->spec[i - is].error = spc->spec[i].error;
      new_spc->spec[i - is].count = spc->spec[i].count;
      new_spc->spec[i - is].weight = spc->spec[i].weight;
      new_spc->spec[i - is].contam = spc->spec[i].contam;
      new_spc->spec[i - is].dq = spc->spec[i].dq;
      
    }
  new_spc->warning = spc->warning;

  //fprintf(stderr,"new length: %d\n",new_spc->spec_len);
  
  return new_spc;
}
