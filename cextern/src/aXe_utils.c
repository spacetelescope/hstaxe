/* 
 * Some utility functions for aXe grism image processing
 */
#ifndef aXe_UTILS_C
  #define aXe_UTILS_C
#endif

#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_utils.h"


#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/**
 * A helper function which populates the error part of an
 * observation structure using simple Poisson noise. Assumes
 * data is in electrons.
 *
 *  @param img a pointer to an exisiting observation structure
 *
 *  @return a pointer to a GDL matrix containing the error array
*/
gsl_matrix     *
simulate_errors_t(gsl_matrix * img, const double exptime,
                  const double rdnoise)
{

  gsl_matrix *err;

  double nc;
  double sqr_rdnoise;

  int i, j;

  err = gsl_matrix_alloc(img->size1, img->size2);

  sqr_rdnoise = rdnoise*rdnoise;

  for (i = 0; i < (int)img->size1; i++) {
    for (j = 0; j < (int)img->size2; j++) {

      nc = fabs(gsl_matrix_get(img, i, j));

      nc = sqrt(nc*exptime + sqr_rdnoise)/exptime;

      gsl_matrix_set(err, i, j, nc);
    }
  }
  return err;
}

/**
 * This function uses the DQ information contained in an observation
 * structure together with a bitmask (represented as an integer) to
 * determine if a pixel is to be flagged as bad or not. Flagged
 * pixel are replaces by GSL_NAN values.
 *
 * @param obs a pointer to an exisiting observation structure
 * @param dqmask an integer representation of a bit mask.
 */
void
apply_DQ(observation * obs, int dqmask)
{
  int             i, j;
  fprintf(stdout, "DQMASK: %4i\n", dqmask);

  if ((obs->dq->size1 != obs->grism->size1) || (obs->dq->size2 != obs->grism->size2)) {
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "DQ array and DATA arrays have different sizes!\n");
  }
  if (obs->dq != NULL) {
    for (i = 0; i < (int)obs->dq->size1; i++) {
      for (j = 0; j < (int)obs->dq->size2; j++) {
        if ((int)gsl_matrix_get(obs->dq, i, j) & dqmask) {
          gsl_matrix_set(obs->grism, i, j, GSL_NAN);
        }
      }
    }

  }
}


/**
 * Function to load a multiple extension FITS file into an
 * observation structure
 *
 * @param fname a pointer to a string containing the name of the file to open
 * @param hdunum_data the extension number containing the data array (first=1)
 * @param hdunum_err the extension number containing the err array (first=1)
 * @param hdunum_dq the extension number containing the DQ array (first=1)
 * @param dqmask the bit mask (represented as an int) to apply using the
 *               DQ info to flag a pixel as good or bad
 * @param exptime exposure time (only used to determine count noise)
 * @param gain gain factor (only used to determine count noise)
 * @return observation a pointer to a new allocated, populated, DQ
 *         masked observation structure
 */
observation    *
load_image_t(const char *const fname, int hdunum_data, int hdunum_err,
             int hdunum_dq, int dqmask, const double exptime,
             const double rdnoise)
{
  observation    *obs;
  FILE           *infile;

  // check whether the image file
  // exists or not.
  if (!(infile = fopen(fname, "r")))
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "Could not find %s\n", fname);
  else
    fclose(infile);

  // allocate memory
  obs = malloc(sizeof(observation));

  // load the science extension
  fprintf(stdout,"Loading DATA from: %s...",fname);
  obs->grism = FITSimage_to_gsl(fname, hdunum_data, 1);
  fprintf(stdout,". Done.\n");


  // check whether there
  // is a dedicated error extension
  if (hdunum_err != -1)
    {
      // load the error extention
      fprintf(stdout, "Loading ERR...");
      obs->pixerrs = FITSimage_to_gsl(fname, hdunum_err, 0);
      fprintf(stdout,". Done.\n");
    }
  else
    {
      // simulate the error extension
      fprintf(stdout, "Simulating ERR with EXPTIME=%fs and RDNOISE= %fe ...", exptime, rdnoise);
      obs->pixerrs = simulate_errors_t(obs->grism, exptime, rdnoise);
      fprintf(stdout,". Done.\n");
    }


  // check whether there is a dedicated
  // dq extention
  if (hdunum_dq != -1)
    {
      // load the dq extention
      fprintf(stdout, "Loading DQ...");
      obs->dq = FITSimage_to_gsl(fname, hdunum_dq, 0);

      // apply the dq filter
      if (dqmask>0)
        {
          fprintf(stdout,"Applying MASK and seting DATA to NaN...");
          apply_DQ(obs, dqmask);
          fprintf(stdout,". Done.\n");
        }
    fprintf(stdout,". Done.\n");
    }
  else
    {
      fprintf(stdout, " Not using DQ. Done.\n");
      obs->dq = NULL;
    }
  fprintf(stdout,"\n");

  // return the structure created
  return obs;
}



/**
    Function to load a multiple extension FITS file into an observation structure

    @param fname a pointer to a string containing the name of the file to open
    @param hdunum_data the extension number containing the data array (first=1)
    @param hdunum_err the extension number containing the err array (first=1)
    @param hdunum_dq the extension number containing the DQ array (first=1)
    @param dqmask the bit mask (represented as an int) to apply using the DQ info to flag a pixel as
        good or bad
        @param exptime exposure time (only used to determine count noise)
        @param gain gain factor (only used to determine count noise)
    @return observation a pointer to a new allocated, populated, DQ masked observation structure
*/
observation    *
load_sci_image(const char *const fname, int hdunum_data)
{
        observation    *obs;
        FILE           *infile;

        if (!(infile = fopen(fname, "r"))) {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                            "Could not find %s\n", fname);
        } else {
                fclose(infile);
        }
        obs = malloc(sizeof(observation));

        fprintf(stdout,"Loading DATA from: %s...",fname);
        obs->grism = FITSimage_to_gsl(fname, hdunum_data, 1);
        fprintf(stdout,". Done.\n");

        obs->pixerrs = NULL;
        obs->dq = NULL;

        fprintf(stdout,"\n");

        return obs;
}

/**
    Read a FITS file header and return the value of the derived
    keyword. If the passed keyword string contains an actual number,
    then this number is converted to a float and returned.
    If keyword=None, then NULL is returned.

    @param filename a pointer to a char array containing the name of the FITS file
    @param hdunum the number of the extension to access
    @param keyword a pointer to a char array containing the name of the keyword, or a number

    @return the floating point value read from the FITS file or converted from the passed string
*/
float get_float_from_keyword(char *filename, int hdunum, char *keyword)
{
  fitsfile *input;
  int f_status = 0, hdutype;
  float val=-1.0;
  char comment[FLEN_COMMENT];

  /* If keyword contains an actual number, convert it to a float and
     return the value */
  //     fprintf(stdout, "opening: %s\n", filename);
  if (isnum2(keyword)) {
    val =  atof(keyword);
    return val;
  }

  /* If keyword == None return NULL */
  if (!strcmp(keyword,"None")) return GSL_NAN;

  //  Open the file for creating/appending
  //   fprintf(stdout, "opening: %s\n", filename);
  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_float_from_keyword: " " Could not open file: %s",
                   filename);
    }

  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_float_from_keyword: "
                   "Could not read extention %d from file: %s", hdunum,
                   filename);
    }


  fits_read_key_flt (input, keyword, &val, comment, &f_status);
  if (f_status)
    {
      //      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
      //                   "get_float_from_keyword: Could not read keyword %s from file: %s \n",
      //                   keyword,filename);
      f_status=0;
      fits_close_file (input, &f_status);
      if (f_status)
        {
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "get_float_from_keyword: " "Error closing file: %s ",
                       filename);
        }
      return GSL_NAN;
    }


  fits_close_file (input, &f_status);
  if (f_status)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_float_from_keyword: " "Error closing file: %s ",
                   filename);
    }

  return val;
}

/**
    Return the number of extensions in a FITS image

  @param fname File name of the FITS image
  @return the number of extensions
*/
int FITSextnum(const char *const fname)
{
     fitsfile *input;
     int f_status = 0;
     int hdunum;

     fits_open_file (&input, fname, READONLY, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Could not open file: %s", fname);
       }
     fits_get_num_hdus (input, &hdunum, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Could get number of extension in file: %s", fname);
       }

     fits_close_file (input, &f_status);

     return hdunum;
}

/**
  Read a FITS image, DO not read data but return an empty gsl array
  This function is meant to be a faster way to get information about
  an observation rather than the actual data.

  @param fname File name of the FITS image
  @return a gsl_matrix (overriden in aXe_grism.h) containing the
  image data or NULL
*/
gsl_matrix *
FITSnaxes_to_gsl (const char *const fname, int hdunum)
{
     fitsfile *input;
     int f_status = 0;
     int bitpix, naxis;
     long naxes[2];
     gsl_matrix *im;
     int hdutype;

     fits_open_file (&input, fname, READONLY, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Could not open" " file:", fname);
       }
     fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Could not read extention %d" " from file: %s",
                         hdunum, fname);
       }
     if (hdutype != 0)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Extension %i of %s is not " " an image type",
                         hdunum, fname);
       }
     ffgipr(input, 2, &bitpix, &naxis, naxes, &f_status);

     if (naxis != 2)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "%s is not a 2D" " FITS image.", fname);
       }

     im = gsl_matrix_alloc (naxes[0], naxes[1]);

     gsl_matrix_set_all (im, 0.0);


     return im;
}


/**
  Read a FITS image into a gsl array

  @param fname File name of the FITS image
  @param hdunum the number of the extension to read in (first primary = 1)
  @param fatal if non zero an error loading the data will trigger an exception
  @return a gsl_matrix (overriden in aXe_grism.h) containing the
  image data or NULL
*/
gsl_matrix *
FITSimage_to_gsl (const char *const fname, int hdunum, int fatal)
{
  fitsfile *input;
  int f_status = 0;
  float nulval = 0;
  int anynul;
  int bitpix, naxis;
  long naxes[2];
  long zero[2] = { 1, 1 };
  gsl_matrix *im;
  PIXEL_T *storage, *dp;
  int x, y;
  int hdutype;

  fits_open_file (&input, fname, READONLY, &f_status);
  if (f_status)
    {
      if (fatal)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not open file %s:", fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not open file %s:", fname);
          return NULL;
        }
    }
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      if (fatal)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not read extention %d"
                       " from file: %s", hdunum, fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not read extention %d"
                       " from file: %s", hdunum, fname);
          f_status = 0;
          fits_close_file (input, &f_status);
          return NULL;
        }
    }
  if (hdutype != 0)
    {
      if (fatal)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Extension %i of %s is not "
                         " an image type", hdunum, fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Extension %i of %s is not "
                       " an image type", hdunum, fname);
          return NULL;
        }
    }
  ffgipr (input, 2, &bitpix, &naxis, naxes, &f_status);

  if (naxis != 2)
    {
      if (fatal)
        {
            fits_report_error (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "FITSimage_to_gsl: " "%s is not a 2D FITS image.",
                         fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "%s is not a 2D FITS image.",
                       fname);
          return NULL;
        }
    }

  if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
    {
      fits_report_error (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "FITSimage_to_gsl: " "Out of memory");
    }

  ffgpxv (input, TFLOAT, zero, naxes[0] * naxes[1], &nulval,
          storage, &anynul, &f_status);
  if (f_status)
    {
      if (fatal)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not read data from extention %d"
                       " from file: %s", hdunum, fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not read data from extention %d"
                         " from file: %s", hdunum, fname);
          f_status = 0;
          fits_close_file (input, &f_status);
          return NULL;
        }
    }
  fits_close_file (input, &f_status);
  if (f_status)
    {
      if (fatal)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not close file: %s", fname);
        } else {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "FITSimage_to_gsl: " "Could not close  file: %s", fname);
          f_status = 0;
          fits_close_file (input, &f_status);
          return NULL;
        }
    }
  im = gsl_matrix_alloc (naxes[0], naxes[1]);
  dp = storage;
  for (y = 0; y < naxes[1]; y++)
    {
      for (x = 0; x < naxes[0]; x++)
        {
          gsl_matrix_set (im, x, y, *dp++);

        }
    }

  free (storage);
  storage = NULL;
  return im;
}
/**
  Read a FITS image into a gsl array

  @param fname File name of the FITS image
  @param hdunum the number of the extension to read in (first primary = 1)
  @param fatal if non zero an error loading the data will trigger an exception
  @return a gsl_matrix (overriden in aXe_grism.h) containing the
  image data or NULL
*/
px_point
get_npixels(const char *const fname, int hdunum)
{
     fitsfile *input;
     int f_status = 0;
     //float nulval = 0;
     //int anynul;
     int bitpix, naxis;
     long naxes[2];
     //long zero[2] = { 1, 1 };
     //int x;
     //int y;
     int hdutype;
     px_point numpix;

     fits_open_file (&input, fname, READONLY, &f_status);
     if (f_status)
       {
         fits_report_error (stderr, f_status);
         aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                      "FITSimage_to_gsl: " "Could not open file %s:", fname);
       }
     fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
     if (f_status)
       {
         fits_report_error (stderr, f_status);
         aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                      "FITSimage_to_gsl: " "Could not read extention %d"
                      " from file: %s", hdunum, fname);
       }
     if (hdutype != 0)
       {
            fits_report_error (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "FITSimage_to_gsl: " "Extension %i of %s is not "
                         " an image type", hdunum, fname);
       }
     ffgipr (input, 2, &bitpix, &naxis, naxes, &f_status);
     if (naxis != 2)
       {
         fits_report_error (stderr, f_status);
         aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                      "FITSimage_to_gsl: " "%s is not a 2D FITS image.",
                      fname);
       }

     numpix.x = (int)naxes[0];
     numpix.y = (int)naxes[1];
     fits_close_file (input, &f_status);
     if (f_status)
       {
         fits_report_error (stderr, f_status);
         aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                      "FITSimage_to_gsl: " "Could not close file: %s", fname);
       }
     return numpix;
}

/**
    Function to create a new FITS image file with an empty first extension.
    If overwrite is non zero, any existing file is deleted.

    @param filename a pointer to an array of character containing the filename
    @param overwrite if set to 1 exisiting file is deleted first
    @return a fitsifle pointer to a newly created FITS file
*/
fitsfile *
create_FITSimage_opened (char filename[], int overwrite)
{
     fitsfile *output;
     int f_status = 0;
     int hdunum,hdutype;

     // Try to open the file
     {
          FILE *in_file;
          in_file = fopen (filename, "r");

          if ((overwrite == 1) && (in_file != NULL))
            {
                 //aXe_message (aXe_M_WARN3, __FILE__, __LINE__,
                //            "create_FITSimage_opened: File %s "
                //            "exits. Overwriting it.", filename);
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
                    "get_ID_index_to_MEPET:" " Could not open file: %s",
                    filename);
            }
            fits_get_num_hdus (output, &hdunum, &f_status);
            if (f_status)
            {
                ffrprt (stderr, f_status);
                aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                    "gsl_to_FITSimage: Could not get"
                    " number of HDU from: %s", filename);
            }

            /* Move to last HDU */
            fits_movabs_hdu (output, hdunum, &hdutype, &f_status);
            if (f_status)
            {
                ffrprt (stderr, f_status);
                aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                    "gsl_to_FITSimage: Could not mov"
                    " to HDU number %d in file: %s", hdunum, filename);
            }
            return output;
            }
     }

     fits_create_file (&output, filename, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "create_FITSimage_opened: Could not open" " file: %s", filename);
       }

     // Create empty HDU
     {
          int naxis = 0;
          long naxes[2];
          ffiimg (output, 16, naxis, naxes, &f_status);
          if (f_status)
            {
                 aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                              "create_FITSimage_opened: Error creating "
                              " empty first HDU in: %s", filename);
            }
     }
    return output;
}

/**
    Function to create a new FITS image file with an empty first extension.
    If overwrite is non zero, any existing file is deleted.

    @param filename a pointer to an array of character containing the filename
    @param overwrite if set to 1 exisiting file is deleted first
*/
void
create_FITSimage (char filename[], int overwrite)
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
                //            "create_FITSimage: File %s "
                //            "exits. Overwriting it.", filename);
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
                         "create_FITSimage: Could not open" " file: %s", filename);
       }

     // Create empty HDU
     {
          int naxis = 0;
          long naxes[2];
          ffiimg (output, 16, naxis, naxes, &f_status);
          if (f_status)
            {
                 aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                              "create_FITSimage: Error creating "
                              " empty first HDU in: %s", filename);
            }
     }

     fits_close_file (output, &f_status);
     if (f_status)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "create_FITSimage: Could not" " close file: %s", filename);
       }
}

/**
    Store the content of a gsl_matrix into a FITS file

    @param data a pointer to an existing gsl_matrix structure
    @param fname a pointer to a string containing the name of the FITS file to write to
    @param overwrite if set to 1 then existing file is deleted, otherwise an extension is added
    @param ID an ID string which will be used to write an EXTNAME keyword in the extension
    @param HDUfile
    @param HDUfilenum

    @return the number of the HDU written to
*/
void
gsl_to_FITSimage_opened (gsl_matrix * data, fitsfile *output, int overwrite, char ID[])
{
     long naxes[2];
     int f_status = 0;
     PIXEL_T *storage, *dp;
     int x, y;

     if (data!=NULL) {
         naxes[0] = data->size1;
         naxes[1] = data->size2;
    } else {
        naxes[0] = 0;
        naxes[1] = 0;
    }

     //create_FITSimage(filename,overwrite);


     /* Moving pixels around into a "fitsio.h" friendly array */
     if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
       }
     dp = storage;
     for (y = 0; y < naxes[1]; y++)
       {
            for (x = 0; x < naxes[0]; x++)
              {
                   *dp = gsl_matrix_get (data, x, y);
                   dp++;
              }
       }

     fits_create_img (output, -32, 2, naxes, &f_status);

     fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                     &f_status);

     /* Add an EXTNAME to this extension */
     if (ID!=NULL) {
         char comment[FLEN_COMMENT];
         char keyname[FLEN_KEYWORD];
        //fits_write_key_lng (input, ID, (long) hdunum, comment, &f_status);
        sprintf(keyname,"EXTNAME");
        sprintf(comment,"name of this extension");
        fits_update_key (output, TSTRING, keyname, ID, comment, &f_status);
     }
     fits_write_date (output, &f_status);


     free (storage);
     storage = NULL;

}

/**
    Store the content of a gsl_matrix into a FITS file

    @param data a pointer to an existing gsl_matrix structure
    @param fname a pointer to a string containing the name of the FITS file to write to
    @param overwrite if set to 1 then existing file is deleted, otherwise an extension is added
    @param ID an ID string which will be used to write an EXTNAME keyword in the extension. If NULL, no extension name is written.
    @param HDUfile
    @param HDUfilenum

    @return the number of the HDU written to
*/
int
gsl_to_FITSimage (gsl_matrix * data, char filename[], int overwrite, char ID[])
{
  fitsfile *output;
  long naxes[2];
  int f_status = 0;
  PIXEL_T *storage, *dp;
  int x, y;
  int hdunum = 1, hdutype;

  if (data!=NULL) {
    naxes[0] = data->size1;
    naxes[1] = data->size2;
  } else {
    naxes[0] = 0;
    naxes[1] = 0;
  }

  create_FITSimage(filename,overwrite);


  /* Moving pixels around into a "fitsio.h" friendly array */
  if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  dp = storage;
  for (y = 0; y < naxes[1]; y++)
    {
      for (x = 0; x < naxes[0]; x++)
        {
          *dp = gsl_matrix_get (data, x, y);
          dp++;
        }
    }

  //  Open the file for creating/appending
  fits_open_file (&output, filename, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: " "Could not open file: %s",
                   filename);
    }


  fits_get_num_hdus (output, &hdunum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could not get"
                   " number of HDU from: %s", filename);
    }

  /* Move to last HDU */
  fits_movabs_hdu (output, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could not mov"
                   " to HDU number %d in file: %s", hdunum, filename);
    }
  /* Get current HDU number */
  fits_get_hdu_num (output, &hdunum);

  fits_create_img (output, -32, 2, naxes, &f_status);

  fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                  &f_status);

  /* Add an EXTNAME to this extension */
  if (ID!=NULL) {
    char comment[FLEN_COMMENT];
    char keyname[FLEN_KEYWORD];
    //fits_write_key_lng (input, ID, (long) hdunum, comment, &f_status);
    sprintf(keyname,"EXTNAME");
    sprintf(comment,"name of this extension");
    fits_update_key (output, TSTRING, keyname, ID, comment, &f_status);
  }
  fits_write_date (output, &f_status);

  fits_get_hdu_num(output,&hdunum);

  fits_close_file (output, &f_status);

  free (storage);
  storage = NULL;

  return hdunum;
}
/**
    Store the content of a gsl_matrix into a FITS file

    @param data a pointer to an existing gsl_matrix structure
    @param fname a pointer to a string containing the name of the FITS file to write to
    @param overwrite if set to 1 then existing file is deleted, otherwise an extension is added
    @param ID an ID string which will be used to write an EXTNAME keyword in the extension.
           If NULL, no extension name is written.
    @param hdunum, an integer pinting to the HDU the image should be written to
    @param HDUfile
    @param HDUfilenum

    @return the number of the HDU written to
*/
int
gsl_to_FITSimageHDU (gsl_matrix * data, char filename[], int overwrite, char ID[], int hdu_num)
{
     fitsfile *output;
     long naxes[2];
     int f_status = 0;
     PIXEL_T *storage, *dp;
     int x, y;
     int hdunum = 1, hdutype;

     if (data!=NULL) {
         naxes[0] = data->size1;
         naxes[1] = data->size2;
    } else {
        naxes[0] = 0;
        naxes[1] = 0;
    }

     create_FITSimage(filename,overwrite);


     /* Moving pixels around into a "fitsio.h" friendly array */
     if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
       }
     dp = storage;
     for (y = 0; y < naxes[1]; y++)
       {
            for (x = 0; x < naxes[0]; x++)
              {
                   *dp = gsl_matrix_get (data, x, y);
                   dp++;
              }
       }

     //  Open the file for creating/appending
     fits_open_file (&output, filename, READWRITE, &f_status);
     if (f_status)
       {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "gsl_to_FITSimage: " "Could not open file: %s",
                         filename);
       }


     //     fits_get_num_hdus (output, &hdunum, &f_status);
     //     if (f_status)
     //       {
     //     ffrprt (stderr, f_status);
     //     aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
     //                  "gsl_to_FITSimage: Could not get"
     //                  " number of HDU from: %s", filename);
     //       }

     /* Move to last HDU */
     fits_movabs_hdu (output, hdu_num, &hdutype, &f_status);
     if (f_status)
       {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "gsl_to_FITSimage: Could not mov"
                         " to HDU number %d in file: %s", hdunum, filename);
       }
     /* Get current HDU number */
     //     fits_get_hdu_num (output, &hdunum);

     //     fits_create_img (output, -32, 2, naxes, &f_status);

     fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                     &f_status);

     /* Add an EXTNAME to this extension */
     if (ID!=NULL) {
         char comment[FLEN_COMMENT];
         char keyname[FLEN_KEYWORD];
        //fits_write_key_lng (input, ID, (long) hdunum, comment, &f_status);
        sprintf(keyname,"EXTNAME");
        sprintf(comment,"name of this extension");
        fits_update_key (output, TSTRING, keyname, ID, comment, &f_status);
     }
     fits_write_date (output, &f_status);

     fits_get_hdu_num(output,&hdunum);

     fits_close_file (output, &f_status);

     free (storage);
     storage = NULL;

     return hdunum;
}

/**
  computes the bbox of a general quadrangle.

  @param corners the four corners of the quadrangle
  @param ll a px_point to hold the lower left corner
  @param ur a px_point to hold the upper right corner
*/
void
quad_to_bbox (const px_point * const corners, px_point * const ll,
              px_point * const ur)
{
     ll->x =
          MIN (MIN (MIN (corners[0].x, corners[1].x), corners[2].x),
               corners[3].x);
     ll->y =
          MIN (MIN (MIN (corners[0].y, corners[1].y), corners[2].y),
               corners[3].y);
     ur->x =
          MAX (MAX (MAX (corners[0].x, corners[1].x), corners[2].x),
               corners[3].x);
     ur->y =
          MAX (MAX (MAX (corners[0].y, corners[1].y), corners[2].y),
               corners[3].y);
}


/**
  Saves a gsl matrix to a pgm image file, properly scaled and all.

  @param fname name of pgm file to write
  @param m the matrix
  @return 0 for success, -1 for failure
*/
int
gsl_matrix_to_pgm (const char *const fname, const gsl_matrix * const m)
{
     FILE *f;
     double max = gsl_matrix_float_max (m);
     double min = gsl_matrix_float_min (m);
     int x, y;
     int szx = m->size1;
     int szy = m->size2;

     f = fopen (fname, "w");
     if (!f)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "gsl_matrix_to_pgm: Could not open file: %s", fname);
       }
     fprintf (f, "P5\n%d %d\n255\n", szx, szy);
     for (y = 0; y < szy; y++)
       {
            for (x = 0; x < szx; x++)
              {
                   fprintf (f, "%c",
                            (int)
                            floor ((gsl_matrix_get (m, x, szy - 1 - y) /
                                    (max - min) - min) * 255));
              }
       }
     fclose (f);
     return 0;
}

/**
    Returns a NULL terminated array of *char which each contain a card from a given
    hdu from a given FITS file. NAXES, NAXIS, BITPIX, SIMPLE, EXTEND, ID*,
    DATE cards are ignored.

    @param filename a pointer to a string containing the name of the FITS file to write to
    @param hdu the number of the HDU in the input FITS file to access

    @return a NULL terminated array of char arrays

*/
FITScards *get_FITS_cards (char filename[], int hdu)
{
    int f_status=0;
    int ninc, nexc;
    fitsfile *input;
    FITScards *cards=NULL;
    int hdutype,i,n;
    char card[FLEN_CARD];
    char *inclist[1] = {
    "*"
    };

    char *exclist[7] = {
    "NAXES",
    "NAXIS*",
    "BITPIX",
    "SIMPLE",
    "EXTEND",
    "ID*",
    "DATE"
    };

    ninc = 1;
    nexc = 7;


    //  Open the file for creating/appending
    fits_open_file (&input, filename, READONLY, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Could not open file: %s",
            filename);
    }

    // Move to the desired HDU
    fits_movabs_hdu (input, hdu, &hdutype, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Could not move to extension %d in file: %s",
            hdu,filename);
    }

    // Find the number of cards in this HDU
    n = 0;
    do {
        fits_find_nextkey (input, inclist, ninc, exclist,
         nexc, card, &f_status);
        n++;
    }
    while (f_status==0);
    if (f_status!= KEY_NO_EXIST) {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Problem reading keys in extension %d of file: %s",
            hdu,filename);
    }
    f_status = 0;


    // Allocate enough room (n+1 for NULL termination) for all the cards
    //cards = (char **) malloc(sizeof(char *)*(n+1));
    //for(i=0;i<(n+1);i++) cards[i] = malloc(sizeof(char)*FLEN_CARD);
    cards = allocate_FITScards(n);

    // Move back tot he top of the HDU
    fits_read_record (input, 0, card, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Could not move to first card in extension %d in file: %s",
            hdu,filename);
    }

    // Read all cards
    i=0;
    do {
        fits_find_nextkey (input, inclist, ninc, exclist,
         nexc, cards->cards[i], &f_status);
        //fprintf(stderr,"%s\n",cards->cards[i]);
        i++;
    }
    while (f_status==0);
    if (f_status!= KEY_NO_EXIST) {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Problem reading keys in extension %d of file: %s",
            hdu,filename);
    }
    f_status = 0;

    // NULL terminate the card array
    //cards[n] = NULL;

    fits_close_file (input, &f_status);
    return cards;
}


FITScards *get_FITS_cards_opened (fitsfile *input)
{
    int f_status=0;
    int ninc, nexc;
    FITScards *cards=NULL;
    //int hdutype;
    int i,n;
    char card[FLEN_CARD];
    char *inclist[1] = {
    "*"
    };

    char *exclist[15] = {
    "NAXES",
    "NAXIS*",
    "BITPIX",
    "SIMPLE",
    "EXTEND",
    "ID*",
    "DATE",
    "EXTNAME",
    "TTYPE*",
    "TFORM*",
    "TUNIT*",
    "TFIELDS",
    "XTENSION",
    "PCOUNT",
    "GCOUNT"
    };

    ninc = 1;
    nexc = 15;


    // Move back tot he top of the HDU
    fits_read_record (input, 0, card, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Could not move to first card in FITS extension");
    }
    // Find the number of cards in this HDU
    n = 0;
    do {
        fits_find_nextkey (input, inclist, ninc, exclist,
         nexc, card, &f_status);
        n++;
    }
    while (f_status==0);
    if (f_status!= KEY_NO_EXIST) {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Problem reading keys in FITS extension");
    }
    f_status = 0;

    // Allocate enough room for all the cards
    cards = allocate_FITScards(n);
    // Move back tot he top of the HDU
    fits_read_record (input, 0, card, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Could not move to first card in FITS extension");
    }

    // Read all cards
    i=0;
    do {
        fits_find_nextkey (input, inclist, ninc, exclist,
         nexc, cards->cards[i], &f_status);
        i++;
    } while (f_status==0);
    if (f_status!= KEY_NO_EXIST) {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "get_FITS_cards: " "Problem reading keys in FITS extension");
    }
    f_status = 0;

    return cards;
}

/**
    routine to write a set of FITS header cards into a FITS header of a given
    extension. Input cards are automatically formatted propely if they are not.

    @param filename a pointer to a char array containing the name fo the existing FITS file
    @param hdu the number of the extension to write to (-1 for last one)
    @param cards a NULL terminated array of char arrays containing properly formatted cards

*/
void put_FITS_cards (char filename[], int hdu, FITScards *cards)
{
    fitsfile *output;
    int hdutype,f_status = 0;
    int i, keytype;
    char card[FLEN_CARD];

    // Do nothing if there is nothing to do
    if (cards->n == 0) return;

    //  Open the file for creating/appending
    fits_open_file (&output, filename, READWRITE, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "put_FITS_cards: " "Could not open file: %s",
            filename);
    }

    // Move to the desired HDU
    // If -1 go to the last one
    if (hdu==-1) {
        fits_get_num_hdus (output, &hdu, &f_status);
        if (f_status)
        {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "put_FITS_cards: Could not get"
                " number of HDU from:", filename);
        }
    }
    fits_movabs_hdu (output, hdu, &hdutype, &f_status);
    if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "put_FITS_cards: " "Could not move to extension %d in file: %s",
            hdu,filename);
    }

    for (i=0; i<cards->n;i++)
    {
        fits_parse_template (cards->cards[i], card, &keytype, &f_status);
//fprintf(stderr,"%s\n",cards->cards[i]);
        if (f_status)
        {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "put_FITS_cards: " "Could not reformet card:\n%s",cards[i]);
        }
        fits_write_record (output, card, &f_status);
        if (f_status)
        {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "put_FITS_cards: " "Could not write card in extension %d in file: %s",
                hdu,filename);
        }
    }
    fits_write_date (output, &f_status);
    fits_close_file (output, &f_status);

}

/**
    routine to write a set of FITS header cards into a FITS header of a given
    extension. Input cards are automatically formatted propely if they are not.

    @param input a fitsfile pointer pointing to the extension to write the cards to.
    @param hdu the number of the extension to write to (-1 for last one)
    @param cards a NULL terminated array of char arrays containing properly formatted cards

*/
void put_FITS_cards_opened (fitsfile *output, FITScards *cards)
{
    int f_status = 0;
    int i, keytype;
    char card[FLEN_CARD];

    // Do nothing if there is nothing to do
    if (cards->n == 0) return;

    for (i=0; i<cards->n;i++)
    {
        fits_parse_template (cards->cards[i], card, &keytype, &f_status);
        if (f_status)
        {
            ffrprt (stderr, f_status);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "put_FITS_cards_opened: " "Could not reformat card:\n%s",cards[i]);
        }
        fits_write_record (output, card, &f_status);
        if (f_status)
        {
            ffrprt (stderr, f_status);
            fprintf(stderr, "FITScard:   %s\n", cards->cards[i]);
            fprintf(stderr, "FITSrecord: %s\n", card);
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "put_FITS_cards_opened: %i, %s" "Could not write cards.", i, card);
        }
    }
    fits_write_date (output, &f_status);
}

/**
    A helper function which takes a filename and returns the real filename
    and the number of the extension if the input filename was formatted with the format
    file[hdu].

    @param file a pointer to a string containing a filename, possibly formatted as file[hdu]
    @param filename a pointer to an allocated char array to receive the real file name
    @param hdu an integer pointer to receive the hdu number

*/
void
get_filename_and_hdunum (char file[], char *filename, int *hdu)
{
     char *find_strings[] = { "[" };
     char *here;
     int here_len, file_len;

     /* Look for "[" */
     init_search (find_strings[0]);
     here = strsearch (file);

     if (here)
       {
            /* If non-null */
            here_len = strlen (here);
            file_len = strlen (file);

            /* copy file into filename */
            strcpy (filename, file);
            /* truncate filename */
            if (((file_len - here_len) >= 0)
                && ((file_len - here_len) <= file_len))
                 filename[file_len - here_len] = 0;
            /* get the hdu number */
            sscanf (here, "[%d]", hdu);
       }
     else
       {
            /* if null ("[" not found) */
            /* copy input file to output file */
            strcpy (filename, file);
            /* default hdu to invalid hdu */
            *hdu = -1;
       }

}

/**
    A helpder function which given the name of an environment variable and a string
    containing a filename, this function
    returns the complete path to the file

    @param envpathname a pointer to a cahr string containing the name of the
    environment variable to use for the path to the file
    @param inputfile a pointer to a string containing a filename, possibly formatted as file[hdu].
    This is replaced with the name of the file without the extension appended to it.
    @param filename a pointer to an allocated char array to receive the real file name
*/
void
build_path (char envpathname[], const char inputfile[], char filename[MAXCHAR])
{
     char *envvar;
     //int envpathname_len;

     /*
       Changed for aXe-1.4 such that the pathnames are OPTIONAL only.
      */
     if (!(envvar = getenv (envpathname))){
       strcpy (filename, "./");
     } else {
       strcpy (filename, envvar);
       strcat (filename, "/");
     }
     strcat (filename, inputfile);

     //   old version:
     //     if (!(envvar = getenv (envpathname)))
     //   aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
     //                "Environment variable %s not found", envpathname);
     //
     //     envpathname_len = strlen (envpathname);
     //
     //     strcpy (filename, envvar);
     //
     //     strcat (filename, "/");
     //
     //     strcat (filename, inputfile);

}

/**
 * Function: alloc_char_arr
 * The function allocates memory for a character array.
 * The dimensions of the character array are specified
 * with the number of columns per row and the number
 * of rows in the input. The pointer to the allocated
 * memory is returned.
 *
 * Parameters:
 * @param nrows   - number of rows to allocate
 * @param ncols   - number of columns per row to allocate
 *
 * Returns:
 * @return char_arr - pointer to the allocated array
 */
char **
alloc_char_arr(int nrows, int ncols)
{
  char **char_arr;
  int i;

  // allocate the memory for the pointers to the individual rows
  char_arr = (char **) malloc(sizeof(char *)*nrows);

  // for each row:
  for (i=0; i < nrows; i++)
    // allocate the columns
    char_arr[i] = malloc(sizeof(char)*ncols);

  // return the result
  return char_arr;
}

/**
 * Function: free_char_arr
 * The function releases the memory allocated
 * in a memory array.
 *
 * Parameters:
 * @param char_arr - pointer to the character array
 * @param nrows    - the number of rows in the character array
 *
 */
void
free_char_arr(char **char_arr, int nrows)
{
  int i;

  // for each row
  for(i=0; i < nrows; i++)
    // release the memory
    free(char_arr[i]);

  // release the rest
  free(char_arr);
}

/**
 * Function: build_config_files
 * The function parses a long character string and separates
 * individual tokens along the character ",". The tokens are
 * stored in a character array. A check is performed such
 * that not more tokens are copied than there exist rows
 * in the character array.
 *
 * Parameters:
 * @param list_file  - the character string to parse
 * @param list_files - the character array
 * @param maxnrows   - the number of rows in the character array
 *
 * Returns:
 * @return nlist   - the number of tokens copied
 */
int
build_config_files(char list_file[], char **list_files, int maxnrows)
{

  int nlist=0;
  int ii = 0;
  int i=0;
  int mod_maxrows;

  // create a modified maximum number fo an easy check
  mod_maxrows = maxnrows -1;

  // go over the input string
  for (i=0; i<(int)strlen(list_file); i++)
    {

      // check whether the actual character
      // indicates the beginning of the next token
      if (list_file[i] != ',' && list_file[i] != '\0')
        {
          // check whether ther is space in the character array
          if (nlist > mod_maxrows)
            // give an error in case the char array is full
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "Too many configuration files! Maximum number: %i !\n", maxnrows);
          else
            // copy the actual character into the token
            list_files[nlist][ii++] = list_file[i];
        }
    else
      {
        // write an end to the actual token
        list_files[nlist][ii] = '\0';

        // enhance the token counter
        nlist++;

        // reset the position within the token
        ii = 0;
      }
    }

  // return the number of tokens
  return ++nlist;
}


void
replace_file_extension (char infile[], char outfile[], char from_ext[],
                        char to_ext[], int hdu)
{
     char chdu[MAXCHAR];
     if (!infile)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "replace_file_extension: Input filename invalid");
       }
     if (!from_ext)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "replace_file_extension: Input extension invalid");
       }
     if (!to_ext)
       {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "replace_file_extension: output extension invalid");
       }

     strcpy (outfile, infile);
     outfile[strlen (outfile) - strlen (from_ext)] = 0;
     if (hdu > 0)
       {
            sprintf (chdu, "_%d", hdu);
            strcat (outfile, chdu);
       }
     strcat (outfile, to_ext);


}

/**
    Free a NULL terminated array of objects

    @param sobjs, a NULL terminated array containing pointers to objects

*/
void
free_oblist (object ** oblist)
{
  int i, j, nobjs = 0;

  /* Find the number of objects in sobjs */
  while (oblist[nobjs])
    nobjs++;

  for (i = 0; i < nobjs; i++)
    {
      for (j=0;j < oblist[i]->nbeams; j++){
        free(oblist[i]->beams[j].spec_trace->data);
        free(oblist[i]->beams[j].spec_trace);

        if (oblist[i]->beams[j].flux != NULL){
          gsl_vector_free (oblist[i]->beams[j].flux);
          //      fprintf(stderr, "FLux was freed\n");
        }
        else{
          //      fprintf(stderr, "FLux is NULL\n");
        }
        //free(oblist[i]->*beams[j]);
      }
      free (oblist[i]);
      oblist[i] = NULL;
    }
  free (oblist);
  oblist = NULL;
}


/**
    output the content of an object structure

*/
void
fprintf_object (FILE * output, object * o)
{
     int i, j;

     fprintf (output, "object ID: %d\n", o->ID);
     fprintf (output, "object nbeams: %d\n", o->nbeams);

     for (i = 0; i < o->nbeams; i++)
       {
            fprintf (output, "object beam ID: #%d ID\n", o->beams[i].ID);
            fprintf (output, "object beam #%d refpixel: %f %f\n",
                     o->beams[i].ID, o->beams[i].refpoint.x,
                     o->beams[i].refpoint.y);

            for (j = 0; j < 4; j++)
              {
                   fprintf (output, "object beam #%d corner[%d]: %d %d\n",
                            o->beams[i].ID, j, o->beams[i].corners[j].x,
                            o->beams[i].corners[j].y);
              }

            for (j = 0; j < 2; j++)
              {
                   fprintf (output, "object beam #%d bbox[%d]: %d %d\n",
                            o->beams[i].ID, j, o->beams[i].bbox[j].x,
                            o->beams[i].bbox[j].y);
              }

            fprintf (output, "object beam #%d width: %f\n", o->beams[i].ID,
                     o->beams[i].width);
            fprintf (output, "object beam #%d orient: %f\n", o->beams[i].ID,
                     o->beams[i].orient);


       }



}


/**
   Function to parse the online arguments and look for an option
   (starting with a "-") with a certain name. If found, a pointer
   to a string pointing to the value of the option (following a "=")
   is returned. If no "=" is found, then the name of the option
   is returned. If the names is not found, then NULL is returned.

   @param option_name the name of the option to look for (without the "-")
   @param argc the number of online parameter
   @param *argv a pointer to an array of strings containing all the
   online parameters.

   @return a pointer to an array containing the value of the option, the
   name of the option if it is found but not value was found (i.e. toggle),
   and NULL if the option was not found at all.
*/
char *
get_online_option (char option_name[], int argc, char *argv[])
{
     char *parname;
     char *val;
     static char wrkarg[MAXCHAR];
     int i;

     for (i = 1; i < argc; i++)
       {
            if (argv[i][0] == '-')
              {
                   strcpy (wrkarg, argv[i] + 1);
                   parname = strtok (wrkarg, "=");
                   if (!strcmp (parname, option_name))
                     {
                          val = strtok (NULL, "=");
                          if (val == NULL)
                               return parname;
                          return val;
                     }
              }
       }
     return NULL;

}
int
get_online_option2 (char option_name[], char option_value[], int argc, char *argv[])
{
     char *parname;
     char *val;
     static char wrkarg[MAXCHAR];
     int i;

     for (i = 1; i < argc; i++)
       {
            if (argv[i][0] == '-')
              {
                   strcpy (wrkarg, argv[i] + 1);
                   parname = strtok (wrkarg, "=");
                   if (!strcmp (parname, option_name))
                     {
                          val = strtok (NULL, "=");
                          if (val == NULL)
                            {
                              strcpy(option_value,parname);
                              return 1;}

                          strcpy(option_value,val);
                          return 1;
                     }
              }
       }
     return 0;

}

/**
    Function to free an observation struture

    @param obs an observation pointer
*/
void
free_observation (observation * obs)
{
  if (obs->grism != NULL)
    gsl_matrix_free (obs->grism);
  if (obs->pixerrs != NULL)
    gsl_matrix_free (obs->pixerrs);
  if (obs->dq != NULL)
    gsl_matrix_free (obs->dq);
  free (obs);
}



/**
   A helper function which returns the extension number in a FITS
   file which has a given name and also optionally containing a
   given keyword/key combination. The leyword/key pair is meant to
   diferentiate between different extensions os the same types. For
   example between two SCI extension containing CCDNUM=1 or CCDNUM=2
   keys.

   @param filename a pointer to a string containing the name
   of an existing FITS file
   @param extname a pointer to a string containing the name
   of the extension to look for
   @param keyword a pointer to a string containing the name
   of an extra keyword to lolok for. Can be set to NULL
   @param a pointer to a string  containing the value of the optional key
   we are looking for
*/
int get_hdunum_from_hduname(char filename[], char extname[],
                           char keyword[], char keyval[], int extver)
{
  fitsfile *input;
  int f_status = 0;
  char comment[FLEN_COMMENT];
  char sval[FLEN_VALUE];
  int i,hdunum,hduext;

  /* If the extname is a number then we simply return that number */
  if (isnum2(extname)) {
    hdunum = (long) atoi(extname);
    return (int) hdunum;
  }

  /* If it not then we first look for the names extension */
  //  fprintf(stderr,"Extension keyword %s %i kkkk \n",filename, f_status);
  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      return -1;
    }
  i=0;
  while (!f_status) {
    /* Keep looking for extensions matching the required */
    /* extname */
    i++;
    fits_movabs_hdu (input, i, NULL, &f_status);
    if (f_status) {
      /* No more HDU found. If we are here then we failed to find */
      /* a matching extension and return -1 */
      f_status = 0;
      fits_close_file (input, &f_status);
      return -1;
    }
    /* Read extname and try to match it */
    fits_read_key_str (input, "EXTNAME", sval, comment, &f_status);
    if (f_status) {
      /* Did not manage to get EXTNAME, keep looking */
      f_status = 0;
      continue;
    } else {
      /* found EXTNAME, check if it is the one we want, or continue */
      //      fprintf(stdout,"Extension name %s kkkk \n",sval);
      if (strcmp(extname,sval)) {
        f_status = 0;
        continue;
      }
    }
    /* If the passed keyword is non NULL then we attempt to match it */
    /* otherwise, we are done */
    if (!strcmp(keyword,"None")) {
      fits_get_hdu_num(input,&hdunum);
      fits_close_file (input, &f_status);
      return hdunum;
    }
    fits_read_key_str (input, keyword, sval, comment, &f_status);
    if (f_status) {
/* keyword was not found, try whether the EXTVER keyword matches */
      if (extver > -1){
        f_status = 0;
        fits_read_key_str (input, "EXTVER", sval, comment, &f_status);
        if (f_status) {
          f_status = 0;
          continue;
        }
        else{
          hduext = atoi(sval);
          if (hduext == extver){
            fits_get_hdu_num(input,&hdunum);
            fits_close_file (input, &f_status);
            return hdunum;
          }
          else{
          f_status = 0;
          continue;
          }
        }
      }
      else {
        f_status = 0;
        continue;
      }
    } else {
      if (strcmp(sval,keyval)) {
        /* keyword was found but does not match the passed value */
        /* skip to next HDU */
        continue;
      } else {
        /* keyword was found and the key matches, return the current */
        /* HDU number */
        fits_get_hdu_num(input,&hdunum);
        fits_close_file (input, &f_status);
        return hdunum;
      }
    }
  }
  return -1;
}
int drzprep_fitstrans(char filename[], int hdunum, fitsfile *PET_ptr){

  char kname[FLEN_KEYWORD];
  char kvalue[FLEN_VALUE];
  char kcomment[FLEN_COMMENT];
  fitsfile *input;
  int      f_status=0, hdutype=0;
  //long tmp;

  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "drzprep_fitstrans: " "Could not open file: %s",
            filename);
    }
  // Move to the desired HDU
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not move to extension %d in file: %s",
                   hdunum,filename);
    }

  strcpy(kname, "OBJECTID");
  //  fits_read_key_lng (input, kname, &tmp, kcomment, &f_status);
  fits_read_key_str(input, kname, kvalue, kcomment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not get keyword: %s",
                   kname);
    }

  fits_write_key_str(PET_ptr, kname, kvalue, kcomment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not write keyword: %s",
                   kname);
    }
  strcpy(kname, "BEAMID");
  //  fits_read_key_lng (input, kname, &tmp, kcomment, &f_status);
  fits_read_key_str(input, kname, kvalue, kcomment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not get keyword: %s",
                   kname);
    }

  fits_write_key_str(PET_ptr, kname, kvalue, kcomment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not write keyword: %s",
                   kname);
    }
  fits_close_file (input, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "drzprep_fitstrans: " "Could not close fitsfile.",
                   kname);
    }
  return 0;
}

/**
 * Function: get_npixel
 * The function extracts the number of pixels in x an y
 * or the science image associated to an
 * 'observation'-structure.
 *
 * Parameters:
 * @param  obs     - the observation with the science image
 *
 * Returns:
 * @return npixels - the struct to store the pixel numbers in x and y
 */
px_point
get_npixel(const observation *obs)
{

  px_point npixels;

  npixels.x = obs->grism->size1;
  npixels.y = obs->grism->size2;

  return npixels;
}
