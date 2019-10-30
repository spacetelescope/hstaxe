/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"
#include "wcs.h"

struct WorldCoor *get_wcs_l(char filename[], int hdunum);
char      *get_fits_header_l(char filename[], int hdunum);

int main(int argc, char *argv[])
{
  //char *opt;
  char direct_image[30];
  //char direct_image_path[30];
  //int direct_hdunum;

  char grism_image[30];
  //char grism_image_path[30];
  //int grism_hdunum;

  struct WorldCoor *dirim_wcs;
  struct WorldCoor *grisi_wcs;

  double x_in   = 1000.0;
  double y_in   = 1000.0;
  double x_out  =    0.0;
  double y_out  =    0.0;
  double ra_in  =    0.0;
  double dec_in =    0.0;

  int offscl=0;

  if (argc<3)
    {
      fprintf (stdout, "ST-ECF European Coordinating Facility\nCopyright (C) 2002 Sofisof\naXe_TEST Version 0.01: \n\nUsage:\n"
                "aXe_TEST direct_image.fits grism_image.fits\n\n");
      exit (1);
    }
  fprintf (stdout, "aXe_TEST: Starting...\n");

  // copy the input over
  strcpy (direct_image, argv[1]);
  strcpy (grism_image,  argv[2]);

  // give feedback on the input
  fprintf (stdout, "aXe_TEST: direct image: %s\n", direct_image);
  fprintf (stdout, "aXe_TEST: grism  image: %s\n", grism_image);

  // load the wcs of both images
  dirim_wcs = get_wcs_l (direct_image, 2);
  grisi_wcs = get_wcs_l (grism_image, 2);

  // convert a point onto the sky
  pix2wcs (dirim_wcs, x_in, y_in, &ra_in, &dec_in);

  fprintf (stdout, "aXe_TEST: (%.1f,%.1f) --> (%e,%e)\n", x_in, y_in, ra_in, dec_in);

  wcs2pix (grisi_wcs, ra_in, dec_in, &x_out, &y_out, &offscl);
  fprintf (stdout, "aXe_TEST: (%e,%e) --> (%.1f,%.1f)\n", ra_in, dec_in, x_out, y_out);

  /*
  // load the wcs of both images
  dirim_wcs = get_wcs (direct_image, 3);
  grisi_wcs = get_wcs (grism_image, 3);

  fprintf(stdout, "\n\n");

  // convert a point onto the sky
  pix2wcs (dirim_wcs, x_in, y_in, &ra_in, &dec_in);

  fprintf (stdout, "aXe_TEST: (%.1f,%.1f) --> (%e,%e)\n", x_in, y_in, ra_in, dec_in);

  wcs2pix (grisi_wcs, ra_in, dec_in, &x_out, &y_out, &offscl);
  fprintf (stdout, "aXe_TEST: (%e,%e) --> (%.1f,%.1f)\n", ra_in, dec_in, x_out, y_out);
*/

  // print a last status message and exit
  fprintf(stdout, "aXe_TEST: Done.\n");
  exit(0);
}


/**
 * Function: get_fits_header
 * Returns the header of extension hdunum from a fits file filename
 *
 * Parameters:
 * @param filename The name of the FITS file
 * @param hdunum the number of the HDU to access
 *
 * Returns:
 * @return a pointer to a newly allocated string containing the entire header
 */
char *
get_fits_header_l (char filename[], int hdunum)
{

  fitsfile *input;
  int f_status = 0, hdutype;
  int nkeys, i;
  char card[FLEN_CARD];

  char *header, *p;

  //  Open the file for creating/appending
  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      fprintf (stdout, "ERROR: get_fits_header: Could not open file: %s", filename);
    }
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      fprintf (stdout, "ERROR: get_fits_header: Could not read extention %d from file: %s", hdunum, filename);
    }



  fits_get_hdrspace (input, &nkeys, NULL, &f_status);   /* get # of keywords */

  header = (char *) malloc ((nkeys + 1) * FLEN_CARD * sizeof (char));
  p = header;
  for (i = 0; i < nkeys; i++)
    {                   /* Read and print each keywords */

      if (fits_read_record (input, i + 1, card, &f_status))
        break;

      sprintf(p,"%-80s", card);
      p = p + 80;
    }
  /* Add END keyword */

  sprintf (card, "END");
  sprintf(p,"%-80s", card);

  return header;
}

/**
 * Function: get_wcs
 * Returns a pointer to a WoldCoor structure filled using the extension
 * hdunum of file filename.
 * Sets things to use linear WCS and J2000 equinox.
 *
 * Parameters:
 * @param filename The name of the FITS file
 * @param hdunum the number of the HDU to access
 *
 * Return:
 * @return a pointer to a WoldCoor structure. NULL if no WCS was found
 */
struct WorldCoor * get_wcs_l (char filename[], int hdunum)
{
  char *header;
  struct WorldCoor *wcs;

  // Read the FITS or IRAF image file header
  header = get_fits_header_l (filename, hdunum);

  // extract the WCS from the header
  wcs = wcsinit (header);

  // release memory
  free (header);

  // return the result
  return wcs;
}
