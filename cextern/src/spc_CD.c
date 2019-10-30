/**
 * WCS related routines
 *
 */
#include "spc_CD.h"

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
get_fits_header (char filename[], int hdunum)
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
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_fits_header: Could not open" " file:", filename);
    }
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_fits_header: "
		   "Could not read extention %d from file: %s", hdunum,
		   filename);
    }
  


  fits_get_hdrspace (input, &nkeys, NULL, &f_status);	/* get # of keywords */
  
  header = (char *) malloc ((nkeys + 1) * FLEN_CARD * sizeof (char));
  p = header;
  for (i = 0; i < nkeys; i++)
    {			/* Read and print each keywords */
      
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
struct WorldCoor *
get_wcs (char filename[], int hdunum)
{
  char *header;
  struct WorldCoor *wcs;
  
  // Read the FITS or IRAF image file header
  header = get_fits_header (filename, hdunum);
  
  // extract the WCS from the header
  wcs = wcsinit (header);
  
  // release memory
  free (header);

  // return the result
  return wcs;
}



/**
 * Function: ij_to_radec
 * Returns the world coordinates corresponding to the 
 * point ij in an image which as a WCS wcs.
 *
 * Parameter:
 * @param wcs, a pointer to a WorldCoor WCS structure
 * @param ij, a point inthe image
 *
 * Return:
 * @return a pointer to a d_point structure containing the world
 *  coordinates at point ij in degrees.
 */
sky_coord *
ij_to_radec (struct WorldCoor * wcs, d_point ij)
{
  sky_coord *pos;

  pos = malloc (sizeof (d_point));
  pix2wcs (wcs, ij.x, ij.y, &pos->ra, &pos->dec);
  
  return pos;
}


/**
 * Function: radec_to_ij
 * Returns the detector pixel coordinates corresponding to the 
 * point radec in the sky sgiven the WCS wcs.
 *
 * Parameter:
 * @param wcs a pointer to a WorldCoor WCS structure
 * @param radec a point with celestial coordinates
 *
 * Return:
 * @return a pointer to a d_point structure containing the image
 *  pixel coordinates.
 */
d_point *
radec_to_ij (struct WorldCoor * wcs, sky_coord radec)
{
  d_point *pos;
  int offscl;
  
  pos = malloc (sizeof (d_point));
  wcs2pix (wcs, radec.ra, radec.dec, &pos->x, &pos->y, &offscl);
  
  return pos;
}
