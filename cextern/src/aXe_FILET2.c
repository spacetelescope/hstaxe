/*
 */
#include <stdio.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
  fitsfile *infptr, *outfptr1, *outfptr2, *outfptr3;   /* FITS file pointers defined in fitsio.h */
  int status = 0;       /* status must always be initialized = 0  */
  int hdu_start=0;

  hdu_start = atoi(argv[2]);

  /* Open the input file */
  if ( !fits_open_file(&infptr, argv[1], READONLY, &status) )
    {
      /* Create the output file */
      fits_create_file(&outfptr1, argv[3], &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_movabs_hdu(infptr, hdu_start, NULL, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_copy_hdu(infptr, outfptr1, 0, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_close_file(outfptr1,  &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}

      fits_create_file(&outfptr2, argv[4], &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_movrel_hdu (infptr, 1, NULL, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_copy_hdu(infptr, outfptr2, 0, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_close_file(outfptr2,  &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}

      fits_create_file(&outfptr3, argv[5], &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_movrel_hdu (infptr, 1, NULL, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_copy_hdu(infptr, outfptr3, 0, &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}
      fits_close_file(outfptr3,  &status);
      if (status){
	fits_report_error(stderr, status);
	return(status);}

      fits_close_file(infptr, &status);
    }

  /* if error occured, print out error message */
  if (status) fits_report_error(stderr, status);
  return(status);
}
