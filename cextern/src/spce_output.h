/**
 * Interface definitions for the postscript output routines
 */

#ifndef _SPCE_OUTPUT_H
#define _SPCE_OUTPUT_H
#include "aper_conf.h"
#include <gsl/gsl_matrix.h>

/*
 * Structure to store all information on
 * a stamp image which is going to be
 * created. This information is steps in x
 * start coos in x,y and image dimensions.
 */
typedef struct
{
  double resolution;

  double xstart;
  double ystart;

  long   xsize;
  long   ysize;
}
drzstamp_dim;


/*
 * Quadrangle structure to be filled with the
 * corners of a pixel in the drizzled coo-system.
 * Those corners are then passed to the boxer
 * subroutine to determine the overlap for a
 * particular pixel
 */
typedef struct
{
  double xmax;  // the maximum value in x[]
  double xmin;  // the minimum value in x[]
  double ymax;  // the maximum value in y[]
  double ymin;  // the minimum value in y[]

                // (x[0],y[0]), (x[1],y[1])
  double x[4];  // (x[2],y[2]) and (x[3],y[3])
  double y[4];  // are the corner points of a pixel
                // in the coo-system of the drizzled image
}
quadrangle;

typedef struct
{
  gsl_matrix *counts;
  gsl_matrix *weight;
}
drzstamp;


/*
 * Structure for all stamp images
 * of a beam in 'drzprep'
 */
typedef struct
{
  gsl_matrix *counts;
  gsl_matrix *error;
  gsl_matrix *cont;
  gsl_matrix *model;
  gsl_matrix *vari;
}
drzprep;


extern gsl_vector_int *
get_trace_inds (const ap_pixel * const ap_p);

extern gsl_matrix *
stamp_img (const ap_pixel * const ap_p, float width, d_point *stp_min);

extern drzprep *
stamp_img_drzprep (const int opt_extr, const ap_pixel * const ap_p, const ap_pixel * const se_p,
		   float width, float nullval, int usemode,
		   drzstamp_dim dimension, gsl_matrix *drzcoeffs,
		   double exptime, double sky_cps, double rdnoise,
		   const int bckmode);

extern void
free_drzprep(drzprep *drzprep_stamps);

extern gsl_matrix *
rectified_stamp_img (const ap_pixel * const ap_p, float width, d_point *stp_min);

extern d_point
get_minxy_from_PET(const ap_pixel * const ap_p);

extern d_point
get_maxxy_from_PET(const ap_pixel * const ap_p);

extern drzstamp *
drizzled_stamp_img (const ap_pixel * const ap_p, double width,
		    const double orient, const drzstamp_dim dimension);

extern drzstamp_dim
get_drzprep_dim(const ap_pixel *const ap_p, float width,
		int boxwidth, int boxheight);

extern drzstamp_dim
get_stamp_dim(const ap_pixel * const ap_p, float width,
	      aperture_conf *conf, const int beamID, d_point *stp_min);

extern quadrangle
get_quad_from_pixel(const ap_pixel *cur_p, const double orient, const drzstamp_dim dimension);

extern gsl_matrix *
drizzled_stamp_img_orig (const ap_pixel * const ap_p, float width,
		    aperture_conf *conf);

extern void
interpolate_over_NaN (gsl_matrix *data);

extern void
free_stamp_img(gsl_matrix *stamp);

extern void
free_drzstamp(drzstamp *stamp);
#endif
