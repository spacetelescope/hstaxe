/**
 */
#ifndef _CROSSDISP_UTILS_H
#define _CROSSDISP_UTILS_H


#include "spc_trace_functions.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include <math.h>
#include "trace_conf.h"

/**
 * Structure: drz_pars
 *  A description of an ellipse with a,b, theta.
 */
typedef struct
{
  gsl_matrix *coeffs;
  d_point     offset;
  px_point    npixels;
}
drz_pars;

int
drizzle_distort(const gsl_vector * x, void *params, gsl_vector * f);


d_point
undistort_point(gsl_matrix *coeffs, const px_point pixmax,
                  d_point xy_image);

d_point
distort_point(gsl_matrix *coeffs, const px_point pixmax, d_point xy_image);

extern gsl_matrix *
get_crossdisp_matrix(char * filename, int sci_numext);

extern d_point
get_axis_scales(beam actbeam, gsl_matrix * drzcoeffs, px_point pixmax);

extern double 
get_crossdisp_scale(trace_func  *trace, d_point refpnt, gsl_matrix * drzcoeffs, px_point pixmax);

extern d_point
get_drz_position(d_point in_point, gsl_matrix * drzcoeffs, px_point pixmax);

extern d_point
get_drz_position_free(d_point in_point, gsl_matrix * drzcoeffs, px_point pix_in, px_point pix_out);

extern trace_func *
get_tracefunc_at(char *filename, d_point p);

extern double 
evaln(double x, double y, gsl_matrix * drzcoeffs, int row);

extern double 
devalndx(double x, double y, gsl_matrix * drzcoeffs, int row);

extern double 
devalndy(double x, double y, gsl_matrix * drzcoeffs, int row);

extern gsl_matrix *
get_jacobian(int i, int j, gsl_matrix * drzcoeffs, int width, int height);

extern double 
get_det_jacobian(int i, int j, gsl_matrix * drzcoeffs, int width, int height);

#endif
