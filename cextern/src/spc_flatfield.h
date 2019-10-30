/**
 * Interface to flatfielding
 */

#ifndef _SPC_FLATFIELD_H
#define _SPC_FLATFIELD_H 1

#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "aXe_grism.h"
#include "aXe_utils.h"


typedef struct
{
  int poly_order;      /* order of the polynomial */
  gsl_matrix **coeffs; /* the list of gsl_matrix'es containing */
                       /* the polynomial coefficients of the FF */
  double wmin;         /* Minumum wavelength value for the coeff. normalization */
  double wmax;         /* Maximum wavelength value for the coeff. normalization */
}
poly_cube_flatfield;

/**
  Descriptor of flatfield with a single flatfield image and a
  polynomial flatfield dependence.
*/
typedef struct
{
  gsl_matrix *flatfield; /* The flatfield image, giving the relative
                            efficiency for each pixel in the aperture */
  gsl_matrix *errors;	 /* relative errors for flatfield */
  int poly_order;	 /* order of the polynomial */
  double *poly_coeffs;	 /* coefficients for the polynom */
}
polynom_flatfield;


/**
   Descriptor of flatfield with multiple flatfields at different
   wave lengths
*/
typedef struct
{
  int num_flats;                /* number of (lambda,flatfield) pairs */
  gsl_matrix **flatfields;      /* flatfield images, giving the relative
				   efficiency for each pixel in the aperture */
  gsl_matrix **errors;	        /* relative errors for flatfields */
  double *lambdas;              /* wave lengths for entries in flatfields */
}
multi_flatfield;




/**
  Descriptor of a flatfielding function, including the function
  to compute the flat field value of a given pixel, the bounding box
  we handle, and auxillary data (usually a flatfield image)
*/
typedef struct s_flatfield
{
  void (*func) (const double lambda, const int x, const int y,
		PIXEL_T * const val, PIXEL_T * const err,
		const struct s_flatfield * const self);
  int ll_x, ll_y, w, h;	/* bounding box we handle */
  union
  {
    polynom_flatfield *poly;
    multi_flatfield *multi;
  }
  data;
}
flatfield_d;


poly_cube_flatfield * load_flat_poly_cube(char *fname);
void free_flat_poly_cube(poly_cube_flatfield *poly_cube);

flatfield_d *make_poly_flatfield (const beam * const curbeam,
				  const char *const flat_name,
				  const int order,
				  const double *const coeffs);
void free_poly_flatfield (flatfield_d * const flat);
flatfield_d *make_multi_flatfield (const beam * const curbeam,
				   const int num_flats, const double lambda,
				   const char *const flat_name, ...);
void free_multi_flatfield (flatfield_d * const flat);

void apply_flatfield (ap_pixel * const ap_p, const flatfield_d * const flat);
double poly_cube_flatfield_lambda (const double lambda, const int x, const int y,
				   poly_cube_flatfield *poly_cube);

#endif
