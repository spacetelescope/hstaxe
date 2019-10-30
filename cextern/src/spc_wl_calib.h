#ifndef _SPC_WL_CALIB_H
#define _SPC_WL_CALIB_H 1
#include "disp_conf.h"

/**
 * A calibration function, mostly containing the coefficients of
 * the polynomial.  To be filled in with create_calib
 */
typedef struct
{
  int order;	       /* The order of the polynomial */
  gsl_vector *pr_range; /* the valid range for prisms */
  double *coeffs; /* a pointer order+1 doubles, 0th coeff at coeffs[0] */
  double (*func) (const double xi, const int order,
		  const double *const coeffs);
}
calib_function;


calib_function *
create_calib_from_gsl_vector (const int for_grism, const gsl_vector * a);

extern void
free_calib (calib_function * const calib);

extern void
wl_calib (ap_pixel * ap_p, const calib_function * const calib);

extern void
pwise_wl_calib(const global_disp *gdisp, const d_point pixel,
	       const beam actbeam, const int for_grism, ap_pixel *ap_p, const calib_function *old_calib);

extern ap_pixel *
prange_cut(ap_pixel * ap_p, const gsl_vector * lambda_range,
	   const calib_function * wl_calibration);

#endif
