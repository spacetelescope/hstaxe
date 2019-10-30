/**
 * The interface to the binning and weight routines
 */
#ifndef _SPCE_BINNING_H
#define _SPCE_BINNING_H 1

#include <math.h>
#include <gsl/gsl_vector.h>

#include "spce_output.h"
#include "spc_trace_functions.h"
#include "aXe_grism.h"
#include "spc_spc.h"


/**
  defines a pixel for weighting purposes.  Here, p0, p1, p2, p3 are
  the coordinates of the points where the square first contribute,
  where its contribution becomes constant, where it slopes down again,
  and where it doesn't contribute any more.  For example, a pixel at (0,0)
  viewed at from a side has p0=p1=0, p2=p3=1, whereas the same pixel
  viewed from a corner has p0=-1, p1=p2=0, p3=1;  slope is the
  ascent on the non-constant parts, fmax the maximum contribution.
  angle is the angle viewing.  Its tangent tana or its cotangent
  cota may be undefined for certain viewing angles.  The function
  weight_function is selected such that this is unimportatnt.
  x0, y0, and size are the coordinates of the lower left corner and
  the size of the square.
*/
typedef struct w_pixel_s
{
     double p0, p1, p2, p3;
     double angle;
     double tana, cota;
     double x0, y0, size;
     double fmax, slope;
     double (*weight_function) (const double x1, const double y1,
				const double x2, const double y2,
				const struct w_pixel_s * const pix);
}
w_pixel;

extern spectrum *
bin_naive (const ap_pixel * const ap_p, const double ob_width,
	   const double ob_orient, const int quant_cont);

extern spectrum *
bin_weighted (const ap_pixel * const ap_p, const double ob_orient,
	      const trace_func * const trace, const int n_sub,
	      const int flags);

extern spectrum *
bin_optimal (const ap_pixel * const ap_p, const beam curbeam,
	     const int quant_cont, const gsl_matrix *weights,
	     const drzstamp_dim dimension, gsl_matrix *coverage);

#endif
