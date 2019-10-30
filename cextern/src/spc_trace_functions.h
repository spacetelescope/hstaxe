/**
\ * External definitions to parametrize spectrum traces
 * in aXe grism exposures
 */

#ifndef _SPC_TRACE_FUNCTIONS_H

#define _SPC_TRACE_FUNCTIONS_H

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>

#include "aXe_errors.h"

/* public */

/** A function for parametrizing spectrum traces; deriv and/or path_len
   (i.e. its derivative or the path length from parameter 0) may be
   NULL if they are not known analytically.
   func should return the the y offset for a given pixel column x of
   the spectrum trace 
*/
typedef struct
{
  int type;			/* A unique type ID for the trace function */
  double (*func) (const double, const void *const pars);
  double (*deriv) (const double, const void *const pars);
  double (*path_len) (const double, const void *const pars);
  void *data;		/* private data for trace function, should
			   be passed as a second parameter */
}
trace_func;


/* The create_xxx functions return a pointer to an allocated trace_func
	 or NULL if the allocation failed.  */
extern trace_func *
create_poly2 (const double a0, const double a1, const double a2);

extern trace_func *
create_polyN (gsl_vector *v);

extern trace_func *
vector_to_trace_poly2 (gsl_vector *v);

extern trace_func *
vector_to_trace_polyN(gsl_vector *v);

extern double
polyN_ds(double x, void *pars);

extern void
free_poly2(trace_func * func);

extern void
free_polyN(trace_func * func);

#endif /* ! _SPC_TRACE_FUNCTIONS_H */
