/**
 * Interface of spce_pathlength.h
 */

#ifndef _SPCE_PATHLENGTH_H

#define _SPCE_PATHLENGTH_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "spc_trace_functions.h"

extern int
abscissa_to_pathlength (const trace_func * const func, gsl_vector * const data);

#endif /* !_SPCE_PATHLENGTH_H */
