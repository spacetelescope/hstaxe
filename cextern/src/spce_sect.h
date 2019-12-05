/*
 * File: spce_sect.h
 * External definitions for spce_sect.c
 *
 *
 */

#ifndef spce_sect_H
#define spce_sect_H

#include <math.h>
#include <gsl/gsl_roots.h>
#include "aXe_grism.h"
#include "spc_trace_functions.h"


/**
  @package spce_sect
*/

/** 
    A structure to support old gsl 0.9 gsl_interval
*/
typedef struct
{
     double lower;
     double upper;
}
gsl_interval;

/** 
  A structure describing the function that has a zero at the section point,
  including the measuring point and the slope of the line through it, 
  the function describing the spectrum trace, functions related to
  gsl and a flag to initiate special handling of vertically oriented
  objects.

  @see fill_in_sectionfun
  @see free_sectionfun
*/
typedef struct
{
  int vertical;	     /* Special handling if orientation is close to vertical */
  double m;	     /* slope of line through x0, y0 that is to */
  double x0, y0;     /* intersect the trace */
  trace_func *func;  /* Parametrization of the trace */

     /* the GSL stuff has to be kept in here to avoid excessive re-allocing
        of the solver for each pixel */
  gsl_interval *interv;
  gsl_function *gslfun;
  gsl_root_fsolver *solver;
}
sectionfun;

/* public */

int fill_in_sectionfun (sectionfun * const sf, const double inclination,
			const beam * const b);
int find_section_point (sectionfun *sf, const double x, const double y,
			double *const res);
void free_sectionfun (sectionfun * const sf);
#endif
