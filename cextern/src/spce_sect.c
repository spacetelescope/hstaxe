/*
 * File: scpe_sect.c
 * A function find the abscissa of the section between a curve and
 *   a line of given inclination through a number of points
 *   (and its helpers)
 *
 * Usage:
 *  sectionfun sf;
 * double res;
 *  if (fill_in_sectionfun(&sf, init_lo, init_hi, inclination, 
 *    curve_func, curve_data))
 *    # init_hi and init_lo are the starting values for the bracketing 
 *    Error
 *  for (i=0; i<20; i++)
 *    if (!find_section_point(&sf,i,i,&res))
 *      printf("%f\n",res);
 *  free_sectionfun(&sf);
 *
 *  The coordinates passed to find_section_point should be relative
 *  to the reference point within the framework of spc_extract
 *
 *
 */



#include "spce_sect.h"

/**
  computes the value of the section function.

  @param x abscissa at which to evaluate (should be relative to the beam's
  refpoint).
  @param data a sectionfun struct (declared void to please gsl).
  @return the value of the section function.
*/
static double
compute_sectionfun (double x, void *data)
{
  const sectionfun *const sf = data;
  double trace_y, sect_y;
  
  trace_y = sf->func->func (x, sf->func->data);
  sect_y = sf->m * (x - sf->x0) + sf->y0;
  return sect_y - trace_y;
}


/**
  finds the section point between a line of given slope through a given 
  point and a not too unreasonable function.

  @param sf the section function descriptor, created by fill_in_sectionfun.
  @param x abscissa of the given point.
  @param x ordinate of the given point.
  @param res abscissa of section point.
  @returns 0 if section point could be found, -1 if not.
*/
int
find_section_point (sectionfun *sf, const double x, const double y,
		    double *const res)
{
  int iter_count = 0;
  double x_lower, x_upper;
  
  if (sf->vertical)
    {
      *res = x;
      return 0;
    }
  
  sf->x0 = x;
  sf->y0 = y;
  
  if (compute_sectionfun (sf->interv->lower, sf) *
      compute_sectionfun (sf->interv->upper, sf) >= 0)
    {
      return -1;
    }
  gsl_root_fsolver_set (sf->solver, sf->gslfun, sf->interv->lower,
			sf->interv->upper);

  while (1)
    {
      gsl_root_fsolver_iterate (sf->solver);
      x_lower = gsl_root_fsolver_x_lower (sf->solver);
      x_upper = gsl_root_fsolver_x_upper (sf->solver);
      if (gsl_root_test_interval (x_lower, x_upper, 0.0001, 0.000001) ==
	  GSL_SUCCESS)
	break;
      
      if (iter_count++ > 100)
	{
	  return -1;
	}
    }

  *res = gsl_root_fsolver_root (sf->solver);
  return 0;
}


/**
  prepares a sectionfun structure

  @param sf a pointer to the structure to fill
  @param inclination the slope of the line (orientation of the object)
  @param b the current beam 
  @return 0 if successful, -1 otherwise
*/
int
fill_in_sectionfun (sectionfun * const sf, const double inclination,
		    const beam * const b)
{
  gsl_function *gf;
  gsl_interval *interv;
  int i;
  double delta_x, delta_y, slope, x0, y0, x, xr, yr;

  if (fabs (inclination) < 0.001)
    {
      sf->vertical = 1;
      return 0;
    }
  else
    {
      sf->vertical = 0;
    }
  
  sf->m = tan (inclination);
  sf->func = b->spec_trace;
  
  if (!(gf = malloc (sizeof (gsl_function))))
    {
      return -1;
    }
  gf->function = &compute_sectionfun;
  gf->params = sf;
  
  if (!(interv = malloc (sizeof (gsl_interval))))
    {
      free (gf);
      gf = NULL;
      return -1;
    }

  interv->lower = b->bbox[0].x - b->refpoint.x;
  interv->upper = b->bbox[1].x - b->refpoint.x;
  
  /* find out slope of object w.r.t. the spectrum trace */
  xr = b->refpoint.x;
  yr = sf->func->func (0, sf->func->data) + b->refpoint.y;
  delta_y =
    b->spec_trace->func (interv->lower,
			 b->spec_trace->data) -
    b->spec_trace->func (interv->upper, b->spec_trace->data);
  delta_x = interv->lower - interv->upper;
  if (fabs (delta_x) < 1e-7)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Vanishing bounding box");
    }
  slope = delta_y / delta_x;
  if (fabs (sf->m - slope) < 1e-7)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Spectrum trace and object"
		   " orientation are almost parallel");
    }

  /* extend interval such that the section point should be found at
     least when everything is linear */
  for (i = 0; i < 4; i++)
    {
      x0 = b->corners[i].x - xr;
      y0 = b->corners[i].y - yr;
      x = (y0 - x0 * sf->m) / (slope - sf->m);
      if (x < interv->lower)
	interv->lower = x;
      if (x > interv->upper)
	interv->upper = x;
    }
  
  /* choose something harmless to initialize the solver */
  sf->x0 = 1;
  sf->y0 = 1;
  if (compute_sectionfun (interv->lower, sf) *
      compute_sectionfun (interv->upper, sf) >= 0)
    {
      //aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "Invalid spectrum trace "
      //  "function");
    }
  sf->gslfun = gf;
  sf->interv = interv;
  sf->solver = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  return 0;
}


/**
  releases the memory allocated in the fields of a sectionfun structure
  (but not the structure itself)

  @param sf the struct the children of which should be freed
*/

void
free_sectionfun (sectionfun * const sf)
{
  if (sf->vertical)
    {
      return;
    }
  gsl_root_fsolver_free (sf->solver);
  free (sf->interv);
  sf->interv = NULL;
  free (sf->gslfun);
  sf->gslfun = NULL;
}
