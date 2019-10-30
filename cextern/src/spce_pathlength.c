/**
 * Functions to replace an x coordinate with the path length of a given
 *   function of that x coordinate.
 *
 * Sample usage:
 *  trace_func *f;
 *  gsl_vector *data;
 *  
 *  abscissa_to_pathlength(f, data);
 *  
 *  So, this is a little trivial...
 *
 */

#include "spce_pathlength.h"


/**
  transforms a list of abscissas to a list of path lengths if the path length
  is known analytically.

  @param func the spectrum trace function to use in the transformation
  @param data the table of abscissas relative to the reference point
*/
static int
absc_to_pathl_path_len (const trace_func * const func,
			gsl_vector * const data)
{
  int i;

  for (i = 0; i < (int)data->size; i++)
    {
      gsl_vector_set (data, i,
		      func->path_len (gsl_vector_get (data, i),
				      func->data));

    }
  return 0;
}

/* computes sqrt(1+f'(x)^2), i.e. the function to integrate to obtain
 *   the path length.
 * 
 * @param x the abscissa to compute the function for
 * @pars a pointer to the trace_func struct, declared void to please gsl
 */
static double
pathlength_fun (double x, void *pars)
{
  const trace_func *const func = pars;
  double y = func->deriv (x, func->data);
  
  return sqrt (1 + y * y);
}


/* Compute the path length of (x,func->func(x)) from x0 to x1 (helper
 *   function for absc_to_pathl_deriv)
 * 
 * @param func the spectrum trace
 * @param x0 start of parameter range
 * @param x1 end of parameter range
 * 
 */
static double
pathlength_integral (const trace_func * const func, const double x0,
		     const double x1)
{
  double res, abserr;
  size_t neval;
  gsl_function f;
  
  f.function = &pathlength_fun;
  f.params = (void *) func;
  
  gsl_integration_qng (&f, x0, x1, 1e-6, 1e-7, &res, &abserr, &neval);
  return res;
}


#define NUM_INTERPOINTS 50	/* number of interpolation points */
/**
 *  transforms a list of abscissas to a list of path lengths if the derivative
 * but not the path length is known analytically.  The strategy is to compute
 * NUM_INTERPOINTS actual path lengths and the interpolate between them.
 * This should work with sufficient accuracy if the trace function is
 * reasonably smooth.
 *
 * @param func the spectrum trace function to use in the transformation
 * @param data the table of abscissas relative to the reference point
 */
static int
absc_to_pathl_deriv (const trace_func * const func, gsl_vector * const data)
{
  int i;
  double min_absc = gsl_vector_min (data), max_absc =
    gsl_vector_max (data);
  double step;
  double xvals[NUM_INTERPOINTS], yvals[NUM_INTERPOINTS];
  double res;
  gsl_spline *spline;
  gsl_interp_accel *acc;
  
  //gsl_interp_factory factory = gsl_interp_factory_cspline_natural;
  //gsl_interp_obj *interpolator;
  //gsl_interp_accel *accelerator;



  /* thresh in a little space left and right to avoid interpolation
     disasters at the edges of the interval */
  min_absc = min_absc - (max_absc - min_absc) / NUM_INTERPOINTS * 2;
  max_absc = max_absc + (max_absc - min_absc) / NUM_INTERPOINTS * 2;
  step = (max_absc - min_absc) / NUM_INTERPOINTS;
  xvals[0] = min_absc;
  yvals[0] = pathlength_integral (func, 0, xvals[0]);
  for (i = 1; i < NUM_INTERPOINTS; i++)
    {
      xvals[i] = min_absc + i * step;
      yvals[i] =
	yvals[i - 1] + pathlength_integral (func, xvals[i - 1],
					    xvals[i]);
    }
  
  //interpolator = factory.create(xvals, yvals, NUM_INTERPOINTS);
  //accelerator = gsl_interp_accel_new();
  
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, NUM_INTERPOINTS);
  gsl_spline_init (spline, xvals, yvals, NUM_INTERPOINTS);

  for (i = 0; i < (int)data->size; i++)
    {
      //gsl_interp_eval_impl(interpolator, xvals, yvals, 
      //  gsl_vector_get(data, i), accelerator, &res);
      //res = gsl_spline_eval (spline, gsl_vector_get (data, i), acc);
      res = spline->interp->type->eval(spline->interp->state, spline->x, spline->y, spline->interp->size, gsl_vector_get (data, i), acc, &yvals[i]);
      gsl_vector_set (data, i, res);
    }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}


/**
 * transforms a list of abscissas to a list of path lengths if not even
 * the derivative is known.
 *
 * @param func the spectrum trace function to use in the transformation
 * @param data the table of abscissas relative to the reference point
 * @todo pathlength for traces without derivatibes is not implemented
 */
static int
absc_to_pathl_val_only (const trace_func * const func,
			gsl_vector * const data)
{
  fprintf (stderr, "Not implemented yet\n");
  abort ();
  return 0;
}


/**
 * transforms a list of abscissas to a list of path lengths.  This function
 * decides what of the absc_to_pathl_len, absc_to_pathl_deriv, or
 * absc_to_pathl_val_only to use.
 *
 * @param func the spectrum trace function to use in the transformation
 * @param data the table of abscissas relative to the reference point
 */
int
abscissa_to_pathlength (const trace_func * const func,
			gsl_vector * const data)
{
  if (func->path_len)
    {
      return absc_to_pathl_path_len (func, data);
    }
  if (func->deriv)
    {
      return absc_to_pathl_deriv (func, data);
    }
  return absc_to_pathl_val_only (func, data);
}
