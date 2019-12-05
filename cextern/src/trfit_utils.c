/**
 * File: trfit_utils.c
 *
 */
#include "trfit_utils.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_sect.h"
#include "spce_fitting.h"
#include "spce_is_in.h"
#include "spc_back.h"
#include "spce_pathlength.h"
#include "lmmin.h"
#include "lm_eval.h"

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

int gauss_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  int ndata      = ((struct function_data *)params)->n;
  double *x_data = ((struct function_data *)params)->x;
  double *y_data = ((struct function_data *)params)->y;
  //double *e_data = ((struct function_data *)params)->sig;

  double a  = gsl_vector_get (x, 0);
  double b  = gsl_vector_get (x, 1);
  double c  = gsl_vector_get (x, 2);
  double x0 = gsl_vector_get (x, 3);

  double Yi = 0.0;

  int i;


  for (i = 0; i < ndata; i++)
    {
      /* Model Yi = a * exp(-b*(Xi-X0) + c */
      Yi = a * exp (-b * (x_data[i]-x0)*(x_data[i]-x0)) + c;
      //      gsl_vector_set (f, i, (Yi - y_data[i])/e_data[i]);
      gsl_vector_set (f, i, Yi - y_data[i]);
      //      fprintf(stderr, "Yi: %e, f: %e\n", Yi, (Yi - y_data[i])/e_data[i]);
    }

  return GSL_SUCCESS;
}



int
gauss_df(const gsl_vector * x, void *params,
	 gsl_matrix * J)
{
  int ndata      = ((struct function_data *)params)->n;
  double *x_data = ((struct function_data *)params)->x;
  //double *y_data = ((struct function_data *)params)->y;
  //double *e_data = ((struct function_data *)params)->sig;

  double a  = gsl_vector_get (x, 0);
  double b  = gsl_vector_get (x, 1);
  //double c  = gsl_vector_get (x, 2);
  double x0 = gsl_vector_get (x, 3);

  int i;
  double value;
  double xdiff;

  for (i = 0; i < ndata; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /* Model Yi = a * exp(-b*(Xi-X0) + c   */
      /* and the xj are the parameters (a,b,c,x0) */

      xdiff = x_data[i]-x0;

      value = exp(-b*(xdiff)*(xdiff));

      gsl_matrix_set (J, i, 0, value);
      gsl_matrix_set (J, i, 1, -a*(xdiff)*(xdiff)*value);
      gsl_matrix_set (J, i, 2, 1.0);
      gsl_matrix_set (J, i, 3,a*b*value*2*xdiff);
    }
  return GSL_SUCCESS;
}

int
gauss_fdf (const gsl_vector * x, void *params, gsl_vector *f, gsl_matrix *J)
{
  gauss_f(x, params, f);
  gauss_df(x, params, J);

  return GSL_SUCCESS;
}

/**
 * Function: fit_beamtrace
 * The function determines a new trace solution for a given beam. A gaussian
 * is fitted to the pixel values around the nominal trace position.
 * A kappa-sigma-klipping algorithm eliminates fitted trace centers
 * which are off due to cosmics or similar.
 * Finally the a new trace solution is fitted to the set of trace centers
 * The new slope is directly written to the beam structure, the
 * new y-value of the reference point is returned.
 *
 * Parameters:
 * @param conf     - the aXe configuration structure
 * @param obs      - the current observation
 * @param beamID   - the ID of the current beam
 * @param act_beam - the beam to fit
 *
 * Returns:
 * @return yf - the y-value of the reference point for the fitted trace
 */

gsl_vector *
fit_beamtrace(const aperture_conf  *conf, observation *obs,
	      const int beamID, beam act_beam)
{
  //gsl_vector_int *yvec;
  gsl_vector     *fit_params;
  gsl_vector     *lin_fit;
  gsl_vector     *fit_result;
  trace_func     *tracefun;

  px_point        xborder;
  px_point        tpoint;

  fit_data       *f_data;

  double *x;
  double *y;
  double *w;

  double *trace;

  double yf, yf_err;

  int i;
  int index;

  // This is the return result
  fit_result = gsl_vector_alloc(5);
  gsl_vector_set(fit_result, 0, 0);
  gsl_vector_set(fit_result, 1, 0);
  gsl_vector_set(fit_result, 2, 0);
  gsl_vector_set(fit_result, 3, 0);
  gsl_vector_set(fit_result, 4, 0);

  // define the beam and the trace function
  tracefun = act_beam.spec_trace;
  trace = (double *)(tracefun->data);

  // If this beam's ignore flag is set to 1 then do nothing
  if (act_beam.ignore == 1)
    return fit_result;

  // determine the start and end point in x
  //xborder = get_fit_xrange(conf, obs, beamID, act_beam);
  xborder = get_fit_xrange(conf, obs, act_beam);

  // allocate arrays for the fitting
  x = (double*)malloc((xborder.y-xborder.x+1)*sizeof(double));
  y = (double*)malloc((xborder.y-xborder.x+1)*sizeof(double));
  w = (double*)malloc((xborder.y-xborder.x+1)*sizeof(double));

  // Loop over all columns
  index = 0;
  for (i = xborder.x; i < xborder.y; i++)
    //for (i=130; i < 150; i++)
    {
      // determine the pixel which is
      // closest to the nominal trace position
      tpoint.x = i;
      tpoint.y = (int)floor(tracefun->func((double)i-act_beam.refpoint.x,
					   tracefun->data)
			    + act_beam.refpoint.y+0.5);


      // get the data to be fitted
      // in the current column
      f_data = get_fitdata(obs, N_COLUMNS_FIT, tpoint);
      /*
      if (tpoint.x > 130 && tpoint.x < 150)
	{
	  fprintf(stdout, "x: %i\n", tpoint.x);
	  print_fit_data(f_data);
	  fprintf(stdout, "\n\n");
	}
      if (tpoint.x == 135)
      {*/
      //fprintf(stdout, "x: %i\n", tpoint.x);
      //print_fit_data(f_data);
      //fprintf(stdout, "\n\n");
	  // fit the gaussian to the data
	  fit_params = fit_wuttke(f_data);
	  /*	}
      else{
	fit_params = fit_wuttke_talk(f_data);
	}*/
      // check whether the fit was successful
      if (gsl_vector_get(fit_params, 4) == 1)
	{
	  // transfer the results of the fit
	  // to this column in the vector
	  // which contains all fit results
	  x[index] = (double)tpoint.x;
	  y[index] = gsl_vector_get(fit_params, 3);
	  w[index] = 1.0;
	  index++;
	}

      // freem memory
      free_fit_data(f_data);
      gsl_vector_free(fit_params);
    }

  //for (i=0; i < index; i++)
  //  fprintf(stdout, "%f %f %f\n", x[i], y[i], w[i]);

  // do a kappa-sigma klipping on
  // the array with the trace locations
  // to improve the robustness of the result
  comp_kappasigma_interp(x, y, w, index, 1, N_KAPPASIG_ITER,
			 N_KAPPASIG_SIG, obs, 0);

  // make a linear fit to the trace positions
  lin_fit = det_vector_linear(x, y, w, index, 0);
  //for (i=0; i < index; i++)
  //fprintf(stdout, "%f %f %f\n", x[i], y[i], w[i]);

  // compute the y-position of the beam reference position
  // according to the new trace solution
  gsl_fit_linear_est (act_beam.refpoint.x-gsl_vector_get(lin_fit, 0), gsl_vector_get(lin_fit, 1),
		      gsl_vector_get(lin_fit, 2), gsl_vector_get(lin_fit, 3),
		      gsl_vector_get(lin_fit, 4), gsl_vector_get(lin_fit, 5),
		      &yf, &yf_err);

  gsl_vector_set(fit_result, 0, yf);
  gsl_vector_set(fit_result, 1, yf_err);
  gsl_vector_set(fit_result, 2, gsl_vector_get(lin_fit, 2));
  gsl_vector_set(fit_result, 3, sqrt(gsl_vector_get(lin_fit, 5)));
  gsl_vector_set(fit_result, 4, gsl_vector_get(lin_fit, 6));


  fprintf(stdout, "\nFitting result XXX: xpos = %f , ypos = %f +- %f , slope = %f +- %f , chi^2 = %f, x_pivot = %f, y_abs = %f\n",
	  act_beam.refpoint.x, gsl_vector_get(fit_result, 0), gsl_vector_get(fit_result, 1),
	  gsl_vector_get(fit_result, 2), gsl_vector_get(fit_result, 3), gsl_vector_get(lin_fit, 6),
	  gsl_vector_get(lin_fit, 0), gsl_vector_get(lin_fit, 1));


  // copy the slope of the new trace solution
  // to the according position in the beam description
  trace[1] = gsl_vector_get(lin_fit, 2);
  trace[2] = sqrt(gsl_vector_get(lin_fit, 5));

  // free the memory
  gsl_vector_free(lin_fit);
  free(x);
  free(y);
  free(w);

  // return the new y-position of
  // the reference point
  //  return yf;
  return fit_result;
}

/**
 * Function: gagauss
 * The function computes and returns the value
 * of a gauss function y = a * exp(-b*(x-c)^2) - d
 * at a position x.
 *
 * Parameters:
 * @param x      - the x-poition to evaluate the function
 * @param params - the function parameters
 *
 * Returns:
 * @return - the function value
 */
double
gagauss(const double x, const double *params)
{
  // just compute and return the result
  return (params[0] * exp (-params[1] * (x-params[3])*(x-params[3])) + params[2]);
}

/**
 * Function: comp_shift
 *
 * Returns:
 * @return - the function value
 */
double
comp_shift(const double x, const double *params, const double *fpars)
{
  double dy;
  double tmp;

  // Parameters in the original
  // notation fo W. Freudling:
  // fpars[0] = 'contrast'
  // fpars[1] = 'exp'
  // fpars[2] = 'tan_ang'
  // gpars[0] = 'shift'

  // compute the shift for the
  // x-value
  dy = x * fpars[2] + params[0];

  // compute the sin()-expression
  tmp = (sin((2.0*dy-0.5)*M_PI) + fpars[1])/2.0;

  // compute and return the whole expression
  return tmp * tmp * fpars[0] + (1.0 - fpars[0]);
}


/**
 * Function: fit_wuttke
 * Fits a gaussian to the data given in the input structure.
 * Uses a library for nn-linear fits (made by Wuttke)
 *
 * Parameters:
 * @param f_data - the fit data
 *
 * Returns:
 * @return params - the fit results
 */
gsl_vector *
fit_wuttke(const fit_data *f_data)
{
  gsl_vector *params;

  int m_dat = f_data->n_data;
  int n_p =  4;

  double p[4];

  // auxiliary settings:
  lm_control_type control;
  lm_data_type data;

  // allocate the parameter vector
  params = gsl_vector_alloc(5);

  // give a first guess
  // for the gauss parameters
  p[0] = f_data->y_values[f_data->n_data / 2];
  p[1] = 1.0;
  p[2] = (f_data->y_values[f_data->n_data-1]-f_data->y_values[0])/2.0;
  p[3] = f_data->x_values[0] + (f_data->x_values[f_data->n_data-1]-f_data->x_values[0])/2.0;

  // initialize the fit
  lm_initialize_control( &control );

  // define the parameters for the fitting
  data.user_func = gagauss;
  data.user_t = f_data->x_values;
  data.user_y = f_data->y_values;

  // perform the fit:
  lm_minimize( m_dat, n_p, p, lm_evaluate_default, lm_print_nothing,
	       &data, &control );

  //fprintf(stdout, "Result: %i, amp: %f, width: %f, const:  %f, pos: %f\n", control.info, p[0], p[1], p[2], p[3]);

  // copy the fit results to the output structure
  gsl_vector_set(params, 0, p[0]);
  gsl_vector_set(params, 1, p[1]);
  gsl_vector_set(params, 2, p[2]);
  gsl_vector_set(params, 3, p[3]);
  gsl_vector_set(params, 4, (double)control.info);

  // return the output structure
  return params;
}
gsl_vector *
fit_wuttke_talk(const fit_data *f_data)
{
  gsl_vector *params;

  int m_dat = f_data->n_data;
  int n_p =  4;

  double p[4];

  // auxiliary settings:
  lm_control_type control;
  lm_data_type data;

  // allocate the parameter vector
  params = gsl_vector_alloc(5);

  // give a first guess
  // for the gauss parameters
  p[0] = 1.0;
  p[1] = 1.0;
  p[2] = 0.0;
  p[3] = f_data->x_values[0] + (f_data->x_values[f_data->n_data-1]-f_data->x_values[0])/2.0;

  // initialize the fit
  lm_initialize_control( &control );

  // define the parameters for the fitting
  data.user_func = gagauss;
  data.user_t = f_data->x_values;
  data.user_y = f_data->y_values;

  // perform the fit:
  lm_minimize( m_dat, n_p, p, lm_evaluate_default, lm_print_default,
	       &data, &control );

  // copy the fit results to the output structure
  gsl_vector_set(params, 0, p[0]);
  gsl_vector_set(params, 1, p[1]);
  gsl_vector_set(params, 2, p[2]);
  gsl_vector_set(params, 3, p[3]);
  gsl_vector_set(params, 4, (double)control.info);

  // return the output structure
  return params;
}


d_point
find_grav_center(const fit_data *f_data)
{
  d_point ret={0.0,0.0};

  int i=0;

  double sum    = 0.0;
  double w_sum    = 0.0;
  double weight = 0.0;

  for (i=0; i < f_data->n_data; i++)
    {
      weight = fabs(f_data->y_values[i]/f_data->e_values[i]);
      sum   += weight * f_data->x_values[i];
      w_sum += weight;
    }

  ret.x = sum/w_sum;
  ret.y = weight;

  return ret;
}

/*
 * Function: get_fit_xrange
 * The function determines and returns the mininmum
 * and the maximum x-values for a beam.
 *
 * Parameters:
 * @param conf     - the aXe configuration structure
 * @param obs      - the current observation
 * @param beamID   - the ID of the current beam
 * @param act_beam - the beam to fit
 *
 * Returns:
 * @return ret - x_min and x_max of the beam
 */
px_point
get_fit_xrange(const aperture_conf  *conf, const observation *obs,
	       beam act_beam)
  {
    px_point ret;

    // compute one extremum
    ret.x = (int)floor(act_beam.refpoint.x
		       + conf->beam[act_beam.ID].offset.dx0 +.5);

    // compute the secund extremum
    ret.y = (int)floor(act_beam.refpoint.x
		       + conf->beam[act_beam.ID].offset.dx1 +.5);

    // limit the start and end column
    // to the image size
    ret.x = MAX(ret.x, 0);
    ret.y = MIN(ret.y, (int)obs->grism->size1);

    // return the min and the max value
    return ret;
  }

gsl_vector *
fit_function(const fit_data *f_data)
{
  gsl_vector *fit_params = gsl_vector_alloc(4);

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  gsl_multifit_function_fdf f;
  gsl_vector *x_init = gsl_vector_alloc(4);

  int status;
  //size_t i, 
  size_t iter = 0;

  const size_t n = f_data->n_data;
  const size_t p = 4;

  //gsl_matrix *covar = gsl_matrix_alloc (p, p);

  struct function_data d = {f_data->n_data, f_data->x_values,
			    f_data->y_values, f_data->e_values};

  fit_params = gsl_vector_alloc(p);

  gsl_vector_set(x_init, 0, 1.0);
  gsl_vector_set(x_init, 1, 1.0);
  gsl_vector_set(x_init, 2, 0.0);
  gsl_vector_set(x_init, 3,
		 f_data->x_values[0]
		 +(f_data->x_values[f_data->n_data-1]-f_data->x_values[0])/2.0);

  f.f      = gauss_f;
  f.df     = gauss_df;
  f.fdf    = gauss_fdf;
  f.n      = n;
  f.p      = p;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, x_init);

  print_fit_state (iter, s);


  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status %i = %s\n", status, gsl_strerror (status));

      print_fit_state (iter, s);


      if (status != GSL_CONTINUE)
       	break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-12, 1e-12);
      printf ("status  a = %s\n", gsl_strerror (status));
      fprintf(stdout, "stat: %i, GSL_CONT: %i\n", status, GSL_CONTINUE);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_fdfsolver_free (s);
  return fit_params;
}

/*
 * Function: get_fitdata
 * The function collects the pixel values of at most np pixels
 * on both sides of the nominal trace position and stores
 * these values in an apppropriate structure.
 *
 * Parameters:
 * @param obs      - the current observation
 * @param np       - the number of pixels to sample around the nominal position
 * @param tr_point - the nominal traxe pixel
 *
 * Returns:
 * @return f_data  - the pixel data around the nominal trace position
 *
 */
fit_data *
get_fitdata(observation *obs, int np, px_point tr_point)
{
  gsl_vector_int *tmp;

  fit_data *f_data;

  int np_act=0;

  //int y_low;
  //int y_upp;

  int l_space=1;
  int u_space=1;

  int l_act;
  int u_act;

  int l_np=0;
  int u_np=0;

  int ncols=obs->grism->size2;

  int ii;
  int index;

  // allocate the vector
  tmp = gsl_vector_int_alloc(2*np+2);
  gsl_vector_int_set_all(tmp, -1);

  // limit the starting ppoint of the search
  // to values within the image dimension
  tr_point.y = MAX(0,tr_point.y);
  tr_point.y = MIN((int)obs->grism->size2,tr_point.y);

  // initialize the row numbers
  // to search up- and downwards
  l_act = tr_point.y - 1;
  u_act = tr_point.y;

  // as long as interpolation points are missing
  // and one direction, either up or down,
  // is 'open', continue searching
  while (np_act < 2*np && (l_space || u_space))
    {

      // check whether the direction
      // downwards is still open
      if (l_space)
	{
	  // if you are at the end of the frame
	  if (l_act < 0)
	    {
	      // close the direction downwards
	      l_space=0;
	    }
	  else
	    {
	      if (!isnan(gsl_matrix_get(obs->grism, tr_point.x, l_act)))
		{
		  // else store the interpolation point,
		  // do the various increments
		  gsl_vector_int_set(tmp, np_act, l_act);
		  //		  fprintf(stderr, "put in %i\n", l_act);
		  np_act++;
		  l_np++;
		}
	      l_act--;
	    }
	}

      // check whether the direction
      // upwards is still open
      if (u_space)
	{
	  // if you are at the end of the frame
	  if (u_act >= ncols)
	    {
	      // close the direction upwards
	      u_space=0;
	    }
	  else
	    {
	      if (!isnan(gsl_matrix_get(obs->grism, tr_point.x, u_act)))
		{
		  // else store the interpolation point,
		  // do the various increments
		  gsl_vector_int_set(tmp, np_act, u_act);
		  np_act++;
		  u_np++;
		}
	      u_act++;
	    }
	}
    }

  // sort the pixel indices
  gsl_sort_vector_int(tmp);

  // allocate a data structure
  f_data = alloc_fit_data(np_act);

  // fill the pixel values into the data structure
  index = 0;
  for (ii=((int)tmp->size - np_act); ii < (int)tmp->size; ii++)
    {
      f_data->x_values[index] = (double)gsl_vector_int_get(tmp, ii);
      f_data->y_values[index] = gsl_matrix_get(obs->grism, tr_point.x,
					       gsl_vector_int_get(tmp, ii));
      f_data->e_values[index] = gsl_matrix_get(obs->pixerrs, tr_point.x,
					       gsl_vector_int_get(tmp, ii));
      index++;
    }

  // release memory
  gsl_vector_int_free(tmp);

  //  return the result
  return f_data;
}

/*
 * Function: get_ipc_coldata
 * The function checks for valid pixel data on both side of the trace.
 * The maximum distance away from the trace is given as a parameter.
 * Valid data is collected and returned in a vector.
 *
 * Parameters:
 * @param img_data - the current observation
 * @param np       - the number of pixels to sample around the nominal position
 * @param tr_point - the nominal traxe pixel
 *
 * Returns:
 * @return ipc_row  - the pixel data around the nominal trace position
 *
 */
gsl_vector *
get_ipc_coldata(gsl_matrix *img_data, int np, px_point tr_point)
{
  gsl_vector *tmp;

  gsl_vector *ipc_row = NULL;

  int np_act=0;

  int y_act;

  int ncols = img_data->size2;

  int ii;

  // check whether the central tracepoint
  // is in the image
  if (tr_point.y < 0 || tr_point.y > ncols-1)
    // give back the NULL array
    return ipc_row;

  // allocate the vector
  tmp = gsl_vector_alloc(2*np+1);
  gsl_vector_set_all(tmp, -1.0);

  // set the current y-position
  y_act = tr_point.y;

  // check whether there is good data
  // at the current position
  if (!isnan(gsl_matrix_get(img_data, tr_point.x, y_act)))
    {
      // store the data point and enhance the counter
      gsl_vector_set(tmp, np_act, gsl_matrix_get(img_data, tr_point.x, y_act));
      np_act++;
    }

  // check on both sides
  // of the central tracepoint
  for (ii=1; ii <= np; ii++)
    {
      // go down
      y_act = tr_point.y - ii;

      // check whether we are still inside the image and
      // whether we have valid data
      if (y_act > -1 && !isnan(gsl_matrix_get(img_data, tr_point.x, y_act)))
	{
	  // store the data point and enhance the counter
	  gsl_vector_set(tmp, np_act, gsl_matrix_get(img_data, tr_point.x, y_act));
	  np_act++;
	}

      // go up
      y_act = tr_point.y + ii;

      // check whether we are still inside the image and
      // whether we have valid data
      if (y_act < ncols && !isnan(gsl_matrix_get(img_data, tr_point.x, y_act)))
	{
	  // store the data point and enhance the counter
	  gsl_vector_set(tmp, np_act, gsl_matrix_get(img_data, tr_point.x, y_act));
	  np_act++;
	}
    }

  // check whether good data was found
  if (np_act > 0)
    {
      // allocate the return vector
      ipc_row = gsl_vector_alloc(np_act);

      // transfer the data to the
      // return vector
      for (ii=0; ii < np_act; ii++)
	gsl_vector_set(ipc_row, ii, gsl_vector_get(tmp, ii));
    }

  // release memory
  gsl_vector_free(tmp);

  //  return the vector
  return ipc_row;
}

/*
 * Function: alloc_fit_data
 * The function allocates memory for a fit-data structure
 *
 * Parameters:
 * @param n_data - the number of points to allocate
 *
 * Returns:
 * @return f_data - the allocated structure
 */
fit_data *
alloc_fit_data(const int n_data)
{
  fit_data *f_data;

  f_data = (fit_data *)malloc(sizeof(fit_data));

  f_data->x_values = (double *)malloc(n_data * sizeof(double));
  f_data->y_values = (double *)malloc(n_data * sizeof(double));
  f_data->e_values = (double *)malloc(n_data * sizeof(double));

  f_data->n_data = n_data;

  return f_data;
}

/*
 * Function: print_fit_data
 * The function prints the content of a fit-data structure
 * onto the screen (error-out).
 *
 * Parameters:
 * @param f_data - the allocated structure
 *
 * Returns:
 * @return -
 */
void
print_fit_data(fit_data *f_data)
{
  int i;

  for (i=0; i < f_data->n_data; i++)
    {
      fprintf(stdout, "%f %f %f\n", f_data->x_values[i],
	      f_data->y_values[i], f_data->e_values[i]);
      //      fprintf(stderr, "x: %f, y: %f, err: %f\n", f_data->x_values[i],
      //	      f_data->y_values[i], f_data->e_values[i]);
    }
}

/*
 * Function: print_fit_data
 * The function releases the memory allocated in
 * a fit-data structure
 *
 * Parameters:
 * @param f_data - the allocated structure
 *
 * Returns:
 * @return -
 */
void
free_fit_data(fit_data *f_data)
{
  free(f_data->x_values);
  free(f_data->y_values);
  free(f_data->e_values);

  free(f_data);

  f_data = NULL;
}

void
print_fit_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3zu x = %e %e %e %e "
          "|f(x)| = %g |dx| = %g\n\n", iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_blas_dnrm2 (s->f),
          gsl_blas_dnrm2 (s->dx));
}

/**
 * Function: comp_intrel_max
 *
 * Parameters:
 * @param f_data  - the input data
 *
 * Returns:
 * @return rel12 - "second brightes row pixel divided by brighest row pixel"
 */
double
comp_intrel_max(const fit_data *f_data)
{

  int i;
  double rel12=-100000.0;

  gsl_vector *val_vector;

  // check whether there is enough
  // data to make the relation
  if (f_data->n_data < 2)
    return rel12;

  // allocate a vector
  val_vector = gsl_vector_alloc(f_data->n_data);

  // fill in the data
  for (i=0; i < f_data->n_data; i++)
    gsl_vector_set(val_vector, i, f_data->y_values[i]);

  // sort the values
  gsl_sort_vector(val_vector);

  // prevent division by zero
  if (gsl_vector_get(val_vector, 0))
    // compute the quantity second brightes pixel
    // divided by brighest pixel
    rel12 = gsl_vector_get(val_vector, f_data->n_data-2) / gsl_vector_get(val_vector, f_data->n_data-1);

  // free the memory
  gsl_vector_free(val_vector);

  // return the result
  return rel12;
}

/**
 * Function: bcksub_observation
 *
 * Parameters:
 * @param obs           - the grism image
 * @param bck           - the background image
 *
 * Returns:
 * @return ret - background subtracted data
 */
gsl_matrix *
bcksub_observation(observation *obs, observation *bck)
{
  gsl_matrix *ret;

  int ii;
  int jj;

  // initialize the return
  ret = gsl_matrix_alloc(obs->grism->size1, obs->grism->size2);

  // set all values to NAN
  gsl_matrix_set_all(ret, GSL_NAN);

  // go over all cols
  for (ii=0; ii < (int)obs->grism->size1; ii++)
    // go over all rows
    for (jj=0; jj < (int)obs->grism->size2; jj++)
      {
	// check whether there is valid data in
	// object and background image
	if (!isnan(gsl_matrix_get(obs->grism, ii, jj))
	    && !isnan(gsl_matrix_get(bck->grism, ii, jj))
	    && gsl_matrix_get(bck->grism, ii, jj))

	  // compute and set the pixel
	  gsl_matrix_set(ret, ii, jj, gsl_matrix_get(obs->grism, ii, jj) - gsl_matrix_get(bck->grism, ii, jj));
      }

  // return matrix
  return ret;
}

/**
 * Function: evaluate_ipc_fit
 *
 * Parameters:
 * @param f_data  - the input for the independent values
 * @param gpars   - the fited parameters
 * @param fpars   - the fixed parameters
 *
 * Returns:
 * @return y_values - the function values
 */
double *
evaluate_ipc_fit(fit_data *f_data, double *gpars, double *fpars)
{
  int index;

  //double yvalue;

  double *y_values;

  // allocate memory for the function values
  y_values = (double *)malloc(f_data->n_data * sizeof(double));

  // go over all independent data
  for (index=0; index < f_data->n_data; index++)
      // compute and store a function value
      y_values[index] = comp_shift(f_data->x_values[index], gpars, fpars);

  // return the function values
  return y_values;
}

/**
 * Function: write_ipc_fit
 *
 * Parameters:
 * @param f_data   - the input data
 * @param y_values - the model values
 * @param ipc_file - name of the file to write to
 *
 * Returns:
 * @return -
 */
void
write_ipc_fit(fit_data *f_data, double *y_values, const char *ipc_file, const double qual, double *yfract)
{
  FILE *fout;
  char Buffer[MAXCHAR];

  int index;

  // open the file for the values
  fout = fopen (ipc_file, "w");

  // write norm (= quality of fit) to comment line
  sprintf (Buffer, "# norm of the residue vector: %e\n", qual);
  fputs (Buffer, fout);
  // write norm (= quality of fit) to comment line
  sprintf (Buffer, "#%e\n", qual);
  fputs (Buffer, fout);
  sprintf (Buffer, "# Wanderer kommst Du nach Sparta: %e\n", qual);
  fputs (Buffer, fout);
  sprintf (Buffer, "# berichtige dorten: %e\n", qual);
  fputs (Buffer, fout);
  sprintf (Buffer, "# Du habest uns hier liegen gesehen: %e\n", qual);
  fputs (Buffer, fout);
  sprintf (Buffer, "# wie es die Geschichte befahl. %e\n", qual);
  fputs (Buffer, fout);
  sprintf (Buffer, "#%e\n", qual);
  fputs (Buffer, fout);

  // go over the data
  for (index=0; index < f_data->n_data; index++)
    {
      // write the data point to the buffer
      sprintf (Buffer, "%e %e %e %e %e\n", f_data->x_values[index], f_data->y_values[index], y_values[index],  f_data->e_values[index], yfract[index]);
      // write the buffer to the file
      fputs (Buffer, fout);
    }

  // close the file
  fclose(fout);
}

/**
 * Function: fit_ipc_data
 *
 * Parameters:
 * @param f_data  - the input data
 * @param gpars   - the guess or fit parameters
 * @param n_gpars - number of guess or fit parameters
 * @param fpars   - the fixed parameters
 * @param n_fpars - number of fixed parameters
 *
 * Returns:
 * @return control.fnorm - goodness of the fit value
 */
double
fit_ipc_data(fit_data *f_data, double *gpars, int n_gpars,
	     double *fpars, int n_fpars)
{
  int m_dat = f_data->n_data;

  // auxiliary settings:
  lm_control_type control;
  lm_data_fpar_type data;


  // initialize the fit
  lm_initialize_control( &control );

  // define the parameters for the fitting
  data.user_func = comp_shift;
  data.user_t = f_data->x_values;
  data.user_y = f_data->y_values;
  data.fpars  = fpars;

  // print the initial
  // value
  fprintf(stdout, "input: %e, ", gpars[0]);

  // perform the fit:
  lm_minimize( m_dat, n_gpars, gpars, lm_evaluate_fpar, lm_print_nothing,
  	       &data, &control );

  // print the new value plus additional information
  fprintf(stdout, " output: %e, info: %i, niter: %i norm: %e\n",
	  gpars[0], control.info, control.nfev, control.fnorm);

  // return the variable
  // describing the goodness
  // of the fit
  return control.fnorm;
}

/**
 * Function: get_ipc_fdata
 *
 * Parameters:
 * @param conf          - the aXe configuratuion structure
 * @param obs           - the grism image
 * @param data_matrix   - background subtracted grism image
 * @param actbeam       - the beam to examine
 *
 * Returns:
 * @return f_final - the filled data structure
 */
fit_data *
get_ipc_fdata(const aperture_conf  *conf,  observation *obs,
	      gsl_matrix *data_vals, beam act_beam)
{
  trace_func     *tracefun;

  px_point        xborder;
  px_point        tpoint;

  gsl_vector      *ipc_data;

  double rel12;

  int i;
  int index;

  fit_data *fdata  = NULL;
  fit_data *f_final= NULL;

  // define the beam and the trace function
  tracefun = act_beam.spec_trace;

  // If this beam's ignore flag is set to 1 then do nothing
  if (act_beam.ignore == 1)
    return fdata;

  // determine the start and end point in x
  xborder = get_fit_xrange(conf, obs, act_beam);

  // an assertion to prevent
  // catastrophic failures
  if  (xborder.x > xborder.y)
    return fdata;
  // introduced after making release 1 reduction


  // allocate memory for the fitting structure
  fdata = alloc_fit_data(xborder.y-xborder.x);

  // Loop over all columns
  index = 0;
  for (i = xborder.x; i < xborder.y; i++)
    {
      // determine the pixel which is
      // closest to the nominal trace position
      tpoint.x = i;
      tpoint.y = (int)floor(tracefun->func((double)i-act_beam.refpoint.x,
					   tracefun->data)
			    + act_beam.refpoint.y+0.5);


      // collect pixel data
      ipc_data = get_ipc_coldata(data_vals, N_COLUMNS_FIT, tpoint);

      // go to the next
      // row if no data is available
      if (!ipc_data)
	continue;

      // sort the values
      gsl_sort_vector(ipc_data);

      // at least two values are needed;
      // prevent division by zero
      if (ipc_data->size > 1 && gsl_vector_get(ipc_data, ipc_data->size-1))
	{
	  // compute the relation
	  rel12 = gsl_vector_get(ipc_data, ipc_data->size-2) / gsl_vector_get(ipc_data, ipc_data->size-1);

	  fdata->x_values[index] = (double)i - act_beam.refpoint.x;
	  fdata->y_values[index] = rel12;
	  fdata->e_values[index] = 1.0;
	  index++;
	}

      // free the allocated memory
      gsl_vector_free(ipc_data);
    }


  // check whether the ipc data covers the minimum length
  if (fabs(fdata->x_values[fdata->n_data-1] - fdata->x_values[0]) > MIN_IPC_LENGTH)
    // reduce the data to what is necessary
    f_final = strip_fitdata(fdata, index);
  else
    fprintf(stdout, "beam tooo short\n");

  // release the original structure
  free_fit_data(fdata);

  // return the fit data
  return f_final;
}


/**
 * Function: get_wf_special
 *
 * Parameters:
 * @param conf          - the aXe configuratuion structure
 * @param obs           - the grism image
 * @param data_matrix   - background subtracted grism image
 * @param actbeam       - the beam to examine
 *
 * Returns:
 * @return f_final - the filled data structure
 */
double *
get_wf_special(const aperture_conf  *conf,  observation *obs,
	       gsl_matrix *data_vals, beam act_beam, double yoffs)
{
  trace_func     *tracefun;

  px_point        xborder;
  px_point        tpoint;

  gsl_vector      *ipc_data;

  double rel12;

  int i;
  int index;

  fit_data *fdata  = NULL;
  //fit_data *f_final= NULL;

  double yyy;

  double *yfract=NULL;
  double *yfinal = 0;


  // define the beam and the trace function
  tracefun = act_beam.spec_trace;

  // If this beam's ignore flag is set to 1 then do nothing
  if (act_beam.ignore == 1)
    return yfinal;

  // determine the start and end point in x
  xborder = get_fit_xrange(conf, obs, act_beam);

  // allocate memory for the fitting structure
  fdata = alloc_fit_data(xborder.y-xborder.x);

  // allocate memory
  yfract = (double *)malloc((xborder.y-xborder.x)*sizeof(double));

  // Loop over all columns
  index = 0;
  for (i = xborder.x; i < xborder.y; i++)
    {
      // determine the pixel which is
      // closest to the nominal trace position
      tpoint.x = i;
      tpoint.y = (int)floor(tracefun->func((double)i-act_beam.refpoint.x,
					   tracefun->data)
			    + act_beam.refpoint.y+0.5);

      yyy           = tracefun->func((double)i-act_beam.refpoint.x, tracefun->data) + act_beam.refpoint.y + yoffs;
      yfract[index] = yyy - floor(yyy);

      // collect pixel data
      ipc_data = get_ipc_coldata(data_vals, N_COLUMNS_FIT, tpoint);

      // go to the next
      // row if no data is available
      if (!ipc_data)
	continue;

      // sort the values
      gsl_sort_vector(ipc_data);

      // at least two values are needed;
      // prevent division by zero
      if (ipc_data->size > 1 && gsl_vector_get(ipc_data, ipc_data->size-1))
	{
	  // compute the relation
	  rel12 = gsl_vector_get(ipc_data, ipc_data->size-2) / gsl_vector_get(ipc_data, ipc_data->size-1);

	  fdata->x_values[index] = (double)i - act_beam.refpoint.x;
	  fdata->y_values[index] = rel12;
	  fdata->e_values[index] = 1.0;
	  index++;
	}

      // free the allocated memory
      gsl_vector_free(ipc_data);
    }


  yfinal = strip_wf_special(fdata, yfract, index);

  // release the original structure
  free_fit_data(fdata);
  if (yfract)
    free(yfract);

  // return the fit data
  return yfinal;
}



/**
 * Function: strip_wf_special
 *
 * Parameters:
 * @param old_fdata - the experimental data
 * @param i_max     - the maximum index to check for valid data
 *
 * Returns:
 * @return fdata - the new, completely filled data structure
 */
double *
strip_wf_special(fit_data *old_fdata, double *y_fract, int i_max)
{
  double *y_new=NULL;

  int ii;
  int i_count;

  // check that SOME values
  // are accepted
  if (i_max < 1)
    // return NULL if not
    return y_new;

  // count the number of
  // values with weight
  i_count=0;
  for (ii=0; ii < i_max; ii++)
    if (old_fdata->e_values[ii])
      i_count++;

  // check whether there
  // is data at all
  if (i_count < 1)
    // return NULL if not
    return y_new;

  // allocate a new structure
  y_new = (double *)malloc(i_count * sizeof(double));

  // set the counter
  i_count=0;

  // go over the old data
  // to the maximum index
  for (ii=0; ii < i_max; ii++)
    {
      // check whether the
      // data is valid
      if (old_fdata->e_values[ii])
	{
	  // transfer valid data
	  y_new[i_count] = y_fract[ii];

	  // enhance the counter
	  i_count++;
	}
    }

  // return the new structure
  return y_new;
}


/**
 * Function: strip_fitdata
 *
 * Parameters:
 * @param old_fdata - the experimental data
 * @param i_max     - the maximum index to check for valid data
 *
 * Returns:
 * @return fdata - the new, completely filled data structure
 */
fit_data *
strip_fitdata(fit_data *old_fdata, int i_max)
{
  fit_data *fdata = NULL;

  int ii;
  int i_count;

  // check that SOME values
  // are accepted
  if (i_max < 1)
    // return NULL if not
    return fdata;

  // count the number of
  // values with weight
  i_count=0;
  for (ii=0; ii < i_max; ii++)
    if (old_fdata->e_values[ii])
      i_count++;

  // check whether there
  // is data at all
  if (i_count < 1)
    // return NULL if not
    return fdata;

  // allocate a new structure
  fdata = alloc_fit_data(i_count);

  // set the counter
  i_count=0;

  // go over the old data
  // to the maximum index
  for (ii=0; ii < i_max; ii++)
    {
      // check whether the
      // data is valid
      if (old_fdata->e_values[ii])
	{
	  // transfer valid data
	  fdata->x_values[i_count] = old_fdata->x_values[ii];
	  fdata->y_values[i_count] = old_fdata->y_values[ii];
	  fdata->e_values[i_count] = old_fdata->e_values[ii];

	  // enhance the counter
	  i_count++;
	}
    }

  // return the new structure
  return fdata;
}

/**
 * Function: kappa_sigma_klipp_ipc
 *
 * Parameters:
 * @param conf          - the aXe configuratuion structure
 * @param obs           - the grism image
 * @param data_matrix   - background subtracted grism image
 * @param actbeam       - the beam to examine
 * @param ipc_file_path - file name for the phase data
 *
 * Returns:
 * @return ret - the change in y-shift and the goodness of the fit
 */
d_point
kappa_sigma_klipp_ipc(const aperture_conf  *conf, observation *obs,
		      gsl_matrix *data_matrix, beam act_beam,
		      const char *ipc_file)
{
  fit_data *f_data   = NULL;
  fit_data *f_work   = NULL;

  double   *y_values = NULL;
  //double   *test = NULL;
  double    sigma;

  d_point   ret;

  trace_func *tracefun  = act_beam.spec_trace;
  double     *tracedata =  (double *)(tracefun->data);

  double p[1];
  double f[3] = {0.6, 0.9, tracedata[2]};

  int index=0;
  int nold_data;

  // REMOVEME ASAP !!!!!
  double *yfract=NULL;

  // compute a first guess for the shift
  p[0] = act_beam.refpoint.y - floor(act_beam.refpoint.y + tracedata[1]);

  // remember the shift
  ret.x = p[0];

  // give infinity for
  // the goodness of the fit
  ret.y = 1.0e+20;

  // get the data for the fit
  f_data = get_ipc_fdata(conf, obs, data_matrix, act_beam);

  // check whether
  // there is data
  if (f_data)
    {

      // fix the number of initial data points
      nold_data= f_data->n_data;

      // make the frist fit
      ret.y = fit_ipc_data(f_data, p, 1, f, 3);

      // evaluate the fit; comute the sigma
      y_values = evaluate_ipc_fit(f_data, p, f);
      sigma = calculate_sigma(f_data, y_values);


      // start the iteration
      while (index < N_ITER_IPC)
	{

	  // reject points according to the sigma
	  // generate a new struct to work with
	  reject_kappasigma(f_data, y_values, KAPPA_IPC*sigma);
	  // free memory in the old work structure

	  if (f_work)
	    free_fit_data(f_work);
	  f_work = strip_fitdata(f_data, f_data->n_data);

	  // check whether there
	  // is still data
	  if (!f_work)
	    break;

	  // check whether data was rejected
	  // brek if 'yes' transfer the new number if 'no'
	  if (f_work->n_data == nold_data)
	    break;
	  else
	    nold_data = f_work->n_data;

	  // make a new fit
	  ret.y = fit_ipc_data(f_work, p, 1, f, 3);

	  // evaluate fit values; comute the sigma
	  if (y_values)
	    free(y_values);
	  y_values = evaluate_ipc_fit(f_data, p, f);
	  sigma = calculate_sigma(f_data, y_values);

	  // enhance
	  // the kappa counter
	  index++;
	}

      // evaluate the fit values
      if (y_values)
	free(y_values);
      y_values = evaluate_ipc_fit(f_data, p, f);

      // REMOVEME ASAP
      yfract = get_wf_special(conf, obs, data_matrix, act_beam, p[0] - ret.x);

      // CHANGEME back ASAP
      // write the phase shift file
      write_ipc_fit(f_data, y_values, ipc_file, ret.y, yfract);

      // release memory
      if (f_data)
	free_fit_data(f_data);
      if (f_work)
	free_fit_data(f_work);
      if (y_values)
	free(y_values);
      if (yfract)
	free(yfract);
    }

  // compute the
  // difference in y-shift
  ret.x = p[0] - ret.x;

  // return the result
  return ret;
}


/**
 * Function: reject_kappasigma
 *
 * Parameters:
 * @param f_data   - the experimental data
 * @param y_values - the model data
 * @param max_diff - the maximum allowed difference
 *                   between data and model
 *
 * Returns:
 * @return -
 */
void
reject_kappasigma(fit_data *f_data, double *y_values, double max_diff)
{
  int index;

  // go over all data
  for (index=0; index < f_data->n_data; index++)
    // check whether the data point is valid and OVER the maximum allowed difference
    if (f_data->e_values[index] && fabs(f_data->y_values[index] - y_values[index]) > max_diff)
      // make the data point INVALID
      f_data->e_values[index] = 0.0;
}

/**
 * Function: calculate_sigma
 *
 * Parameters:
 * @param f_data   - the experimental data
 * @param y_values - the model data
 *
 * Returns:
 * @return sigma - the standard deviation
 */
double
calculate_sigma(fit_data *f_data, double *y_values)
{
  double *tmp;

  double sigma;
  double mean=0.0;
  double diff;

  int n_data=0;
  int index;

  // count the number of valid data points
  for (index=0; index < f_data->n_data; index++)
    if (f_data->e_values[index])
      n_data++;

  // allocate a temporary vector for
  // the valid data
  tmp = (double *)malloc(n_data * sizeof(double));

  // re-set the counter;
  // go over all data points
  n_data=0;
  for (index=0; index < f_data->n_data; index++)
    {
      // check whether the data point
      // is valid
      if (f_data->e_values[index])
	{
	  // compute the difference between data and model
	  diff =  f_data->y_values[index] - y_values[index];

	  // store the difference
	  tmp[n_data] = diff;

	  // add to the mean
	  mean += diff;

	  // enhance the counter
	  n_data++;
	}
    }

  // comopute the mean
  mean /= (double)n_data;

  // compute standard deviation
  sigma = gsl_stats_sd_m(tmp, 1, n_data, mean);

  // relese the memory
  free(tmp);

  // return the sigma
  return sigma;
}
