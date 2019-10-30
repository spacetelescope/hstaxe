/**
* Subroutines for fitting functions to data.
 * Most of the functions are used for the
 * background determination.
 */

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include "spce_fitting.h"


/**
 * Function: det_vector_polyN
 * A function using GSL fitting functions to fit a set x,y,w of n elements
 * such that an Nth degree polynomial is fitted to the x and y vectors 
 * using the weights w. This function avoids NaN values in ys AND elements
 * with an associated weight that is zero.
 *
 * Parameters:
 * @param m   - order of the fit 
 * @param xs  - double vector containing the x values
 * @param ys  - double vector containing the y values
 * @param ws  - double vector containing the weights associated with ys
 * @param n   - number of points in xs,ys, and ws (must be greater than m!)
 * @param c   - the vector with the fitted coefficients
 * @param cov - the covariance matrix
 *
 * Returns:
 * @return interp - vector with additional fitting information
 */
gsl_vector *
det_vector_polyN (int m, const double *const xs, double *const ys,
		   double *const ws, const int n, gsl_vector *c,
		   gsl_matrix *cov)
{
  int i, j, nn, ii;
  double chisq;
  gsl_matrix *X;
  gsl_vector *y, *w;
  
  double *tmp, *xmp;
  double median, xean;

  gsl_vector *interp;

  // allocate the return vector
  interp = gsl_vector_alloc(5);

  // allocate temporary vectors
  tmp = (double *) malloc (n * sizeof (double));
  xmp = (double *) malloc (n * sizeof (double));


  // fill weights and independent values
  // for the background pixels into the arrays
  nn = 0;
  for (i = 0; i < n; i++)
    {
      if (ws[i] > 0.0)
	{
	  //	  tmp[nn] = ws[i];
	  tmp[nn] = 1.0/(ws[i]*ws[i]);
	  xmp[nn] = xs[i];
	  nn++;
	}
    }

  // compute the median weight
  // to be used for all pixels
  gsl_sort( tmp, 1, m);
  median = gsl_stats_median_from_sorted_data(tmp, 1, nn);

  // determine the anchor point in the
  // independent variable
  xean = gsl_stats_mean (xmp, 1, nn);

  // check whether the desired
  // intepolation degree is feasible
  if (nn < m)
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		   "Not enough points (%d)  found to perform the %d order fit. Doing %d order instead.\n",
		   nn, m, nn);
      m = nn;
    }

  // allocate more temporary space
  X = gsl_matrix_alloc (nn, m); 
  y = gsl_vector_alloc (nn);
  w = gsl_vector_alloc (nn);
  //  c = gsl_vector_alloc (m);
  //  cov = gsl_matrix_alloc (m, m);
  
  // transfer independent/dependent and weight
  // values from the background pixels to the
  // temporary arrays
  ii = 0;
  for (i = 0; i < n; i++)
    {
      if (ws[i] > 0.0)
	{
	  // shift the x-values
	  // and store them in the matrix
	  for (j = 0; j < m; j++)
	    {
	      gsl_matrix_set (X, ii, j, pow (xs[i] - xean, j));
	    }
	  gsl_vector_set (y, ii, ys[i]);
	  gsl_vector_set (w, ii, median);
	  ii++;
	}
      
    }

  // allocate space and do the fit; release the space
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (nn, m);
  gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free (work);


  // release memory
  gsl_matrix_free (X); 
  gsl_vector_free (y);
  gsl_vector_free (w);

  // free the tmp-arrays
  free(tmp);
  free(xmp);

  // fill the mean x-value
  // and the order into the return vector
  gsl_vector_set(interp, 0, xean);
  gsl_vector_set(interp, 1, (double)m);
  

  return interp;
}

/**
 * Function: fill_polyN_interp
 * Evaluates a polynomial of any order and fills the
 * result plus the error, compute using the covariance
 * matrix, in vectors.
 *
 * Parameters:
 * @param xs     - double vector containing the x values
 * @param ys     - double vector containing the y values
 * @param yi     - double vector containing the intermediate y-values
 * @param ws     - double vector containing the weights associated with ys
 * @param n      - number of points in xs,ys, and ws (must be greater than m!)
 * @param c      - the vector with the fitted coefficients
 * @param cov    - the covariance matrix
 * @param interp - additional fitting ifnormation
 * @param final  - flagg to indicate final fit and
 *                 therefore a replacement
 */
void
fill_polyN_interp(const double *const xs, double *const ys,
		  double *const ws, double *yi, const int n,
		  gsl_vector *coeffs, gsl_matrix *cov,
		  gsl_vector *interp, const int final)
{
  int i, j, m;
  double xf, yf, sq_yf_err, yf_err;
  double xean;

  // get the mean x-value and
  // polynomial order from the vector
  xean = gsl_vector_get(interp, 0);
  m    = (int)gsl_vector_get(interp, 1);

  // put the interpolated values
  // for object pixels back into the
  // inut arrays
  for (i = 0; i < n; i++)
    {

      // compute and store interpolated
      // values an corresponding error
      xf = xs[i]-xean;
      yf = 0.0;
      yf_err = 0.0;
      sq_yf_err = 0.0;
      
      // compute the polynimials
      for (j = 0; j < m; j++)
	{
	  // first compute the y-value
	  yf += gsl_vector_get(coeffs,j)*pow (xf, j);

	  // compute the associated error
	  sq_yf_err += gsl_matrix_get(cov,j,j) * pow (xf, j) * pow (xf, j);
	}

      // put the intepolated value in the intermediate vector
      yi[i] = yf;

      // if requested, fill the interpolated values
      // into the origina; arrays
      if (final && ws[i] == 0.0)
	{
	  ys[i] = yf;
	  ws[i] = sqrt(sq_yf_err);
	}
    }
}




/**
 * Function: comp_vector_average
 * The function provides the computation of the average of
 * the y-values and store that average in temporary
 * or the original y-vector.
 *
 * @param xs    - absissa
 * @param ys    - value at xs[]
 * @param ws    - weight in ys[]
 * @param yi    - temprary y-values
 * @param n     - number of points in xs, ys, ws
 * @param final - flagg to store to temp (0) or final (1) vectors
 */
void
comp_vector_average (const double *xs, double *ys,
		     double *ws, double *yi, const int n, const int final)
{
  double ave, std;
  //int i, m;

  det_vector_average (xs, ys, ws, n, &ave, &std);
  fill_const_value(ys, ws, yi, n, ave, std, final);
}


/**
 * Function: comp_vector_median
 * The function provides the computation of the median of
 * the y-values and store that median in temporary
 * or the original y-vector
 *
 * @param xs - absissa
 * @param ys - value at xs[]
 * @param ws - weight in ys[]
 * @param yi    - temprary y-values
 * @param n     - number of points in xs, ys, ws
 * @param final - flagg to store to temp (0) or final (1) vectors
 */
void
comp_vector_median (const double *xs, double *ys,
		    double *ws, double *yi, const int n, const int final)
{
  double med, std;

  det_vector_median (xs, ys, ws, n, &med, &std);
  fill_const_value(ys, ws, yi, n, med, std, final);

}

/**
 * Function: comp_vector_linear
 * The function fits a linear function to the y-values
 * and fills the linear inteprolated values in a temporary
 * vector and, if requested, in the original vector as well.
 *
 * @param xs    - absissa
 * @param ys    - value at xs[]
 * @param ws    - weight in ys[]
 * @param yi    - temprary y-values
 * @param n     - number of points in xs, ys, ws 
 * @param final - flagg to store to temp (0) or final (1) vectors
 */
void
comp_vector_linear (const double *xs, double *ys,
		    double *ws, double *yi, const int n, const int final)
{
  gsl_vector *interp;

  // make a weighted linear fit
  interp = det_vector_linear(xs, ys, ws, n, 1);

  // make the interpolations,
  // using the linear fit
  fill_linear_interp(xs, ys, ws, yi, n, interp, final);

  gsl_vector_free(interp);
}

/**
 * Function: comp_vector_polyN
 * The function fits a polynomial function to the y-values
 * and fills the evaluated polynomial values in a temporary
 * vector and, if requested, in the original vector as well.
 *
 * @param xs    - absissa
 * @param ys    - value at xs[]
 * @param ws    - weight in ys[]
 * @param yi    - temprary y-values
 * @param n     - number of points in xs, ys, ws
 * @param final - flagg to store to temp (0) or final (1) vectors
 */
void
comp_vector_polyN (const int m, const double *xs, double *ys,
		    double *ws, double *yi, const int n, const int final)
{
  gsl_vector *interp, *coeffs;
  gsl_matrix *cov;

  coeffs = gsl_vector_alloc(m);
  cov = gsl_matrix_alloc(m,m);

  interp = det_vector_polyN (m, xs, ys, ws, n, coeffs, cov);
  fill_polyN_interp(xs, ys, ws, yi, n, coeffs, cov, interp, final);

  gsl_vector_free(interp);
  gsl_vector_free(coeffs);
  gsl_matrix_free(cov);
}

/**
 * Function: det_vector_linear
 * Performs a linear fit of the xs, ys, and ws arrays.
 * Ignores NaN values. Fitted values are returned in
 * ys and errors in ws.
 *
 * Parameters:
 * @param xs     - absissa
 * @param ys     - value at xs[]
 * @param ws     - error in ys[]
 * @param n      - number of points in xs, ys, ws
 * @param weight - make weighted fit (=1) or not (=0)
 *
 * Returns:
 * @return ret - the covariance values
 */
gsl_vector *
det_vector_linear(const double *xs, double *ys, double *ws,
		  const int n, const int weight)
{
  int i, m;
  double c0, c1, cov00, cov01, cov11, chisq;
  double *xx, *yy, *ww, *tmp, *xmp;
  double median, xean;
  //double xf, yf, yf_err;

  gsl_vector *ret;

  ret = gsl_vector_alloc(10);

  // allocate space for temporary vectors 
  tmp = (double *) malloc (n * sizeof (double));
  xmp = (double *) malloc (n * sizeof (double));

  // fill the temporary vectors 
  // with background values and weights
  m = 0;
  for (i = 0; i < n; i++)
    {
      if (ws[i] > 0.0)
	{
	  //	  tmp[m] = ws[i];
	  tmp[m] = 1.0/(ws[i]*ws[i]);
	  xmp[m] = xs[i];
	  m++;
	}
    }


  // determine the median weight which
  // is used for all data
  gsl_sort( tmp, 1, m);
  median = gsl_stats_median_from_sorted_data(tmp, 1, m);

  // determine the mean in the independent variable
  xean = gsl_stats_mean (xmp, 1, m);

  // allocate more temporary vectors
  xx = (double *) malloc (m * sizeof (double));
  yy = (double *) malloc (m * sizeof (double));
  ww = (double *) malloc (m * sizeof (double));

  // fill the temporary vectors with
  // background independent/dependent variables
  // plus the weight
  m = 0;
  for (i = 0; i < n; i++)
    {
      if (ws[i] > 0.0)
	{
	 xx[m] = xs[i]-xean;
	 yy[m] = ys[i];
	 ww[m] = median;
	 m++;
	}
    }

  // check for weighted fit
  if (weight)
    // doe the linear fit
    gsl_fit_wlinear (xx, 1, ww, 1, yy, 1, m, &c0, &c1, &cov00,
		     &cov01, &cov11, &chisq);
  else
    // does the linear fit
    gsl_fit_linear (xx, 1, yy, 1, m, &c0, &c1, &cov00,
		    &cov01, &cov11, &chisq);


  gsl_vector_set(ret, 0, xean);
  gsl_vector_set(ret, 1, c0);
  gsl_vector_set(ret, 2, c1);
  gsl_vector_set(ret, 3, cov00);
  gsl_vector_set(ret, 4, cov01);
  gsl_vector_set(ret, 5, cov11);
  gsl_vector_set(ret, 6, chisq);


  // free the arrays
  free(tmp);
  free(xmp);
  free(xx);
  free(yy);
  free(ww);

  return ret;
}

/**
 * Function: fill_linear_interp
 * The function computes and stores linear interpolated values
 * in one or several output vectors.
 *
 * Parameters:
 * @param xs     - double vector containing the x values
 * @param ys     - double vector containing the y values
 * @param ws     - double vector containing the weights associated with ys
 * @param yi     - double vector containing the intermediate y-values
 * @param n      - number of points in xs,ys, and ws (must be greater than m!)
 * @param interp - vector with the fitted values
 * @param final  - flagg to indicate final fit and
 *                 therefore a replacement
 */
void
fill_linear_interp(const double *const xs, double *const ys, double *const ws,
		   double *yi, const int n, gsl_vector *interp,
		   const int final)
{
  double c0, c1, cov00, cov01, cov11, xean;
  double xf, yf, yf_err;
  int i;

  xean  = gsl_vector_get(interp, 0);
  c0    = gsl_vector_get(interp, 1);
  c1    = gsl_vector_get(interp, 2);
  cov00 = gsl_vector_get(interp, 3);
  cov01 = gsl_vector_get(interp, 4);
  cov11 = gsl_vector_get(interp, 5);

  // determine the interpolated values
  // and store them in the input arrays
  for (i = 0; i < n; i++)
    {

      xf = xs[i]-xean;
      gsl_fit_linear_est (xf, c0, c1, cov00, cov01, cov11, &yf, &yf_err);
      yi[i] = yf;

      // for object pixels,
      // fill in the interpolated values
      if (final && ws[i] == 0.0)
	{
	  ys[i] = yf;
	  ws[i] = yf_err;
	}
    }
}


/**
 * Function: det_vector_average
 * The function computes the average and the standard
 * deviation from the values in a vector
 *
 * @param xs  - absissa
 * @param ys  - value at xs[]
 * @param ws  - weight in ys[]
 * @param n   - number of points in xs, ys, ws
 * @param avg - the average
 * @param std - the standard deviation
 */
void
det_vector_average (const double *xs, double *ys,
		    double *ws, const int n, double *avg, double *std)
{
  double *tmp, sum = 0.0;
  int nn = 0;
  int j;

  // count the number of
  // interpolation points
  nn = 0;
  for (j = 0; j < n; j++)
    {
      if (ws[j] != 0.0)
	nn++;
    }

  // allocate a tmp-array
  tmp = malloc (nn * sizeof (double));

  // sum up all values in the
  // interpolation points,
  // and transfer the values to the 
  // tmp-array
  nn = 0;
  for (j = 0; j < n; j++)
    {
      if (ws[j] != 0.0)
	{
	  sum += ys[j];
	  tmp[nn] = ys[j];
	  nn++;
	}
    }

  // compute the average
  if (nn > 0)
    *avg = sum / (double)nn;

  // compute the standard deviation
  *std = gsl_stats_sd (tmp, 1, nn);
}

/**
 * Function: det_vector_median
 * The function computes the median and the standard
 * deviation from the values in a vector
 *
 * @param xs  - absissa
 * @param ys  - value at xs[]
 * @param ws  - weight in ys[]
 * @param n   - number of points in xs, ys, ws
 * @param med - the median
 * @param std - the standard deviation
 */
void
det_vector_median (const double *xs, double *ys,
		   double *ws, const int n, double *med, double *std)
{
  double *tmp;
  int nn = 0;
  int j;

  // count the number of
  // interpolation points
  for (j = 0; j < n; j++)
    {
      if (ws[j] != 0.0)
	  nn++;
    }

  // allocate a tmp-array
  tmp = malloc (nn * sizeof (double));

  // transfer the values of the interpolation
  // points to the tmp-array
  nn = 0;
  for (j = 0; j < n; j++)
    {
      if (ws[j] != 0.0)
	{
	  tmp[nn] = ys[j];
	  nn++;
	}
    }

  // sort the tmp array
  gsl_sort (tmp, 1, nn);

  // compute the median
  *med = gsl_stats_median_from_sorted_data (tmp, 1, nn);

  // compute the standard deviation
  *std = gsl_stats_sd (tmp, 1, nn);
}

/**
 * Function:fill_const_value
 * The function fills two constant values into
 * up to three vectors.
 *
 * @param ys    - value at xs[]
 * @param ws    - weight in ys[]
 * @param yi    - value at xs[]
 * @param n     - number of points in xs, ys, ws
 * @param cval  - number of points in xs, ys, ws
 * @param stdev - number of points in xs, ys, ws
 * @param final - flagg which indicates which vectors to fill
 */
void
fill_const_value(double *ys, double *ws, double *yi, const int n,
		 double cval, double stdev, const int final)
{

  int j;
  int m=0;

  // fill the constant value and the stdev
  // in the whole value and error
  // arrays, respectively
  for (j = 0; j < n; j++)
    {
      yi[j] = cval;

      if (final && ws[j] == 0.0)
	{
	  ys[j] = cval;
	  ws[j] = stdev;
	  m++;
	}
    }
}
