/**
 * Set of functions to handle object slitless background subtraction
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include "fitsio.h"
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_sect.h"
#include "spce_fitting.h"
#include "spce_is_in.h"
#include "spc_back.h"

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))



/**
 * Function: is_pt_in_a_beam
 * The function checks whether a particular image pixel is part
 * of a beam or not. The pixel is checked against an array
 * of beams it could be part of.
 *
 * Parameters:
 * @param apoint  - the point to check
 * @param iids    - an array of is_in_descriptors for all the beams that
 *                  should be checked
 * @param tnbeams - number of beams in iids
 *
 * Returns:
 * @return 1/0   - fixed values
 */
int
is_pt_in_a_beam (const px_point * const apoint,
                 const is_in_descriptor * const iids, const int tnbeams)
{
  int i;

  for (i = 0; i < tnbeams; i++)
    {
      if (apoint->x < (iids + i)->mini)
        continue;
      if (apoint->x > (iids + i)->maxi)
        continue;
      if (apoint->y < (iids + i)->minj)
        continue;
      if (apoint->y > (iids + i)->maxj)
        continue;
      if (is_in (apoint->x, apoint->y, iids + i))
        {
          return 1;
        }
    }
  return 0;
}

/**
 * Function: get_window_points
 * The subroutine searches in a window around the tracepoint
 * for pixels which are suited for the background determination.
 * The extend of the window is specified in the beam structured,
 * and pixel which are part of any beam can not be used.
 *
 * Parameters:
 * @param  obs      - the observation
 * @param  bck_mask - the background mask
 * @param  actbeam  - the beam to search interp points for
 * @param  tr_point - the integer trace point
 *
 * Returns:
 * @return ret      - the vector with the row numbers
 */
/*
gsl_vector_int *
get_window_points(observation *obs, gsl_matrix *bck_mask, beam actbeam,
                  px_point tr_point)
{
  gsl_vector_int *tmp;
  gsl_vector_int *ret;

  int np_act=0;
  int l_act;
  int u_act;
  int ii;
  int ncols;
  int onbeam=0;

  double umax, lmax;

  // extract the lower and upper window extend
  umax = actbeam.backwindow.x;
  lmax = actbeam.backwindow.y;

  // get the number of columns
  ncols=obs->grism->size2;

  // allocate the vector
  tmp = gsl_vector_int_alloc((int)ceil(lmax) + (int)ceil(umax) + 1);

  // limit the starting point of the search
  // to values within the image dimension
  tr_point.y = MAX(0,tr_point.y);
  tr_point.y = MIN((int)obs->grism->size2,tr_point.y);

  // initialize the row numbers
  // to search up- and downwards
  l_act = tr_point.y -1;
  u_act = tr_point.y;

  if (!gsl_matrix_get(bck_mask,tr_point.x,tr_point.y))
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,"The closest pixel to the",
                 "trace at %i, %i must be part of an object!",tr_point.x,
                 tr_point.y );


  // search downwards from the trace,
  // assume to start on the beeam
  onbeam=1;
  while (l_act > -1 && fabs((double)l_act - (double)tr_point.y) < lmax)
    {
      // if the pixel is not part of any beam and is not NAN
      if (gsl_matrix_get(bck_mask,tr_point.x,l_act) == 0 &&
          !isnan(gsl_matrix_get(obs->grism, tr_point.x, l_act)))
        {
          // store the pixel index
          gsl_vector_int_set(tmp, np_act, l_act);

          if (onbeam)
            onbeam = 0;

          // enhance the counter
          np_act++;
        }
      // if you are still on the beam
      else if (onbeam && gsl_matrix_get(bck_mask,tr_point.x,l_act))
        {
          // store the pixel index
          gsl_vector_int_set(tmp, np_act, l_act);

          // enhance the counter
          np_act++;
        }
      // decrease the search index
      l_act--;
    }

  // search upwards from the trace,
  // start always searching on the beam
  onbeam=1;
  while (u_act < ncols && fabs((double)u_act - (double)tr_point.y) < umax)
    {
      // if the pixel is not part of any beam and is not NAN
      if (gsl_matrix_get(bck_mask,tr_point.x,u_act) == 0 &&
          !isnan(gsl_matrix_get(obs->grism, tr_point.x, u_act)))
        {
          // store the pixel index
          gsl_vector_int_set(tmp, np_act, u_act);

          if (onbeam)
            onbeam = 0;

          // enhance the counter
          np_act++;
        }
      // if you are still on the beam
      else if (onbeam && gsl_matrix_get(bck_mask,tr_point.x,u_act))
        {
          // store the pixel index
          gsl_vector_int_set(tmp, np_act, u_act);

          // enhance the counter
          np_act++;
        }

      // enhance the search index
      u_act++;
    }

  // check whether the beam extends over
  // the image. Add an artificial start or
  // end point if necessary
  if (l_act < 0)
    gsl_vector_int_set(tmp, np_act++, -1);
  if (u_act > ncols-1)
    gsl_vector_int_set(tmp, np_act++, ncols);

  if (np_act >0)
    {
      // transfer the row numbers to a
      // vector of the right size
      ret = gsl_vector_int_alloc(np_act);
      for (ii=0; ii < np_act; ii++)
        gsl_vector_int_set(ret, ii, gsl_vector_int_get(tmp, ii));

      // sort the row numbers
      gsl_sort_vector_int(ret);

    }
  else
    {
      ret = NULL;
    }

  gsl_vector_int_free(tmp);

  //  return the result
  return ret;
}
*/

/**
 * Function: get_interp_points
 * The function searches for 2*n interpolation points above and
 * below the trace. The search is iteratively from the trace
 * in both directions. In case that the fram borders are met,
 * the respective direction is "closed". The subroutine
 * then tries to get more background pixels on the other
 * size to reach the desired number.
 *
 * Parameters:
 * @param  obs      - the observation
 * @param  bck_mask - the background mask
 * @param  np       - the desired number of interp. points
 * @param  tr_point - the integer trace point
 *
 * Returns:
 * @return ret      - the vector with the row numbers
 */
gsl_vector_int *
get_interp_points(observation *obs, gsl_matrix *bck_mask,
                  int np, px_point tr_point)
{
  gsl_vector_int *tmp;
  gsl_vector_int *ret;
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

  // allocate the vector
  tmp = gsl_vector_int_alloc(2*np+2);

  // limit the starting ppoint of the search
  // to values within the image dimension
  tr_point.y = MAX(0,tr_point.y);
  tr_point.y = MIN((int)obs->grism->size2,tr_point.y);

  // initialize the row numbers
  // to search up- and downwards
  l_act = tr_point.y -1;
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
          // move downward until you either find
          // and interp. point or the end of the frame
          while (l_act > -1 && (gsl_matrix_get(bck_mask,tr_point.x,l_act)
                                !=0 ||
                                isnan(gsl_matrix_get(obs->grism,
                                                     tr_point.x, l_act))))
            {
              l_act--;
            }

          // if you are at the end of the frame
          if (l_act < 0)
            {
              // close the direction downwards
              l_space=0;
            }
          else
            {
              // else store the interpolation point,
              // do the various increments
              gsl_vector_int_set(tmp, np_act, l_act);
              np_act++;
              l_np++;
              l_act--;
            }
        }

      // check whether the direction
      // upwards is still open
      if (u_space)
        {
          // move upward until you either find
          // and interp. point or the end of the frame
          while (u_act < ncols && (gsl_matrix_get(bck_mask,tr_point.x,u_act)
                                   !=0 ||
                                   isnan(gsl_matrix_get(obs->grism,
                                                        tr_point.x, u_act))))
            {
              u_act++;
            }

          // if you are at the end of the frame
          if (u_act >= ncols)
            {
              // close the direction upwards
              u_space=0;
            }
          else
            {
              // else store the interpolation point,
              // do the various increments
              gsl_vector_int_set(tmp, np_act, u_act);
              np_act++;
              u_np++;
              u_act++;
            }
        }
    }

  // check whether the beam extends over
  // the image. Add an artificial start or
  // end point if necessary
  if (!l_np)
    gsl_vector_int_set(tmp, np_act++, -1);
  if (!u_np)
    gsl_vector_int_set(tmp, np_act++, ncols);

  // transfer the row numbers to a
  // vector of the right size
  ret = gsl_vector_int_alloc(np_act);
  for (ii=0; ii < np_act; ii++)
    gsl_vector_int_set(ret, ii, gsl_vector_int_get(tmp, ii));

  // sort the row numbers
  gsl_sort_vector_int(ret);

  // release memory
  gsl_vector_int_free(tmp);

  //  return the result
  return ret;
}

/**
 * Function: compute_background
 * The function extracts possible background pixels for a given
 * beam using a specified interpolation functions.
 * After possibly rejecting cosmics, the background is
 * interpolated on the areas which are masked out. The interpolated
 * values are filled into a background structure. If kappa-sigma klipping
 * is applied, the clipped pixels are flagged in the dq-array
 * of the background image.
 *
 * Parameters:
 * @param  obs         - the object list
 * @param  actbeam     - the beam to compute the background for
 * @param  bck_mask    - the background mask
 * @param  fib         - the baground structure
 * @param  npoints     - the number of interpolation points
 * @param  interporder - the interpolation order
 * @param  niter_med   - the number of iterations around the median
 * @param  niter_fit   - the number of iterations around the fit
 * @param  kappa       - the kappa value
 */
void
compute_background (observation *obs, beam actbeam, gsl_matrix *bck_mask,
                    fullimg_background *fib, int npoints,
                    int interporder, const int niter_med,
                    const int niter_fit, const double kappa)
{
  px_point xborder;
  gsl_vector_int *yvec;
  trace_func *tracefun;
  px_point tpoint;
  double *ys, *fs, *ws, *yi;
  //double var;

  int i, ii;
  int j, n;
  int min_y, max_y;


  // define the beam and the trace function
  tracefun = actbeam.spec_trace;

  /* If this beam's ignore flag is set to 1 then do nothing */
  if (actbeam.ignore == 1)
    return;

  // determine the start and end point in x
  xborder = get_xrange(obs, actbeam);


  // Loop over all columns
  for (i = xborder.x; i < xborder.y; i++)
    {

      // determine the pixel closest to the trace
      tpoint.x = i;
      tpoint.y = (int)floor(tracefun->func((double)i-actbeam.refpoint.x,
                                           tracefun->data)
                            + actbeam.refpoint.y+0.5);

      // determine the interpolation points
      // around the trace
      yvec = get_interp_points(obs, bck_mask, npoints, tpoint);

      //-----------------------------------------------------------------
      //   some code for FORS2 MXU
      //      if (actbeam.backwindow.x < 0.0)
      //        yvec = get_interp_points(obs, bck_mask, npoints, tpoint);
      //      else
      //        yvec = get_window_points(obs, bck_mask, actbeam, tpoint);


      // give a warning and go to the next
      // column if there are no background points
      if (!yvec)
        {
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "No backgound points could be found for beam %C. "
                       "at collumn %d %d",
                       BEAM(actbeam.ID), i, tpoint.y );
          continue;
        }

      // extract maximum and minimum 0f the y-values
      min_y = gsl_vector_int_get(yvec, 0);
      max_y = gsl_vector_int_get(yvec, yvec->size-1);
      //      if (actbeam.backwindow.x == 23.3 && actbeam.backwindow.y == 4.0)
      //        {
          //    fprintf(stdout, "%i <--> %i; ", min_y, max_y);
      //      if (min_y > tpoint.y - actbeam.width
      //          || max_y < tpoint.y + actbeam.width)
      //        fprintf(stdout, "%i %f <--> %f %i;  ", min_y, tpoint.y - actbeam.width,
      //                tpoint.y + actbeam.width, max_y);
              //        }

      // determine the size of the double vectors
      // to make the background determination
      // allocater the space and initialize
      // all values to default
      n = max_y - min_y + 1;
      ys = (double *) malloc (n * sizeof (double));
      fs = (double *) malloc (n * sizeof (double));
      ws = (double *) malloc (n * sizeof (double));
      yi = (double *) malloc (n * sizeof (double));
      for (ii = 0; ii < n; ii++)
        {
          ys[ii] = min_y + ii;
          fs[ii] = 0.0;
          ws[ii] = 0.0;
        }


      // transfer the values from the image column
      // to the double vectors, set the weight
      for (ii = 0; ii < (int)yvec->size; ii++)
        {
          // extract the row number
          j = gsl_vector_int_get (yvec, ii);

          // check whether the row is inside the imag
          // and whether the is no object on the pixel
          if ((j != -1) && (j != (int)obs->grism->size2)
              && !gsl_matrix_get(bck_mask, tpoint.x, j))
            {
              // set the values and the weight
              fs[j - min_y] = gsl_matrix_get (obs->grism, tpoint.x, j);
              ws[j - min_y] = gsl_matrix_get (obs->pixerrs, tpoint.x, j);
            }
        }

      // iterate on the background points to
      // reject e.g. cosmics
      if (niter_med > 0 || niter_fit > 0)
        {
          // iterate on the median
          for (j=0; j < niter_med; j++)
            kappa_sigma_clipp(ys, fs,ws, n,kappa, obs, tpoint.x);

          // iterate on the fit
          if (niter_fit > 0)
            comp_kappasigma_interp( ys, fs, ws, n, interporder,
                                    niter_fit, kappa, obs, tpoint.x);
        }

      // do the final background determination
      comp_vector_interp( ys, fs, ws, yi, n, interporder, 1);

      // copy the intepolated values
      // to the background matrix
      for (j = 0; j < n; j++)
        {
          if ((ys[j] < 0) || (ys[j] >= obs->grism->size2))
            continue;
          gsl_matrix_set (fib->bck, tpoint.x, (int) floor (ys[j]), fs[j]);
          gsl_matrix_set (fib->err, tpoint.x, (int) floor (ys[j]), ws[j]);
        }

      // release memory
      free (ys);
      ys = NULL;
      free (fs);
      fs = NULL;
      free (yi);
      yi = NULL;
      free (ws);
      ws = NULL;

      // release memory
      gsl_vector_int_free(yvec);
    }
}

/**
 * Function: get_xrange
 * The subroutine determines the extend of a beam in x-direction.
 * The minimum and maximum value in x of pixels which are part
 * of the particular beam are determined and returned.
 *
 * Parameters:
 * @param  obs     - the object list
 * @param  actbeam - the beam to determine the extent for
 *
 * Returns:
 * @return ret     - Min/Max values of the beam in x
 */
px_point
get_xrange(observation *obs, beam actbeam)
{

  px_point ret;

  // Find the object starting and ending column
  // for the beam of interest actbeam.corners
  ret.x = MIN (MIN (MIN (actbeam.corners[0].x, actbeam.corners[1].x),
                    actbeam.corners[2].x), actbeam.corners[3].x);
  ret.y  = MAX (MAX (MAX (actbeam.corners[0].x, actbeam.corners[1].x),
                     actbeam.corners[2].x), actbeam.corners[3].x);

  // limit the start and end column
  // to the image size
  ret.x = MAX(ret.x, 0);
  ret.y = MIN(ret.y, (int)obs->grism->size1);

  // warn if the beam is completely outside
  // the grism image
  if (ret.x > (int)obs->grism->size1 || ret.y < 0)
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                   "Object  is not in the image start_i:%d end_i:%d",
                   ret.x, ret.y);
    }

  // return the result
  return ret;
}

/**
 * Function: comp_kappasigma_interp
 * The function performs a kappa-sigma clipping rejection
 * on a set of data given in vectors for indipendent, dipendent
 * and weight values. The differences to base the clipping on
 * is "true value minus background value", where background
 * value is determined using the identical algorithm as the
 * final background determination.
 *
 * Parameters:
 * @param xs     - array for independent values
 * @param ys     - array for dependent values
 * @param ws     - weight array
 * @param n      - number of pixels
 * @param interp - number indicating the interpolation scheme
 * @param niter  - the number of iterations
 * @param kappa  - the kappa value for rejection
 * @param obs    - the object list
 * @param colnum - the column number
 */
void
comp_kappasigma_interp(const double *const xs, double *const ys,
                       double *const ws, const int n,
                       const int interp, const int niter, const double kappa,
                       observation *obs, int colnum)
{
  double *ys_tmp;
  double *yi_tmp;
  double *y_diff;

  double stdev;

  int *iindex;


  int i, m, j;

  // allocate temporary vectors
  ys_tmp = (double *) malloc (n * sizeof (double));
  yi_tmp = (double *) malloc (n * sizeof (double));
  y_diff = (double *) malloc (n * sizeof (double));
  iindex = (int *) malloc (n * sizeof (int));

  // transfer the dependent value
  // to the tmp vector
  for (i=0; i < n; i++)
    ys_tmp[i] = ys[i];


  // do niter times
  for (j=0; j < niter; j++)
    {

      // make the background determination
      comp_vector_interp(xs, ys_tmp, ws, yi_tmp, n, interp, 0);

      // calculate for all background
      // pixels the differences between
      // the original and background value
      m=0;
      for (i=0; i < n; i++)
        {
          if (ws[i] != 0.0)
            {
              y_diff[m] = ys[i] - yi_tmp[i];
              iindex[m] = i;
              m++;
            }
        }

      // compute the standard deviation
      // on the differences
      stdev = gsl_stats_sd (y_diff, 1, m);

      // do the clipping
      for (i=0; i < m; i++)
        {
          // check for pixels to exclude
          if (fabs(y_diff[i]) > kappa*stdev)
            {
              // set its weight to 0.0
              ws[iindex[i]] = 0.0;
              // mark the pixel in the dq-array
              if ( xs[iindex[i]] >= 0 &&  xs[iindex[i]] < obs->dq->size2)
                gsl_matrix_set(obs->dq, colnum, xs[iindex[i]], DQ_KAPPA_SIGMA);
            }
        }
    }


  for (i=0; i < n; i++)
    ys[i] = ys_tmp[i];

  // free memory
  free(ys_tmp);
  free(yi_tmp);
  free(y_diff);
  free(iindex);
}


/**
 * Function: kappa_sigma_clipp
 * The subroutine performs one kappa-sigma step on the data
 * given in various vectors for independent, dipendent and
 * weight data. The difference to apply the clipping on
 * is "original value <minus> median of the data set".
 *
 * Parameters:
 * @param xs     - array for independent values
 * @param ys     - array for dependent values
 * @param ws     - weight array
 * @param n      - number of pixels
 * @param kappa  - the kappa value for rejection
 * @param obs    - the object list
 * @param colnum - the column number
 */
void
kappa_sigma_clipp(const double *const xs, double *const ys, double *const ws,
                  const int n, const double kappa, observation *obs,
                  int colnum)
{
  double *ys_tmp, *ys_med;

  int *iindex;

  int ii, m=0, npixel=0;

  double median=0.0, stdev=0.0;

  // allocate memory for temporay vectors
  ys_tmp = (double *) malloc (n * sizeof (double));
  ys_med = (double *) malloc (n * sizeof (double));
  iindex = (int *) malloc (n * sizeof (int));

  // store the background values
  // in the temporary vectors
  for (ii=0; ii < n; ii++)
    {
      if (ws[ii] != 0.0)
        {
          ys_tmp[m] = ys[ii];
          ys_med[m] = ys[ii];
          iindex[m] = ii;
          m++;
        }
    }

  // derive median and standard deviation
  gsl_sort (ys_tmp, 1, m);
  median = gsl_stats_median_from_sorted_data(ys_tmp, 1, m);
  stdev = gsl_stats_sd_m (ys_med, 1, m, median);

  // apply the clipping
  npixel=0;
  for (ii=0; ii < m; ii++)
    {
      // check whether the pixel should
      // be clipped
      if (fabs(ys_med[ii]-median) > kappa*stdev)
        {
          // set the weight of a clipped pixel to 0.0
          ws[iindex[ii]] = 0.0;
          // mark the clipped pixel in the dq-array
          gsl_matrix_set(obs->dq, colnum, xs[iindex[ii]], DQ_KAPPA_SIGMA);
        }
      else
        {
          npixel++;
        }

    }

  // release memory
  free(ys_med);
  free(ys_tmp);
  free(iindex);
}

/**
 * Function: comp_vector_interp
 * The function passes the interpolation data to the
 * desired interpolating function.
 *
 * Parameters:
 * @param xs     - double vector containing the x values
 * @param ys     - double vector containing the y values
 * @param ws     - double vector containing the weights associated with ys
 * @param n      - number of points in xs,ys, and ws (must be greater than m!)
 * @param interp - desired interpolation type
 * @param final  - indicates a final interpolation
 */
void
comp_vector_interp(const double *const xs, double *const ys,
                   double *const ws, double *const yi, const int n,
                   const int interp, const int final)
{
  /* Median the background */
  if (interp == -1)
    {
      comp_vector_median(xs, ys, ws, yi, n, final);
    }
  /* Straight average of the background */
  else if (interp == 0)
    {
      comp_vector_average(xs, ys, ws, yi, n, final);
    }
  /* Linear interpolation of the background */
  else if (interp == 1)
    {
      comp_vector_linear(xs, ys, ws, yi, n, final);
    }
  /* n(>1) order interpolation of the background */
  else if (interp > 1)
    {
      comp_vector_polyN (interp + 1, xs, ys, ws, yi, n, final);
    }
  else
    {
      aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                  "Do not know what to do with interpolation: %i %s\n", interp);
    }
}


void
compute_global_background (object **oblist, const int obj_index,
                           gsl_matrix *bck_mask, fullimg_background * fib,
                           int interporder)
{
  int i, j, n;
  double *ys, *fs, *ws, *yi;
  //double *ws0;
  observation *grism = oblist[obj_index]->grism_obs;
  //long ma;
  //double var;

  for (i = 0; i < (int)grism->grism->size1; i++)
    {
      /* Loop over the columns of interest */

      n = grism->grism->size2;

      ys = (double *) malloc (n * sizeof (double));
      fs = (double *) malloc (n * sizeof (double));
      ws = (double *) malloc (n * sizeof (double));
      yi = (double *) malloc (n * sizeof (double));

      for (j = 0; j < n; j++)
        {
          ys[j] = j;
          fs[j] = 0.0;
          ws[j] = 0.0;
          if ((gsl_matrix_get(bck_mask,i,j)==0) &&
              (!(isnan(gsl_matrix_get(grism->grism, i,j)))))
            {
              fs[j]=gsl_matrix_get (grism->grism,i,j);
              ws[j]=gsl_matrix_get (grism->grism,i,j);
            }
        }


      /* Median the background */
      if (interporder == -1)
        {
          double *tmp, med, std;
          int nn = 0;
          for (j = 0; j < n; j++)
            {
              if (ws[j] != 0.0)
                {
                  nn++;
                }
            }
          tmp = malloc (nn * sizeof (double));
          nn = 0;
          for (j = 0; j < n; j++)
            {
              if (ws[j] != 0.0)
                {
                  tmp[nn] = fs[j];
                  nn++;
                }
            }
          gsl_sort (tmp, 1, nn);
          med = gsl_stats_median_from_sorted_data (tmp, 1, nn);
          std = gsl_stats_sd (tmp, 1, nn);
          //fprintf(stderr,"med: %g\n",med);
          for (j = 0; j < n; j++)
            {
              if (ws[j] == 0.0)
                {
                  fs[j] = med;
                  ws[j] = std;
                }
              //              else
              //                {
              //                  ws[j] = 1.0/sqrt(ws[j]);
              //                }
            }
          free (tmp);
          tmp = NULL;
        }

      /* Straight average of the background */
      if (interporder == 0)
        {
          double *tmp, sum = 0.0, avg = 0.0, std = 0.0;
          int nn = 0;

          nn = 0;
          for (j = 0; j < n; j++)
            {
              if (ws[j] != 0.0)
                {
                  nn++;
                }
            }
          tmp = malloc (nn * sizeof (double));
          nn = 0;
          for (j = 0; j < n; j++)
            {
              if (ws[j] != 0.0)
                {
                  sum += fs[j];
                  tmp[nn] = fs[j];
                  nn++;
                }
            }
          if (nn > 0)
            avg = sum / nn;
          std = gsl_stats_sd (tmp, 1, nn);
          for (j = 0; j < n; j++)
            {
              if (ws[j] == 0.0)
                {
                  fs[j] = avg;
                  ws[j] = std;
                }
              //              else
              //                {
              //                  ws[j] = 1.0/sqrt(ws[j]);
              //                }
            }
          free (tmp);
          tmp = NULL;
        }

      /* Linear interpolation of the background */
      if (interporder == 1)
        {
          comp_vector_linear (ys, fs, ws, yi, n, 1);
          //      fit_vector_linear_t (ys, fs, ws, n);
        }

      /* n(>1) order interpolation of the background */
      if (interporder > 1)
        {
          comp_vector_polyN (interporder + 1, ys,fs, ws, yi, n, 1);
      //          fit_vector_poly_N_t (interporder + 1, ys, fs, ws, n);
        }

      for (j = 0; j < n; j++)
        {
          if ((ys[j] < 0) || (ys[j] >= grism->grism->size2))
            continue;
            gsl_matrix_set (fib->bck, i, j,fs[j]);
            gsl_matrix_set (fib->err, i, j,ws[j]);
        }
      free (ys);
      ys = NULL;
      free (fs);
      fs = NULL;
      free (ws);
      ws = NULL;
    }
}



/**
 * Function: fullimg_background_function
 * Returns the background level for the point x, y if a complete
 * gsl_matrix of the background is available.  This function is
 * exposed to the outside through a pointer in the observation
 * structure.
 *
 * Parameters:
 * @param x ad nauseam
 * @param y ad nauseam
 * @param val pointer to a double to leave the background value in
 * @param err pointer to a double to leave the absolute error of val in
 * @param pars a gsl_matrix containing the background
 */
void
fullimg_background_function (const int x, const int y, PIXEL_T * const val,
                             PIXEL_T * const err,
                             const background * const back)
{
  const fullimg_background *const fib = back->pars;

  *val = gsl_matrix_get (fib->bck, (int) rint (x), (int) rint (y));
  if (fib->err)
    {
      *err = gsl_matrix_get (fib->err, (int) rint (x), (int) rint (y));
    }
  else
    {
      *err = 0;
    }
}



/**
 * Function: compute_fullimg_background
 * Computes a background image. All beams with ignore=1 are
 * completely ignored. Anything else is not.
 *
 * Parameters:
 * @param obs         - a pointer to the observation structure to fill out
 * @param oblist      - a list of all objects of the observation
 * @param npoints     - number of points to examine when fitting the background
 * @param interporder - order of the polynomial to fit to the background
 * @param niter_med   - order of the polynomial to fit to the background
 * @param niter_fit   - order of the polynomial to fit to the background
 * @param kappa       - order of the polynomial to fit to the background
 *
 * Returns:
 * @return background - an allocated background structure
 *
 */
background *
compute_fullimg_background (observation *obs, object **oblist,
                            int npoints, int interporder, const int niter_med,
                            const int niter_fit, const double kappa,
                            int nor_flag, const int sm_length, const double fwhm)
{
  fullimg_background *fib;
  background *backg;
  gsl_matrix *bck_mask;
  int i, j;
  //object *const *obp;

  // allocate space for the backgrounds
  fib   = (fullimg_background *)malloc (sizeof (fullimg_background));
  fib->bck = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);
  gsl_matrix_set_all (fib->bck, 0.);
  fib->err = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);
  gsl_matrix_set_all (fib->err, 0.);

  if (nor_flag)
    {
      for (i = 0; i < (int)fib->bck->size1; i++)
        {
          for (j = 0; j < (int)fib->bck->size2; j++)
            {
              gsl_matrix_set (fib->bck, i, j,
                              gsl_matrix_get (obs->grism, i, j));
              gsl_matrix_set (fib->err, i, j,
                              gsl_matrix_get (obs->pixerrs, i, j));
            }
        }
    }

  // allocate memory
  backg = (background *)malloc (sizeof (background));

  // create the mask image
  bck_mask = aperture_mask(obs,oblist);

  // in case that pixels may get dq values,
  // make sure that there will be a dq-array
  //  if (niter_med > 0 || niter_fit > 0)
  //    {
  if (!obs->dq)
    obs->dq = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);

  // initialize the dq-array
  //      gsl_matrix_set_all (obs->dq, 0.0);
  //    }

  if (oblist != NULL)
    {
      // Now compute background for each beam, one after the other
      // go over the whole object list
      i=0;
      while (oblist[i] != NULL) {
        // go over each beam
        for (j = 0; j < oblist[i]->nbeams; j++)
          {
            // check for beams to be neglected
            if (oblist[i]->beams[j].ignore == 1)
              {
                continue;
              }
            else
              {
                // start the background interpolation
                // for a specific beam
                fprintf(stdout,"Computing background of BEAM %d%c.",
                        oblist[i]->ID,BEAM(oblist[i]->beams[j].ID));
                compute_background(obs, oblist[i]->beams[j], bck_mask,
                                   fib, npoints, interporder, niter_med,
                                   niter_fit, kappa);
                fprintf(stdout," Done.\n");
              }
          }
        // increment the counter
        i++;
      }
    }
  // replace all NAN's with values 0.0
  for (i = 0; i < (int)fib->bck->size1; i++)
    {
      for (j = 0; j < (int)fib->bck->size2; j++)
        {
          if(isnan(gsl_matrix_get (fib->bck, i, j)))
            gsl_matrix_set (fib->bck, i, j,0.0);
          if(isnan(gsl_matrix_get (fib->err, i, j)))
            gsl_matrix_set (fib->err, i, j,0.0);
        }
    }

  if (sm_length && fwhm)
    gsmooth_background (bck_mask, sm_length, fwhm, fib);

  // put together the result
  backg->pars = fib;
  backg->bck_func = &fullimg_background_function;

  // free space
  gsl_matrix_free(bck_mask);

  // return the result
  return backg;
}

/**
 * Function: compute_backsub_mask
 * Computes a mask image for the background subtraction. All pixels
 * covered by at least one beam are set to a vale -100000 to distinguish
 * them from pixels which are part of the background.
 *
 * Parameters:
 * @param obs a pointer to the observation structure to fill out
 * @param oblist a list of all objects of the observation
 * @param npoints number of points to examine when fitting the background
 * @param interporder order of the polynomial to fit to the background
 *
 * Returns:
 * @return an allocated background structure
 *
 */
background *
compute_backsub_mask (observation *obs, object **oblist)
{
  fullimg_background *fib;// = malloc (sizeof (fullimg_background));
  background *backg; // = malloc (sizeof (background));
  gsl_matrix *bck_mask;
  int i, j;
  //object *const *obp;
  //int tnbeams;

  // make the mask image
  bck_mask = aperture_mask(obs,oblist);

  // allocate memory for the return structure
  // and fill in some dummy values
  fib = malloc (sizeof (fullimg_background));
  fib->bck = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);
  gsl_matrix_set_all (fib->bck, 0.);
  fib->err = NULL;

  // allocate memory for the return structure
  backg = malloc (sizeof (background));

  // transfer the pixel values from the grism image
  // to the mask image
  for (i = 0; i < (int)fib->bck->size1; i++)
    {
      for (j = 0; j < (int)fib->bck->size2; j++)
        {
          gsl_matrix_set (fib->bck, i, j,
                          gsl_matrix_get (obs->grism, i, j));
        }
    }

  // set pixels occupied by beams
  // to the value -10000000
  for (i = 0; i < (int)fib->bck->size1; i++)
    {
      for (j = 0; j < (int)fib->bck->size2; j++)
        {

          if (gsl_matrix_get (bck_mask, i, j))
            gsl_matrix_set (fib->bck, i, j, -1000000.0);
        }
    }

  // release the dq-structure,
  // otherwise it is saved to
  // the mask image
  if (obs->dq != NULL)
    {
      gsl_matrix_free (obs->dq);
      obs->dq=NULL;
    }

  // release memory
  gsl_matrix_free(bck_mask);

  // assemble the return structure
  backg->pars = fib;
  backg->bck_func = &fullimg_background_function;

  // return the result
  return backg;
}


/**
 * Function: compute_fullimg_global_background
 * The subroutine computes the background image based on all pixels
 * in a column which are not part of an object.
 *
 * Parameters:
 * @param obs         - a pointer to the observation structure to fill out
 * @param oblist      - a list of all objects of the observation
 * @param interporder - order of the polynomial to fit to the background
 *
 * Returns:
 * @return background - an allocated background structure
 *
 */
background *
compute_fullimg_global_background(observation *obs, object **oblist,
                                  int interporder, const int sm_length, const double fwhm)
{
  fullimg_background *fib;
  background *backg;
  gsl_matrix *bck_mask;
  int i, j;
  object *const *obp;
  int tnbeams;

  // allocate memory
  fib = (fullimg_background *)malloc (sizeof (fullimg_background));
  backg = (background *)malloc (sizeof (background));

  // allocate memory
  bck_mask = aperture_mask(obs,oblist);
  fib->bck = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);
  fib->err = gsl_matrix_alloc (obs->grism->size1, obs->grism->size2);

  // initialize the new arrays
  gsl_matrix_set_all (fib->bck, 0.);
  gsl_matrix_set_all (fib->err, 0.);

  /* Count beams in observation */
  if (oblist!=NULL)
    {
      tnbeams = 0;
      for (obp = oblist; *obp; obp++)
        {
          for (i = 0; i < (*obp)->nbeams; i++)
            {
              // if (((*obp)->beams[i]).ignore!=1)
              tnbeams++; /* Count ALL the beams */
            }
        }
    }


  /* Now compute global background using first beam grism info */
  if (oblist != NULL) {
    compute_global_background (oblist, 0, bck_mask,
                               fib, interporder);
  }


  // replace all NAN's with values 0.0
  for (i = 0; i < (int)fib->bck->size1; i++)
    {
      for (j = 0; j < (int)fib->bck->size2; j++)
        {
          if(isnan(gsl_matrix_get (fib->bck, i, j)))
            {
              gsl_matrix_set (fib->bck, i, j,0.0);
            }
          if(isnan(gsl_matrix_get (fib->err, i, j)))
            {
              gsl_matrix_set (fib->err, i, j,0.0);
            }
        }
    }

  // make a Gaussian smoothing
  // if requested
  if (sm_length && fwhm)
    gsmooth_background (bck_mask, sm_length, fwhm, fib);

  // compose the result structure
  backg->pars = fib;
  backg->bck_func = &fullimg_background_function;

  // return the result
  return backg;
}

/**
 * Function: free_fullimg_background
 * Frees a background structure allocated using make_fullimg_background
 *
 * Parameters:
 * @param backg - a pointer to the background structure to free
 */
void
free_fullimg_background (background * backg)
{
  fullimg_background *fib = backg->pars;

  gsl_matrix_free (fib->bck);
  if (fib->err)
    {
      gsl_matrix_free (fib->err);
    }
  free (backg);
  backg = NULL;
}



/**
 * Function: aperture_mask
 * This functions returns an image mask where pixels that are within
 * an aperture are set to the number of aperture they appear.
 * Skips any beams which are set to be ignored (ignore=1)
 *
 * Parameters:
 * @param  obs    - a pointer to an observation structure
 * @param  oblist - a pointer to an object list
 *
 * Returns:
 * @return bck   - a pointer to a gsl_matrix
 */
gsl_matrix *
aperture_mask (observation * const obs, object **oblist)
{

  int x, y,i=0,j;

  double xrel, yrel, width;

  is_in_descriptor iid;
  px_point ll, ur;
  gsl_matrix *bck;
  sectionfun sf;

  // create the image,
  // set all values to zero
  bck = gsl_matrix_alloc(obs->grism->size1,obs->grism->size2);
  gsl_matrix_set_all(bck,0.0);

  // Return a zero image when there is no
  // beam in the list
  if (oblist==NULL)
    return bck;

  // go over each object
  while(oblist[i]!=NULL)
    {
      // go over each beam
      for(j=0;j<oblist[i]->nbeams;j++)
        {
          // continue if the beam is to be ignored
          if ((oblist[i]->beams[j]).ignore==1)
            continue;

          // check whther the trace is second order or
          // even higher
          if ((oblist[i]->beams[j]).spec_trace->type > 1)
            {
              // create the section function for second order
              // traces such as FORSII
              fill_in_sectionfun (&sf, (oblist[i]->beams[j]).orient,
                                  &(oblist[i]->beams[j]));
            }
          else
            {
              // create the descriptor to set up the quadrangle routines
              fill_is_in_descriptor (&iid, (oblist[i]->beams[j]).corners);
            }

          // check for the corners of the quadrangle.
          // go over each point in x and y.
          quad_to_bbox ((oblist[i]->beams[j]).corners, &ll, &ur);
          for (x = ll.x; x <= ur.x; x++)
            {
              for (y = ll.y; y <= ur.y; y++)
                {

                  // neglect the pixel if it is a]outside the
                  // image area
                  if ((x < 0) || (y < 0) || (x >= (int)obs->grism->size1)
                      || (y >= (int)obs->grism->size2))
                    continue;

                  // check the order of the trace polynomial
                  if ((oblist[i]->beams[j]).spec_trace->type > 1)
                    {

                      // determine the trace distance and
                      // enhance the image value if the pixel
                      // satisfies the distance criterium
                      xrel = x-(oblist[i]->beams[j]).refpoint.x;
                      yrel = y-(oblist[i]->beams[j]).refpoint.y;
                      width = (oblist[i]->beams[j]).width+0.5;
                      if (tracedist_criteria(xrel, yrel, &sf, (oblist[i]->beams[j]).spec_trace, width))
                        gsl_matrix_set (bck,x,y,gsl_matrix_get (bck,x,y)+1);
                    }
                  else
                    {
                      //                      if (oblist[i]->ID == 11 && x == 59)
                      //                        fprintf(stdout, "xx: %i, yy: %i: %i\n", x, y, is_in (x, y, &iid));

                      // check whether the pixel is inside the quadrangle
                      // which defines the beam are, and enhance the image
                      // value if yes.
                      if (is_in (x, y, &iid))
                        gsl_matrix_set (bck,x,y,gsl_matrix_get (bck,x,y)+1);
                    }
                }
            }
        }

      // enhance the object counter
      i++;
    }

  // return the resulting image
  return bck;
}


/**
 * Function: background_to_FITSimage
 * Function to write the data and error content of a backgound
 * observation pars component (if it is of type fullimg_background only)
 * into the main HDU (data) and first extension (error) of a FITS file.
 * A GQ array can be appended if an observation with a non NULL DQ is
 * passed to this function
 *
 * Parameters:
 * @param filename - name of the image
 * @param bck      - filled background structure
 * @param obs      - the observation the background is based upon
 */
void
background_to_FITSimage (char filename[], background * bck, observation *obs)
{
  fitsfile *output;
  long naxes[2];
  int f_status = 0;
  PIXEL_T *storage, *dp;
  int x, y;
  fullimg_background *pars;
  int hdunum,hdutype;

  pars = (fullimg_background *) bck->pars;

  // Open the file for creating/appending
  create_FITSimage (filename, 1);
  fits_open_file (&output, filename, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: " "Could not open file: %s",
                   filename);
    }

  // count the number of extentions
  fits_get_num_hdus (output, &hdunum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could not get"
                   " number of HDU from: %s", filename);
    }

  // Move to last HDU
  fits_movabs_hdu (output, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could not mov"
                   " to HDU number %d in file: %s", hdunum, filename);
    }
  /* Get current HDU number */
  fits_get_hdu_num (output, &hdunum);

  // Deal with the SCI part of the background
  naxes[0] = pars->bck->size1;
  naxes[1] = pars->bck->size2;

  // Allocate storage room
  if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");

  // Populate the storage array
  dp = storage;
  for (y = 0; y < naxes[1]; y++)
    {
      for (x = 0; x < naxes[0]; x++)
        {
          *dp = gsl_matrix_get (pars->bck, x, y);
          dp++;
        }
    }
  /* create HDU extname */
  fits_create_img (output, -32, 2, naxes, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could create SCI HDU in file: %s", filename);
    }
  fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                  &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "gsl_to_FITSimage: Could write SCI HDU in file: %s", filename);
    }

  /* write the HDU EXTNAME */
  {
    char comment[FLEN_COMMENT];
    char str[FLEN_KEYWORD];
    strcpy (str, "SCI");
    strcpy (comment, "Extension name");
    fits_write_key_str (output, "EXTNAME", str, comment, &f_status);
  }
  free (storage);
  storage = NULL;

  if (pars->err)
    {
      /* Deal with the error part of the background */
      naxes[0] = pars->err->size1;
      naxes[1] = pars->err->size2;
      /* Allocate storage room */
      if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
        {
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
        }
      /* Populate the storage array */
      dp = storage;
      for (y = 0; y < naxes[1]; y++)
        {
          for (x = 0; x < naxes[0]; x++)
            {
              *dp = gsl_matrix_get (pars->err, x, y);
              dp++;
            }
        }
      /* create HDU extname */
      /* Get current HDU number */
      fits_get_hdu_num (output, &hdunum);

      fits_create_img (output, -32, 2, naxes, &f_status);
      fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                      &f_status);

      /* Get current HDU number */
      fits_get_hdu_num (output, &hdunum);
      /* write the HDU EXTNAME */
      {
        char comment[FLEN_COMMENT];
        char str[FLEN_KEYWORD];
        strcpy (str, "ERR");
        strcpy (comment, "Extension name");
        fits_write_key_str (output, "EXTNAME", str, comment, &f_status);
      }
      free (storage);
      storage = NULL;
    }

  /* Deal with the DQ part of the background */
  if (obs->dq != NULL)
    {
      naxes[0] = obs->dq->size1;
      naxes[1] = obs->dq->size2;
      /* Allocate storage room */
      if (!(storage = malloc (naxes[0] * naxes[1] * sizeof (double))))
        {
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
        }
      /* Populate the storage array */
      dp = storage;
      for (y = 0; y < naxes[1]; y++)
        {
          for (x = 0; x < naxes[0]; x++)
            {
              *dp = gsl_matrix_get (obs->dq, x, y);
              dp++;
            }
        }
      /* create HDU extname */
      fits_create_img (output, 16, 2, naxes, &f_status);
      fits_write_img (output, TFLOAT, 1, naxes[0] * naxes[1], storage,
                      &f_status);
      /* write the HDU EXTNAME */
      {
        char comment[FLEN_COMMENT];
        char str[FLEN_KEYWORD];
        strcpy (str, "DQ");
        strcpy (comment, "Extension name");
        fits_write_key_str (output, "EXTNAME", str, comment, &f_status);
      }
      free (storage);
      storage = NULL;
    }

  /* close file */
  fits_close_file (output, &f_status);
}


/**
 * Function: gsmooth_background
 * Smooth all interpolated pixel in the background using a Gaussian
 * function. The smoothing is done exclusively towards the x-values.
 * The smoothing should help to reduce the noise from the limited number
 * of background pixels.
 *
 * Parameters:
 * @param  bck_mask     - the background mask
 * @param  smoot_length - number of pixels on either side to use for smoothing
 * @param  fwhm         - fwhm of the Gaussian
 * @param  fib          - the background structure
 */
void
gsmooth_background (const gsl_matrix *bck_mask, const int smooth_length,
                    const double fwhm, fullimg_background *fib)
{
  double efactor;

  int ix, iy;

  gsl_vector *pixvalues;
  gsl_vector *weights;
  gsl_vector *pmask;

  gsl_matrix *new_bck;

  // allocate memory for the vectors
  pixvalues = gsl_vector_alloc(2 * smooth_length + 1);
  weights   = gsl_vector_alloc(2 * smooth_length + 1);
  pmask     = gsl_vector_alloc(2 * smooth_length + 1);

  // allocate the new background
  // and initialize it
  new_bck = gsl_matrix_alloc(bck_mask->size1, bck_mask->size2);
  gsl_matrix_set_all(new_bck, 0.0);

  // prepare the Gaussian
  efactor = compute_efactor(fwhm);

  // fill the weights
  for (ix=0; ix < (int)weights->size; ix++)
    gsl_vector_set(weights, ix, compute_gvalue((double)ix - (double)smooth_length, efactor));


  // go over all rows
  for (iy=0; iy < (int)fib->bck->size2; iy++)
    {
      // go over all columns
      for (ix=0; ix < (int)fib->bck->size1; ix++)
        {
          // check whether the pixel IS part of a beam
          if (!(gsl_matrix_get(bck_mask, ix, iy) != 0.0
                && gsl_matrix_get(fib->bck, ix, iy) !=0.0))
            {
              // transfer the old background value
              gsl_matrix_set(new_bck, ix, iy, gsl_matrix_get(fib->bck, ix, iy));
            }
          else
            {
              // fill the pixel and mask values
              fill_pixvalues(bck_mask, smooth_length, fib, ix, iy, pixvalues, pmask);

              // fill in the weighted mean
              gsl_matrix_set(new_bck, ix, iy, get_weighted_mean(pixvalues, weights, pmask));
            }
        }
    }

  // release the allocated dspace
  gsl_vector_free(pixvalues);
  gsl_vector_free(weights);
  gsl_vector_free(pmask);

  // release the memory of the
  // old background
  gsl_matrix_free(fib->bck);

  // transfer the new background
  // to the background structure
  fib->bck = new_bck;
}

/**
 * Function: get_weighted_mean
 * The function computes the weighted mean of values stored in a value vector
 * and a weight vector. As mask vector marks values pixels not to be considered.
 * Using a separate mask vectorhas the advantage the weights can be kept
 * constant and do not have to be re-calculated in repeated runs.
 *
 * Parameters:
 * @param pixvalues - vector with pixel values
 * @param weights   - vector with weights
 * @param pmask     - mask vector
 *
 * Returns:
 * @return sum/www  - the weighted mean
 */
double
get_weighted_mean(const gsl_vector *pixvalues, const gsl_vector *weights,
                  const gsl_vector *pmask)
{
  int index;

  // initialize the total
  // sum and weight
  double sum=0.0;
  double www=0.0;

  for (index=0; index < (int)pixvalues->size; index ++)
    {

      if (gsl_vector_get(pmask, index))
        {
          // enhance the total sum
          sum += gsl_vector_get(pixvalues, index) * gsl_vector_get(weights, index);

          // enhance the total weight
          www += gsl_vector_get(weights, index);
        }
    }

  // return the total sum,
  // divided by the total weight
  return sum / www;
}


/**
 * Function: fill_pixvalues
 * The function provides the essential information for Gaussian smoothing
 * for a single pixel. It fills a vector with the values of all pixels
 * within the smoothing length. Not interpolated pixels are excluded.
 * A mask vector provides the information on which position is filled
 * with pixel values.
 *
 * Parameters:
 * @param  bck_mask     - the background mask
 * @param  smoot_length - number of pixels on either side to use for smoothing
 * @param  fib          - the background structure
 * @param  bck_mask     - the background mask
 * @param  ix           - x-value of central pixel
 * @param  iy           - y-value of central pixel
 * @param  pixvalues    - vector for pixel values
 * @param  pmask        - vector for pixel mask
 */
void
fill_pixvalues(const gsl_matrix *bck_mask, const int smooth_length,
               const fullimg_background *fib, const int ix, const int iy,
               gsl_vector *pixvalues, gsl_vector *pmask)
{
  int iact;
  int index;

  // initialize the pixel values
  // and the mask
  gsl_vector_set_all(pixvalues, 0.0);
  gsl_vector_set_all(pmask, 0.0);

  // iterate over the x-direction
  index=0;
  for (iact=ix-smooth_length; iact<=ix+smooth_length; iact++)
    {
      // check if you are within the chip
      if (iact > -1 && iact < (int)bck_mask->size1)
        {
          // check whether the pixel was interpolated
          // and whether the background is non-zero
          if (gsl_matrix_get(bck_mask, iact, iy) != 0.0
               && gsl_matrix_get(fib->bck, iact, iy) != 0.0)
            {
              // if yes, get the pixel value and set the mask
              gsl_vector_set(pixvalues, index, gsl_matrix_get(fib->bck, iact, iy));
              gsl_vector_set(pmask,     index, 1.0);
            }
        }
      // enhance the vector index
      index++;
    }
}


/**
 * Function: compute_gvalue
 * The function computes the values of a Gauss function
 * [exp(factor * xdiff^2)]. NO normalization factor is applied.
 *
 * Parameters:
 * @param xdiff   - the value [x-x_0]
 * @param efactor - the factor for the exponent
 *
 * Returns:
 * @return value - the Gaussian value
 */
double
compute_gvalue(const double xdiff, const double efactor)
{
  double value;

  // just compose the exp-function
  value = exp(efactor * xdiff * xdiff);

  // return the value
  return value;
}

/**
 * Function: compute_efactor
 * The function computes the appropriate factor of the Gaussian for any
 * FWHM given as input. This speeds up any later computation of the
 * Gaussian.
 *
 * Parameters:
 * @param fwhm - the input FWHM
 *
 * Returns:
 * @return expfactor - the factor for the Gaussian
 */
double
compute_efactor(const double fwhm)
{
  double sigma;
  double expfactor;

  // get the sigma value [sig = fwhm / (2 sqrt(2 ln(2)))]
  sigma = fwhm / 2.3458;

  // compute the factor for the Gaussian
  expfactor =  -1.0 / (2.0* sigma * sigma);

  // return the factor
  return expfactor;
}
