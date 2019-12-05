/**
 * File: nicback_utils.c
 * Subroutine for the NICMOS background subtraction
 *
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_FITScards.h"
#include "spce_sect.h"
#include "spce_fitting.h"
#include "spce_is_in.h"
#include "spc_back.h"
#include "spce_pathlength.h"
#include "trfit_utils.h"
#include "nicback_utils.h"

/**
 * Function: make_nicmos_back
 * The subroutine computes and stores the background image for a specific
 * NICMOS grism image. The non-masked pixels in the grism image are
 * correlated with their corresponding values in the master background image.
 * A linear relation is fitted to the value pairs. Deviating points are
 * excluded using kappa-sigma clipping.
 * Using the slope of the final fit, the background image is computed
 * and stored as a fits image.
 *
 * Parameters:
 * @param obs         - the observation for the grism image
 * @param msk_name    - the full name of the mask file
 * @param master_bck  - the full name of the master sky background
 * @param bck_name    - the full name of the background image
 *
 * Returns:
 * @return -
 */
void
make_nicmos_back(const observation * const obs, const char *msk_name,
		 const char *master_bck, char bck_name[], char plist_name[],
		 const char corr_bck[])
{
  px_point npixels;

  double skypix_frac;

  gsl_matrix *msk_img;
  gsl_matrix *mbck_img;
  gsl_matrix *corr_img=NULL;
  gsl_matrix *grism_bck;

  gsl_vector *lfit;

  fitbck_data *fbck_data;
  //int f_status=0;

  FITScards *cards;

  char  ID[MAXCHAR];

  // fix the extension name
  sprintf (ID, "BCK");

  // load the master background image
  fprintf(stdout,"Loading DATA from: %s...", master_bck);
  mbck_img = FITSimage_to_gsl(master_bck, 1, 1);
  fprintf(stdout,". Done.\n");

  if (strlen(corr_bck) > 0)
    {
      // load the correction background image
      fprintf(stdout,"Loading DATA from: %s...", corr_bck);
      corr_img = FITSimage_to_gsl(corr_bck, 1, 1);
      fprintf(stdout,". Done.\n");
    }

  // load the mask image
  fprintf(stdout,"Loading DATA from: %s...", msk_name);
  msk_img = FITSimage_to_gsl(msk_name, 2, 1);
  fprintf(stdout,". Done.\n");

  // get the dimension of
  // the grism image
  npixels = get_npixel(obs);

  // allocate memory for the fit data
  fbck_data = alloc_fitbck_data(npixels.x * npixels.y);

  // fill the fit data into the structure
  fill_fbck_data(obs, msk_img, mbck_img, fbck_data, plist_name, corr_img);

  // make the fit, including
  // kappa-sigma iterations
  lfit = make_ksig_linfit(fbck_data);

  // report the result of the fit onto the screen
  fprintf(stdout, "\nFitting result: c0 = %f +- %f , c1 = %f +- %f , chi^2 = %f\n\n",
	  gsl_vector_get(lfit, 1), gsl_vector_get(lfit, 3),
	  gsl_vector_get(lfit, 2), gsl_vector_get(lfit, 4), gsl_vector_get(lfit, 6));

  // compute the background for the grism frame
  grism_bck = compute_nicmos_back(mbck_img, gsl_vector_get(lfit, 2),
				  gsl_vector_get(lfit, 1), corr_img);

  // write the background image to a file
  fprintf(stdout,"Writing data to: %s...", bck_name);
  gsl_to_FITSimage (grism_bck, bck_name, 1, ID);
  fprintf(stdout,". Done.\n");

  // determine the fraction of pixels
  // used in the background determination
  skypix_frac = get_skypix_frac(fbck_data);

  // create the fits keywords
  cards = nicbck_info_to_FITScards(skypix_frac, gsl_vector_get(lfit, 2),
				   gsl_vector_get(lfit, 1));

  // transfer the fits keywords
  // to the background image
  put_FITS_cards (bck_name, 2, cards);

  // release memory
  free_fitbck_data(fbck_data);
  gsl_matrix_free(mbck_img);
  if (corr_img != NULL)
    gsl_matrix_free(corr_img);
  gsl_matrix_free(msk_img);
  gsl_matrix_free(grism_bck);
  gsl_vector_free(lfit);
  free_FITScards(cards);
}


/**
 * Function: alloc_fitbck_data
 */
fitbck_data *
alloc_fitbck_data(const int n_data)
{
  fitbck_data *fbck_data;

  fbck_data = (fitbck_data *)malloc(sizeof(fitbck_data));

  fbck_data->x_values = (double *)malloc(n_data * sizeof(double));
  fbck_data->y_values = (double *)malloc(n_data * sizeof(double));
  fbck_data->e_values = (double *)malloc(n_data * sizeof(double));
  fbck_data->x_pos = (int *)malloc(n_data * sizeof(int));
  fbck_data->y_pos = (int *)malloc(n_data * sizeof(int));

  fbck_data->n_data = n_data;

  return fbck_data;
}

/**
 * Function: get_skypix_frac
 * The subroutine computes the fraction of image pixels which were used for
 * determining the scaling factor of the master background and the computation
 * of the bacground image. Pixels not used were masked out, since they are part
 * of an object beam or rejected in the kappa-sigma clipping.
 *
 * Parameters:
 * @param fbck_data - the data used for the fit
 *
 * Returns:
 * @return skypix_frac - fraction of pixels used for background determination
 */
double
get_skypix_frac(const fitbck_data *fbck_data)
{
  int index;
  int nvalues=0;

  double skypix_frac=0.0;

  // go over all data values
  for (index=0; index < fbck_data->n_data; index++)
    {
      // check whether the weight is not NULL
      if (fbck_data->e_values[index])
	// enhance the number of used pixels
	nvalues++;
    }

  // compute the fration of used pixels
  skypix_frac = (double)nvalues / (double)fbck_data->n_data;

  // return the fraction
  return skypix_frac;
}

/**
 * Function: compute_nicmos_back
 * The subroutine applies a scaling factor to a master background
 * in order to get the background image for a specific grism image.
 *
 * Parameters:
 * @param mbck_img - the master background pixel
 * @param factor   - the scaling factor
 *
 * Returns:
 * @return grism_bck - the background image as gsl-matrix
 */
gsl_matrix *
compute_nicmos_back(const gsl_matrix *mbck_img, const double factor,
		    const double offset, const gsl_matrix *corr_img)
{
  gsl_matrix *grism_bck;

  int ii, jj;

  // allocate the space for the background image
  grism_bck = gsl_matrix_alloc(mbck_img->size1, mbck_img->size2);

  // go over each row
  for (ii=0; ii < (int)mbck_img->size1; ii++)
    {
      // go over each column
      for (jj=0; jj < (int)mbck_img->size2; jj++)
	{
	  // compute and set the value in the background image
	  if (corr_img != NULL)
	    // using the pedestal image
	    gsl_matrix_set(grism_bck, ii, jj,
	        (gsl_matrix_get(mbck_img, ii, jj) + gsl_matrix_get(corr_img, ii, jj)) * factor);
	  else
	    // using the fit values only
	    gsl_matrix_set(grism_bck, ii, jj, gsl_matrix_get(mbck_img, ii, jj) * factor + offset);
	}
    }

  // return the background image
  return grism_bck;
}

/**
 * Function: compute_nicmos_back
 * The subroutine iteratively fits a linear relation to the pairs of
 * grism image - master background pixel values. After each fit,
 * pairs with large deviations are determined and excluded.
 * The result of the final fit is returned.
 *
 * Parameters:
 * @param fbck_data - the data used for the fits
 *
 * Returns:
 * @return lfit - the results of the final fit
 */
gsl_vector *
make_ksig_linfit(fitbck_data *fbck_data)
{
  gsl_vector *lfit=NULL;

  int index=0;

  int clipped=1;

  // check whether iterations still
  // must be done and whether
  // the data was changed from clipping
  while (index < N_KSIGBCK_ITER && clipped)
    {
      if (check_fbck_data(fbck_data))
	{
	  // make a non-weighted linear fit
	  lfit = det_vector_linear(fbck_data->x_values, fbck_data->y_values, fbck_data->e_values,
				   fbck_data->n_data, 0);

	  // make a clipping iteration
	  clipped = clipp_fbck_data(fbck_data, lfit, N_KSIGBCK_KAPPA);

	  // inhance the clipping counter
	  index++;
	}
      else
	{
	  // release the space
	  // for the old vector
	  if (lfit)
	    gsl_vector_free(lfit);

	  // get a new dummy vector
	  lfit = get_fbck_defaults();

	  break;
	}
    }

  // check whether the scale is within the boundaries
  if (gsl_vector_get(lfit, 2) < BCK_SCALE_MIN || gsl_vector_get(lfit, 2) > BCK_SCALE_MAX)
    {

      fprintf(stdout, "aXe_NICBACK: Background scale is out of bounds. It is re-adjusted.\n");
      // release the space
      // for the old vector
      gsl_vector_free(lfit);

      // get a new dummy vector
      lfit = get_fbck_defaults();
    }

  // return the fit result
  return lfit;
}

int
check_fbck_data(fitbck_data *fbck_data)
{
  int check=1;
  int ix_min=60000;
  int ix_max=0;
  int iy_min=60000;
  int iy_max=0;
  int index;
  double frac;
  double x_ext, y_ext;

  // get the fraction of background pixel
  frac =  get_skypix_frac(fbck_data);

  if (frac < FRAC_BPIX_MIN)
    {
      check = 0;
      fprintf(stdout, "aXe_NICBACK: Not enough background pixels. Fraction: %e\n", frac);
    }
  else
    {
      // go over all data values
      for (index=0; index < fbck_data->n_data; index++)
	{
	  // check whether the weight is not NULL
	  if (fbck_data->e_values[index])
	    {
	      if (fbck_data->x_pos[index] > ix_max)
		ix_max = fbck_data->x_pos[index];
	      if (fbck_data->x_pos[index] < ix_min)
		ix_min = fbck_data->x_pos[index];
	      if (fbck_data->y_pos[index] > iy_max)
		iy_max = fbck_data->y_pos[index];
	      if (fbck_data->y_pos[index] < iy_min)
		iy_min = fbck_data->y_pos[index];
	    }
	}
      x_ext = (float)(ix_max-ix_min) / NPIX_DIM1;
      y_ext = (float)(iy_max-iy_min) / NPIX_DIM2;

      if (x_ext < AREA_EXTEND_MIN || y_ext < AREA_EXTEND_MIN)
	{
	  check = 0;
	  fprintf(stdout, "aXe_NICBACK: Not enough coverage. xcover,ycover: %e,%e\n", x_ext, y_ext);
	}
    }

  return check;
}

gsl_vector*
get_fbck_defaults()
{
  gsl_vector *ret;

  ret = gsl_vector_alloc(10);

  /*
  gsl_vector_set(ret, 0, xean);
  gsl_vector_set(ret, 1, c0);
  gsl_vector_set(ret, 2, c1);
  gsl_vector_set(ret, 3, cov00);
  gsl_vector_set(ret, 4, cov01);
  gsl_vector_set(ret, 5, cov11);
  gsl_vector_set(ret, 6, chisq);
  */

  gsl_vector_set(ret, 0, 0.0);
  gsl_vector_set(ret, 1, 0.0);
  gsl_vector_set(ret, 2, 1.0);
  gsl_vector_set(ret, 3, 1.0);
  gsl_vector_set(ret, 4, 0.0);
  gsl_vector_set(ret, 5, 1.0);
  gsl_vector_set(ret, 6, 0.0);

  return ret;
}

/**
 * Function: clipp_fbck_data
 * The subroutine applies kappa-sigma clipping to a set of data points.
 * The difference between the measured values and the expected values
 * according to linear regression are determined.
 * Data which is outside of the accepted interval is excluded by
 * setting the weight to '0.0'. If at least one value was clipped,
 * '1' is returned, otherwise '0' (if all values are within the
 * accepted interval.
 *
 * Parameters:
 * @param fbck_data - the data used for the fits
 * @param lfit      - the result of the fit
 * @param kappa     - the kappa value
 *
 * Returns:
 * @return clipped - boolean to indicate clipped values
 */
int
clipp_fbck_data(fitbck_data *fbck_data, const gsl_vector *lfit, const float kappa)
{
  int clipped=0;
  int nclipps=0;

  int ndata=0;
  int index;

  double dy_mean=0.0;
  double stdev=0.0;
  double absdev;
  double dy_act;

  double *y_diff;

  // allocate memory for the tmp-vector
  y_diff = (double *) malloc(fbck_data->n_data * sizeof(double));

  // go over all data values
  for (index=0; index < fbck_data->n_data; index++)
    {
      // check whether the weight is not NULL
      if (fbck_data->e_values[index])
	{
	  // compute the difference from the data point
	  // to the prediction from the fit
	  dy_act = (fbck_data->x_values[index] - gsl_vector_get(lfit, 0)) * gsl_vector_get(lfit, 2)
	    + gsl_vector_get(lfit, 1)-fbck_data->y_values[index];

	  // fill the tmp vector
	  y_diff[ndata] = dy_act;

	  // pre-comute the mean
	  dy_mean += dy_act;

	  // enhance the counter
	  // of non-NULL values
	  ndata++;
	}
    }

  // compute the mean offset
  dy_mean = dy_mean / (double)ndata;

  // compute the standard deviation
  stdev = comp_stdev_from_array(y_diff, ndata, dy_mean);

  // get the maximum allowed
  // deviation
  absdev = kappa*stdev;

  // go over the data
  for (index=0; index < fbck_data->n_data; index++)
    {
      // check whether the data has weigth
      if (fbck_data->e_values[index])
	{
	  // compute the absolute difference from the data point
	  // to the prediction from the fit
	  dy_act = fabs((fbck_data->x_values[index] - gsl_vector_get(lfit, 0)) * gsl_vector_get(lfit, 2)
			+ gsl_vector_get(lfit, 1)-fbck_data->y_values[index]);

	  // check whether the difference
	  // is larger than allowed
	  if (dy_act > absdev)
	    {
	      // make the weight NULL
	      fbck_data->e_values[index] = 0.0;

	      // set the clipped flagg
	      clipped=1;
	      nclipps++;
	    }
	}
    }

  fprintf(stdout, "Number of clipped pixels: %i\n", nclipps);

  // release memory
  free(y_diff);

  // return the integer
  // indicating clipped values
  return clipped;
}

/**
 * Function: comp_stdev_from_array
 * The subroutine computes the standard deviation
 * of a data vector, using the known mean value.
 *
 * Parameters:
 * @param data  - the data vector
 * @param ndata - number of data points to use
 * @param mean  - mean value of the data points
 *
 * Returns:
 * @return clipped - boolean to indicate clipped values
 */
double
comp_stdev_from_array(const double *data, const int ndata, const double mean)
{
  int index=0;

  double stdev=0.0;

  // go over the data array
  for (index=0; index < ndata; index++)
    {
      // compute and add the square difference to the mean
      stdev += (data[index] - mean) * (data[index] - mean);
    }

  // compute the average square difference
  if (ndata > 1)
    stdev /= ((double)ndata-1);

  // compute the final
  // standard deviation
  stdev = sqrt(stdev);

  // return the stdev
  return stdev;
}


/**
 * Function: make_nicmos_back
 * The subroutine creates and fills a structure with pairs of correpsonding
 * from a grism image and a master sky background. Dq-masked pixels values
 * and pixels which are part of an object beam are used, but marked with
 * a zero weight.
 *
 * Parameters:
 * @param obs       - the observation for the grism image
 * @param msk_image - the maks image values
 * @param mbck_img  - the master background values
 * @param fbck_data    - structure with the data pairs
 *
 * Returns:
 * @return -
 */
void
fill_fbck_data(const observation * const obs, const gsl_matrix *msk_img,
	       const gsl_matrix *mbck_img, fitbck_data *fbck_data, char plist_name[],
	       const gsl_matrix *corr_img)
{
  int ix=0;
  int iy=0;
  int index=0;
  char Buffer[10240];
  FILE *fout;

  int px_min=20;
  int px_max=180;
  int py_min=50;
  int py_max=245;
  //int px_min=0;
  //int px_max=255;
  //int py_min=0;
  //int py_max=255;
  //10:240,15:245[10:180,20:245][20:180,50:245]

  // intitialize the index for
  // the fit_data structure
  index = 0;

  double diffval;

  fout = fopen(plist_name, "w");

  // go over all pixels
  // in the grism image
  for (ix=0; ix < (int)obs->grism->size1; ix++)
    {
      for (iy=0; iy < (int)obs->grism->size2; iy++)
	{
	  // set the master background pixel as independent variable
	  fbck_data->x_values[index] = gsl_matrix_get(mbck_img, ix, iy);

	  if (corr_img != NULL)
	    // subtract the correction value from the background value
	    diffval = gsl_matrix_get(obs->grism, ix, iy) - gsl_matrix_get(corr_img, ix, iy);
	  else
	    // use the naked background value
	    diffval = gsl_matrix_get(obs->grism, ix, iy);

	  // set the grism image pixel as dependent variable
	  //fbck_data->y_values[index] = gsl_matrix_get(obs->grism, ix, iy);
	  // set the difference as dependent variable
	  fbck_data->y_values[index] = diffval;

	  // set the x- and y-positions
	  fbck_data->x_pos[index] = ix;
	  fbck_data->y_pos[index] = iy;

	  // check whether the pixel was masked out
	  // or is part of an object
	  if (isnan(gsl_matrix_get(obs->grism, ix, iy))
		    || gsl_matrix_get(msk_img, ix, iy) < -9.0e+05)
	    // mask the pixel in the array
	    fbck_data->e_values[index] = 0.0;
	  else
	    {
	      // give the pixel full weight
	      //x,y,value_back,value_grism_imag
	      fbck_data->e_values[index] = 1.0;
	      sprintf (Buffer, "%i %i %e %e\n",ix, iy, gsl_matrix_get(obs->grism, ix, iy), gsl_matrix_get(mbck_img, ix, iy));
	      fputs (Buffer, fout);
	    }

	  // check whether the pixel is in the
	  // allowed are
	  if (ix < px_min || ix > px_max
	      || iy < py_min || iy > py_max)
	    // mask the pixel in the array
	    fbck_data->e_values[index] = 0.0;
	  // enhance the counter
	  index++;
	}
    }

  fclose (fout);
}

/*
 * Function: free_fitbck_data
 * The function releases the memory allocated in
 * a fit-data structure
 *
 * Parameters:
 * @param fbck_data - the allocated structure
 *
 * Returns:
 * @return -
 */
void
free_fitbck_data(fitbck_data *fbck_data)
{
  free(fbck_data->x_values);
  free(fbck_data->y_values);
  free(fbck_data->e_values);
  free(fbck_data->x_pos);
  free(fbck_data->y_pos);

  free(fbck_data);

  fbck_data = NULL;
}
