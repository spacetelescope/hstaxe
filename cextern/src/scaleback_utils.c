/**
 * File: scaleback_utils.c
 * Subroutine for the WFC3 background scaling
 *
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
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
#include "scaleback_utils.h"

/**
 * Function: fit_to_FITScards
 * Fills the results of the fitting into a set of fits cards.
 *
 * Parameters:
 * @param bck_vals   - vector with the fit results
 * @param npixels    - dimension of all images
 *
 * Returns:
 * @return cards - the list of fits header cards
 */
FITScards *fit_to_FITScards(const gsl_vector* bck_vals, const px_point npixels)
{
    char templt[FLEN_CARD];
    int i=0,keytype, f_status=0;
    int npix_tot;

    FITScards *cards;

    // allocate the cards
    cards = allocate_FITScards(7);

    // compute the number of pixels
    npix_tot = npixels.x * npixels.y;

    i=0;
    sprintf(templt,"SCALVAL = %e / computed scale value", gsl_vector_get(bck_vals, 0));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"SCALERR = %e / error for scale value", gsl_vector_get(bck_vals, 1));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"NPIXINI = %e / initial fill value", gsl_vector_get(bck_vals, 2));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"NPIXFIN = %e / final fill value", gsl_vector_get(bck_vals, 3));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"FRACINI = %e / initial fill value", gsl_vector_get(bck_vals, 2) / (float)npix_tot);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"FRACFIN = %e / final fill value", gsl_vector_get(bck_vals, 3) / (float)npix_tot);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"NITER = %d / number of iterations", (int)gsl_vector_get(bck_vals, 4));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

    // return the filled cards
    return cards;
}

/**
 * Function: make_scale_back
 * Loads the various image data inputs and performs a fit using a kappa-sigma
 * clipping iteration. Then a scaled version of one of the input images
 * is produced and written to disk. Optionally, an ASCII list with the
 * pixel values before the fitting is produced.
 *
 * Parameters:
 * @param grism_image - full pathname to the grism image
 * @param grism_mask  - full pathname to the mask image
 * @param conf_file   - full pathname to the configuration file
 * @param scale_image - full pathname to the scaling (=master sky) image
 * @param bck_image   - full pathname to the scaled (=background) image
 * @param plist_name  - full pathname to the pixel list
 * @param scale_to_master - integer/boolean to specify what should be scaled as output
 * @param make_plis   - integer/boolean to request a pixel list as output
 *
 * Returns:
 * @return -
 */
void
make_scale_back(char grism_image[], const char grism_mask[], char conf_file[],
		const char scale_image[], char bck_image[], const char plist_name[],
		const int scale_to_master, const int make_plis)
{
  px_point npixels;

  gsl_matrix *gr_img;
  gsl_matrix *gr_dqval;
  gsl_matrix *gr_mask;
  gsl_matrix *sc_img;
  gsl_matrix *bck_img;

  gsl_vector *bck_vals;

  fitbck_data *fbck_data;
  //int f_status=0;

  FITScards *cards;

  aperture_conf  *conf;

  char  ID[MAXCHAR];

  // fix the extension name
  sprintf (ID, "BCK");

  // read the configuration file
  // determine all extensions
  conf = get_aperture_descriptor(conf_file);
  get_extension_numbers(grism_image, conf, conf->optkey1, conf->optval1);

  // load the scale image
  fprintf(stdout,"Loading DATA from: %s...", scale_image);
  sc_img = FITSimage_to_gsl(scale_image, 1, 1);
  fprintf(stdout,". Done.\n");

  // load the grism image
  fprintf(stdout,"Loading DATA from: %s...", grism_image);
  gr_img = FITSimage_to_gsl(grism_image, conf->science_numext, 1);
  fprintf(stdout,". Done.\n");

  if (conf->dq_numext < 0)
    {
    // allocate the space for the dq-image
    // set all values to 0.0
    gr_dqval = gsl_matrix_alloc(gr_img->size1, gr_img->size2);
    gsl_matrix_set_all (gr_dqval, 0.0);
    }
  else
    {
    // load the scale dq values
    fprintf(stdout,"Loading DQ from: %s...", grism_image);
    gr_dqval = FITSimage_to_gsl(grism_image, conf->dq_numext, 1);
    fprintf(stdout,". Done.\n");
    }

  // load the grism mask
  fprintf(stdout,"Loading DATA from: %s...", grism_mask);
  gr_mask = FITSimage_to_gsl(grism_mask, 2, 1);
  fprintf(stdout,". Done.\n");

  // get the number of pixels
  npixels.x = sc_img->size1;
  npixels.y = sc_img->size2;

  // report the number of pixels
  //fprintf(stdout,"Loading DATA from: %i pix\n", npixels.x * npixels.y);

  // allocate memory for the fit data
  fbck_data = alloc_fitbck_data(npixels.x * npixels.y);

  // fill the data into the structure
  if (scale_to_master)
    fill_cont_data(gr_img, gr_dqval, gr_mask, sc_img, fbck_data, conf->dqmask);
  else
    fill_mask_data(gr_img, gr_dqval, gr_mask, sc_img, fbck_data, conf->dqmask);


  // make a pixel list if desired
  if (make_plis)
    print_plis(fbck_data, plist_name);

  // make the fit
  bck_vals = make_ksig_scalefit(fbck_data);

  // report the result of the fit onto the screen
  fprintf(stdout, "\nScale result image %s : c0 = %f +- %f", grism_image, gsl_vector_get(bck_vals, 0), gsl_vector_get(bck_vals, 1));
  fprintf(stdout, "\nInitial fill factor: %.1f%%, final: %.1f%%", 100.0 * gsl_vector_get(bck_vals, 2) / ((float)npixels.x * (float)npixels.y), 100.0 *gsl_vector_get(bck_vals, 3) / ((float)npixels.x * (float)npixels.y));
  fprintf(stdout, "\nNumber of iterations: % 2i\n\n", (int)gsl_vector_get(bck_vals, 4));

  if (scale_to_master)
    bck_img = compute_scale_grism(gr_img, gr_dqval, gr_mask, conf->dqmask, bck_vals);
  else
    bck_img = compute_scale_master(sc_img, bck_vals);

  // write the background image to a file
  fprintf(stdout,"Writing data to: %s...", bck_image);
  gsl_to_FITSimage (bck_img, bck_image, 1, ID);
  cards = fit_to_FITScards(bck_vals, npixels);
  put_FITS_cards(bck_image, 1, cards);
  fprintf(stdout,". Done.\n");

  // release memory
  free_fitbck_data(fbck_data);
  free_FITScards(cards);
  gsl_matrix_free(sc_img);
  gsl_matrix_free(gr_img);
  gsl_matrix_free(gr_dqval);
  gsl_matrix_free(gr_mask);
  gsl_matrix_free(bck_img);
  gsl_vector_free(bck_vals);
  free_aperture_conf(conf);
}

/**
 * Function: print_plis
 * Prints the content of a pixel list structure
 * to an ASCII file.
 *
 * Parameters:
 * @param fbck_data  - the list of pixel values
 * @param plist_name - full pathname to the pixel list
 *
 * Returns:
 * @return -
 */
void
print_plis(const fitbck_data *fbck_data, const char plist_name[])
{
    int index=0;
    char Buffer[10240];
    FILE *fout;

    // open the pixel list
    fout = fopen(plist_name, "w");

    // go over all data values
    for (index=0; index < fbck_data->n_data; index++)
      {
      // put the values into the buffer
      sprintf (Buffer, "%i %i %e %e\n",fbck_data->x_pos[index], fbck_data->y_pos[index],
          fbck_data->x_values[index], fbck_data->y_values[index]);

      // push the buffer to the  file
      fputs (Buffer, fout);
      }

    // close the pixel file
    fclose(fout);
}


/**
 * Function: compute_scale_grism
 * Computes a scaled version of the input grism images, using the
 * scale given as input. Pixels masked in other images are given
 * a fixed value.
 *
 * Parameters:
 * @param gr_img  - the list of pixel values
 * @param gr_dqval - full pathname to the pixel list
 * @param gr_mask - full pathname to the pixel list
 * @param bck_vals - full pathname to the pixel list
 *
 * Returns:
 * @return bck_img - the scaled grism image
 */
gsl_matrix *
compute_scale_grism(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
		int dqmask, gsl_vector *bck_vals)
{
    gsl_matrix *bck_img;

    int ii, jj;
    double scale=0.0;

    // get the scale
    scale = gsl_vector_get(bck_vals, 0);

    // allocate the space for the background image
    bck_img = gsl_matrix_alloc(gr_img->size1, gr_img->size2);

    // go over each row
    for (ii=0; ii < (int)gr_img->size1; ii++)
      // go over each column
      for (jj=0; jj < (int)gr_img->size2; jj++)
        // compute and set the value in the background image
        if (gsl_matrix_get(gr_mask, ii, jj) > 0.0 ||
            (int)gsl_matrix_get(gr_dqval, ii, jj) & dqmask)
          gsl_matrix_set(bck_img, ii, jj, MASK_VALUE);
        else
          // using the fit values only
          gsl_matrix_set(bck_img, ii, jj, gsl_matrix_get(gr_img, ii, jj) * scale);

    // return the background image
    return bck_img;
}

/**
 * Function: compute_scale_master
 * Computes a scaled version of the scaling (=master sky) image,
 * using the scale given as input.
 *
 * Parameters:
 * @param sc_img  - the list of pixel valuessc_img
 * @param bck_vals - full pathname to the pixel list
 *
 * Returns:
 * @return bck_img - the scaled grism image
 */
gsl_matrix *
compute_scale_master(const gsl_matrix *sc_img, const gsl_vector *bck_vals)
{
    gsl_matrix *bck_img;

    int ii, jj;
    double scale=0.0;

    // get the scale
    scale = gsl_vector_get(bck_vals, 0);

    // allocate the space for the background image
    bck_img = gsl_matrix_alloc(sc_img->size1, sc_img->size2);

    // go over each row
    for (ii=0; ii < (int)sc_img->size1; ii++)
      // go over each column
      for (jj=0; jj < (int)sc_img->size2; jj++)
        // scale the master sky
        gsl_matrix_set(bck_img, ii, jj, gsl_matrix_get(sc_img, ii, jj) * scale);

    // return the background image
    return bck_img;
}

/**
 * Function: fill_mask_data
 * Fills the pixel value structure with the data for the x/y-positions,
 * the grism image and the scaling image values and a weigth. dq-flagged
 * pixels and masked pixel are neglected. The mask is supposed to be a background
 * mask where pixel with a value < 900000 shall be neglected.
 * The weights are initially all set to 1.
 *
 * Parameters:
 * @param gr_img    - the grism image array
 * @param gr_dqval  - the dq-value array
 * @param gr_mask   - the grism mask array
 * @param sc_img    - the scaling image array
 * @param fbck_data - the pixel list structure
 * @param dqmask    - dq value from the configuration file
 *
 * Returns:
 * @return -
 */
void
fill_mask_data(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
		gsl_matrix *sc_img, fitbck_data *fbck_data, int dqmask)
{
    int ix=0;
    int iy=0;
    int index=0;

    // intitialize the index for
    // the fit_data structure
    index = 0;

    // go over all pixels
    // in the grism image
    for (ix=0; ix < (int)sc_img->size1; ix++)
      {
      for (iy=0; iy < (int)sc_img->size2; iy++)
        {
        // set the scale image pixel as independent variable
        fbck_data->x_values[index] = gsl_matrix_get(gr_img, ix, iy);

        // set the grism image pixel as dependent variable
        fbck_data->y_values[index] = gsl_matrix_get(sc_img, ix, iy);;

        // set the x- and y-positions
        fbck_data->x_pos[index] = ix;
        fbck_data->y_pos[index] = iy;

        // check whether the pixel was masked out
        // or is part of an object
        if (gsl_matrix_get(sc_img, ix, iy) < 0.0 ||
            gsl_matrix_get(gr_mask, ix, iy) < -900000.0 ||
            ((int)gsl_matrix_get(gr_dqval, ix, iy) & dqmask))
          {
          continue;
          }
        else
          {
          // give the pixel full weight
          //x,y,value_back,value_grism_imag
          fbck_data->e_values[index] = 1.0;

          // enhance the counter
          index++;
          }
        }
      }

    fbck_data->n_data = index;
    fprintf(stdout, "\nNumber of pixels in structure: %i\n\n", fbck_data->n_data);

    if (fbck_data->n_data < 1)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
          "aXe_SCALEBCK: no data left to determine the scale!\n");
}

/**
 * Function: fill_cont_data
 * Fills the pixel value structure with the data for the x/y-positions,
 * the grism image and the scaling image values and a weigth. dq-flagged
 * pixels and masked pixel are neglected. The mask is supposed to be a geometric
 * contamination where pixel with a value > 0 shall be neglected.
 * The weights are initially all set to 1.
 *
 * Parameters:
 * @param gr_img    - the grism image array
 * @param gr_dqval  - the dq-value array
 * @param gr_mask   - the grism mask array
 * @param sc_img    - the scaling image array
 * @param fbck_data - the pixel list structure
 *
 * Returns:
 * @return -
 */
void
fill_cont_data(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
                gsl_matrix *sc_img, fitbck_data *fbck_data, int dqmask)
{
    int ix=0;
    int iy=0;
    int index=0;

    // intitialize the index for
    // the fit_data structure
    index = 0;

    // go over all pixels
    // in the grism image
    for (ix=0; ix < (int)sc_img->size1; ix++)
      {
      for (iy=0; iy < (int)sc_img->size2; iy++)
        {
        // set the scale image pixel as independent variable
        fbck_data->x_values[index] = gsl_matrix_get(sc_img, ix, iy);

        // set the grism image pixel as dependent variable
        fbck_data->y_values[index] = gsl_matrix_get(gr_img, ix, iy);;

        // set the x- and y-positions
        fbck_data->x_pos[index] = ix;
        fbck_data->y_pos[index] = iy;

        // check whether the pixel was masked out
        // or is part of an object
        if (gsl_matrix_get(sc_img, ix, iy) < 0.0 ||
            gsl_matrix_get(gr_mask, ix, iy) > 0.0 ||
            ((int)gsl_matrix_get(gr_dqval, ix, iy) & dqmask))
          continue;
        else
          {
          // give the pixel full weight
          //x,y,value_back,value_grism_imag
          fbck_data->e_values[index] = 1.0;

          // enhance the counter
          index++;
          }
        }
      }

    fbck_data->n_data = index;
    fprintf(stdout, "\nNumber of pixels in structure: %i\n\n", fbck_data->n_data);

    if (fbck_data->n_data < 1)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
          "aXe_SCALEBCK: no data left to determine the scale!\n");
}

/**
 * Function: make_ksig_scalefit
 * Determines the scale from the values in the pixel value structure.
 * The results are determined iteratively while rejecting pixels via
 * kappa-sigma clipping to get robust measurements. The results are
 * the scale value, the standard deviation, the final number of
 * not rejected points and the number of iterations.
 *
 * Parameters:
 * @param fbck_data - the pixel value structure
 *
 * Returns:
 * @return bck_vals - a vector with the fit results
 */
gsl_vector *
make_ksig_scalefit(fitbck_data *fbck_data)
{
  gsl_vector *bck_vals=NULL;

  float old_scale;
  int index=0;
  int clipped=1;

  bck_vals = gsl_vector_alloc(6);

  // check whether iterations still
  // must be done and whether
  // the data was changed from clipping
  old_scale= 1.0e+06;
  while (index < N_BCKSCALE_ITER && clipped)
    {

    // make a non-weighted linear fit
    get_bck_scale(fbck_data->x_values, fbck_data->y_values, fbck_data->e_values,
        fbck_data->n_data, bck_vals);

    // make a clipping iteration
    clipped = clipp_scale_data(fbck_data, bck_vals, N_BCKSCALE_KAPPA);

    // break the iteration if the scale change is below a certain threshold
    if (fabs((gsl_vector_get(bck_vals, 0)-old_scale)/gsl_vector_get(bck_vals, 0)) < N_BCKSCALE_ACCUR)
      clipped = 0;

    // store the new scale
    old_scale = gsl_vector_get(bck_vals, 0);

    // inhance the clipping counter
    index++;
    }

  // set the number of iterations
  gsl_vector_set(bck_vals, 4, (double)index);

  // return the fit result
  return bck_vals;
}

/**
 * Function: get_bck_scale
 * Computes the scale value from the pixel values in two arrays.
 * The scale value is the median value of all input with weight.
 * Also the standard deviation, the total number of input and
 * the number of input with weight is returned.
 *
 * Parameters:
 * @param xs       - pixel values from one image
 * @param ys       - pixel values from second image
 * @param ws       - weight values
 * @param n_elem   - number of elements
 * @param bck_vals - vector for the result
 *
 * Returns:
 * @return -
 */
void
get_bck_scale(const double *xs, double *ys, double *ws,
		  const int n_elem, gsl_vector *bck_vals)
{
  int i, m;
  double *tmp;
  double median;
  double stdev;

  // allocate space for temporary vectors
  tmp = (double *) malloc (n_elem * sizeof (double));

  // initialize the array
  for (i = 0; i < n_elem; i++)
    tmp[i] = 0.0;

  // fill the temporary vectors
  // with scale values and weights
  m = 0;
  for (i = 0; i < n_elem; i++)
  {
    if (ws[i] > 0.0 && ys[i] != 0.0)
      {
      tmp[m] = xs[i] / ys[i];
      m++;
      }
  }

  // sort the vector
  gsl_sort( tmp, 1, m);

  // just confirm the sorting
  for (i = 1; i < m; i++)
    if (tmp[i] < tmp[i-1])
      fprintf(stdout, "Wrong: %f <--> %f\n", tmp[i],tmp[i-1]);

  // take the home made median
  median = tmp[(int)m/2];

  // get the standard deviation
  stdev = comp_stdev_from_array(tmp, m, median);

  // put the results in the vector
  gsl_vector_set(bck_vals, 0, median);
  gsl_vector_set(bck_vals, 1, stdev);
  gsl_vector_set(bck_vals, 2, (double)n_elem);
  gsl_vector_set(bck_vals, 3, (double)m);

  // release memory
  free(tmp);
}

/**
 * Function: clipp_scale_data
 * Set the weight in the pixel value structure such as to clip
 * point that deviate from the mean/characteristic value by
 * more than the allowed amount.
 *
 * Parameters:
 * @param fbck_data - pixel values from one image
 * @param bck_vals  - vector for the result
 * @param kappa     - vector for the result
 *
 * Returns:
 * @return nclip - integer/boolean indicating new clipping
 */
int
clipp_scale_data(fitbck_data *fbck_data, const gsl_vector *bck_vals,
		const float kappa)
{
    int index;
    int nclip=0;
    int newclip=0;
    double mean;
    //double stdev;
    double absdev;

    // get the mean value
    mean = gsl_vector_get(bck_vals, 0);

    // compute the maximum allowed deviation
    absdev = kappa * gsl_vector_get(bck_vals, 1);

    // go over all data values
    for (index=0; index < fbck_data->n_data; index++)
      {
      // avoid zero values
      if (fbck_data->y_values[index] != 0.0)
        {
        // check whether the actual value is outside the allowed range
        if (fabs(mean - (fbck_data->x_values[index]/fbck_data->y_values[index])) > absdev)
          {
          // check and mark clip changes
          if (fbck_data->e_values[index])
            newclip = 1;

          // set the weight, change the counter
          fbck_data->e_values[index] = 0.0;
          nclip++;
          }
        else
          {
          // set the weight
          fbck_data->e_values[index] = 1.0;
          }
        }
      else
        {
        // set the weight, change the counter
        fbck_data->e_values[index] = 0.0;
        nclip++;
        }
      }

    // return the newclip indicator
    return newclip;
}
