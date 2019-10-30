/**
 */
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sort_double.h>
#include "inout_aper.h"
#include "aXe_grism.h"
#include "spce_sect.h"
#include "spc_back.h"
#include "spce_PET.h"
#include "spc_wl_calib.h"
#include "aXe_errors.h"
#include "fringe_conf.h"
#include "spc_resp.h"
#include "spce_pathlength.h"
#include "aper_conf.h"
#include "specmodel_utils.h"
#include "model_utils.h"
#include "spc_fluxcube.h"
#include "fringe_conf.h"
#include "dirimage_model.h"


#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define SQR(x) ((x)*(x))

int
compute_dirimage_model(char dirim_file[], char conf_file[], char tpass_file[],
		       char specmod_file[], char objmod_file[], char aper_file[],
		       const double model_scale, const double tel_area, const double lambda_psf,
		       observation *obs, char map_file[])
{
  object          **oblist;
  dirobject       **dirlist;
  spectral_models *spec_mod;
  object_models   *obj_mod;

  interpolator    *tpass;

  gsl_matrix *dirimage_matrix=NULL;

  px_point npixels;

  // load the object list
  fprintf (stdout, "aXe_DIRIMAGE: Loading object aperture list...");
  oblist = file_to_object_list_seq (aper_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // check whether highres models
  // are given
  if (strlen(specmod_file) > 0)
    // load the spectral models
    spec_mod = load_spectral_models(specmod_file);
  else
    // or set the struct to NULL
    spec_mod = NULL;

  // check whether direct emission models
  // are given
  if (strlen(objmod_file) > 0)
    obj_mod = load_object_models(objmod_file);
  else
    obj_mod = NULL;

  // get the sensitivity curve of the total passband
  tpass =  get_filter_sensitivity(tpass_file, tel_area);

  // get the image dimensions
  npixels = get_npixel(obs);

  // create the list of direct objects to be simulated
  // the interpolation type is fixed to linear
  dirlist = oblist_to_dirlist2(dirim_file, conf_file, npixels, oblist,
			       spec_mod, obj_mod, model_scale, 1);

  // determine the XOFF and YOFF values
  // for the various beams
  fill_xy_offsets(dirlist, conf_file);

  // create and fill the pixel matrix for the direct images
  dirimage_matrix = make_dirimage(oblist, dirlist, npixels, lambda_psf, tpass);

  // store the contamination image
  gsl_to_FITSimage (dirimage_matrix, map_file, 1, NULL);

  // check if memory must be released;
  // do it if necessary
  if (spec_mod != NULL)
    free_spectral_models(spec_mod);
  if (obj_mod != NULL)
    free_object_models(obj_mod);
  if (oblist !=NULL)
    free_oblist (oblist);

  // always release these
  // objects
  free_dirlist(dirlist);
  gsl_matrix_free(dirimage_matrix);
  free_interp(tpass);

  // return always '1'
  return 1;
}


gsl_matrix *
make_dirimage(object **oblist, dirobject **dirlist, const px_point npixels,
	      const double lambda_psf, interpolator *tpass)
{
  dirobject       *actdir;
  beam             actbeam;

  gsl_matrix      *dirimage_matrix;

  double cps        = 0.0;
  double sval       = 0.0;
  double value      = 0.0;

  int ii=0;
  int nx=0;
  int ny=0;

  d_point dpixel;

   // allocate memory for the image matrix
  dirimage_matrix = gsl_matrix_alloc(npixels.x, npixels.y);
  gsl_matrix_set_all(dirimage_matrix,0.0);

  // go over each beam model
  ii = 0;
  while (oblist[ii] != NULL)
    {

      // get the direct object for the actual model spectrum
      actdir = get_dirobject_from_list(dirlist, oblist[ii]->ID);

      // get the beam for the actual model spectrum
      actbeam = oblist[ii]->beams[0];

      // report onto the screen
      fprintf(stdout, "aXe_DIRIMAGE: modelling object %i ...", oblist[ii]->ID);

      // get the integrated direct image intensity in cps
      cps = get_cps_for_dirobject(tpass, actdir);

      // correct the direct image positions for
      // the offset values introduced by building
      // the direct image objects around the reference
      // position, which is shifted from the true
      // direct image position by XOFF and YOFF as
      // given in the aXe configuration file
      actdir->ix_min -= (int)floor(actdir->xy_off[actbeam.ID].x + 0.5);
      actdir->ix_max -= (int)floor(actdir->xy_off[actbeam.ID].x + 0.5);
      actdir->iy_min -= (int)floor(actdir->xy_off[actbeam.ID].y + 0.5);
      actdir->iy_max -= (int)floor(actdir->xy_off[actbeam.ID].y + 0.5);
      actbeam.refpoint.x -= actdir->xy_off[actbeam.ID].x;
      actbeam.refpoint.y -= actdir->xy_off[actbeam.ID].y;

      // go over each pixel in the direct object area
      for (nx=actdir->ix_min; nx<=actdir->ix_max; nx++)
	{
	  // make sure to be inside the image
	  if (nx < 0 || nx >=  (int)dirimage_matrix->size1)
	    // skip column if not
	    continue;
	  for (ny=actdir->iy_min; ny<=actdir->iy_max; ny++)
	    {
	      // make sure to be inside the image
	      if (ny < 0 || ny >=  (int)dirimage_matrix->size2)
		// skip element if not
		continue;

	      // fill the dpixel structure
	      dpixel.x = (double)nx;
	      dpixel.y = (double)ny;

	      // check whether there is
	      // an image template
	      if (actdir->dirim)
		// get the pixel intensity from the image template
		sval = get_diremission_value(actdir->dirim, dpixel.x - actbeam.refpoint.x, dpixel.y - actbeam.refpoint.y);
	      else
		// get the Gaussian pixel intensity;
		//
		// do a subsampling over the pixel
		// to get a more appropriate value for the
		// emission val
		sval = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);

	      // compute and set the new pixel values
	      value = gsl_matrix_get(dirimage_matrix, nx, ny) + sval*cps;
	      gsl_matrix_set(dirimage_matrix, nx, ny, value);
	    }
	}

      // release the space for the various structures
      fprintf(stdout, " Done\n");
      ii++;
    }

  // return the image matrix
  return dirimage_matrix;
}



double
get_cps_for_dirobject(interpolator *tpass, dirobject *actdir)
{
  interpolator *combine;

  double cps=0.0;

  // create aa new multiplicator
  // interpolator
  combine = combine_tpass_SED(tpass, actdir);

  // integrate over the interpolator
  cps = integrate_interpolator(combine);

  // release memory
  free_interp(combine);

  // return the cps value
  return cps;
}


double
integrate_interpolator(interpolator *combine)
{
  double integral=0.0;

  int index;

  for (index=1; index < combine->nvals-1; index++)
    // add the increment
    integral +=  combine->yvals[index] * (combine->xvals[index+1] - combine->xvals[index-1]);

  // dont forget the start and end piece
  integral +=  combine->yvals[0]                * (combine->xvals[1]                - combine->xvals[0]);
  integral +=  combine->yvals[combine->nvals-1] * (combine->xvals[combine->nvals-1] - combine->xvals[combine->nvals-2]);

  // then return the half
  return integral/2.0;
}

interpolator *
combine_tpass_SED(interpolator *tpass, dirobject *actdir)
{
  interpolator *combine;

  int index=0;
  int n_additional=0;
  int act_index;

  gsl_vector     *indep_data;

  // go over all SED data points
  for (index = 0; index < actdir->SED->npoints; index++)
    {
      // transform to Angstrom
      actdir->SED->wavelength[index] *= 10.0;

      // check whether the independent data point
      // falls in the area of the sensitivity
      if (actdir->SED->wavelength[index] < tpass->xmax
	  && actdir->SED->wavelength[index] >  tpass->xmin)
	// enhance the number
	// of additional data points
	n_additional += 1;
    }

  // allocate memory for the new independent value list
  indep_data   = gsl_vector_alloc(tpass->nvals + n_additional);

  // go over all sensitivity data points and fill in
  // the independent values
  for (index = 0; index < tpass->nvals; index++)
    gsl_vector_set(indep_data, index, tpass->xvals[index]);

  // set the counter to the next free index
  act_index = tpass->nvals;

  // go over all SED data points
  for (index = 0; index < actdir->SED->npoints; index++)
    // check whether the independent data point
    // falls in the area of the sensitivity
    if (actdir->SED->wavelength[index] < tpass->xmax
	&& actdir->SED->wavelength[index] >  tpass->xmin)
      {
	// fill in the additional independent data point
	// from the SED
	gsl_vector_set(indep_data, act_index, actdir->SED->wavelength[index]);

	// enhance the counter
        act_index += 1;
      }

  // sort the vector
  gsl_sort_vector(indep_data);


  // compose a new interpolator from the independent values
  combine = get_combined_tpass_SED(indep_data, tpass, actdir);

  // release memory in the weight vector
  gsl_vector_free(indep_data);

  // go over all SED data points
  for (index = 0; index < actdir->SED->npoints; index++)
    // transform to nm
    actdir->SED->wavelength[index] /= 10.0;

  // return the interpolator
  return combine;
}

interpolator *
get_combined_tpass_SED(gsl_vector *indep_data, interpolator *tpass, dirobject *actdir)
{
  interpolator *combine;

  gsl_vector_int *indep_weight;

  int n_new=0;
  int index;
  int act_index;

  double *xvals_new;
  double *yvals_new;

  // allocat memory for the weight vector
  indep_weight = gsl_vector_int_alloc(indep_data->size);

  // set all the weights
  gsl_vector_int_set_all(indep_weight, 1);

  // go over the independent data
  for (index = 1; index < (int)indep_data->size; index++)
    // check whether the current entry is equal the previous one
    if (gsl_vector_get(indep_data, index) ==  gsl_vector_get(indep_data, index-1))
      // set the weight of the current entry to zero
      gsl_vector_int_set(indep_weight,index, 0);


  // initialize the
  // final number
  n_new = 0;

  // go over the weight array
  for (index=0; index < (int)indep_weight->size; index++)
    // just count the weights
    n_new += gsl_vector_int_get(indep_weight, index);

  // allocate memory for the interpolator data
  xvals_new = (double *)malloc(n_new * sizeof(double));
  yvals_new = (double *)malloc(n_new * sizeof(double));

  // go over the weight array
  act_index=0;
  for (index=0; index < (int)indep_weight->size; index++)
    // check for weight
    if (gsl_vector_int_get(indep_weight, index))
      {
	// take the independent value from the vector
	xvals_new[act_index] = gsl_vector_get(indep_data, index);

	// compute the new dependent value bye multiplying
	// the brightness with the sensitivity
	yvals_new[act_index] = eval_interp(tpass, xvals_new[act_index])
	  * get_flux_from_SED(actdir->SED, xvals_new[act_index]);

	// enhance the counter
	act_index++;
      }

  // create a new interpolator
  combine = create_interp(n_new, TPASS_INTERP_TYPE, xvals_new, yvals_new);


  // release memory in the weight vector
  gsl_vector_int_free(indep_weight);

  // return the interpolator
  return combine;
}


/*
 * Function: get_filter_sensitivity
 * The function reads in a total passband file from fits format
 * and stroes the data as an interpolation function. The dependent
 * data values are transformed from throughput to sensitivity,
 * and the interpolation function is returned.
 *
 * Parameters:
 * @param tpass_file  - the full name of the total passband
 * @param tel_are     - the collecting area of the telescope
 *
 * Returns:
 * @return tpass    - the filter sensitivity
 */
interpolator *
get_filter_sensitivity(const char tpass_file[], const double tel_area)
{
  interpolator    *tpass;

  double h_erg  = 6.6260693E-27;
  double c_cm   = 2.99792458E+10;

  double factor = 0.0;

  int index=0;

  fprintf(stdout, "Load Total bandpass table :%s\n", tpass_file);
  tpass =  create_interp_ftable(tpass_file, 2,"WAVELENGTH","THROUGHPUT", TPASS_INTERP_TYPE);

  for (index=0; index < tpass->nvals; index++)
    {
      // compute the conversion factor
      factor = tel_area / (h_erg * c_cm / (1.0E-08 * tpass->xvals[index]));

      // apply the conversion factor
      tpass->yvals[index] = tpass->yvals[index] * factor;
    }

  // return the filter sensitivity
  return tpass;
}
