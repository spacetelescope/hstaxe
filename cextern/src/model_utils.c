/**
 * File: model_utils.c
 * Subroutines to calculate the
 * various contamination models
 */
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

#include "model_utils.h"
#include "fitsio.h"
#include "inout_aper.h"
#include "aXe_grism.h"
#include "spce_PET.h"
#include "spc_wl_calib.h"
#include "aXe_errors.h"
#include "fringe_conf.h"
#include "spc_resp.h"
#include "spce_pathlength.h"
#include "aper_conf.h"
#include "crossdisp_utils.h"
#include "spc_FITScards.h"


#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define SQR(x) ((x)*(x))


dirim_emission *
model_gauss_dirim(dirobject *actdir, beam actbeam, aperture_conf *conf, double psf_offset)
{
  int nx, ny;
  int ii, jj;

  double sval=0.0;

  d_point dpixel;

  dirim_emission *gauss_dirim = NULL;

  // check whether something can be done
  if (actdir->dirim || (conf->psfcoeffs && conf->psfrange) || psf_offset)
	  // return NULL if not
	  return NULL;

  // dimension the new matrix;
  // leave one border row more
  nx = actdir->ix_max - actdir->ix_min + 3;
  ny = actdir->iy_max - actdir->iy_min + 3;

  // allocate memory for the structure
  gauss_dirim = (dirim_emission *) malloc(sizeof(dirim_emission));

  // load the image in the gsl
  gauss_dirim->modimage = gsl_matrix_alloc(nx, ny);

  // transfer the image dimension
  gauss_dirim->dim_x = (int)gauss_dirim->modimage->size1;
  gauss_dirim->dim_y = (int)gauss_dirim->modimage->size2;

  // compute and store the mean image coordinates
  gauss_dirim->xmean = (float)(gauss_dirim->dim_x-1) / 2.0;
  gauss_dirim->ymean = (float)(gauss_dirim->dim_y-1) / 2.0;

  // go over all pixels in the area
  for (ii=0; ii < gauss_dirim->dim_x; ii++)
    {
	  for (jj=0; jj < gauss_dirim->dim_y; jj++)
	    {
		  // fill the dpixel structure with the position
		  // RELATIVE to the reverence position of the beam
		  dpixel.x = (double)ii - gauss_dirim->xmean + actbeam.refpoint.x;
		  dpixel.y = (double)jj - gauss_dirim->ymean + actbeam.refpoint.y;

		  // do a subsampling over the pixel
          // to get a more appropriate value for the
          // emission value
          sval = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);

          // set the emission value in the matrix
          gsl_matrix_set(gauss_dirim->modimage, ii, jj, sval);

	    }
    }

  // return the matrix
  return gauss_dirim;
}

/**
 * Function: get_calib_function
 * The function extracts the necessary data for the wavelength
 * calibration from the configuration file. The wavelength
 * calibration is assembled for a particular beam of a particular
 * object.
 *
 * Parameters:
 * @param actspec   - the model spectrum this is done for
 * @param actdir    - the dirobject this is done for
 * @param CONF_file - the full filename of the configuration file
 * @param conf      - the configuration structure
 *
 * Returns:
 * @return wl_calibration - the wavelength calibration
 */
calib_function *
get_calib_function(beamspec *actspec, dirobject *actdir, char CONF_file[],
                   const aperture_conf * conf)
{
  calib_function *wl_calibration;
  dispstruct     *disp;
  d_point pixel;
  int for_grism=0;

  // get the reference point right
  pixel.x = actdir->refpoint.x - conf->refx;
  pixel.y = actdir->refpoint.y - conf->refy;

  // look whether we are for grisms or prisms
  for_grism = check_for_grism (CONF_file,actspec->beamID);

  // get the dispersion structure
  disp = get_dispstruct_at_pos(CONF_file, for_grism,
                               actspec->beamID,pixel);


  // transform the dispersion structure into the
  // wavelength calibration
  wl_calibration = create_calib_from_gsl_vector(for_grism, disp->pol);
  if (!for_grism)
    wl_calibration->pr_range = get_prange (CONF_file, actspec->beamID);

  free_dispstruct(disp);

  // return the wavelength calibration
  return wl_calibration;
}

/**
 * Function: get_throughput_spec
 * The function determines the throughput file for a certain
 * model spectrum. The filename is extracted  from the configuration
 * file. Then the throughput file is loaded into a spectrum structure.
 * The spectrum structure is returned.
 *
 * Parameters:
 * @param actspec   - the model spectrum this is done for
 * @param CONF_file - the full filename of the configuration file
 *
 * Returns:
 * @return resp     - the response function as a spectrum structure
 */
spectrum *
get_throughput_spec(beamspec *actspec, char CONF_file[])
{
  spectrum *resp;
  char through_file[MAXCHAR];
  char through_file_path[MAXCHAR];

  // determine the filename of the throughput file from the configuration
  get_troughput_table_name(CONF_file, actspec->beamID, through_file);

  // build up the full filename
  build_path (AXE_CONFIG_PATH, through_file, through_file_path);

  // load the throughput into the spectrum struct
  resp=get_response_function_from_FITS(through_file_path,2);

  // return the spectrum struct
  return resp;
}

/**
 * Function: compute_tracedata
 * The function creates a tracedata structure for a specific
 * beam model. A tracedata structure consists of all relevant
 * information (position, wavelength, dispersion) for the pixels
 * along the trace of a specific beam. Since the reference point
 * remains constant during the modelling of a beam, it is
 * better to store the relevant data for the use in each pixel.
 *
 * Parameters:
 * @param  actbeam        - the beam of the model spectrum
 * @param  actdir         - the direct object of the model spectrum
 * @param  wl_calibration - the wavelength calibration of the model spectrum
 * @param  actspec        - the model spectrum
 *
 * Returns:
 * @return acttrace       - the tracedata structure
 */
tracedata *
compute_tracedata(const beam actbeam, const dirobject *actdir,
                  const calib_function *wl_calibration,
                  const beamspec *actspec)
{

  tracedata  *acttrace;
  trace_func *tracefun;

  gsl_vector *dx;
  gsl_vector *dy;
  gsl_vector *xi;
  gsl_vector *lambda;
  gsl_vector *dlambda;
  gsl_vector *flux;
  gsl_vector *gvalue;

  double tmp1, tmp2;

  int dx_min, dx_max;
  int npoints;
  int i;
  int nentries;

  // allocate space for the tracedata
  acttrace = (tracedata *)malloc(sizeof(tracedata));
  if (acttrace == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "compute_tracedata:" " Could not allocate"
                 " memory for a tracedata object!");

  tracefun =  actbeam.spec_trace;

  // use the x-range of the direct obect area to
  // compute the dx-range for this beam
  dx_min = MIN(actspec->model_ref.x - actdir->ix_min, actspec->model_ref.x - actdir->ix_max)-5;
  dx_max = MAX(actspec->model_ref.x - actdir->ix_min, actspec->model_ref.x - actdir->ix_max)
    + (int)actspec->model->size1+5;

  // compute the number of points in the dx-range
  npoints = dx_max - dx_min + 1;

  // allocate space for the tracedata
  dx      = gsl_vector_alloc(npoints);
  dy      = gsl_vector_alloc(npoints);
  xi      = gsl_vector_alloc(npoints);
  lambda  = gsl_vector_alloc(npoints);
  dlambda = gsl_vector_alloc(npoints);
  flux    = gsl_vector_alloc(npoints);
  gvalue  = gsl_vector_alloc(npoints);


  // go over each point in the dx-range
  for (i=0; i < npoints; i++)
    {

      // compute the dx- and dy-values
      tmp1 = (double)(dx_min + i);
      tmp2 = tracefun->func (tmp1, tracefun->data);

      // store dx and dy, and xi (=dx)
      gsl_vector_set(dx, i, tmp1);
      gsl_vector_set(dy, i, tmp2);
      gsl_vector_set(xi, i, tmp1);
      gsl_vector_set(gvalue, i, 1.0);
      gsl_vector_set(flux, i, 1.0);
    }

  // compute the true xi values
  abscissa_to_pathlength (tracefun, xi);

  // go over each point in the dx-range
  for (i=0; i < npoints; i++)
    {

      // compute  and store lambda
      gsl_vector_set(lambda, i, wl_calibration->func(gsl_vector_get(xi, i), wl_calibration->order, wl_calibration->coeffs));

      // compute  and store dlambda
      tmp1 = wl_calibration->func(gsl_vector_get(xi, i)-0.5, wl_calibration->order, wl_calibration->coeffs);
      tmp2 = wl_calibration->func(gsl_vector_get(xi, i)+0.5, wl_calibration->order, wl_calibration->coeffs);
      gsl_vector_set(dlambda, i, fabs(tmp2-tmp1));
    }

  for (i=0; i < npoints; i++)
    {
      if (i == 0)
        {
          tmp1 = fabs(gsl_vector_get(lambda,i+1)-gsl_vector_get(lambda,i));
        }
      else if (i == npoints-1)
        {
          tmp1 = fabs(gsl_vector_get(lambda,i)-gsl_vector_get(lambda,i-1));
        }
      else
        {
          tmp1 = fabs(gsl_vector_get(lambda,i+1) - gsl_vector_get(lambda,i-1)) / 2.0;
        }
      gsl_vector_set(dlambda, i, tmp1);
    }

  // transfer the quantities to the structure
  acttrace->npoints = npoints;
  acttrace->dx_start= dx_min;

  acttrace->dx      = dx;
  acttrace->dy      = dy;
  acttrace->xi      = xi;
  acttrace->lambda  = lambda;
  acttrace->dlambda = dlambda;
  acttrace->flux    = flux;
  acttrace->gvalue  = gvalue;

  // for prism data: constrain the tracedata
  if (wl_calibration->pr_range != NULL)
    {
     nentries = get_valid_tracedata(acttrace, wl_calibration);
     select_tracedata(acttrace, wl_calibration,nentries);
    }
  // return the tracestructure
  return acttrace;
}
/**
 * Function: compute_short_tracedata
 *
 * Parameters:
 * @param  actbeam        - the beam of the model spectrum
 * @param  actdir         - the direct object of the model spectrum
 * @param  wl_calibration - the wavelength calibration of the model spectrum
 * @param  actspec        - the model spectrum
 *
 * Returns:
 * @return acttrace       - the tracedata structure
 */
tracedata *
compute_short_tracedata(const aperture_conf *conf, const beam actbeam,
		                const dirobject *actdir,
                        const calib_function *wl_calibration,
                        const beamspec *actspec)
{

  tracedata  *acttrace;
  trace_func *tracefun;

  gsl_vector *dx;
  gsl_vector *dy;
  gsl_vector *xi;
  gsl_vector *lambda;
  gsl_vector *dlambda;
  gsl_vector *flux;
  gsl_vector *gvalue;

  double tmp1, tmp2;

  int dx_min, dx_max;
  int npoints;
  int i;
  int nentries;
  int dx0, dx1;

  // allocate space for the tracedata
  acttrace = (tracedata *)malloc(sizeof(tracedata));
  if (acttrace == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "compute_tracedata:" " Could not allocate"
                 " memory for a tracedata object!");

  tracefun =  actbeam.spec_trace;

  dx0 = (double)conf->beam[actspec->beamID].offset.dx0;
  dx1 = (double)conf->beam[actspec->beamID].offset.dx1;
  dx_min = MIN(dx0, dx1)-1;
  dx_max = MAX(dx0, dx1)+1;

  // use the x-range of the direct obect area to
  // compute the dx-range for this beam
  //dx_min = MIN(actspec->model_ref.x - actdir->ix_min, actspec->model_ref.x - actdir->ix_max)-5;
  //dx_max = MAX(actspec->model_ref.x - actdir->ix_min, actspec->model_ref.x - actdir->ix_max)
  //  + (int)actspec->model->size1+5;

  // compute the number of points in the dx-range
  npoints = dx_max - dx_min + 1;

  // allocate space for the tracedata
  dx      = gsl_vector_alloc(npoints);
  dy      = gsl_vector_alloc(npoints);
  xi      = gsl_vector_alloc(npoints);
  lambda  = gsl_vector_alloc(npoints);
  dlambda = gsl_vector_alloc(npoints);
  flux    = gsl_vector_alloc(npoints);
  gvalue  = gsl_vector_alloc(npoints);


  // go over each point in the dx-range
  for (i=0; i < npoints; i++)
    {

      // compute the dx- and dy-values
      tmp1 = (double)(dx_min + i);
      tmp2 = tracefun->func (tmp1, tracefun->data);

      // store dx and dy, and xi (=dx)
      gsl_vector_set(dx, i, tmp1);
      gsl_vector_set(dy, i, tmp2);
      gsl_vector_set(xi, i, tmp1);
      gsl_vector_set(gvalue, i, 1.0);
      gsl_vector_set(flux, i, 1.0);
    }

  // compute the true xi values
  abscissa_to_pathlength (tracefun, xi);

  // go over each point in the dx-range
  for (i=0; i < npoints; i++)
    {

      // compute  and store lambda
      gsl_vector_set(lambda, i, wl_calibration->func(gsl_vector_get(xi, i), wl_calibration->order, wl_calibration->coeffs));

      // compute  and store dlambda
      tmp1 = wl_calibration->func(gsl_vector_get(xi, i)-0.5, wl_calibration->order, wl_calibration->coeffs);
      tmp2 = wl_calibration->func(gsl_vector_get(xi, i)+0.5, wl_calibration->order, wl_calibration->coeffs);
      gsl_vector_set(dlambda, i, fabs(tmp2-tmp1));
    }

  for (i=0; i < npoints; i++)
    {
      if (i == 0)
        {
          tmp1 = fabs(gsl_vector_get(lambda,i+1)-gsl_vector_get(lambda,i));
        }
      else if (i == npoints-1)
        {
          tmp1 = fabs(gsl_vector_get(lambda,i)-gsl_vector_get(lambda,i-1));
        }
      else
        {
          tmp1 = fabs(gsl_vector_get(lambda,i+1) - gsl_vector_get(lambda,i-1)) / 2.0;
        }
      gsl_vector_set(dlambda, i, tmp1);
    }

  // transfer the quantities to the structure
  acttrace->npoints = npoints;
  acttrace->dx_start= dx_min;

  acttrace->dx      = dx;
  acttrace->dy      = dy;
  acttrace->xi      = xi;
  acttrace->lambda  = lambda;
  acttrace->dlambda = dlambda;
  acttrace->flux    = flux;
  acttrace->gvalue  = gvalue;

  // for prism data: constrain the tracedata
  if (wl_calibration->pr_range != NULL)
    {
     nentries = get_valid_tracedata(acttrace, wl_calibration);
     select_tracedata(acttrace, wl_calibration,nentries);
    }
  // return the tracestructure
  return acttrace;
}

/*
 * Function: get_valid_tracedata
 * The derives the number of valid tracedata points for prism data.
 * For prism data, only tracepoints with a certain offset
 * from the singularity at tr_length = a_0 is accepted.
 *
 * Parameters:
 * @param acttrace       - the tracedata
 * @param wl_calibration - the calibration structure
 *
 * Returns:
 * @return nentries      - the number of valid tracedata points
 */
double
get_valid_tracedata(tracedata *acttrace, const calib_function *wl_calibration)
{

  double lower, upper, a_0, d_xi;

  int nentries=0;
  int i=0;

  // get the lower and upper boundaries of the accepted range
  lower = gsl_vector_get(wl_calibration->pr_range, 0);
  upper = gsl_vector_get(wl_calibration->pr_range, 1);

  // get the a0 coefficient
  a_0 = wl_calibration->coeffs[0];

  // go over all point int he tracedata
  for (i=0; i < acttrace->npoints; i++)
    {
      // compute the offset from the singularity
      d_xi = gsl_vector_get(acttrace->xi, i) -  a_0;

      // enhance a counter if the offset is within the accepted range
      if (d_xi >= lower && d_xi <= upper)
        nentries++;
    }

  // return the number of valid data
  return nentries;
}


/*
 * Function: select_tracedata
 * The function applies additional constraints for prism data.
 * In prisms not all trace positions can be accepted, since
 * the odd calibration function results into strange wavelengths
 * and values close and beyond the singularity.
 * The function creates a new tracedata structure with values
 * only from the accepted trace range.
 *
 * Parameters:
 * @param acttrace       - the tracedata
 * @param wl_calibration - the calibration structure
 * @param nentries       - the number of valid entries
 */
void select_tracedata(tracedata *acttrace,
                      const calib_function *wl_calibration, const int nentries)
{
  gsl_vector *dx;
  gsl_vector *dy;
  gsl_vector *xi;
  gsl_vector *lambda;
  gsl_vector *dlambda;
  gsl_vector *flux;
  gsl_vector *gvalue;

  double lower, upper, a_0, d_xi;

  int i=0, iact;

  // get the lower and upper boundaries of the accepted range
  lower = gsl_vector_get(wl_calibration->pr_range, 0);
  upper = gsl_vector_get(wl_calibration->pr_range, 1);

  // get the a0 coefficient
  a_0 = wl_calibration->coeffs[0];

  // check wehter thereis valid data at all
  if (nentries < 1)
    {
      // if no valid data, set all vectors to NULL
      dx      = NULL;
      dy      = NULL;
      xi      = NULL;
      lambda  = NULL;
      dlambda = NULL;
      flux    = NULL;
      gvalue    = NULL;

      acttrace->dx_start = 0.0;
    }
  else
    {
      // if there is valid data:
      // allocate space for the new tracedata structure
      dx      = gsl_vector_alloc(nentries);
      dy      = gsl_vector_alloc(nentries);
      xi      = gsl_vector_alloc(nentries);
      lambda  = gsl_vector_alloc(nentries);
      dlambda = gsl_vector_alloc(nentries);
      flux    = gsl_vector_alloc(nentries);
      gvalue  = gsl_vector_alloc(nentries);

      // go over all tracedata points
      iact=0;
      for (i=0; i < acttrace->npoints; i++)
        {

          // compute the offset from the singularity
          d_xi = gsl_vector_get(acttrace->xi, i) -  a_0;

          // if the offside is in the accepted range
          if (d_xi >= lower && d_xi <= upper)
            {

              // transfer the data from the old to the new structure
              gsl_vector_set(dx     , iact, gsl_vector_get(acttrace->dx, i));
              gsl_vector_set(dy     , iact, gsl_vector_get(acttrace->dy, i));
              gsl_vector_set(xi     , iact, gsl_vector_get(acttrace->xi, i));
              gsl_vector_set(lambda , iact, gsl_vector_get(acttrace->lambda, i));
              gsl_vector_set(dlambda, iact, gsl_vector_get(acttrace->dlambda, i));
              gsl_vector_set(gvalue , iact, 1.0);

              // enhance the counter
              iact++;
            }
        }

      // set the intitial value
      acttrace->dx_start = gsl_vector_get(dx, 0);
    }

  // set the number of points
  acttrace->npoints = nentries;

  // release the old vectors in the structure
  gsl_vector_free(acttrace->dx);
  gsl_vector_free(acttrace->dy);
  gsl_vector_free(acttrace->xi);
  gsl_vector_free(acttrace->lambda);
  gsl_vector_free(acttrace->dlambda);
  gsl_vector_free(acttrace->flux);
  gsl_vector_free(acttrace->gvalue);

  // fill the new vectors into the tracedata structure
  acttrace->dx      = dx;
  acttrace->dy      = dy;
  acttrace->xi      = xi;
  acttrace->lambda  = lambda;
  acttrace->dlambda = dlambda;
  acttrace->flux    = flux;
  acttrace->gvalue  = gvalue;
}

/*
 * Function: fill_fluxfrom_SED
 * This function fills flux values at different wavelengths
 * into the tracedata structure. The flux values are computed
 * via the SED in the direct objects.
 *
 * Parameters:
 * @param actdir   - the direct object
 * @param acttrace - the tracedata structure
 */
void
fill_fluxfrom_SED(const dirobject *actdir, tracedata *acttrace)
{
  int i=0;

  // go over all tracedata points
  for (i=0; i<acttrace->npoints; i++)
    // compute and fill in the flux values at each wavelength
    //gsl_vector_set(acttrace->flux, i, get_flux_from_SED(actdir->SED, gsl_vector_get(acttrace->lambda, i)/10.0));
    gsl_vector_set(acttrace->flux, i,
                   get_aveflux_from_SED(actdir->SED,
                                        gsl_vector_get(acttrace->lambda, i)/10.0, gsl_vector_get(acttrace->dlambda, i)/10.0));
}

/*
 * Function: print_tracedata
 * This function prints some tracedata values
 * onto the screen. Used mainly for debugging
 * reasons.
 *
 * Parameters:
 * @param acttrace - the tracedata structure
 */
void
print_tracedata(tracedata *acttrace)
{
  int i=0;

  fprintf(stdout, "Trace starting at: %f\n", acttrace->dx_start);
  for (i=0; i < acttrace->npoints; i++)
    {
      fprintf(stdout, "dx: %f, dy: %f, xi: %f, lambda: %f, dlambda: %f, gdata: %f, flux: %g\n",
              gsl_vector_get(acttrace->dx, i), gsl_vector_get(acttrace->dy, i),gsl_vector_get(acttrace->xi, i),
              gsl_vector_get(acttrace->lambda, i), gsl_vector_get(acttrace->dlambda, i), gsl_vector_get(acttrace->gvalue, i),
              gsl_vector_get(acttrace->flux, i));
    }
}

/**
 * Function: oblist_to_dirlist
 * Transforms a list of "object"'s into a list of
 * "dirobject"'s. The list of "dirobject"'s is terminated
 * with a NULL-object at the end.
 *
 * Parameters:
 * @param  grism_file  - the full path to the grism image
 * @param  CONF_file   - the full path to the config file
 * @param  npixels     - the dimension of the grism image
 * @param  oblist      - the object list to start from
 * @param  spec_mod    - the model spectra
 * @param  model_scale - the scale for the size of the direct object area
 * @param  int_type    - interpolation type
 *
 * Returns:
 * @return dirlist     - the dirobject list created
 */
dirobject **
oblist_to_dirlist(char grism_file[], char CONF_file[], const  px_point npixels,
                  object  **oblist, spectral_models *spec_mod,
                  const double model_scale, const int int_type)
{

  dirobject **dirlist;
  aperture_conf *conf;
  gsl_matrix *drzcoeffs;

  int nobjects=0;
  int i=0;
  int j=0;
  int beamID;
  int max_offs;

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  // load the extension numbers
  get_extension_numbers(grism_file, conf,conf->optkey1,conf->optval1);

  // get the  matrix with the drizzle coefficients
  drzcoeffs = get_crossdisp_matrix(grism_file, conf->science_numext);
  if (drzcoeffs->size1 < 2 || !drzcoeffs->size2)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "oblist_to_dirlist:" " Could not get"
                 " the drizzle coefficients in file: %s\n", grism_file);

  // determine an offset from the PSF_OFFSET
  max_offs = (int)ceil(get_max_offset(conf));

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // allocate space for the dirobject list
  dirlist = (dirobject **) malloc((nobjects+1) * sizeof(dirobject *));
  if (dirlist == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "oblist_to_dirlist:" " Could not allocate"
                 " memory for pointers to %i dirobject objects", nobjects+1);

  // loop over all objectsconfig->camera
  for (i = 0; i < nobjects; i++)
    {
      // make sure that there are beams in the object
      if (oblist[i]->nbeams > 0)
        {
          // create a dirobject for each object
          dirlist[j] = fill_dirobject(oblist[i], npixels, drzcoeffs, model_scale, max_offs);
          fill_spectrum(oblist[i], dirlist[j], spec_mod, int_type);
          //      fill_dpsf_function(conf, dirlist[j]);

          // fill the xoffset and yoffset values.
          // within the gaussian models they are dummys.
          for (beamID=0; beamID < conf->nbeams; beamID++)
            {
              dirlist[j]->xy_off[beamID].x = 0.0;
              dirlist[j]->xy_off[beamID].y = 0.0;
            }
          j++;
        }
    }

  // terminate the dirobject list with NULL
  dirlist[j] = NULL;

  // release the memory for the drizzle matrix
  gsl_matrix_free(drzcoeffs);

  free_aperture_conf(conf);

  return dirlist;
}

/**
 * Function: oblist_to_dirlist2
 * Transforms a list of "object"'s into a list of
 * "dirobject"'s. The list of "dirobject"'s is terminated
 * with a NULL-object at the end.
 *
 * Parameters:
 * @param  grism_file  - the full path to the grism image
 * @param  CONF_file   - the full path to the config file
 * @param  npixels     - the dimension of the grism image
 * @param  oblist      - the object list to start from
 * @param  spec_mod    - the model spectra
 * @param  obj_mod     - te direct emission models
 * @param  model_scale - the scale for the size of the direct object area
 * @param  int_type    - interpolation type
 *
 * Returns:
 * @return dirlist     - the dirobject list created
 */
dirobject **
oblist_to_dirlist2(char grism_file[], char CONF_file[], const  px_point npixels,
                   object  **oblist, spectral_models *spec_mod, object_models *obj_mod,
                   const double model_scale, const int int_type)
{

  dirobject     **dirlist;
  aperture_conf  *conf;
  gsl_matrix     *drzcoeffs;

  int nobjects=0;
  int i=0;
  int j=0;
  int beamID;
  int max_offs;

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  // load the extension numbers
  get_extension_numbers(grism_file, conf,conf->optkey1,conf->optval1);

  // get the  matrix with the drizzle coefficients
  drzcoeffs = get_crossdisp_matrix(grism_file, conf->science_numext);
  if (drzcoeffs->size1 < 2 || !drzcoeffs->size2)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "oblist_to_dirlist2:" " Could not get"
                 " the drizzle coefficients in file: %s\n", grism_file);

  // determine an offset from the PSF_OFFSET
  max_offs = (int)ceil(get_max_offset(conf));

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // allocate space for the dirobject list
  dirlist = (dirobject **) malloc((nobjects+1) * sizeof(dirobject *));
  if (dirlist == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "oblist_to_dirlist:" " Could not allocate"
                 " memory for pointers to %i dirobject objects", nobjects+1);

  // loop over all objectsconfig->camera
  for (i = 0; i < nobjects; i++)
    {
      // make sure that there are beams in the object
      if (oblist[i]->nbeams > 0)
        {

          if (has_aperture_dirim(obj_mod, oblist[i]))
            dirlist[j] = fill_dirobj_fromdirim(oblist[i], obj_mod);
          else
            // create a dirobject for each object
            dirlist[j] = fill_dirobject(oblist[i], npixels, drzcoeffs, model_scale, max_offs);

          fill_spectrum(oblist[i], dirlist[j], spec_mod, int_type);
          //      fill_dpsf_function(conf, dirlist[j]);

          // fill the xoffset and yoffset values.
          // within the gaussian models they are dummys.
          for (beamID=0; beamID < conf->nbeams; beamID++)
            {
              dirlist[j]->xy_off[beamID].x = 0.0;
              dirlist[j]->xy_off[beamID].y = 0.0;
            }
          j++;
        }
    }

  // terminate the dirobject list with NULL
  dirlist[j] = NULL;

  // release the memory for the drizzle matrix
  gsl_matrix_free(drzcoeffs);

  free_aperture_conf(conf);

  return dirlist;
}

/**
 * Function: fill_dirobject
 * Transforms one object into one dirobject.
 * Either transfers or computes the content
 * of the dirobject from the object.
 *
 * Parameters:
 * @param  actobject   - the object to be transformed
 * @param  npixels     - the dimension of the grism image
 * @param  drzcoeffs   - the drizzle coefficients
 * @param  model_scale - the scale for the size of the direct object area
 *
 * Returns:
 * @return actdir    - the dirobject created
 */
dirobject *
fill_dirobject(const object *actobject, const  px_point npixels,
               gsl_matrix *drzcoeffs, const double model_scale, const int max_offset)
{

  dirobject *actdir;
  d_point dirmod[4];
  double delx_a, dely_a;
  double delx_b, dely_b;
  gsl_vector *extention;

  // allocate space for the dirobject
  actdir = (dirobject *) malloc (sizeof (dirobject));
  if (actdir == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "fill_dirobject:" " Could not allocate"
                 " memory for a dirobject object");

  // compute the vector of maximum elongation along the large half axis
  delx_a = model_scale * actobject->beams[0].awidth * cos(actobject->beams[0].aorient);
  dely_a = model_scale * actobject->beams[0].awidth * sin(actobject->beams[0].aorient);

  // compute the vector of maximum elongation along the small half axis
  delx_b = model_scale * actobject->beams[0].bwidth * cos(actobject->beams[0].aorient + M_PI/2.0);
  dely_b = model_scale * actobject->beams[0].bwidth * sin(actobject->beams[0].aorient + M_PI/2.0);

  // corner refpoint + large vector + small vector
  dirmod[0].x = actobject->beams[0].refpoint.x + delx_a + delx_b;
  dirmod[0].y = actobject->beams[0].refpoint.y + dely_a + dely_b;

  // corner refpoint + large vector - small vector
  dirmod[1].x = actobject->beams[0].refpoint.x + delx_a - delx_b;
  dirmod[1].y = actobject->beams[0].refpoint.y + dely_a - dely_b;

  // corner refpoint - large vector + small vector
  dirmod[2].x = actobject->beams[0].refpoint.x - delx_a + delx_b;
  dirmod[2].y = actobject->beams[0].refpoint.y - dely_a + dely_b;

  // corner refpoint - large vector - small vector
  dirmod[3].x = actobject->beams[0].refpoint.x - delx_a - delx_b;
  dirmod[3].y = actobject->beams[0].refpoint.y - dely_a - dely_b;

  // transfer refpoint and ID
  actdir->ID = actobject->ID;
  actdir->refpoint.x = actobject->beams[0].refpoint.x;
  actdir->refpoint.y = actobject->beams[0].refpoint.y;

  // get the extentions on the reference points
  extention = get_refpoint_ranges(actobject);

  // derive and store min/max in x/y for the corners
  actdir->ix_min = (int)floor(MIN(MIN(dirmod[0].x,dirmod[1].x),MIN(dirmod[2].x,dirmod[3].x))+ gsl_vector_get(extention, 0) + 0.5) - max_offset;
  actdir->ix_max = (int)floor(MAX(MAX(dirmod[0].x,dirmod[1].x),MAX(dirmod[2].x,dirmod[3].x))+ gsl_vector_get(extention, 1) + 0.5) + max_offset;
  actdir->iy_min = (int)floor(MIN(MIN(dirmod[0].y,dirmod[1].y),MIN(dirmod[2].y,dirmod[3].y))+ gsl_vector_get(extention, 2) + 0.5) - max_offset;
  actdir->iy_max = (int)floor(MAX(MAX(dirmod[0].y,dirmod[1].y),MAX(dirmod[2].y,dirmod[3].y))+ gsl_vector_get(extention, 3) + 0.5) + max_offset;

  // derive the distortion scales along the major an minor axis
  actdir->drzscale = get_axis_scales(actobject->beams[0], drzcoeffs, npixels);

  actdir->dirim = NULL;

  // free the memory for the vector
  gsl_vector_free(extention);

  // return the dirobject
  return actdir;
}


/**
 * Function: fill_dirobj_fromdirim
 * Transforms an object which has a direct emission model associated
 * to into a direct object. The direct object created is returned
 *
 * Parameters:
 * @param  actobject     - the object to be transformed
 * @param  object_models - the structure with the direct emission models
 *
 * Returns:
 * @return actdir    - the dirobject created
 */
dirobject *
fill_dirobj_fromdirim(const object *actobject, object_models *objmodels)
{

  dirobject      *actdir;
  gsl_vector     *extention;
  dirim_emission *dirim;

  // allocate space for the dirobject
  actdir = (dirobject *) malloc (sizeof (dirobject));
  if (actdir == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "fill_dirobj_fromdirim:" " Could not allocate"
                 " memory for a dirobject object");

  // assign a direct emission model to a local variable
  dirim = get_dirim_emission(objmodels, actobject->beams[0].modimage);

  // transfer refpoint and ID
  actdir->ID = actobject->ID;
  actdir->refpoint.x = actobject->beams[0].refpoint.x;
  actdir->refpoint.y = actobject->beams[0].refpoint.y;

  // get the extentions on the reference points
  extention = get_refpoint_ranges(actobject);

  // derive and store min/max in x/y for the corners
  actdir->ix_min = (int)floor(actobject->beams[0].refpoint.x - dirim->xmean + gsl_vector_get(extention, 0) + 0.5);
  actdir->ix_max = (int)floor(actobject->beams[0].refpoint.x + dirim->xmean + gsl_vector_get(extention, 1) + 0.5);
  actdir->iy_min = (int)floor(actobject->beams[0].refpoint.y - dirim->ymean + gsl_vector_get(extention, 2) + 0.5);
  actdir->iy_max = (int)floor(actobject->beams[0].refpoint.y + dirim->ymean + gsl_vector_get(extention, 3) + 0.5);

  // fix the drizzle scale, since it has no
  // business in this context
  actdir->drzscale.x = 1.0;
  actdir->drzscale.y = 1.0;

  // store the direct emission
  // object in the direct object
  actdir->dirim = dirim;

  // free the memory for the vector
  gsl_vector_free(extention);

  // return the dirobject
  return actdir;
}

/**
 *
 * Function: get_refpoint_ranges
 *
 * Determines the differences of the reference point positions within
 * the beams of an object. The min/max values in x/y with respect to
 * the reference point of the first beam are determined and returned
 * in a vector.
 *
 * Parameters:
 * @param  actobject   - the object to be transformed
 *
 * Returns:
 * @return ret    - the dirobject created
 */
gsl_vector *
get_refpoint_ranges(const object *actobject)
{
  gsl_vector *ret;
  int j;
  //d_point reference;

  // allocate memory
  ret = gsl_vector_alloc(4);

  // initialize the vector with the default value
  gsl_vector_set(ret, 0, actobject->beams[0].refpoint.x);
  gsl_vector_set(ret, 1, actobject->beams[0].refpoint.x);
  gsl_vector_set(ret, 2, actobject->beams[0].refpoint.y);
  gsl_vector_set(ret, 3, actobject->beams[0].refpoint.y);


  // go over all beams in the object
  for (j=1; j < actobject->nbeams; j++)
    {
      // get the the new absolute mins and maxs in the reference points
      gsl_vector_set(ret, 0, MIN(gsl_vector_get(ret, 0) , actobject->beams[j].refpoint.x));
      gsl_vector_set(ret, 1, MAX(gsl_vector_get(ret, 1) , actobject->beams[j].refpoint.x));
      gsl_vector_set(ret, 2, MIN(gsl_vector_get(ret, 2) , actobject->beams[j].refpoint.y));
      gsl_vector_set(ret, 3, MAX(gsl_vector_get(ret, 3) , actobject->beams[j].refpoint.y));
    }

  // transform the absolute ranges into relative ones
  gsl_vector_set(ret, 0, gsl_vector_get(ret, 0) - actobject->beams[0].refpoint.x);
  gsl_vector_set(ret, 1, gsl_vector_get(ret, 1) - actobject->beams[0].refpoint.x);
  gsl_vector_set(ret, 2, gsl_vector_get(ret, 2) - actobject->beams[0].refpoint.y);
  gsl_vector_set(ret, 3, gsl_vector_get(ret, 3) - actobject->beams[0].refpoint.y);

  // return the result
  return ret;
}

/**
 * Function: fill_spectrum
 * Transfers the wavelength information from an object
 * to a dirobject. Store the minimum and maximum wavelength
 * as well as associated flux values in the structure.
 *
 * Parameters:
 * @param actobject - the object to be transformed
 * @param actdir    - the dirobject to store the flux values in
 * @param spec_mod  - the spectral model list
 * @param int_type  - the interpolation type to be used
 */
void
fill_spectrum(const object *actobject, dirobject *actdir,
              spectral_models *spec_mod, const int int_type)
{
  energy_distrib *sed=NULL;
  //  double *sed_wavs;
  //  double *sed_flux;

  beam onebeam;

  //int i_type=0;
  //  int nwavs=0;
  //int i = 0;
  int j = 0;

  // the data are derived from the 1st non-zero beam in the object;
  // go along the beams util you find a non-zero beam
  while (actobject->beams[j].flux == NULL && j < actobject->nbeams)
    j++;
  // store this beam
  onebeam = actobject->beams[j];

  // do something only if the first beam has flux values
  if (onebeam.flux != NULL && j < actobject->nbeams)
    {

      if (onebeam.modspec > 0 && spec_mod)
        {
          // get the model spectrum
          sed = get_model_sed(spec_mod, onebeam.modspec);

          // mark that the SED does NOT come
          // from broad band colours
          actdir->bb_sed=0;
        }
      else
        {

          sed = make_sed_from_beam(onebeam, int_type, actobject->ID);

          /*
          //
          // the code below is a mess!!!!
          // MUST be re-factorized!!
          //

          // determine the number of flux values, allocate space
          nwavs = (onebeam.flux->size)/2;

          // determine the number of flux values, allocate space
          nwavs = (onebeam.flux->size)/2;
          sed = (energy_distrib*) malloc(sizeof(energy_distrib));
          sed_wavs = (double*) malloc(nwavs*sizeof(double));
          sed_flux = (double*) malloc(nwavs*sizeof(double));

          // transfer the data into the 'local' arrays
          for (i=0; i<nwavs; i++)
            {
              sed_wavs[i] = gsl_vector_get(onebeam.flux, 2*i);
              sed_flux[i] = gsl_vector_get(onebeam.flux, 2*i+1);
            }

          // transfer the vector and length information
          // to the SED object
          sed->npoints    = nwavs;
          sed->wavelength = sed_wavs ;
          sed->flux       = sed_flux;

          // check the interpolation type
          i_type = check_interp_type(int_type, nwavs, actobject->ID);
          if (nwavs > 1)
            {
              // if possible, allocate and
              // initialize an interpolation object
              if (i_type == 2)
                sed->interp     = gsl_interp_alloc (gsl_interp_polynomial, (size_t)nwavs);
              else if (i_type == 3)
                sed->interp     = gsl_interp_alloc (gsl_interp_cspline, (size_t)nwavs);
              else
                sed->interp     = gsl_interp_alloc (gsl_interp_linear, (size_t)nwavs);

              sed->accel      = gsl_interp_accel_alloc ();
              gsl_interp_init (sed->interp, sed->wavelength, sed->flux, (size_t)sed->npoints);
            }
          else
            {
              // if no interpolation, set
              // the memeber to NULL
              sed->interp = NULL;
              sed->accel= NULL;
            }
          */
          // mark that the SED DOES come
          // from broad band colours
          actdir->bb_sed=1;
        }
    }
  else
    {
      // make an error if there are no flux values
      aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                  "aXe_PETCONT: " "the OAF does not have flux values\n");
    }

  // transfer the SED object to the direct object
  actdir->SED = sed;
}


/**
 * Function: make_sed_from_beam
 *
 * Parameters:
 * @param onebeam  - the object to be transformed
 * @param int_type - the interpolation type to be used
 * @param objID    - the interpolation type to be used
 */
energy_distrib *
make_sed_from_beam(const beam onebeam, const int int_type, const int objID)
{
  energy_distrib *sed;
  double *sed_wavs;
  double *sed_flux;

  int nwavs;
  int i_type;
  int i;

  // determine the number of flux values, allocate space
  nwavs = (onebeam.flux->size)/2;

  // determine the number of flux values, allocate space
  nwavs = (onebeam.flux->size)/2;
  sed = (energy_distrib*) malloc(sizeof(energy_distrib));
  sed_wavs = (double*) malloc(nwavs*sizeof(double));
  sed_flux = (double*) malloc(nwavs*sizeof(double));

  // transfer the data into the 'local' arrays
  for (i=0; i<nwavs; i++)
    {
      sed_wavs[i] = gsl_vector_get(onebeam.flux, 2*i);
      sed_flux[i] = gsl_vector_get(onebeam.flux, 2*i+1);
    }

  // transfer the vector and length information
  // to the SED object
  sed->npoints    = nwavs;
  sed->wavelength = sed_wavs ;
  sed->flux       = sed_flux;

  // check the interpolation type
  i_type = check_interp_type(int_type, nwavs, objID);
  if (nwavs > 1)
    {
      // if possible, allocate and
      // initialize an interpolation object
      if (i_type == 2)
        sed->interp     = gsl_interp_alloc (gsl_interp_polynomial, (size_t)nwavs);
      else if (i_type == 3)
        sed->interp     = gsl_interp_alloc (gsl_interp_cspline, (size_t)nwavs);
      else
        sed->interp     = gsl_interp_alloc (gsl_interp_linear, (size_t)nwavs);

      sed->accel      = gsl_interp_accel_alloc ();
      gsl_interp_init (sed->interp, sed->wavelength, sed->flux, (size_t)sed->npoints);
    }
  else
    {
      // if no interpolation, set
      // the memeber to NULL
      sed->interp = NULL;
      sed->accel  = NULL;
    }

  // return the sed
  return sed;
}

/**
 * Function: get_dirobject_from_list
 * The function identifies a dirobject
 * in a list of dirobjects. The identification
 * is made on the attribute ID. The identified
 * dirobject is returned. If no dirobject could
 * be identified, the NULL object, which is at
 * the end of each dirobject list, is returned.
 *
 * Parameters:
 * @param  dirlist   - the dirobject list
 * @param  ID        - the ID number to be identified
 *
 * Returns:
 * @return dirobject - the identified dirobject or the 'NULL'-dirobject
 */
dirobject *
get_dirobject_from_list(dirobject ** dirlist, const int ID)
{

  //dirobject * actdir;
  int i, ndirs = 0;

  // count the number of dirobjects in the list
  while (dirlist[ndirs] != NULL)
    ndirs++;

  // loop over all dirobject
  for (i = 0; i < ndirs; i++)
    {
      // try to identify a dirobject,
      // return it in case of a postitive identification
      if (dirlist[i]->ID == ID)
        return dirlist[i];
    }

  // return the NULL-dirobject at the end of the list
  return dirlist[ndirs];
}

/**
 * Function: get_dirobject_meanpos
 * The function computes the average pixel coos of a
 * direct object. This is done in a very simple way, by
 * averaging the maximum and minimum pixel coordinates in
 * both, x and y. Neither the shape nor the
 * the intensity in the individual pixels are taken
 * into account.
 *
 * Parameters:
 * @param  actdir  - the direct object
 *
 * Returns:
 * @return m_point - the means coo's in x and y
 */
d_point
get_dirobject_meanpos(dirobject *actdir)
{
  d_point m_point;

  m_point.x = ((double)actdir->ix_max + (double)actdir->ix_min) / 2.0;
  m_point.y = ((double)actdir->iy_max + (double)actdir->iy_min) / 2.0;

  return m_point;
}
/**
 * Function: get_beam_for_beamspec
 * The function selects for a given beamspec the corresponding
 * beam from an object list. The identification is done
 * via objectID and beamID. An error is thrown in case that
 * no matching beam could be found.
 *
 * Parameters:
 * @param oblist   - the object list to identify a beam from
 * @param nobjects - the number of objects in the object list
 * @param actspec  - the model spectrum to identify a beam for
 *
 * @return actbeam - the identified beam
 */
beam
get_beam_for_beamspec(object **oblist, const int nobjects,
                      const beamspec *actspec)
{

  beam actbeam;

  int i, j;

  // set the beam ID to -1 to identify
  // failed identification
  actbeam.ID = -1;

  // go over all objects in the list
  for (i = 0; i < nobjects; i++)
    {

      // search for a matching object ID
      if (oblist[i]->ID == actspec->objectID)
        {

          // go over all beams in the matchin object
          for (j=0; j < oblist[i]->nbeams; j++)
            {

              // search for a matching beam ID
              if (oblist[i]->beams[j].ID == actspec->beamID)
                actbeam = oblist[i]->beams[j];
            }
        }
    }

  // report an error in case that the identification failed
  if (actbeam.ID == -1)
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "aXe_PETCONT: " "object ID %i, Beam %c not found\n", actspec->objectID, BEAM(actspec->beamID));

  // return the identified beam
  return actbeam;
}

/**
 * Function: get_beamspec_from_list
 * The functions extracts and returns a specific model beam
 * out of the list of model beams. The requested model
 * beam is identified on the basis of the aperture ID and
 * the beam ID. Without positive identification the last modelled
 * spectrum in the ist is returned, which is NULL.
 *
 * Parameters:
 * @param  speclist -
 * @param  aperID   - the object ID of the requested beam
 * @param  beamID   - the beam ID of the requested beam
 *
 * Returns:
 * @return ret      - the identified beam (or the NULL beam at the end of the list)
 */
beamspec *
get_beamspec_from_list(beamspec **speclist, const int aperID, const int beamID)
{
  //beamspec *ret;
  int i=0;

  while (speclist[i] != NULL)
    {
      if (speclist[i]->objectID == aperID && speclist[i]->beamID == beamID)
        {
          break;
        }
      i++;
    }

return speclist[i];
}

/**
 * Function: print_dirobject
 * Prints the data in a dirobject onto the screens.
 * "Natural" format is used, which means printed are first
 * the upper corners, then the refpoint, then
 * the lower corners.
 *
 * Parameters:
 * @param actdir - the dirobject to be printed
 */
void
print_dirobject(const dirobject * actdir){
  int i=0;
  int npoints=0;

  fprintf(stdout, "Object ID: %i\n", actdir->ID);
  fprintf(stdout,"%5i,%5i              %5i,%5i\n", actdir->ix_min, actdir->iy_max, actdir->ix_max, actdir->iy_max);
  fprintf(stdout,"        %7.1f,%7.1f         \n", actdir->refpoint.x, actdir->refpoint.y);
  fprintf(stdout,"%5i,%5i              %5i,%5i\n\n", actdir->ix_min, actdir->iy_min, actdir->ix_max, actdir->iy_min);


  npoints = actdir->SED->npoints;

  fprintf(stdout,"Wavelengths:\n");
  fprintf(stdout,"Minimum: %.1f,%.3e  Maximum: %.1f, %.3e\n", actdir->SED->wavelength[0], actdir->SED->flux[0],
          actdir->SED->wavelength[npoints-1],  actdir->SED->flux[npoints-1]);
  for (i=0; i<npoints;i++)
    fprintf(stdout,"Wavelength: %.1f, Flux: %.5e\n", actdir->SED->wavelength[i], actdir->SED->flux[i]);
  fprintf(stdout,"\n\n");
}


/**
 * Function: check_interp_type
 * The function checks the requested interpolation type
 * against the number of flux values in the models.
 * The interpolation type is adjusted in case that there
 * are not enough flux values for the reuqested type.
 *
 * Parameters:
 * @param inter_type - the requested interpolation
 * @param n_flux     - the number of flux values
 * @param ID         - the object ID
 *
 * Returns:
 * @return type       - the feasable interpolation type
 */
int
check_interp_type(const int inter_type, const int n_flux, const int ID)
{
  int type=0;

  // for interpolation types other than linear at leas
  // three flux values must be there
  if (n_flux < 3 && inter_type > 1)
    {
      // make a warnign message if there are
      // not enough flux value
      if (inter_type > 2)
        aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                     "\naXe_PETCONT: Object: %i: interpolation switched from \"spline\" to \"linear\"!\n", ID);
      else
        aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                     "\naXe_PETCONT: Object: %i: interpolation switched from \"polynomial\" to \"linear\"!\n", ID);

      // adjust the interpolation type to linear
      type = 1;
    }
  else
    {
      // copy the requeted interpolation type
      // if there are enough flux values
      type = inter_type;
    }

  // return the adjusted type
  return type;
}

/**
 *
 * Function: free_dirlist
 * Releases the memory allocated
 * for a list of dirobjects.
 * after releasing memory, all elements
 * are set to NULL/
 *
 * Parameters:
 * @param dirlist  - the list of dirobjects
 * @paran spec_mod - describes what all to free
 */
void
free_dirlist (dirobject ** dirlist)
{
  int i, ndirs = 0;

  // count the number of dirobjects
  while (dirlist[ndirs] != NULL)
    ndirs++;

  // go over each item in the list
  for (i = 0; i < ndirs; i++)
    {

      // check for broad band SED
      if (dirlist[i]->bb_sed)
        // free the energy distribution
        free_enerdist (dirlist[i]->SED);

      // fre the dirobject
      free (dirlist[i]);

      // set the dirobject to NULL
      dirlist[i] = NULL;
    }

  // free the list
  free (dirlist);

  // set the list to NULL
  dirlist = NULL;
}


/**
 * Function: free_enerdist
 * The function releases all the memory
 * of a SED object
 *
 * Parameters:
 * @param sed - the SED object
 */
void
free_enerdist(energy_distrib *sed)
{

  // free the two arrays
  free(sed->wavelength);
  free(sed->flux);

  // free the intepolation
  // structures if defined
  if (sed->interp != NULL)
    gsl_interp_free (sed->interp);
  if (sed->accel != NULL)
    gsl_interp_accel_free (sed->accel);

  // free the object itsel
  free(sed);
}

/**
 * Function: free_speclist
 * The function releases all the memory
 * of a list of beamspec's.
 *
 * Parameters:
 * @param speclist - a list of beamspecs
 *
 */
void
free_speclist(beamspec **speclist)
{

  int i=0;

  // go over each item in the list
  while (speclist[i] != NULL)
    {

      // free the memory in the matrix
      gsl_matrix_free (speclist[i]->model);

      // free the beamspec object itself
      free (speclist[i]);

      // set the beamspec to NULL,
      // increment the counter
      speclist[i++] = NULL;
    }

  // free the memory of the list
  free(speclist);

  // set the list to NULL
  speclist = NULL;
}

/**
 * Function: free_tracedata
 * Releases the memory in
 * a tracedata structure.
 *
 * Parameters:
 * @param acttrace - the structure to be freed
 */
void
free_tracedata(tracedata *acttrace)
{
  if (acttrace->npoints > 0)
    {
      gsl_vector_free(acttrace->dx);
      gsl_vector_free(acttrace->dy);
      gsl_vector_free(acttrace->xi);
      gsl_vector_free(acttrace->lambda);
      gsl_vector_free(acttrace->dlambda);
      gsl_vector_free(acttrace->flux);
      gsl_vector_free(acttrace->gvalue);
    }
  free(acttrace);
}


/**
 * Function: get_flux_from_SED
 *  The functions returns the flux value at a given
 *  wavelength for a certain SED.
 *  In case that the wavelength is beyound the
 *  wavelength interval where the SED is defined,
 *  the closest model value in wavelength is returned.
 *  Otherwise the interpolated value as defined in the
 *  SED object is returned
 *
 * @param  sed     - the SED-object to derive a flux value from
 * @param  in_wave - the wavelength to compute the flux value for
 *
 * @return flux       - the flux value for the requested wavelength
 */
double
get_flux_from_SED(const energy_distrib *sed, double in_wave)
{
  double flux=0.0;

  if (in_wave <= sed->wavelength[0]){
    // give the lowest defined value if
    // the requested wavelength is lower
    flux = sed->flux[0];
}
  else if (in_wave > sed->wavelength[sed->npoints-1]){
    // give the highest defined value if
    // the requested wavelength is higher
    flux = sed->flux[sed->npoints-1];}
  else {
    // derive the interpolated value
    flux = gsl_interp_eval(sed->interp, sed->wavelength, sed->flux, in_wave, sed->accel);}

  // return the interpolated flux value
  return flux;
}


/**
 * Function: get_aveflux_from_SED
 * The function computes and returns the integrated flux of a given SED object
 * over a given interval at a given wavelength position.
 * Beyond the wavelength interval the SED is defined on, the flux is taken
 * as constant.
 *
 * @param  sed         - the SED-object to derive a flux value from
 * @param  in_wave     - the wavelength to average the flux value at
 * @param  wave_interv - the wavelength interval to average over
 *
 * @return flux       - the average for the requested wavelength interval
 */
double
get_aveflux_from_SED(const energy_distrib *sed, double in_wave, double wave_interv)
{
  double wav_min=0.0;
  double wav_max=0.0;

  double range_1=0.0;
  double range_2=0.0;

  double flux=0.0;

  // check whether there is only one
  // flux point
  if (sed->npoints < 2)
    // return the value at the single
    // flux point
    return sed->flux[0];

  // compute the interval coundaries
  wav_min = in_wave - 0.5*wave_interv;
  wav_max = in_wave + 0.5*wave_interv;

  // check whether the upper wavelength
  // interval boundary is just shorter than
  // the SED range
  if (wav_max < sed->wavelength[0]){
    // give the smallest SED flux value
    flux = sed->flux[0];}

  // check whether the lower wavelength
  // interval boundary is just longer than
  // the SED range
  else if (wav_min > sed->wavelength[sed->npoints-1]){
    // give the highest defined value
    flux = sed->flux[sed->npoints-1];}

  // check whether the interval covers the lower
  // end covered by the SED
  else if (wav_min < sed->wavelength[0] && wav_max > sed->wavelength[0])
    {
      // get the lower interval
      range_1 = sed->wavelength[0] - wav_min;

      // check whether the upper boundary is
      // covered by the SED
      if (wav_max <= sed->wavelength[sed->npoints-1])
        {
          //fprintf(stderr,"got here!\n");
          // get the integration interval
          range_2 = wav_max - sed->wavelength[0];

          // average the two integrals
          flux = (range_1*sed->flux[0] + gsl_interp_eval_integ(sed->interp, sed->wavelength, sed->flux, sed->wavelength[0], wav_max, sed->accel)) / wave_interv;
        }
      else
        {
          // get the integration interval
          range_2 = wav_max - sed->wavelength[sed->npoints-1];

          // average the three integrals
          flux = (range_1*sed->flux[0] + range_2*sed->flux[sed->npoints-1]
                  + gsl_interp_eval_integ(sed->interp, sed->wavelength, sed->flux, sed->wavelength[0], sed->wavelength[sed->npoints-1], sed->accel)) / wave_interv;
        }
    }

  // check whether the interval start is beyound
  // the lower
  // end covered by the SED
  else if (wav_min >= sed->wavelength[0] && wav_max >= sed->wavelength[0])
    {
      // check whether the interval is completely
      // covered by the SED
      if (wav_max <= sed->wavelength[sed->npoints-1])
        {
          // compute the integral over the SED and the interval
          flux = gsl_interp_eval_integ(sed->interp, sed->wavelength, sed->flux, wav_min, wav_max, sed->accel) / wave_interv;
        }
      // if the interval is partly out
      // of the SED
      else
        {
          // compute the range which is in
          range_1 = sed->wavelength[sed->npoints-1] - wav_min;
          // compute the range which is out
          range_2 = wav_max - sed->wavelength[sed->npoints-1];

          // compute the integrated flux via weighted summation
          // of the in part and the out part
          flux = (gsl_interp_eval_integ(sed->interp, sed->wavelength, sed->flux, wav_min, sed->wavelength[sed->npoints-1], sed->accel)
                  + range_2*sed->flux[sed->npoints-1]) / wave_interv;
        }
    }
  // return the integrated flux value
  return flux;
}


/**
 *
 * Function: get_flambda_from_magab
 * The subroutine calculates the flambda value for a
 * mag_AB value given with its wvavelength as input
 * parameters.
 *
 * Parameters:
 * @param  mag     - the mag_AB value
 * @param  lambda  - the wavelength for mag_AB
 *
 * Returns:
 * @return flambda - the calculated flambda value
 */
double
get_flambda_from_magab(double mag, double lambda)
{
  double flambda=0.0;
  double fnu=0.0;

  fnu     = pow(10.0, -0.4*(mag+48.6));
  flambda = 1.0e+16*LIGHTVEL*fnu/(lambda*lambda);

  return flambda;
}

/*
 * Function: fill_gaussvalues
 * The function determines the wavelength dependent emission
 * for a given object at a given point in the gaussian emission model.
 * The values are filled into the
 * according vector of the tracedata structure.
 *
 * Parameters:
 * @param dpixel      - the coordinates of the point
 * @param actbeam     - the beam to derive
 * @param actdir      - the direct object
 * @param lambda_ref  - the reference wavelength
 * @param conf        - the configuration structure
 * @param acttrace    - the tracedata structure
 *
 */
void
fill_gaussvalues(const d_point dpixel, const beam actbeam,
                 const dirobject *actdir, const double lambda_ref,
                 const aperture_conf *conf, const double psf_offset,
                 tracedata *acttrace)
{
  beam new_beam;

  double dpsf=0.0;
  //double gval=0.0;
  double lambda=0.0;
  //double factor=0.0;
  //double nom=0.0;
  //double denom=0.0;

  int i=0;
  int minpsf_flagg=0;

  // get the acurate emission value with subgridding
  //gval   = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);
  //fprintf(stdout, "gvalue: %e %f %f", gval, actdir->drzscale.x, actdir->drzscale.y);

  //  0: fastest scenario: fill the vector with
  //     identical values: 34s
  //  gsl_vector_set_all(acttrace->gvalue, i, gval);


  // 1: iterate over the vector and enter the same
  //    value: 35s
  if (conf->psfrange && conf->psfcoeffs)
    {
      // go over all tracedata (that means wavelength) points
      for (i=0; i < acttrace->npoints; i++)
        {
          // derive the wavelength
          lambda = gsl_vector_get(acttrace->lambda, i)/10.0;

          // derive the correction of the psf at the wavelength

          // 3: calculating the psf-offset: 90s
          dpsf   = psf_offset + get_dpsf(lambda_ref, lambda, conf, actbeam);
          //      dpsf=0.0;
          //fprintf(stdout, "dpsf: %f", dpsf);


          // derive a beam with the correct widths at the wavelength

          // 2: creating a new beam with a different size: 50s
          new_beam = get_newbeam(actbeam,dpsf);

          // transport the
          if (new_beam.ID)
            minpsf_flagg=1;

          // derive the correction factor for the emission using the correct widths
          // 4: computing the correction factor for the wavelength: 255s
          //factor = get_emodel_value(dpixel,new_beam,actdir->drzscale)/get_emodel_value(dpixel,actbeam,actdir->drzscale);
          //nom = get_emodel_value(dpixel,new_beam,actdir->drzscale);
          //denom = get_emodel_value(dpixel,actbeam,actdir->drzscale);
          //if (nom != 0.0 && denom != 0.0)
          //  factor = nom / denom;
          //else
          //  factor = 1.0;
          //fprintf(stdout, "factor: %f ", factor);

          // store the emission value, taking into account the correction
          //gsl_vector_set(acttrace->gvalue, i, gval*factor);
          gsl_vector_set_all(acttrace->gvalue, get_sub_emodel_value(dpixel, new_beam, actdir->drzscale));
          //fprintf(stdout, "gvalue * factor: %e ", gval*factor);
          }
    }
  else
    {
      // derive a beam with the correct widths at the wavelength
      new_beam = get_newbeam(actbeam,psf_offset);
      if (new_beam.ID)
        minpsf_flagg=1;

      // derive the correction factor for the emission using the correct widths
      //factor = get_emodel_value(dpixel,new_beam,actdir->drzscale)/get_emodel_value(dpixel,actbeam,actdir->drzscale);
      //nom = get_emodel_value(dpixel,new_beam,actdir->drzscale);
      //denom = get_emodel_value(dpixel,actbeam,actdir->drzscale);
      //if (nom != 0.0 && denom != 0.0)
      //  factor = nom / denom;
      //else
      //  factor = 1.0;

      // store the emission value, taking into account the correction
      //gsl_vector_set_all(acttrace->gvalue, gval*factor);
      gsl_vector_set_all(acttrace->gvalue, get_sub_emodel_value(dpixel, new_beam, actdir->drzscale));

    }
  if (minpsf_flagg)
    aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                 "\naXe_PETCONT: points in PSF of object %i beam %c smaller than PSF_MIN=%f! Set to PSF_MIN\n", actdir->ID, BEAM(actbeam.ID), MINPSF);

}

/*
 * Function: get_newbeam
 * The function creates a new beam with modified widths.
 * The modification of the width is given as the differential
 * value. The new beam is used to compute an emission value
 * for a specific wavelength.
 *
 * Parameters:
 * @param actbeam   - the original beam
 * @param dpsf      - the changes in width
 *
 * Returns:
 * @return new_beam - the new beam
 */
beam
get_newbeam(const beam actbeam, const double dpsf)
{
  beam new_beam;

  new_beam.ID=0;

  // check the sign of the modification
  if (dpsf >0.0)
    {
      // positive sign:
      // add (in quadrature) the psf-modification
      new_beam.awidth   = sqrt(SQR(actbeam.awidth) + SQR(dpsf));
      new_beam.bwidth   = sqrt(SQR(actbeam.bwidth) + SQR(dpsf));
    }
  else
    {

      // check whether awidth is big enough
      if (actbeam.awidth+dpsf < MINPSF)
        {
          // set the minimum sign, and set the psf to MINPSF
          new_beam.ID=1;
          new_beam.awidth = MINPSF;
        }
      else
        {
          // negative sign:
          // subtract (in quadrature) the psf-modification
          new_beam.awidth = sqrt(SQR(actbeam.awidth) - SQR(dpsf));
        }

      // check whether bwidth is big enough
      if (actbeam.bwidth+dpsf < MINPSF)
        {
          // set the minimum sign, and set the psf to MINPSF
          new_beam.ID=1;
          new_beam.bwidth = MINPSF;
        }
      else
        {
          // negative sign:
          // subtract (in quadrature) the psf-modification
          new_beam.bwidth = sqrt(SQR(actbeam.bwidth) - SQR(dpsf));
        }
    }


  // transfer necessary data
  new_beam.aorient = actbeam.aorient;
  new_beam.refpoint = actbeam.refpoint;

  // return the new beam
  return new_beam;
}

/**
 * Function: get_dpsf
 * The function computes the difference of the psf width
 * at two different wavelengths. The dependence of the
 * psf as a function of the wavelength is supposed to be
 * a polynomial
 *
 * Parameters:
 * @param lambda_ref - the reference wavelength
 * @param lambda     - the wavelength to evaluate the difference for
 * @param conf       - the configuration structure
 *
 * Returns:
 * @return dpsf      - the psf difference
 */
double
get_dpsf(const double lambda_ref, const double lambda,
         const aperture_conf *conf, const beam actbeam)
{

  double dpsf=0.0;
  double lmin=0.0;
  double lmax=0.0;

  // get the minimum and maximum wavelength
  // for the dependency
  lmin = gsl_vector_get(conf->psfrange, 0);
  lmax = gsl_vector_get(conf->psfrange, 1);

  // the reference wavelength must be inside of the
  // area where the polynomial is valid
  if (lambda_ref < lmin){
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "aXe_PETCONT:" "Reference wavelength %f"
                 " is smaller than minimum wavelength: %f!\n", lambda_ref, lmin);
  }
  else if (lambda_ref > lmax){
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "aXe_PETCONT:" "Reference wavelength %f"
                 " is larger than maximum wavelength: %f!\n", lambda_ref, lmax);
  }
  else{
    if (lambda < lmin)
      // if the target wavelength is larger,
      // evaluate at the upper border
      dpsf = get_polyN_gsl(lmin, conf->psfcoeffs) - get_polyN_gsl(lambda_ref, conf->psfcoeffs);
    else if (lambda > lmax)
      // if the target wavelength is smaller,
      // evaluate at the lower border
      dpsf = get_polyN_gsl(lmax, conf->psfcoeffs) - get_polyN_gsl(lambda_ref, conf->psfcoeffs);
    else
      // evaluate the exact position
      dpsf = get_polyN_gsl(lambda, conf->psfcoeffs) - get_polyN_gsl(lambda_ref, conf->psfcoeffs);
  }

  // return the result
  return dpsf;
}

/**
 * Function: get_sub_emodel_value
 * The function evaluates the emission of a 2D gauss model
 * with the parameters as given in a beam at a particular position.
 * This function derives the function values on a series of
 * grid positions +-.5pixels in x/y around the requested
 * positions. This avoids rounding problems.
 * The 1D grisdsize is set by the macro NSUB in the
 * header file
 *
 * Parameters:
 * @param dpixel   - the point to evaluate the emission model
 * @param actbeam  - the beam to set up the model for
 * @param drzscale - the relative pixelscale at the model position
 *
 * Returns:
 * @return sval    - the value of the emission model
 */
double
get_sub_emodel_value(const d_point dpixel, const beam actbeam,
                     const d_point drzscale)
{
  d_point dtmp;

  double sval   = 0.0;
  double step   = 0.0;
  double offset = 0.0;

  int irange = 0;
  int kk=0, ll=0;

  // convert the number of steps to a local integer
  irange = (int)NSUB;

  // compute the step size
  step = 1.0/(2.0*(double)NSUB);

  // compute the initial offset
  offset = step/2.0;

  for (kk=-irange; kk < irange; kk++)
    {
      for (ll=-irange; ll < irange; ll++)
        {
          // determine the actual grid position
          dtmp.x = dpixel.x + (double)kk * step + offset;
          dtmp.y = dpixel.y + (double)ll * step + offset;

          // get the value at the grid position
          sval = sval + get_emodel_value(dtmp, actbeam, drzscale);
        }
    }

  // normalize the result
  sval = sval *step * step;

  // return the result
  return sval;
}

/**
 * Function: get_emodel_value
 * The function evaluates the emission of a 2D gauss model
 * with the parameters as given in a beam at a particular position.
 * A different pixelscale at the emission point due to
 * geometrical distortion can be taken into account.
 *
 * Parameters:
 * @param dpixel   - the point to evaluate the emission model
 * @param actbeam  - the beam to set up the model for
 * @param pixscale - the relative pixelscale at the model position
 *
 * Returns:
 * @return evalue  - the value of the emission model
 */
double
get_emodel_value(const d_point dpixel, const beam actbeam,
                 const d_point drzscale)
     //get_emodel_value(const px_point dpixel, const beam actbeam, const d_point drzscale)
{
  double evalue=0.0;

  double amod, bmod;
  double xrel, yrel;
  double angle;

  double arg;

  // apply the correction due to geom. distortion
  amod =  actbeam.awidth / drzscale.x;
  bmod =  actbeam.bwidth / drzscale.y;

  // make some precalculations
  xrel = dpixel.x - actbeam.refpoint.x;
  yrel = dpixel.y - actbeam.refpoint.y;
  angle = actbeam.aorient;

  // determine the argument of the exponent
  arg =
    SQR(( xrel*cos(angle) + yrel*sin(angle)) / amod) +
    SQR((-xrel*sin(angle) + yrel*cos(angle)) / bmod);

  // determine the emission value
  evalue = 0.5/(amod*bmod*M_PI)*exp(-0.5*arg);

  // return the emission value
  return evalue;
}


double
get_psf_WFC(const double lambda)
{
  double psf;

  psf = 1.2994407614072 + 0.11290883734113e-02*lambda - 0.44245634760175e-06*pow(lambda,2.0);

  return psf;
}

double
get_psf_HRC(const double lambda)
{
  double psf = 0.0;

  psf = 8.1986545952891 - 0.82947763250959e-01*lambda          + 0.40134114766915e-03*pow(lambda,2.0) -
                          0.94651169604925e-06*pow(lambda,3.0) + 0.11804479784053e-08*pow(lambda,4.0) -
                          0.74395828991899e-12*pow(lambda,5.0) + 0.18657361078105e-15*pow(lambda,6.0);

  return psf;
}

double
get_psf_SBC(const double lambda)
{
  double psf = 0.0;

  psf = 7.9833718878881 - 0.11283853224345    *lambda          + 0.64235625522476e-03*pow(lambda,2.0)
                        - 0.12471899390220e-05*pow(lambda,3.0);

  return psf;

}

/**
 *
 * Function: get_polyN_gsl
 * The function evaluates a polynomial
 * at a given position.
 *
 * Parameters:
 * @param x      - the point to evaluate the polynomial
 * @param params - gsl-vector with the coefficients of the polynomial
 *
 * Returns:
 * @return p  - the value of the polynomial
 */
double
get_polyN_gsl (const double x, const gsl_vector *params)
{
  int i;

  double p=0.0;

  // go over the gsl-vector and
  // sum up the individual terms
  for (i=0; i<(int)params->size; i++)
    {
      p += gsl_vector_get(params,i)*pow(x,i);
    }

  // return the value
  return p;
}
