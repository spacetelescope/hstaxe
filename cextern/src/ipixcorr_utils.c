/**
 * Subroutines to correct for the intra-pixel sensitivity
 * variation in e.g. NICMOS pixels. There are routines
 * to corect both, an extracted SPC's and a list
 * of PET pixels.
 *
 */
#include "fitsio.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multifit_nlin.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "aXe_errors.h"
#include "spc_spc.h"
#include "spce_pathlength.h"
#include "fringe_conf.h"
#include "spc_wl_calib.h"
#include "trfit_utils.h"
#include "ipixcorr_utils.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

/**
 * Function: get_intpix_corr
 * The function determines and returns a correction factor with
 * values stored in an interpolator at a certain wavelength.
 * To derive the result the dispersion solution is inverted,
 * then the trace position for the wavelength is determined.
 * The the fractional pixel value in y is computed and the
 * corresponding correction factor is returned.
 *
 * Parameters:
 * @param actbeam - beam to correct
 * @param cdisp   - the dispersion coefficients for the beam
 * @param ipcorr  - the correction values stored in an interpolator
 * @param lambda  - the wavelength to correct
 *
 * Returns:
 * @return cfactor - the correction factor
 */
double
get_intpix_corr(beam *actbeam, gsl_vector *cdisp, interpolator *ipcorr,
		float lambda)
{
  double tlength;

  double xvalue;
  double yfract;

  double cfactor=0.0;

  // get the trace length for the dispersion
  // solution and wavelength
  tlength = get_tlength_from_dispersion(cdisp, lambda);

  // compute the x-offset from the reference point
  // for this trace length value
  xvalue = get_xvalue_from_tlength(actbeam, tlength);

  // compute the fractional y-value
  yfract = get_yfract_for_xvalue(actbeam, xvalue);

  // get the correction factor
  cfactor = eval_interp(ipcorr, yfract);

  fprintf(stdout, "xdiff: %e, yfrac: %e, factor: %e ", xvalue, yfract, cfactor);

  //fprintf(stdout, "lambda: %f, tlength: %f, xval: %f, yval: %f, yfract: %f, cfactor: %f\n",
  // 	  lambda, tlength, actbeam->refpoint.x+xvalue, actbeam->refpoint.y+actbeam->spec_trace->func (xvalue, actbeam->spec_trace->data),yfract, cfactor);

  // return the correction factor
  return cfactor;
}


/**
 * Function: get_yfract_for_xvalue
 * The function computes the fractional y-value for a trace
 * position in a beam. The trace position is given as the
 * x-offset position with respect to the reference position
 *
 * Parameters:
 * @param actbeam - beam to correct
 * @param xvalue  - the x-value
 *
 * Returns:
 * @return yfract - the fractional y-value
 */
double
get_yfract_for_xvalue(const beam *actbeam, double xvalue)
{
  double dy;
  double yabs;
  double yfract;

  // compute the y-value relative to the reference point
  dy = actbeam->spec_trace->func (xvalue, actbeam->spec_trace->data);

  // compute the absolute y-value on the chip
  yabs = actbeam->refpoint.y + dy;

  // compute the fractional pixel of this y-value
  yfract = yabs - floor(yabs);

  // return the fractional y-value
  return yfract;
}


/**
 * Function: get_xvalue_from_tlength
 * The function computes the x-offset position from the reference point
 * for a given tracelength in a given beam.
 * The solution is numerically derived and therefore
 * does work for any reasonable tracefunction.
 *
 * Parameters:
 * @param actbeam - the beam to compute the x-offset
 * @param tlength - tracelength
 *
 * Returns:
 * @return xdiff - the x-offset
 */
double
get_xvalue_from_tlength(beam *actbeam, double tlength)
{
  int iter=0;
  int status;

  double xdiff=0.0;
  d_point x_interv;

  // define and initialize the solver
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver            *s = gsl_root_fsolver_alloc (T);

  gsl_function F;
  tlength_pars *tpars;

  // derive the x-interval from the beam boundaries
  x_interv = get_xinterv_from_beam(actbeam);
  //  fprintf(stdout, "xpos: %f, ypos: %f\n", x_interv.x, x_interv.y);

  // allocate and fill the parameters
  tpars = (tlength_pars *) malloc(sizeof(tlength_pars));
  tpars->actbeam = actbeam;
  tpars->tlength = tlength;

  // fille the GSL-function
  F.function = &zero_tlength;
  F.params =   tpars;

  // set the boundaries for the solver
  gsl_root_fsolver_set (s, &F, x_interv.x, x_interv.y);


  // iterate to find the zeropoint
  do
    {
      // increment the iteration counter
      iter++;

      // iterate on the solver
      status = gsl_root_fsolver_iterate (s);

      // get a new guess from the solver
      xdiff = gsl_root_fsolver_root (s);

      // derive and set new boundaries
      x_interv.x = gsl_root_fsolver_x_lower (s);
      x_interv.y = gsl_root_fsolver_x_upper (s);

      // check the accuracy
      status = gsl_root_test_interval (x_interv.x, x_interv.y,
				       0, 0.0001);
      //--------------CAVEAT---------------------------------------
      // until March 28th 08 the code, wriongly was:
      //status = gsl_root_test_interval (x_interv.x, x_interv.x,
      //                                  0, 0.0001);
      // somehow this made no difference.....
      //-----------------------------------------------------------
    }
  // check for the break condition
  while (status == GSL_CONTINUE && iter < MAX_ITER);

  // free the memory
  free(tpars);

  // free the memory
  gsl_root_fsolver_free (s);

  // return the result
  return xdiff;
}

/**
 * Function: zero_tlength
 * Helper function for the 1D-root finding routine.
 * Evaluates a functional value at the position given
 * in the first parameter using the values given
 * as second parameter.
 *
 * Parameters:
 * @param x      - the independent value
 * @param params - ingredients to evaluate the function value
 *
 * Returns:
 * @return tzero - the function value
 */
double
zero_tlength (double x, void *params)
{
  double tzero;

  // extract the components from the parameter-structure
  beam  *actbeam = ((tlength_pars *) params)->actbeam;
  double tlength = ((tlength_pars *) params)->tlength;

  // compute the function value
  tzero = actbeam->spec_trace->path_len (x, actbeam->spec_trace->data) - tlength;

  // return the function value
  return tzero;
}


/**
 * Function: get_xinterv_from_beam
 * Determines the maximum and minimum extension of a given
 * beam along the x-axis. For security the interval is extended
 * at both ends by a fixed amount (quantified in the header-file).
 *
 * Parameters:
 * @param actbeam - the beam
 *
 * Returns:
 * @return x_interv - xmin, xmax covered by the beam
 */
d_point
get_xinterv_from_beam(const beam *actbeam)
{
  d_point x_interv;

  // compute the beam extension
  // towards the negative x-axis
  x_interv.x = (double)MIN(actbeam->corners[0].x,
			   MIN(actbeam->corners[1].x,
			       MIN(actbeam->corners[2].x,
				   actbeam->corners[3].x)))-actbeam->refpoint.x - INTERVEXT;

  // compute the beam extension
  // towards the positive x-axis
  x_interv.y = (double)MAX(actbeam->corners[0].x,
			   MAX(actbeam->corners[1].x,
			       MAX(actbeam->corners[2].x,
				   actbeam->corners[3].x)))-actbeam->refpoint.x + INTERVEXT;

  // return the interval
  return x_interv;
}


/**
 * Function: get_tlength_from_dispersion
 * The function computes and returnes the trace length for
 * a given diseprsion solution and wavelength. The tracelength
 * value is returned. Currently only linear and quadratic
 * dispersio solution are considered.
 *
 * Parameters:
 * @param cdisp  - the dispersion coefficient polynomial
 * @param lambda - the wavelength
 *
 * Returns:
 * @return tlength - the trace length
 */
double
get_tlength_from_dispersion(gsl_vector *cdisp, float lambda)
{
  double tlength=0.0;
  //double tlength_a=0.0;
  double det;

  // check whether the polynomial is linear
  if (cdisp->size == 2)
    {
      // compute the tracelength for linear dispersion
      tlength = (lambda-gsl_vector_get(cdisp, 0))/gsl_vector_get(cdisp, 1);
    }
  // check whether the polynomial is quadratic
  else if (cdisp->size == 3)
    {
      // compute the determinante
      det = gsl_vector_get(cdisp, 1) * gsl_vector_get(cdisp, 1)
		 - 4.0 * gsl_vector_get(cdisp, 0) * gsl_vector_get(cdisp, 2);

      // gove error if det < 0.0
      if (det < 0.0)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "get_tlength_from_dispersion: Can not determine the tracelength since det < 0.0!\n");

      // compute the tracelength
      tlength   = (-gsl_vector_get(cdisp, 1)+sqrt(det))/(2.0*gsl_vector_get(cdisp, 2));
      //    tlength_a = (-gsl_vector_get(cdisp, 1)-sqrt(det))/(2.0*gsl_vector_get(cdisp, 2));
      //    fprintf(stdout, "plus: %f, minus: f\n", tlength, tlength_a);
    }
  else
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_tlength_from_dispersion: The dispersion solution has %i coefficients and can not be inverted!\n", cdisp->size);
    }

  // return the tracelength
  return tlength;
}


/**
 * Function: condense_dispersion
 * The function strips off leading zeroes in the high order terms
 * of the dispersion solution derived from the configuration file.
 * These usually become a problem in the inversion to derive the
 * tracelength.
 *
 * Parameters:
 * @param disp__input - the configuration file dispersion
 *
 * Returns:
 * @return cdisp - the effective dispersion relation
 */
gsl_vector *
condense_dispersion(gsl_vector *disp_input)
{
  int index;
  int i;

  gsl_vector *cdisp;

  // get the size of the incomingg dispersion
  // vector
  index = disp_input->size;

  // determine the true size by subtracting
  // high order zeros
  while (index > -1 && !gsl_vector_get(disp_input, index-1))
    index--;

  // allocate space for the new vector
  cdisp = gsl_vector_alloc(index);

  // tranasfer the relevant values
  // from the old to the new vector
  for (i=0; i < index; i++)
    gsl_vector_set(cdisp, i, gsl_vector_get(disp_input, i));

  // return the new dispersion vector
  return cdisp;
}

/**
 * Function: fitting_ipc_corr
 *
 * Parameters:
 * @param actbeam        - the beam to examine
 * @param conf_file_path - the full pathname too the configuration file
 * @param ipcorr         - the interpolator with the correction factor
 * @param SPC            - the full spectrum to correct
 * @param obs            - the grism image
 * @param data_matrix    - background subtracted grism image
 * @param ipc_file_path  - file name for the phase data
 *
 * Returns:
 * @return icorr - marker whether the correction was applied or not (1/0)
 */
int
fitting_ipc_corr(beam act_beam, char conf_file_path[], interpolator *ipcorr,
		 full_spectr *SPC, observation *obs, gsl_matrix *data_matrix,
		 char ipc_file_path[], beam *beam_ptr)
{
  aperture_conf *conf;

  d_point        fit_qual;

  int            icorr     = 0;

  trace_func    *tracefun  = act_beam.spec_trace;
  double        *tracedata =  (double *)(tracefun->data);

  // load the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // determine the phase data and fix the phase
  // by fitting a function to it
  fit_qual = kappa_sigma_klipp_ipc(conf, obs, data_matrix, act_beam, ipc_file_path);

  // report on the shift and the quality
  fprintf(stdout, "shift-difference: %e, quality: %e\n", fit_qual.x, fit_qual.y);

  // check whether the quality is good enough
  if (fit_qual.y < MAXQUALITY)
    {
      // apply a correction to the
      // reference point to achive the
      // right corrections
      act_beam.refpoint.y += fit_qual.x;
      beam_ptr->refpoint.y += fit_qual.x;
      beam_ptr->ignore = -100;

      // report the new y-shift
      fprintf(stdout, "new y-shift: %e\n", act_beam.refpoint.y - (double)(int)(act_beam.refpoint.y + tracedata[1]));

      // apply the sensitivity correction to the spectrum
      intpix_corr_beam(act_beam, conf_file_path, ipcorr, SPC);

      // set the correction marker
      icorr=1;
    }

  // free allocated memory
  free_aperture_conf(conf);

  // return the marker
  return icorr;
}

/**
 * Function: intpix_corr_beam
 * The functio corrects a full spectrum for the intrapixel sensitivity
 * variations.
 * For every spectral element its fractional y-value is determined,
 * and then the correction factor for this y-value is computed
 * and applied to the appropriate elements of the spectral bin.
 *
 * Parameters:
 * @param actbeam        - the beam to examine
 * @param conf_file_path - the full pathname too the configuration file
 * @param ipcorr         - the interpolator with the correction factor
 * @param SPC            - the full spectrum to correct
 *
 * Returns:
 * @return -
 */
void
intpix_corr_beam(beam actbeam, char conf_file_path[], interpolator *ipcorr,
		 full_spectr *SPC)
{
  int            index=0;
  int            for_grism=1;

  double         cfactor;

  double         lambda;

  gsl_vector     *cdisp;

  d_point        pixel;

  dispstruct     *beam_disp;

  aperture_conf  *conf;

  // load the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // check whether it is grism (for_grism=1)
  // or prism (for_grism=0) data
  // give an error if there is a prism solution
  for_grism = check_for_grism (conf_file_path, actbeam.ID);
  if (!for_grism)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "intpix_corr_beam: Only grism dispersion solution can be corrected.\n");

  // get the wavelength dispersion relation at
  // position "refpoint". conf->refx and conf->refy
  // are used at this point to allow for a non (0,0) centered
  // 2D field dependence.
  pixel.x = actbeam.refpoint.x - conf->refx;
  pixel.y = actbeam.refpoint.y - conf->refy;

  // derive the dispersion at the object position
  beam_disp = get_dispstruct_at_pos(conf_file_path, for_grism,
				    actbeam.ID, pixel);

  // skipp high order zeroes in the dispersion solution
  cdisp = condense_dispersion(beam_disp->pol);

  for (index=0; index < SPC->nelems; index++)
    {
      lambda = SPC->fgr_spec->spec[index].lambda_mean;

      if (!gsl_isnan (lambda) && lambda)
	{
	  cfactor = get_intpix_corr(&actbeam, cdisp, ipcorr, lambda);

    	  fprintf(stdout, "lambda: %e, factor: %e\n", lambda, cfactor);
	  /*
	  SPC->fgr_spec->spec[index].count  = SPC->fgr_spec->spec[index].count / cfactor;
	  SPC->fgr_spec->spec[index].error  = SPC->fgr_spec->spec[index].error / cfactor;
	  SPC->fgr_spec->spec[index].flux   = SPC->fgr_spec->spec[index].flux / cfactor;
	  SPC->fgr_spec->spec[index].ferror = SPC->fgr_spec->spec[index].ferror / cfactor;

	  SPC->bck_spec->spec[index].count  = SPC->bck_spec->spec[index].count / cfactor;
	  SPC->bck_spec->spec[index].error  = SPC->bck_spec->spec[index].error / cfactor;
	  SPC->bck_spec->spec[index].flux   = SPC->bck_spec->spec[index].flux / cfactor;
	  SPC->bck_spec->spec[index].ferror = SPC->bck_spec->spec[index].ferror / cfactor;

	  SPC->obj_spec->spec[index].count  = SPC->obj_spec->spec[index].count / cfactor;
	  SPC->obj_spec->spec[index].error  = SPC->obj_spec->spec[index].error / cfactor;
	  SPC->obj_spec->spec[index].flux   = SPC->obj_spec->spec[index].flux / cfactor;
	  SPC->obj_spec->spec[index].ferror = SPC->obj_spec->spec[index].ferror / cfactor;

	  this is the wqrong version!!
	  SPC->fgr_spec->spec[index].count  = SPC->fgr_spec->spec[index].count * cfactor;
	  SPC->fgr_spec->spec[index].error  = SPC->fgr_spec->spec[index].error * cfactor;
	  SPC->fgr_spec->spec[index].flux   = SPC->fgr_spec->spec[index].flux * cfactor;
	  SPC->fgr_spec->spec[index].ferror = SPC->fgr_spec->spec[index].ferror * cfactor;

	  SPC->bck_spec->spec[index].count  = SPC->bck_spec->spec[index].count * cfactor;
	  SPC->bck_spec->spec[index].error  = SPC->bck_spec->spec[index].error * cfactor;
	  SPC->bck_spec->spec[index].flux   = SPC->bck_spec->spec[index].flux * cfactor;
	  SPC->bck_spec->spec[index].ferror = SPC->bck_spec->spec[index].ferror * cfactor;

	  SPC->obj_spec->spec[index].count  = SPC->obj_spec->spec[index].count * cfactor;
	  SPC->obj_spec->spec[index].error  = SPC->obj_spec->spec[index].error * cfactor;
	  SPC->obj_spec->spec[index].flux   = SPC->obj_spec->spec[index].flux * cfactor;
	  SPC->obj_spec->spec[index].ferror = SPC->obj_spec->spec[index].ferror * cfactor;
	  */
	}
    }

  // free the configuration structure
  free_aperture_conf(conf);

  // free the dispersion struct
  free_dispstruct(beam_disp);

  // free the memory for the dispersion
  gsl_vector_free(cdisp);
}


/**
 * Function: get_ipc_lambdas
 * The function computes the wavelengths for a list of x-offsets from the
 * reference point of a certain beam. The wavelength values are
 * returned as a gsl-vector.
 *
 * Parameters:
 * @param actbeam        - the beam
 * @param conf_file_path - the full path to the aXe configuration file
 * @param xvalues        - the list of x-offsets
 *
 * Returns:
 * @return lambdas - the list of wavelength values
 */
gsl_vector *
get_ipc_lambdas(const beam actbeam, char conf_file_path[], gsl_vector *xvalues)
{
  aperture_conf  *conf;

  dispstruct     *disp;
  calib_function *wl_calib;

  gsl_vector *lambdas;

  d_point pixel;
  int for_grism;

  int i;

  // load the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // check whether it is grism (for_grism=1)
  // or prism (for_grism=0) data
  for_grism = check_for_grism (conf_file_path, actbeam.ID);

  // determine the referencee point position
  pixel.x = actbeam.refpoint.x - conf->refx;
  pixel.y = actbeam.refpoint.y - conf->refy;

  // determine the dispersion at the reference point
  disp = get_dispstruct_at_pos(conf_file_path, for_grism,
			       actbeam.ID,pixel);

  // make a calibration structure from the dispersion
  wl_calib = create_calib_from_gsl_vector(for_grism, disp->pol);

  // convert the x-values to tracelength-values
  abscissa_to_pathlength (actbeam.spec_trace, xvalues);

  // allocate memory for the wavelengths
  lambdas = gsl_vector_alloc(xvalues->size);

  // go over all tracelength values
  for (i=0; i < (int)xvalues->size; i++)
    {
      // determine and store the wavelength for each tracelength
      gsl_vector_set(lambdas, i,
		     wl_calib->func(gsl_vector_get(xvalues, i), wl_calib->order,
				    wl_calib->coeffs));
    }

  // free the configuration structure
  free_aperture_conf(conf);

  // free the memory in the calibration structure
  free_calib(wl_calib);

  // free the dispersion structure
  free_dispstruct(disp);

  // return the wavelengths
  return lambdas;
}


/**
 * Function: get_ipc_cvalues
 * The function computes the intra-pixel correction factors for a certain
 * beam on a set of trace positions specified by their x-offset from
 * the reference point.
 * For every x-offset position the trace positon and its fractional y-pixel
 * (in absolute coordinates) is evaluated. Then the correction factor is
 * determined using the input interpolator.
 *
 * Parameters:
 * @param actbeam - the beam to correct
 * @param ipcorr  - the correction values depending on fractional y-pixel
 * @param xvalues - the list x-offsets from the reference point
 *
 * Returns:
 * @return cvalues - the list of correction values
 */
gsl_vector *
get_ipc_cvalues(const beam actbeam, interpolator *ipcorr, gsl_vector *xvalues)
{
  gsl_vector *cvalues;
  double yfract=0.0;
  int i;

  // allocate memory
  cvalues = gsl_vector_alloc(xvalues->size);

  // go over all x-values
  for (i=0; i < (int)xvalues->size; i++)
    {
      // determine the y-fraction for the x-value
      yfract = get_yfract_for_xvalue(&actbeam, gsl_vector_get(xvalues, i));

      // detyermine and store the correction factor in the array
      gsl_vector_set(cvalues, i, eval_interp(ipcorr, yfract));
    }

  // return the array
  return cvalues;
}


/**
 * Function: get_ipc_xvalues
 * The function computes a list of x-offsets from the reference
 * point of a beam. the regularly space offsets span the range
 * of x-values covered by the pixels of the beam.
 *
 * Parameters:
 * @param actbeam  - the beam to compute the x-offsets for
 *
 * Returns:
 * @return xvalues - the list of x-offsets
 */
gsl_vector *
get_ipc_xvalues(const beam actbeam)
{
  gsl_vector *xvalues;
  d_point xrange;

  int npoints;
  int i;

  // determine the x-range covered by the beam
  xrange =  get_xinterv_from_beam(&actbeam);

  // determine the number of points within the x-range
  npoints = (xrange.y - xrange.x) / XSTEPSIZE + 1;

  // allocate and fill a proper vector
  // with the x-values
  xvalues = gsl_vector_alloc(npoints);
  for (i=0; i < npoints; i++)
    {
      gsl_vector_set(xvalues, i, xrange.x + (double)i * XSTEPSIZE);
    }

  // return the vector with the x-values
  return xvalues;
}


/**
 * Function: get_ipclambda
 * The function creates an interpolator for the intra-pixel correction
 * as a function of wavelength based on this correction as function
 * of fractional y-pixel and the full calibration information
 * on a beam.
 * The interpolator is created after stepping along the beam trace
 * and combining the wavelength values with the correction values
 * at the trace positions.
 *
 * Parameters:
 * @param actbeam        - the beam to find the correction values for
 * @param conf_file_path - the full path to the aXe cofiguration file
 * @param ipcorr         - the correction values as function of y-fraction
 *
 * Returns:
 * @return ipclambda     - the correction values as function of wavelength
 */
interpolator *
get_ipclambda(beam actbeam, char conf_file_path[],
	      interpolator *ipcorr)
{
  gsl_vector *xvalues;
  gsl_vector *cvalues;
  gsl_vector *lambdas;

  double *cv;
  double *lv;

  interpolator *ipclambda;

  int i=0;
  int j=0;

  // get the x values around the reference point
  xvalues = get_ipc_xvalues(actbeam);

  // get the correction factors for these x-values
  cvalues = get_ipc_cvalues(actbeam, ipcorr, xvalues);

  // get the wavelength at the correction factors
  lambdas = get_ipc_lambdas(actbeam, conf_file_path, xvalues);

  // allocate memory for the dependent values
  cv = (double *) malloc(cvalues->size * sizeof(double));
  if (!cv) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // allocate memory for the independent values
  lv = (double *) malloc(cvalues->size * sizeof(double));
  if (!lv) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  if (gsl_vector_get(lambdas, lambdas->size-1) > gsl_vector_get(lambdas, 0))
    {
      // transfer the values from the
      // gsl vectors to the c-vectors
      for (i = 0; i < (int)cvalues->size; i++)
	{
	  cv[i] = gsl_vector_get(cvalues, i);
	  lv[i] = gsl_vector_get(lambdas, i);
	}
    }
  else
    {
      // transfer the values from the
      // gsl vectors to the c-vectors
      // invert the order from the
      // gsl-vectors
      j = cvalues->size - 1;
      for (i = 0; i < (int)cvalues->size; i++)
	{
	  cv[j] = gsl_vector_get(cvalues, i);
	  lv[j] = gsl_vector_get(lambdas, i);
	  j--;
	}
    }

  // create the interpolator
  ipclambda = create_interp(cvalues->size, FILTER_INTERP_TYPE, lv, cv);


  // free the memory allocated to
  // the vectors
  gsl_vector_free(xvalues);
  gsl_vector_free(cvalues);
  gsl_vector_free(lambdas);

  // return the interpolator
  return ipclambda;
}


/**
 * Function: apply_corr_pet
 * The function applies an intra-pixel sensitivity correction to
 * a list of PET pixels. The correction values are given depending
 * on the wavelength of the pixel, and are applied to the PET pixels
 * in place.
 *
 * Parameters:
 * @param ipclambda - the ipc correction factor as function of wavelength
 * @param PET       - the list of PET pixels to correct
 *
 * Returns:
 * -
 */
void
apply_corr_pet(interpolator *ipclambda, ap_pixel *PET)
{
  int j = 0;

  double cvalue=0.0;

  // go along the PET pixels until you meet
  // the last
  while (PET[j].p_x != -1)
    {
      // get the correction factor
      cvalue = eval_interp(ipclambda, PET[j].lambda);

      // apply the corection factor
      // to the counts and the error
      PET[j].count = PET[j].count / cvalue;
      PET[j].error = PET[j].error / cvalue;

      //      fprintf(stdout, "lambda: %f, correction: %f\n", PET[j].lambda, cvalue);
      j++;
    }

}

/**
 * Function: intpix_corr_pet
 * The function applies the intra-pixel sensitivity correction to
 * a list of PET pixels. The values in the pixels are corrected in
 * place using the correction function, gemoetrical parameters and
 * dispersion relation in the various parameters.
 *
 * Parameters:
 * @param actbeam        - the beam to correct
 * @param conf_file_path - full path to the axe configuration file
 * @param ipcorr         - the ipc correction function
 * @param PET            - the list of PET pixels to correct
 *
 * Returns:
 * -
 */
void
intpix_corr_pet(beam actbeam, char conf_file_path[],
		interpolator *ipcorr, ap_pixel *PET)
{
  interpolator *ipclambda;
  aperture_conf  *conf;

  int for_grism=0;

  // load the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // check whether it is grism (for_grism=1)
  // or prism (for_grism=0) data
  // give an error if there is a prism solution
  for_grism = check_for_grism (conf_file_path, actbeam.ID);
  if (!for_grism)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "intpix_corr_beam: Only grism dispersion solution can be corrected.\n");

  // determine the correction factor versus wavelength
  ipclambda = get_ipclambda(actbeam, conf_file_path, ipcorr);

  print_interp(ipclambda);

  // apply the correction to the PET
  apply_corr_pet(ipclambda, PET);

  // free the configuration structure
  free_aperture_conf(conf);

  // free the memory in the interpolator
  free_interp(ipclambda);
}


/**
 * Function: is_pointlike
 * The function evaluates all criteria to decide whether an object
 * is pointlike or not. In the current implementation an object
 * is considered pointlike if a special OAF is given and the
 * the beam is NOT excluded OR if the spatial extension of the
 * object is smaller than the maximal extension.
 *
 * Parameters:
 * @param actbeam   - the beam to examine
 * @param spec_OAF  - indicates a special OAF file
 * @param max_ext   - the maximal allowed extension
 *
 * Returns:
 * @return point_like - pointer to the opened fits file
 */
int
is_pointlike(beam actbeam, int spec_OAF, double max_ext)
{
  int point_like=0;

  int OAF_crit=0;
  int EXT_crit=0;

  // check whether a special OAF is given
  // and the beam should NOT be ignored
  if (spec_OAF && !actbeam.ignore)
    OAF_crit=1;

  // check whether the maximal extension
  // is given and the object extension
  // is smaller;
  // also check whether non-default values
  // for the object width are set
  if (actbeam.awidth > -1.0 && actbeam.bwidth > -1.0 && max_ext && actbeam.awidth < max_ext && actbeam.bwidth < max_ext)
    EXT_crit=1;

  // set to pointlike if at least
  // one of the criteria is set
  if (OAF_crit || EXT_crit)
    point_like=1;

  // return the result
  return point_like;
}


/**
 * Function: get_SPC_opened
 * The function opens an existing SPC file and returns the pointer
 * to it.
 * As of now, the mode of opening it is automatically READWRITE.
 * Later on a differentiation using the free parameter "mode"
 * might be added  to generalize the function.
 *
 * Parameters:
 * @param SPCname - name of the SPC file
 * @param mode    - mode to open it (not yet used)
 *
 * Returns:
 * @return SPC_ptr - pointer to the opened fits file
 *
 */
fitsfile *
get_SPC_opened(char SPCname[], int mode)
{
    fitsfile *SPC_ptr;
    int f_status=0;

  // Open the OPET file for reading/writing
  fits_open_file (&SPC_ptr, SPCname, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_INTPIXCORR: Could not open file: %s\n",
		   SPCname);
    }

  // return the pointer to the fits file
  return SPC_ptr;
}


/**
 * Function: get_ALL_from_next_in_SPC
 * The function creates and fills a full spectrum structure with
 * the content of a SPC table extension. This is only done when,
 * according to the beam ID, this beam should be corrected.
 * For extensions with higher order beams which are not corrected
 * an emply structure is returned.
 *
 * Parameters:
 * @param SPCn_ptr - pointer to the opened SPC file
 * @param aperID   - pointer to aperture identification number
 * @param beamID   - pointer to beam identification number
 *
 * Returns:
 * @return SPC - the full spectrum structure
 */
full_spectr *
get_ALL_from_next_in_SPC(fitsfile *SPC_ptr, int *aperID, int *beamID)
{
  int f_status=0, hdutype;

  long tmp;
  //long nrows=0;
  char comment[FLEN_COMMENT];

  full_spectr *SPC;


  fits_movrel_hdu (SPC_ptr, 1, &hdutype, &f_status);

  if (f_status)
    {
      *aperID = -1;
      *beamID = -1;
      SPC = NULL;
      return SPC;
    }


  // read the beam ID number
  fits_read_key_lng (SPC_ptr, "BEAMID", &tmp, comment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_ALL_from_next_in_SPC: Error getting index keyword OBJECTID");
    }
  *beamID = (int)tmp;

  // check whether this beam should be correct
  if (*beamID > CORRMAX)
    {
      // set it to NULL and return
      SPC = NULL;
      //      fprintf (stdout, "aXe_PETFF: Skipping beam: %c.\n", BEAM(*beamID));
      return SPC;
    }

  // the beam shall be corrected and
  // first must be read in
  SPC = (full_spectr *) malloc (sizeof (full_spectr ));

  // transfer the beam ID
  SPC->beamID = (int)tmp;

  // read the aperture number
  fits_read_key_lng (SPC_ptr, "OBJECTID", &tmp, comment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_ALL_from_next_in_SPC: Error getting index keyword OBJECTID");
    }
  // transfer the aperture ID
  *aperID = (int)tmp;
  SPC->aperID = (int)tmp;


  // Get the number of rows
  fits_get_num_rows (SPC_ptr, &tmp, &f_status);
  if (f_status) {
    ffrprt (stderr, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "get_ALL_from_next_in_SPC: "
		 "Could not determine the number of rows in"
		 " correction function table!");
  }
  SPC->nelems = (int)tmp;

  // load the background subtracted object spectrum
  SPC->obj_spec = get_spectrum_from_SPC(SPC_ptr, "COUNT", "ERROR", SPC->nelems);

  // load the total object spectrum
  SPC->fgr_spec = get_spectrum_from_SPC(SPC_ptr, "TCOUNT", "TERROR", SPC->nelems);

  // load the background spectrum
  SPC->bck_spec = get_spectrum_from_SPC(SPC_ptr, "BCOUNT", "BERROR", SPC->nelems);

  // return the filled structure
  return SPC;
}


/**
 * Function: free_full_spectr
 * The function frees the memory allocated in a full
 * spectrum structure.
 *
 * Parameters:
 * @param SPC - the full spectrum structure
 *
 * Returns:
 * @return -
 */
void
free_full_spectr(full_spectr *SPC)
{
  // free the three spectra in the
  // full spectrum structure
  free_spectrum(SPC->obj_spec);
  free_spectrum(SPC->fgr_spec);
  free_spectrum(SPC->bck_spec);

  // free the full spectrum
  free(SPC);

  // set the structure to NULL
  SPC = NULL;
}


/**
 * Function: get_spectrum_from_SPC
 * The function fills a spectrum structure with data in an SPC extension.
 * There is more data in an SPC extension than can be stored in a
 * spectrum structure. Two parameters in this function specify partly which
 * data should be loaded.
 *
 * Parameters:
 * @param SPC_ptr   - pointer to the SPC extension
 * @param count_col - column to load for 'count'
 * @param error_col - column to load for 'error'
 * @param nelems    - the number ol elements in the spectrum
 *
 * Returns:
 * @return act_spec - the filled spectrum structure
 */
spectrum *
get_spectrum_from_SPC(fitsfile *SPC_ptr, char count_col[],
		      char error_col[], int nelems)
{
  int index=0;

  spectrum *act_spec;

  double *count;
  double *lambda;
  double *error;
  double *flux;
  double *ferror;
  double *weight;
  double *contam;
  long   *dq;

  // allocate the spectrum
  act_spec = allocate_spectrum (nelems);

  // transfer the column entries into a vector
  lambda = get_dcolumn_from_SPC_opened(SPC_ptr, "LAMBDA", nelems);
  count  = get_dcolumn_from_SPC_opened(SPC_ptr, count_col, nelems);
  error  = get_dcolumn_from_SPC_opened(SPC_ptr, error_col, nelems);
  flux   = get_dcolumn_from_SPC_opened(SPC_ptr, "FLUX", nelems);
  ferror = get_dcolumn_from_SPC_opened(SPC_ptr, "FERROR", nelems);
  weight = get_dcolumn_from_SPC_opened(SPC_ptr, "WEIGHT", nelems);
  contam = get_dcolumn_from_SPC_opened(SPC_ptr, "CONTAM", nelems);
  dq     = get_lcolumn_from_SPC_opened(SPC_ptr, "DQ", nelems);

  // transfer the data from the vector into the
  // spectrum structure
  for (index = 0; index < nelems; index++)
    {
      act_spec->spec[index].lambda_mean = lambda[index];
      act_spec->spec[index].count       = count[index];
      act_spec->spec[index].error       = error[index];
      act_spec->spec[index].flux        = flux[index];
      act_spec->spec[index].ferror      = ferror[index];
      act_spec->spec[index].weight      = weight[index];
      act_spec->spec[index].contam      = contam[index];
      act_spec->spec[index].dq          = dq[index];
    }

  // free the arrays
  free(lambda);
  free(count);
  free(error);
  free(flux);
  free(ferror);
  free(weight);
  free(contam);
  free(dq);

  // return the spectrum
  return act_spec;
}


/**
 * Function: get_dcolumn_from_SPC_opened
 * The function reads the values from a double column specified
 * by its name into a array of doubles. The array is returned.
 *
 * Parameters:
 * @param SPC_ptr   - pointer to the SPC extension
 * @param count_col - name of the column to load
 * @param nelems    - number of elements in the column
 *
 * Returns:
 * @return values - pointer to a filled array
 */
double *
get_dcolumn_from_SPC_opened(fitsfile *SPC_ptr, char colname[], int nelems)
{
  int colnum=0;
  int anynul;
  int f_status=0;
  double *values;

  // allocate the return array;
  // give an error if allocation fails
  values = (double *) malloc(nelems * sizeof(double));
  if (!values) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // get the desired column number;
  // give an error if the column name
  // can not be read
  fits_get_colnum (SPC_ptr, CASEINSEN, colname, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_dcolumn_from_SPC_opened: "
		   "Could not determine column %s in "
		   " table!\n", colname);
    }

  // read all data in the column
  fits_read_col (SPC_ptr, TDOUBLE, colnum, 1, 1, nelems, NULL, values,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_dcolumn_from_SPC_opened: "
		   "Could not read column %s"
		   " from BINARY table!",colname);
    }

  // return the filled vector
  return values;
}


/**
 * Function: get_lcolumn_from_SPC_opened
 * The function reads the values from a column with long specified
 * by its name into a array of type long. The array is returned.
 *
 * Parameters:
 * @param SPC_ptr   - pointer to the SPC extension
 * @param count_col - name of the column to load
 * @param nelems    - number of elements in the column
 *
 * Returns:
 * @return values - pointer to a filled array
 */
long *
get_lcolumn_from_SPC_opened(fitsfile *SPC_ptr, char colname[], int nelems)
{
  int colnum=0;
  int anynul;
  int f_status=0;
  long *values;

  // allocate the return array;
  // give an error if allocation fails
  values = (long *) malloc(nelems * sizeof(long));
  if (!values) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // get the desired column number;
  // give an error if the column name
  // can not be read
  fits_get_colnum (SPC_ptr, CASEINSEN, colname, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_dcolumn_from_SPC_opened: "
		   "Could not determine column %s in "
		   " table!\n", colname);
    }

  // read all data in the column
  fits_read_col (SPC_ptr, TLONG, colnum, 1, 1, nelems, NULL, values,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_dcolumn_from_SPC_opened: "
		   "Could not read column %s"
		   " from BINARY table!",colname);
    }

  // return the filled array
  return values;
}


/**
 * Function: create_nlincor
 * The function creates an interpolator for the nonlinearity
 * correction applied to NICMOS data.
 *
 * Parameters:
 *
 * Returns:
 * @return nlincorr - the interpolator created
 */
interpolator *
create_nlincor()
{
  interpolator *nlincorr;
  /*
  double x[14] = {8250.0, 8750.0, 9250.0, 9750.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0};
  double y[14] = {  .069,   .057,   .052,   .050,    .049,    .048,      .041,    .023,    .013,    .008,    .004,    .0,      .0,      .0};
  */
  double *xx;
  double *yy;

  xx = (double *) malloc(14 * sizeof(double));
  yy = (double *) malloc(14 * sizeof(double));

  xx[0]  = 8250.0;
  xx[1]  = 8750.0;
  xx[2]  = 9250.0;
  xx[3]  = 9750.0;
  xx[4]  = 11000.0;
  xx[5]  = 12000.0;
  xx[6]  = 13000.0;
  xx[7]  = 14000.0;
  xx[8]  = 15000.0;
  xx[9]  = 16000.0;
  xx[10] = 17000.0;
  xx[11] = 18000.0;
  xx[12] = 19000.0;
  xx[13] = 20000.0;

  yy[0]  = .069;
  yy[1]  = .057;
  yy[2]  = .052;
  yy[3]  = .050;
  yy[4]  = .049;
  yy[5]  = .048;
  yy[6]  = .041;
  yy[7]  = .023;
  yy[8]  = .013;
  yy[9]  = .008;
  yy[10] = .004;
  yy[11] = .0;
  yy[12] = .0;
  yy[13] = .0;

  // create the interpolator
  nlincorr = create_interp(14, NLINCORR_INTERP_TYPE, xx, yy);

  // return the interpolator
  return nlincorr;
}

/**
 * Function: nlin_corr_beam
 *
 * Parameters:
 * @param nlincorr       - the interpolator with the correction factor
 * @param SPC            - the full spectrum to correct
 *
 * Returns:
 * @return -
 */
void
nlin_corr_beam(interpolator *nlincorr, double adcgain, full_spectr *SPC)
{
  int            index=0;
  //  int            for_grism=1;

  double         cfactor;
  double         cps;

  double         lambda;

  for (index=0; index < SPC->nelems; index++)
    {
      // get the independent spectral values,
      // the wavelenth and the cps value
      lambda = SPC->obj_spec->spec[index].lambda_mean;
      cps    = SPC->obj_spec->spec[index].count;


      // check whether the spectral element
      // is not corrupt
      if (!gsl_isnan (lambda) && lambda)
	{
	  if (cps > 0.0)
	    // compute the correction factor
	    cfactor = get_nlin_corr(nlincorr, lambda, cps/adcgain);
	  else
	    // make a dummy factor
	    cfactor = 1.0;
	  // correct what should be corrected
	  SPC->fgr_spec->spec[index].count  = SPC->fgr_spec->spec[index].count / cfactor;
	  SPC->fgr_spec->spec[index].error  = SPC->fgr_spec->spec[index].error / cfactor;
	  SPC->fgr_spec->spec[index].flux   = SPC->fgr_spec->spec[index].flux / cfactor;
	  SPC->fgr_spec->spec[index].ferror = SPC->fgr_spec->spec[index].ferror / cfactor;

	  SPC->bck_spec->spec[index].count  = SPC->bck_spec->spec[index].count / cfactor;
	  SPC->bck_spec->spec[index].error  = SPC->bck_spec->spec[index].error / cfactor;
	  SPC->bck_spec->spec[index].flux   = SPC->bck_spec->spec[index].flux / cfactor;
	  SPC->bck_spec->spec[index].ferror = SPC->bck_spec->spec[index].ferror / cfactor;

	  SPC->obj_spec->spec[index].count  = SPC->obj_spec->spec[index].count / cfactor;
	  SPC->obj_spec->spec[index].error  = SPC->obj_spec->spec[index].error / cfactor;
	  SPC->obj_spec->spec[index].flux   = SPC->obj_spec->spec[index].flux / cfactor;
	  SPC->obj_spec->spec[index].ferror = SPC->obj_spec->spec[index].ferror / cfactor;
	}
    }

}

/*
 * Function: get_nlin_corr
 *
 * Parameters:
 * @param nlincorr - the interpolator with the parameter
 * @param lamda    - the full spectrum to correct
 * @param cps      - the count rate
 *
 * Returns:
 * @return cfactor - the correction factor
 */
double
get_nlin_corr(interpolator *nlincorr, const double lambda, const double cps)
{
  double cfactor;
  double bbb;

  // evaluate the parameter
  bbb = eval_interp(nlincorr, lambda);

  // compute the correction factor
  cfactor = 1.0 - 2.0 * bbb + bbb * log10(cps);

  // return the correction factor
  return cfactor;
}
