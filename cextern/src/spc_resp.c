/**
 * Set of routines to handle absolute flux calibration of a spectrum
 * structure.
 * The response curves can be in ASCII or (preferably) in FITS
 * binary format (first extension). ASCI format is simply three
 * space separated colums (WAVELENGTH, THROUGHPUT, ERROR). The FITS
 * binary table must contain three columns with names WAVELENGTH,
 * SENSITIVITY, ERROR. The throughput should be DN/s per Erg/cm^2/s/A.
 */


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

#include "aXe_grism.h"
#include "disp_conf.h"
#include "fringe_conf.h"
#include "spc_wl_calib.h"
#include "spc_resp.h"

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/**
 * Function: get_response_function_from_FITS
 * This function reads a FITS table, reads the WAVELENGTH,
 * SENSITIVITY, and ERROR columns and returns a spectrum
 * structure which has been populated with this throughput.
 *
 * Parameters:
 * @param filename - a pointer to a char array containing the
 *                   name of an existing FITS file
 * @param hdunum   - the extension number to read (2 for STSDAS tables)
 *
 * Returns:
 * @return res     - a pointer to a newly allocated spectrum structure
 */
spectrum *get_response_function_from_FITS(char filename[], int hdunum)
{
  fitsfile *input;
  int f_status = 0;
  int hdutype, anynul;
  double *lambda, *resp, *error;
  double tmp;
  long numrows;
  int colnum, i;
  spectrum *res;

  //  Open the file for reading
  fits_open_file (&input, filename, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not open" " file: %s",
                   filename);
    }
  /* Move to first hdu */
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not read extention %d from file: %s",
                   hdunum, filename);
    }

  /* Get number of rows */
  fits_get_num_rows (input, &numrows, &f_status);
  if (f_status) {
    ffrprt (stderr, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "get_response_function_fromFITS: "
                 "Could not determine the number of rows in"
                 " table %s",filename);
  }

  /* Allocate temporary memory space */
  lambda = (double *) malloc(numrows*sizeof(double));
  if (!lambda) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                 "Memory allocation failed");
  }
  resp = (double *) malloc(numrows*sizeof(double));
  if (!resp) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                 "Memory allocation failed");
  }
  error = (double *) malloc(numrows*sizeof(double));
  if (!error) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                 "Memory allocation failed");
  }

  /**************************/
  /* Read the WAVELENGTH column */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, "WAVELENGTH", &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not determine WAVELENGTH column number in "
                   " table %s",filename);
    }

  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, numrows, NULL, lambda,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not read content of WAVELENGTH column "
                   " from BINARY table %s",filename);
    }

  /**************************/
  /* Read the SENSITIVITY column */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, "SENSITIVITY", &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not determine SENSITIVITY column number in "
                   " table %s",filename);
    }
  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, numrows, NULL, resp,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not read content of SENSITIVITY column "
                   " from BINARY table %s",filename);
    }
  /**************************/
  /* Read the ERROR column */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, "ERROR", &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not determine ERROR column number in "
                   " table %s",filename);
    }
  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, numrows, NULL, error,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_response_function_fromFITS: "
                   "Could not read content of ERROR column "
                   " from BINARY table %s",filename);
    }


  /* Allocate memory */
  res = (spectrum *) malloc(sizeof(spectrum));
  if(!res) {
      aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                   "Memory allocation failed");
  }
  res->spec = (spc_entry *) malloc(numrows*sizeof(spc_entry));
  if(!(res->spec)) {
      aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                   "Memory allocation failed");
  }

  for (i=0;i<numrows;i++) {
    res->spec[i].lambda_mean = lambda[i];
    res->spec[i].flux = resp[i];
    res->spec[i].ferror = error[i];
  }

  res->spec_len=i;

  // define the minimum and maximum wavelength
  res->lambdamin=res->spec[0].lambda_mean;
  res->lambdamax=res->spec[numrows-1].lambda_mean;

  // check whether the wavelength is raising.
  // switch minimum and maximum for falling spectrum
  if (res->lambdamin > res->lambdamax)
    {
      tmp = res->lambdamin;
      res->lambdamin = res->lambdamax;
      res->lambdamax = tmp;
    }

  fits_close_file(input,&f_status);
  if (f_status) {
      aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
                   "Could not close %s",filename);
  }

  free(lambda);
  free(resp);
  free(error);

  return res;
}

/**
 * Function: apply_response_function
 * This function applies a thoughput curve to a given spectrum
 * by dividing the count attribute of the spectrum. The flux
 * and ferror attributes of the spectrum are populated.
 * gsl spline interpolation is used to compute the throughput
 * curve at the needed spectral wavelengths.
 *
 * Parameters:
 * @param spec    - a pointer to an existing spectrum array
 * @param resp    - a pointer to a spectrum structure containing
 *                  a throughput curve previously loaded
 * @param exptime - the exposure time
 * @param gain    - the gain
 *
 * Returns
 * @return -
 */
void
apply_response_function(spectrum *spec, spectrum *resp, const int quant_cont)
{
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, resp->spec_len);

  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, resp->spec_len);

  double *x1,*y1, *x2, *y2;
  int i,j;
  double r1, r2, fr, fr1, fr2;
  double min_l, max_l;
  // test
  // int nguess;
  // double tval=0.0;

  if ((int) resp->spec_len <= 0){
    aXe_message(aXe_M_ERROR, __FILE__, __LINE__,
                "apply_response_function: response spec length invalid: %i", resp->spec_len);
      gsl_spline_free(spline1);
      gsl_interp_accel_free(acc1);
      gsl_spline_free(spline2);
      gsl_interp_accel_free(acc2);
      return;
  }
  if ( spec==NULL)
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                   "apply_response_function: spectra empty.");
      gsl_spline_free(spline1);
      gsl_interp_accel_free(acc1);
      gsl_spline_free(spline2);
      gsl_interp_accel_free(acc2);
      return;
    }

  x1 = malloc((double)(resp->spec_len) *sizeof(double));
  y1 = malloc((double)(resp->spec_len) *sizeof(double));
  x2 = malloc((double)(resp->spec_len) *sizeof(double));
  y2 = malloc((double)(resp->spec_len) *sizeof(double));

  min_l = 1e32;
  max_l = -1e32;

  for (i=0;i<resp->spec_len; i++) {
    if (resp->spec[i].lambda_mean > max_l) {
     max_l = resp->spec[i].lambda_mean;
    }
    if (resp->spec[i].lambda_mean < min_l) {
     min_l = resp->spec[i].lambda_mean;
    }
     
    x1[i] = resp->spec[i].lambda_mean;
    y1[i] = resp->spec[i].flux - resp->spec[i].ferror;
    x2[i] = resp->spec[i].lambda_mean;
    y2[i] = resp->spec[i].flux + resp->spec[i].ferror;
  }

  gsl_spline_init (spline1, x1, y1, resp->spec_len);
  gsl_spline_init (spline2, x2, y2, resp->spec_len);
  // test
  // nguess = 0;
  
  for (j=0; j<spec->spec_len; j++) {
    //fprintf(stderr,"%f %f %f\n",spec->spec[j].lambda_mean, min_l, max_l);
    if ((spec->spec[j].lambda_mean>=min_l) && (spec->spec[j].lambda_mean <= max_l)) {
    //     fprintf(stderr,"lambda mean between min and max\n\t%f : %f %f\n", spec->spec[j].lambda_mean, min_l, max_l);
      r1 = gsl_spline_eval(spline1, spec->spec[j].lambda_mean, acc1);
      r2 = gsl_spline_eval(spline2, spec->spec[j].lambda_mean, acc2);
   }
   else {
     r1 = 0;
     r2 = 0;
   }

    /* Need to divide by the pixel width of the
       wavelength calibration at the given lambda */

    if ((r2+r1)!=0) {

      // test
      // spec->spec[j].flux = spec->spec[j].count / tval;
      spec->spec[j].flux = spec->spec[j].count / ((r2+r1)/2.);
      spec->spec[j].flux = spec->spec[j].flux/spec->spec[j].dlambda;
      // fprintf(stdout, "# %f %f %f\n", spec->spec[j].lambda_mean, spec->spec[j].dlambda, spec->spec[j].weight);


      if (quant_cont && (int)spec->spec[j].contam != -1)
        {
          // test
          //      spec->spec[j].contam = spec->spec[j].contam / tval;
          spec->spec[j].contam = spec->spec[j].contam / ((r2+r1)/2.);
          spec->spec[j].contam = spec->spec[j].contam/spec->spec[j].dlambda;
        }
    } else {
      spec->spec[j].flux = GSL_NAN;
    }

    fr1 = spec->spec[j].error / spec->spec[j].count; /* frac. error */
    fr2 = fabs((r2-r1)/((r1+r2)/2.0)); /* frac. error in resp. */
    fr = sqrt(fr1*fr1+fr2*fr2);
    spec->spec[j].ferror = fabs(spec->spec[j].flux) * fr;
  }

  gsl_spline_free(spline1);
  gsl_interp_accel_free(acc1);
  gsl_spline_free(spline2);
  gsl_interp_accel_free(acc2);
  free(x1);
  free(x2);
  free(y1);
  free(y2);

}

/**
 * Function: apply_response_functionII
 * This function applies a thoughput curve to a given spectrum
 * by dividing the count attribute of the spectrum. The flux
 * and ferror attributes of the spectrum are populated.
 *
 * Parameters:
 * @param spec       - a pointer to an existing spectrum array
 * @param resp_func  - a pointer to a spectrum structure containing
 *                     a throughput curve previously loaded
 * @param quant_cont - the exposure time
 *
 * Returns
 * @return -
 */
void
apply_response_functionII(spectrum *spec, response_function *resp_func, const int quant_cont)
{
  long j;
  double fr, fr1, fr2;

  double *resp_vals = malloc( 2 * sizeof(double));

  // immediately return
  // an empty spectrum
  if ( spec==NULL)
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                   "apply_response_function: spectra empty.");
      return;
    }

  // go over all elements
  for (j=0;j<spec->spec_len;j++)
    {
      // get the sensitivity value plus error
      get_response_values(resp_func, spec->spec[j].lambda_mean, resp_vals);

      // make sure you can compute the flux value
      if (resp_vals[0] != 0.0)
        {
          // compute the flux value
          spec->spec[j].flux = spec->spec[j].count / resp_vals[0];
          spec->spec[j].flux = spec->spec[j].flux/spec->spec[j].dlambda;

          // compute the flux error
          fr1 = spec->spec[j].error / spec->spec[j].count; /* frac. error */
          fr2 = fabs(resp_vals[1] / resp_vals[0]); /* frac. error in resp. */
          fr = sqrt(fr1*fr1+fr2*fr2);
          spec->spec[j].ferror = fabs(spec->spec[j].flux) * fr;

          // check whether the contamination value is necessary
          if (quant_cont && (int)spec->spec[j].contam != -1)
            {
              // compute the contamination value
              spec->spec[j].contam = spec->spec[j].contam / resp_vals[0];
              spec->spec[j].contam = spec->spec[j].contam/spec->spec[j].dlambda;
            }
        }
      else
        {
          // set the flux and the error to NaN
          spec->spec[j].flux = GSL_NAN;
          spec->spec[j].ferror = GSL_NAN;
        }
    }

  // free the memory
  free(resp_vals);
}

/**
 * Function: get_smooth_pars_for_beam
 * The function determines the smoothing parameters for an object beam
 * and a given indicator. If smoothed calibration is requested, the source
 * size (in dispersion direction) is compared to the poin-like object size.
 * For non poin-like objects, the adjustment factor is extracted and
 * the smoothing size from the object is derived.
 *
 * Parameters:
 * @param conf        - the configuarion file structure
 * @param smooth_conv - boolean for smooth conversion
 * @param actbeam     - the beam to be calibrated
 *
 * Returns
 * @return smooth_pars - the paramters used in the smoothing
 */
d_point
get_smooth_pars_for_beam(const aperture_conf *conf, const int smooth_conv, beam actbeam)
  {
    d_point smooth_pars;

    // initialize the return
    smooth_pars.x = -1.0;
    smooth_pars.y = -1.0;

    // if there is nothing to
    // do, return the default
    if (!smooth_conv)
      return smooth_pars;

    // check whether the relevant information
    // does exist
    if (smooth_conv && actbeam.slitgeom[2] < 0.0)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
          "Smoothed flux conversion is impossible"
          " since the OAF\ndoes not contain the slit width!\n");

    //fprintf(stdout, "smoothing factor: %f, point-like size: %f\n", conf->smfactor, conf->pobjsize);
    // check whether the objects is smaller than a psf-object
    // in dispersion direction
    if (conf->pobjsize < actbeam.slitgeom[2])
      {
        // compute the smoothing size from the objects;
        // transfer the adjustment factor
        smooth_pars.x = pow((actbeam.slitgeom[2]* actbeam.slitgeom[2]) - (conf->pobjsize*conf->pobjsize), 0.5);
        smooth_pars.y = conf->smfactor;
      }

    // return the
    // smoothing parameters
    return smooth_pars;
  }

/**
 * Function: check_conf_for_smoothing
 * The functions checks whether a smoothed flux conversion
 * is possible or not. In case that keywords in the configuration
 * files are missing, it is NOT possible, and 0 is returned.
 *
 * Parameters:
 * @param conf        - the configuarion file structure
 * @param smooth_conv - boolean for smooth conversion
 *
 * Returns
 * @return is_possible - the paramters used in the smoothing
 */
int
check_conf_for_smoothing(const aperture_conf *conf, const int smooth_conv)
  {
    int is_possible=1;

    // make sure that, if smoothed covnersion is requested,
    // a point-like objects size and an adjustment factor is defined
    if (smooth_conv && (conf->pobjsize < 0.0 || conf->smfactor < 0.0))
      // change the switch
      is_possible = 0;

    // return the pointer
    return is_possible;
  }

/**
 * Function: apply_smoothed_response
 * Convert a spectrum from e/s to flux units in f_lambda. The conversion is done
 * using smoothed sensitivity values that take the object width into account
 * by smoothing the sensitivity function. Gaussian smoothing is applied, taking
 * into account that the original calibration function was derived for a point
 * like object with a finite size.
 *
 * Parameters:
 * @param spec           - a pointer to an existing spectrum array
 * @param resp_func      - a pointer to a spectrum structure containing
 *                         a throughput curve
 * @param quant_cont     - indicator for quantitative contamination
 * @param bwidth         - the B_IMAGE value of the beam
 * @param psf_ext        - the width of a point-like object
 * @param wl_calibration - the wavelength calibration structure
 *
 * Returns
 * @return -
 */
void
apply_smoothed_response(const calib_function *wl_calibration, const int for_grism,
                        const int quant_cont, response_function *resp_func,
                        const d_point smooth_pars, spectrum *spec)
{
  long j;
  double fr, fr1, fr2;
  //double sigma_wav=10.0;

  double *resp_vals = malloc( 2 * sizeof(double));
  //double *tmp_vals  = malloc( 2 * sizeof(double));

  gsl_vector *weights;

  // allocate memory for the vectors
  weights = gsl_vector_alloc(2 * RESP_SMOOTH_LENGTH + 1);

  // immediately return
  // if an empty spectrum
  if ( (spec==NULL) || (spec->spec_len == 0))
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                   "apply_response_function: spectra empty.");
      return;
    }

  // fill the weight vector
  // for smoothing
  fill_weight(weights);

  // output to screen
  fprintf(stdout, "aXe_PET2SPC: Gaussian smoothing in tracelength: %.2f * %.2f = %.2f pix\n", smooth_pars.x, smooth_pars.y, smooth_pars.x * smooth_pars.y);

  // go over all elements
  for (j=0;j<spec->spec_len;j++)
    {
      // get the smoothed sensitivity value plus error
      get_smoothed_response(spec->spec[j].lambda_mean, smooth_pars, wl_calibration,
                            for_grism, weights, resp_func, resp_vals);
      // make sure you can compute the flux value
      if (resp_vals[0] != 0.0)
        {
          // compute the flux value
          spec->spec[j].flux = spec->spec[j].count / resp_vals[0];
          spec->spec[j].flux = spec->spec[j].flux/spec->spec[j].dlambda;

          // compute the flux error
          fr1 = spec->spec[j].error / spec->spec[j].count; /* frac. error */
          fr2 = fabs(resp_vals[1] / resp_vals[0]); /* frac. error in resp. */
          fr = sqrt(fr1*fr1+fr2*fr2);
          spec->spec[j].ferror = fabs(spec->spec[j].flux) * fr;

          // check whether the contamination value is necessary
          if (quant_cont && (int)spec->spec[j].contam != -1)
            {
              // compute the contamination value
              spec->spec[j].contam = spec->spec[j].contam / resp_vals[0];
              spec->spec[j].contam = spec->spec[j].contam/spec->spec[j].dlambda;
            }
        }
      else
        {
          // set the flux and the error to NaN
          spec->spec[j].flux = GSL_NAN;
          spec->spec[j].ferror = GSL_NAN;
        }
    }

  // free the memory
  free(resp_vals);
  gsl_vector_free(weights);
}

/**
 * Function: find_wavelength
 * For a given trace-length value, the function returns the difference
 * between the actual and the targeted wavelength.
 * The function is designed t be called by a gsl root-solving procedure.
 *
 * Parameters:
 * @param x        - guess value for the tracelength
 * @param params   - structure with input parameters
 *
 * Returns
 * @return wav_zero - difference between actual and targeted wavelength
 */
double
find_wavelength(double x, void *params)
  {
    double wav_zero;

    // extract the components from the parameter-structure
    const calib_function *wl_calib = ((trlength_search *) params)->wl_calibration;
    double wavelength        = ((trlength_search *) params)->wavelength;

    // compute the function value
    wav_zero = wl_calib->func (x, wl_calib->order, wl_calib->coeffs) - wavelength;

    // return the function value
    return wav_zero;
  }

/**
 * Function: get_tlength_prismwav
 * The function determines the tracelength value for a given wavelength
 * value and dispersion function. Should be called only for a prism dispersion.
 * The gsl root solver is used.
 *
 * Parameters:
 * @param wavelength     - the wavelength value
 * @param wl_calibration - the dispersion solution
  *
 * Returns
 * @return wav_zero - difference between actual and targeted wavelength
 */
double
get_tlength_prismwav(const double wavelength, const calib_function *wl_calibration)
{
  int iter=0;
  int status;

  d_point x_interv;

  double tr_length=0.0;


  // define and initialize the solver
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver            *s = gsl_root_fsolver_alloc (T);

  gsl_function F;
  trlength_search *tr_pars;

  // allocate and fill the parameters
  tr_pars = (trlength_search *) malloc(sizeof(trlength_search));
  tr_pars->wl_calibration = wl_calibration;
  tr_pars->wavelength     = wavelength;

  // fille the GSL-function
  F.function = &find_wavelength;
  F.params =   tr_pars;

  // initialize the intervall variable
  x_interv.x = gsl_vector_get(wl_calibration->pr_range, 0) + wl_calibration->coeffs[0];
  x_interv.y = gsl_vector_get(wl_calibration->pr_range, 1) + wl_calibration->coeffs[0];

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
      tr_length = gsl_root_fsolver_root (s);

      // derive and set new boundaries
      x_interv.x = gsl_root_fsolver_x_lower (s);
      x_interv.y = gsl_root_fsolver_x_upper (s);

      // check the accuracy
      status = gsl_root_test_interval (x_interv.x, x_interv.y,
                                       0, 0.0001);
    }
  // check for the break condition
  while (status == GSL_CONTINUE && iter < MAX_ITER_TL);

  // free the memory
  free(tr_pars);

  // free the memory
  gsl_root_fsolver_free (s);

  // return the result
  return tr_length;
}

/**
 * Function: get_tlength_prismwav
 * The function determines the tracelength value for a given wavelength
 * value and dispersion function. Should be called only for a prism dispersion.
 * The gsl root solver is used.
 *
 * Parameters:
 * @param wavelength     - the wavelength value
 * @param wl_calibration - the dispersion solution
 * @param for_grism      - boolean to indicate grism/prism dispersion solution
 *
 * Returns
 * @return wav_zero - difference between actual and targeted wavelength
 */
double
get_central_tracelength(const double wavelength, const calib_function *wl_calibration,
                        const int for_grism)
  {
    double tl_central=0.0;
    double a, b, c;

    // checc whether it's a grism
    if (for_grism)
      {
        // solve a linear dispersion
        if (wl_calibration->order == 1)
          {
            tl_central = (wavelength - wl_calibration->coeffs[0]) / wl_calibration->coeffs[1];
          }
        // solve a quadratic dispersion
        else if (wl_calibration->order == 2)
          {
            // make the quantities
            a = wl_calibration->coeffs[2];
            b = wl_calibration->coeffs[1];
            c = wl_calibration->coeffs[0] - wavelength;

            // compute the central tracelength
            //tl_central = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*c);
            tl_central = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
            //fprintf(stdout, "I should be here... %g %g %g %g\n", tl_central, a, b, c);
          }
        else
          {
            // give error message
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                "Smoothed flux conversion is impossible for grism "
                "data with dispersion solution of order: %i", wl_calibration->order);
          }
      }
    else
      {
        // use the root finder to get the tace-length value
        // for the given wavelength
        if (isnan(wavelength))
          tl_central = GSL_NAN;
        else
          tl_central = get_tlength_prismwav(wavelength, wl_calibration);
      }

    // return the central
    // tracelength value
    return tl_central;
  }

/**
 * Function: get_smoothed_response
 * The function determine a smoothe sensitivity value and its error.
 * Gaussian smoothing is applied, but the Gaussian weights are delivered
 * as a function parameter (for perfomrance reasons). The function detemines
 * the response values over the weight window and then combines these values
 * with the weights to get one weighted sensitivity value at the desired
 * wavelength.
 *
 * Parameters:
 * @param wavelength     - the central wavelength
 * @param smooth_pars    - the smoothing parameters
 * @param calib_function - the wavelength calibration function
 * @param for_grism      - boolean to indicate grism/prism solution
 * @param weights        - the weight vector
 * @param resp_func      - the response function
 * @param resp_vals      - vector for the response values
 *
 * Returns
 * @return -
 */
void
get_smoothed_response(const double wavelength, const d_point smooth_pars,
                           const calib_function *wl_calibration, const int for_grism,
                           const gsl_vector *weights, response_function *resp_func, double *resp_vals)
{
  int j;

  double tl_incr;
  double tl_act;
  double tl_central;
  double lambda_act;

  gsl_vector *pixvalues;
  gsl_vector *errvalues;
  gsl_vector *pmask;

  double *tmp_vals = malloc( 2 * sizeof(double));

  // allocate memory for the vectors
  pixvalues = gsl_vector_alloc(2 * RESP_SMOOTH_LENGTH + 1);
  errvalues = gsl_vector_alloc(2 * RESP_SMOOTH_LENGTH + 1);
  pmask     = gsl_vector_alloc(2 * RESP_SMOOTH_LENGTH + 1);

  // set all mask points to zero
  gsl_vector_set_all(pmask, 0.0);

  // determine the tracelength value
  tl_central = get_central_tracelength(wavelength, wl_calibration, for_grism);

  // determine the wavelength increments
  //wav_incr = (double)RESP_SMOOTH_NSIG * sigma_wav / (double)RESP_SMOOTH_LENGTH;
  tl_incr = smooth_pars.x * smooth_pars.y * (double)RESP_SMOOTH_NSIG / (double)RESP_SMOOTH_LENGTH;

  // set the initial wavelength value
  //value = wavelength - sigma_wav * (double)RESP_SMOOTH_NSIG;
  tl_act = tl_central - 1.0 * smooth_pars.x * smooth_pars.y * (double)RESP_SMOOTH_NSIG;

  // print the starting and central wavelength
  //fprintf(stdout, " %.5g < %.5g",
  //    wl_calibration->func(tl_act, wl_calibration->order, wl_calibration->coeffs),
  //    wl_calibration->func(tl_central, wl_calibration->order, wl_calibration->coeffs));

  // go over the array
  for (j=0; j < (int)pixvalues->size; j++)
    {
      // determine the wavelength value
      lambda_act = wl_calibration->func(tl_act, wl_calibration->order, wl_calibration->coeffs);

      // check whether the actual wavelength is inside
      // the defined sensitivity interval
      if  (lambda_act < resp_func->xmin || lambda_act > resp_func->xmax)
        {
          // for outside values:
          // set everything to zero
          gsl_vector_set(pixvalues, j, 0.0);
          gsl_vector_set(errvalues, j, 0.0);
          gsl_vector_set(pmask, j, 0.0);
        }
      else
        {
          // for inside values:
          // compute and set the response
          get_response_values(resp_func, lambda_act, tmp_vals);
          gsl_vector_set(pixvalues, j, tmp_vals[0]);
          gsl_vector_set(errvalues, j, tmp_vals[1]);
          gsl_vector_set(pmask, j, 1.0);
        }

      // increment trace length
      tl_act += tl_incr;
    }
  // print the ending wavelength
  //fprintf(stdout, " > %.5g;", wl_calibration->func(tl_act, wl_calibration->order, wl_calibration->coeffs));

  // get the weighted response values
  get_weighted_sensitivity(pixvalues, errvalues, weights, pmask, resp_vals);

  // free the memory
  free(tmp_vals);
  gsl_vector_free(pixvalues);
  gsl_vector_free(errvalues);
  gsl_vector_free(pmask);
}

/**
 * Function: get_weighted_sensitivity
 * Combines the weigh values, the sensitivity and error at the weight point
 * to a single sensitivity value plus error. A mask array defines at which
 * points values do exist.
 *
 * Parameters:
 * @param pixvalues - the central wavelength
 * @param errvalues - the Gaussian sigma used in smoothing
 * @param weights   - the weight vector
 * @param pmask     - the response function
 * @param resp_vals - vector for the response values
 *
 * Returns
 * @return
 */
void
get_weighted_sensitivity(const gsl_vector *pixvalues, const gsl_vector *errvalues, const gsl_vector *weights,
                         const gsl_vector *pmask, double *resp_vals)
{
  int index;

  // initialize the total
  // sum and weight
  double val_sum = 0.0;
  double err_sum = 0.0;
  double www=0.0;

  // go over all indices of an array
  for (index=0; index < (int)pixvalues->size; index ++)
    {
      // check whether the index is masked
      if (gsl_vector_get(pmask, index) != 0.0)
        {
          // enhance the total sum
          val_sum += gsl_vector_get(pixvalues, index) * gsl_vector_get(weights, index);
          err_sum += gsl_vector_get(errvalues, index) * gsl_vector_get(errvalues, index)
            * gsl_vector_get(weights, index) * gsl_vector_get(weights, index);

          // enhance the total weight
          www += gsl_vector_get(weights, index);
        }
    }

  // check whether there
  // were weights at all
  if (www)
    {
      // return the total sum,
      // divided by the total weight
      resp_vals[0] =  val_sum / www;
      resp_vals[1] =  sqrt(err_sum) / www;
    }
  else
    {
      // set zero sensitivity
      resp_vals[0] = 0.0;
      resp_vals[1] = 0.0;
    }
}

/**
 * Function: fill_weight
 * The function fills an array with Gaussian values.
 * The number of points and the extend are all defined in
 * global variables in the h-file.
 *
 * Parameters:
 * @param weights   - the weight vector
 *
 * Returns
 * @return -
 */
void
fill_weight(gsl_vector *weights)
{
  int j;

  double d_incr;
  double value;

  // determine the increments
  d_incr = (double)RESP_SMOOTH_NSIG / (double)RESP_SMOOTH_LENGTH;

  // compute the initial exponent
  // go over the weights vector
  value = -(double)RESP_SMOOTH_NSIG;
  for (j=0; j < (int)weights->size; j++)
    {
      // compute and set the weight
      gsl_vector_set(weights, j, exp(-1.0*value * value));

      // increment the exponent
      value += d_incr;
    }
}

/**
 * Function: get_troughput_table_name
 * Parses a configuration file and get the name of the SENSITIVITY
 * curve defined for a wanted beamID
 *
 * Parameters:
 * @param filename   - a pointer to a string containing the name of the
 *                     configuration file
 * @param beamID     - the beam ID number (see BEAM())
 * @param table_name - a pointer pointing to a string to contain the name
 *                     of the file.
 *
 * Returns:
 * @return -
 */
void
get_troughput_table_name(char *filename, int beamID, char *table_name)
{
  char beam[MAXCHAR];
  char name[MAXCHAR];
  int i=0;

  struct CfgStrings ThrConfig[] = {
    {NULL, NULL},
    {NULL,NULL}
  };
  ThrConfig[0].name = beam;
  sprintf (beam, "SENSITIVITY_%c", BEAM(beamID));

  CfgRead (filename, ThrConfig);
  if ((ThrConfig[0].name == beam) && (ThrConfig[0].data != NULL))
    {
      strcpy(name,ThrConfig[0].data);
    } else {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_throughput_table_name: %s was not found in "
                   "%s",beam,filename);
    }

  // release memory
  i=0;
  while(ThrConfig[i].name!=NULL)
    free(ThrConfig[i++].data);

  strcpy(table_name,name);
}

/**
 * Function: get_response_value_plus
 * The function determines the sensitivity at a given wavelength
 * in a sensitivity curve. A guess value helps to start searching
 * around the estimated position of the wavelength in the vector
 * of the sensitifity curve. This should considerably speed up the
 * code. the guess value is updated in case that a sensitivity
 * value could be determined. The return value is interpolated
 * between the adjacent tabulated values.
 *
 * Parameters:
 * @param  resp       - the sensitivity spectrum
 * @param  wavelength - the wavelength to determine the sensitivity for
 * @param  nguess     - the guess value for the array to start searching from
 *
 * Returns:
 * @return ret        - the sensitivity value at the input wavelength
 */
double
get_response_value_plus(const spectrum *resp, const double wavelength,
                        int *nguess)
{
  double ret;
  double factor;

  int nact;

  // check whether the wavelength is within the range of
  //the sensitivity curve return GSL_NAN if not
  if (wavelength < resp->lambdamin || wavelength > resp->lambdamax)
    return GSL_NAN;

  // check whether you have to search upwards or downwards
  if (wavelength >= resp->spec[*nguess].lambda_mean)
    {

      // in case that you search upwards, go up
      // the spectrum until you find the right interval
      nact = *nguess + 1;
      while(wavelength > resp->spec[nact].lambda_mean)
        nact++;
    }
  else
    {

      // in case that you search upwards, go up
      // the spectrum until you find the right interval
      nact = *nguess;
      while(wavelength < resp->spec[nact-1].lambda_mean)
        nact--;
    }

  // interpolate within the interval to calculate the
  // sensitivity
  factor = (wavelength - resp->spec[nact-1].lambda_mean) /
    (resp->spec[nact].lambda_mean - resp->spec[nact-1].lambda_mean);
  ret = resp->spec[nact-1].flux + factor * (resp->spec[nact].flux - resp->spec[nact-1].flux);

  // update the guess value
  *nguess = nact;

  // return the sensitivity value
  return ret;
}

/**
 * Function: get_response_values
 * The function evaluates the response function at the desired
 * independent value. The sensitivity and it error are returned
 * in an array parameter.
 *
 * Parameters:
 * @param resp_func  - the response function
 * @param wavelength - the independent value
 * @param resp_vals  - the dependent value and its error
 *
 * Returns:
 * @return -
 */
void
get_response_values(response_function *resp_func, double wavelength, double* resp_vals)
{
  // check whether the independent value is outside the deifned range
  if (wavelength < resp_func->xmin || wavelength > resp_func->xmax)
    {
      // give the default
      // outside values
      resp_vals[0] = 0.0;
      resp_vals[1] = 0.0;
    }
  else
    {
      // get the response value and its error
      resp_vals[0] = eval_interp(resp_func->resp_values, wavelength);
      resp_vals[1] = eval_interp(resp_func->resp_errors, wavelength);
    }
}


/**
 * Function: create_response_function
 * Load the sensitivity file and creates a new response
 * function. Relies heavily on other subroutines
 *
 * Parameters:
 * @param filename - sensitivity file name
 *
 * Returns:
 * @return resp_func - the response function created
 */
response_function *
create_response_function(char *filename)
{
  response_function *resp_func;

  // allocate space for the return structure;
  // complain if this fails
  resp_func = (response_function *)malloc (sizeof (response_function));
  if (resp_func == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "Could not allocate memory for response function");

  // load the two intgerpolators
  resp_func->resp_values = create_interp_ftable(filename, 2, "WAVELENGTH", "SENSITIVITY", RESP_FUNC_INTERP_TYPE);
  resp_func->resp_errors = create_interp_ftable(filename, 2, "WAVELENGTH", "ERROR", RESP_FUNC_INTERP_TYPE);

  // transfer the min and max values from the interpolator
  resp_func->xmin = resp_func->resp_values->xmin;
  resp_func->xmax = resp_func->resp_values->xmax;

  // return the respons function
  return resp_func;
}

/**
 * Function: free_response_function
 * Releases the memory allocated in a response function
 *
 * Parameters:
 * @param resp_func - the response function
 *
 * Returns:
 * @return -
  */
void
free_response_function(response_function *resp_func)
{
  // free the two interpolators
  free_interp(resp_func->resp_values);
  free_interp(resp_func->resp_errors);

  // free the rest
  free(resp_func);

  // set it to NULL
  resp_func = NULL;

}

/**
 * Function: get_calfunc_for_beam
 * The subroutine determines and returns the calibration function
 * for a certain beam and configuration file. The coefficients of
 * the wavelength calibration are evaluated at the position of the
 * beam object.
 *
 * Parameters:
 * @param actbeam   - the current beam
 * @param for_grism - indicator for grism or prism dispersion
 * @param CONF_file - full pathname to configuration file
 * @param conf      - the configuration structure
 *
 * Returns:
 * @return wl_calibration - the calibration function
  */
calib_function *
get_calfunc_for_beam(const beam actbeam, const int for_grism, char CONF_file[],
                     const aperture_conf * conf)
{
  calib_function *wl_calibration;
  dispstruct     *disp;
  d_point         pixel;

  // get the reference point right
  pixel.x = actbeam.refpoint.x - conf->refx;
  pixel.y = actbeam.refpoint.y - conf->refy;

  // get the dispersion structure
  disp = get_dispstruct_at_pos(CONF_file, for_grism, actbeam.ID, pixel);

  // transform the dispersion structure into the
  // wavelength calibration
  wl_calibration = create_calib_from_gsl_vector(for_grism, disp->pol);
  if (!for_grism)
    wl_calibration->pr_range = get_prange (CONF_file, actbeam.ID);

  // free the whole structure
  free_dispstruct(disp);

  // return the wavelength calibration
  return wl_calibration;
}
