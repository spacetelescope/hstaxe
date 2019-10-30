#ifndef _SPC_RESP_H

#define _SPC_RESP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include "fitsio.h"
#include "aXe_grism.h"
#include "spc_cfg.h"
#include "spc_spc.h"

#define RESPBUFFERSIZE 1024

// number of points used
// for smoothing
#define RESP_SMOOTH_LENGTH 50

// the number of sigmas
// to smooth over
#define RESP_SMOOTH_NSIG 3.0

#define MAX_ITER_TL 1000

/*
 * Struct: interpolator
 */
typedef struct
{
  double xmin;
  double xmax;
  interpolator *resp_values;
  interpolator *resp_errors;
}
  response_function;

typedef struct
{
  double wavelength;
  const calib_function *wl_calibration;
}
trlength_search;

extern spectrum *
get_response_function_from_FITS(char filename[], int hdunum);

extern void
apply_response_function(spectrum *spec, spectrum *resp, const int quant_cont);

extern void
apply_response_functionII(spectrum *spec, response_function *resp_func, const int quant_cont);

extern d_point
get_smooth_pars_for_beam(const aperture_conf *conf, const int smooth_conv, beam actbeam);

extern int
check_conf_for_smoothing(const aperture_conf *conf, const int smooth_conv);

extern void
apply_smoothed_response(const calib_function *wl_calibration, const int for_grism,
                        const int quant_cont, response_function *resp_func,
                        const d_point smooth_pars, spectrum *spec);

extern double
find_wavelength(double x, void *params);

extern double
get_tlength_prismwav(const double wavelength, const calib_function *wl_calibration);

extern double
get_central_tracelength(const double wavelength, const calib_function *wl_calibration,
                        const int for_grism);

extern void
get_smoothed_response(const double wavelength, const d_point smooth_pars,
                           const calib_function *wl_calibration, const int for_grism,
                           const gsl_vector *weights, response_function *resp_func,
                           double *resp_vals);
//extern void
//get_smoothed_response(const double wavelength, const double sigma_wav, const gsl_vector *weights,
//                     response_function *resp_func, double *resp_vals);

extern void
get_weighted_sensitivity(const gsl_vector *pixvalues, const gsl_vector *errvalues, const gsl_vector *weights,
                         const gsl_vector *pmask, double *resp_vals);

extern void
fill_weight(gsl_vector *weights);

extern void
get_troughput_table_name(char *filename, int beamID, char *table_name);

extern double
get_response_value_plus(const spectrum *resp, const double wavelength,
                        int *nguess);

extern void
get_response_values(response_function *resp_func, double wavelength, double* resp_vals);

extern  response_function *
create_response_function(char *filename);

extern  void
free_response_function(response_function *resp_func);

extern calib_function *
get_calfunc_for_beam(const beam actbeam, const int for_grism, char CONF_file[], const aperture_conf * conf);

#endif
