/**
 * File: trfit_utils.h
 *
 */
#ifndef _TRFIT_UTILS_H
#define _TRFIT_UTILS_H
#include "aper_conf.h"
#include <gsl/gsl_multifit_nlin.h>

#define N_COLUMNS_FIT 6
#define N_KAPPASIG_ITER 3
#define N_KAPPASIG_SIG  3.0
#define N_ITER_IPC  5
#define KAPPA_IPC  3.0

// the ipc fit to the
// data must be lower than this number
// to initialize a correction of the spectrum
#define MIN_IPC_LENGTH 60.0

typedef struct
{
  int n_data;

  double *x_values;
  double *y_values;
  double *e_values;
}
  fit_data;

struct function_data {
  int      n;
  double * x;
  double * y;
  double * sig;
};


extern gsl_vector *
fit_beamtrace(const aperture_conf  *conf, observation *obs,
	      const int beamID, beam act_beam);

double
gagauss(const double x, const double *params);

extern gsl_vector *
fit_wuttke(const fit_data *f_data);

extern gsl_vector *
fit_wuttke_talk(const fit_data *f_data);

extern d_point
find_grav_center(const fit_data *f_data);

extern px_point
get_fit_xrange(const aperture_conf  *conf, const observation *obs,
	       beam act_beam);

extern gsl_vector *
fit_function(const fit_data *f_data);

extern fit_data *
get_fitdata(observation *obs, int np, px_point tr_point);

extern fit_data *
alloc_fit_data(const int n_data);

void
print_fit_data(fit_data *f_data);

extern void
free_fit_data(fit_data *f_data);

int
gauss_f(const gsl_vector *x, void *params, gsl_vector *f);

int
gauss_df(const gsl_vector *x, void *params, gsl_matrix *J);

int
gauss_fdf (const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);

extern void
print_fit_state (size_t iter, gsl_multifit_fdfsolver * s);

extern double
comp_intrel_max(const fit_data *f_data);

extern gsl_matrix *
bcksub_observation(observation *obs, observation *bck);

extern double *
evaluate_ipc_fit(fit_data *f_data, double *gpars, double *fpars);

extern void
write_ipc_fit(fit_data *f_data, double *y_values, const char *ipc_file, const double qual, double *yfract);

extern double
fit_ipc_data(fit_data *fdata, double *gpars, int n_gpars,
	     double *fpars, int n_fpars);

extern fit_data *
get_ipc_fdata(const aperture_conf  *conf, observation *obs,
	      gsl_matrix *data_vals, beam act_beam);

extern double *
get_wf_special(const aperture_conf  *conf,  observation *obs,
	       gsl_matrix *data_vals, beam act_beam, double yoffs);

extern gsl_vector *
get_ipc_coldata(gsl_matrix *img_data, int np, px_point tr_point);

extern fit_data *
strip_fitdata(fit_data *old_fdata, int i_max);

extern double *
strip_wf_special(fit_data *old_fdata, double *y_fract, int i_max);

extern d_point
kappa_sigma_klipp_ipc(const aperture_conf  *conf, observation *obs,
		      gsl_matrix *data_matrix, beam act_beam,
		      const char *ipc_file);

extern void
reject_kappasigma(fit_data *f_data, double *y_values, double max_diff);

extern double
calculate_sigma(fit_data *f_data, double *y_values);
#endif
