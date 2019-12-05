/**
 * File: scaleback_utils.h
 * header files for scaleback_utils.c
 *
 */
#ifndef _SCALEBACK_UTILS_H
#define _SCALEBACK_UTILS_H

#define N_BCKSCALE_ITER   10
#define N_BCKSCALE_KAPPA  3.0
#define N_BCKSCALE_ACCUR  1.0e-04
#define MASK_VALUE -1.0e+06

FITScards *
fit_to_FITScards(const gsl_vector* bck_vals, const px_point npixels);

void
make_scale_back(char grism_image[], const char grism_mask[], char conf_file[],
		const char scale_image[], char bck_image[], const char plist_name[],
		const int scale_to_master, const int make_plis);

void
print_plis(const fitbck_data *fbck_data, const char plist_name[]);

gsl_matrix *
compute_scale_grism(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
		int dqmask, gsl_vector *bck_vals);

gsl_matrix *
compute_scale_master(const gsl_matrix *sc_img, const gsl_vector *bck_vals);

void
fill_mask_data(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
               gsl_matrix *sc_img, fitbck_data *fbck_data, int dqmask);

void
fill_cont_data(gsl_matrix *gr_img, gsl_matrix *gr_dqval, gsl_matrix *gr_mask,
               gsl_matrix *sc_img, fitbck_data *fbck_data, int dqmask);

gsl_vector *
make_ksig_scalefit(fitbck_data *fbck_data);

void
get_bck_scale(const double *xs, double *ys, double *ws,
		  const int n, gsl_vector *bck_vals);

int
clipp_scale_data(fitbck_data *fbck_data, const gsl_vector *bck_vals,
		const float kappa);

#endif
