/**
 * File: nicback_utils.h
 * header files for nicback_utils.c
 *
 */
#ifndef _NICKBACK_UTILS_H
#define _NICKBACK_UTILS_H

#define N_KSIGBCK_ITER 3
#define N_KSIGBCK_KAPPA  3.0

#define FRAC_BPIX_MIN 0.1
#define NPIX_DIM1 256.0
#define NPIX_DIM2 256.0
#define AREA_EXTEND_MIN 0.5
#define BCK_SCALE_MIN 0.5
#define BCK_SCALE_MAX 1.5


typedef struct
{
  int n_data;

  double *x_values;
  double *y_values;
  double *e_values;
  int *x_pos;
  int *y_pos;
}
  fitbck_data;

extern void
make_nicmos_back(const observation * const obs, const char *msk_name,
		 const char *master_bck, char bck_name[], char plist_name[],
		 const char corr_bck[]);

extern fitbck_data *
alloc_fitbck_data(const int n_data);

extern double
get_skypix_frac(const fitbck_data *fbck_data);

extern  gsl_matrix *
compute_nicmos_back(const gsl_matrix *mbck_img, const double factor,
		    const double offset, const gsl_matrix *corr_img);

extern gsl_vector *
make_ksig_linfit(fitbck_data *fbck_data);

extern double
comp_stdev_from_array(const double *data, const int ndata, const double mean);

extern int
check_fbck_data(fitbck_data *fbck_data);

extern gsl_vector*
get_fbck_defaults(void);

extern int
clipp_fbck_data(fitbck_data *fbck_data, const gsl_vector *lfit, const float kappa);

extern void
fill_fbck_data(const observation * const obs, const gsl_matrix *msk_img,
	       const gsl_matrix *mbck_img, fitbck_data *fbck_data, char plist_name[],
	       const gsl_matrix *corr_img);

extern void
free_fitbck_data(fitbck_data *fbck_data);
#endif
