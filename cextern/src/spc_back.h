/**
 * Header file for the subroutines in spc_back.c
 *
 */
#ifndef _SPC_BACK_H
#define _SPC_BACK_H

#include "spce_is_in.h"

#define DQ_KAPPA_SIGMA -1.0

/**
  A descriptor of a background, consisting of a function to
  compute the background at absolute position x,y and a
  parameter block (e.g., matrices containing the background)
*/
typedef struct s_background
{
  void (*bck_func) (const int x, const int y, PIXEL_T * const val,
                    PIXEL_T * const err, const struct s_background * pars);
  void *pars;
}
background;

/**
  An aggregation of the background and the error for each background pixel,
  used for the fullimg family of background functions.
*/
typedef struct
{
     gsl_matrix *bck;
     gsl_matrix *err;
}
fullimg_background;

extern int
is_pt_in_a_beam (const px_point * const apoint,
                 const is_in_descriptor * const iids, const int tnbeams);

extern void
compute_background2(object **oblist, const int obj_index,
                    const int beamorder, gsl_matrix *bck_mask,
                    const int tnbeams, fullimg_background * fib, int npoints,
                    int interporder);

//gsl_vector_int *
//get_window_points(observation *obs, gsl_matrix *bck_mask, beam actbeam, 
//                px_point tr_point);

extern gsl_vector_int *
get_interp_points(observation *obs, gsl_matrix *bck_mask,
                  int np, px_point tr_point);

extern void
compute_background (observation *obs, beam actbeam, gsl_matrix *bck_mask,
                    fullimg_background *fib, int npoints,
                    int interporder, const int niter_med,
                    const int niter_fit, const double kappa);

extern void
compute_global_background (object **oblist, const int obj_index,
                           gsl_matrix *bck_mask, fullimg_background * fib,
                           int interporder);
extern void
fullimg_background_function (const int x, const int y, PIXEL_T * const val,
                             PIXEL_T * const err,
                             const background * const back);

extern background *
compute_fullimg_background2(observation *obs, object **oblist,
                            int npoints, int interporder);

extern background *
compute_fullimg_background(observation *obs, object **oblist,
                           int npoints, int interporder, const int niter_med,
                           const int niter_fit, const double kappa,
                           int nor_flag, const int sm_length, const double fwhm);

extern background *
compute_backsub_mask (observation *obs, object **oblist);

extern background *
compute_fullimg_global_background (observation *obs, object **oblist,
                                   int interporder, const int sm_length, const double fwhm);

extern void 
free_fullimg_background (background * const fib);

extern void 
background_to_FITSimage (char filename[], background * bck, observation *obs);

extern gsl_matrix *
aperture_mask (observation * const obs, object **oblist);

extern void
comp_vector_interp(const double *const xs, double *const ys,
                   double *const ws, double *const yi, const int n,
                   const int interp, const int final);

extern void
comp_kappasigma_interp(const double *const xs, double *const ys,
                       double *const ws, const int n,
                       const int interp, const int niter, const double kappa,
                       observation *obs, int colnum);
void
kappa_sigma_clipp(const double *const xs, double *const ys, double *const ws, const int n, const double kappa,
                  observation *obs, int colnum);

extern px_point
get_xrange(observation *obs, beam actbeam);

extern void 
gsmooth_background (const gsl_matrix *bck_mask, const int smooth_length,
                    const double fwhm, fullimg_background *fib);

extern double
get_weighted_mean(const gsl_vector *pixvalues, const gsl_vector *weights,
                  const gsl_vector *pmask);

extern void
fill_pixvalues(const gsl_matrix *bck_mask, const int smooth_length,
               const fullimg_background *fib, const int ix, const int iy,
               gsl_vector *pixvalues, gsl_vector *pmask);

extern double
compute_gvalue(const double xdiff, const double efactor);

extern double
compute_efactor(const double fwhm);

#endif
