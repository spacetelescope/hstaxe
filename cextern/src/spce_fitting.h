#ifndef _SPC_FITTING_H
#define _SPC_FITTING_H


#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include "aXe_errors.h"


extern void
comp_vector_average (const double *xs, double *ys,
		       double *ws, double *yi, const int n, const int final);
extern void
comp_vector_median (const double *xs, double *ys,
		      double *ws, double *yi, const int n, const int final);

extern void
comp_vector_linear (const double *xs, double *ys,
		      double *ws, double *yi, const int n, const int final);

extern void
comp_vector_polyN (const int m, const double *xs, double *ys,
		    double *ws, double *yi, const int n, const int final);

extern void
fill_const_value(double *ys, double *ws, double *yi, const int n,
		 double cval, double stdev, const int final);

extern void
det_vector_average (const double *xs, double *ys,
		    double *ws, const int n, double *avg, double *std);
extern void
det_vector_median (const double *xs, double *ys,
		    double *ws, const int n, double *med, double *std);

extern gsl_vector *
det_vector_linear(const double *xs, double *ys, double *ws,
		  const int n, const int weight);

extern gsl_vector *
det_vector_poly_N (int m, const double *const xs, double *const ys,
		   double *const ws, const int n, gsl_vector *c,
		   gsl_matrix *cov);

extern void
fill_linear_interp(const double *const xs, double *const ys,
		   double *const ws, double *yi, const int n,
		   gsl_vector *interp, const int final);

extern void
fill_polyN_interp(const double *const xs, double *const ys,
		  double *const ws, double *yi, const int n,
		  gsl_vector *coeffs, gsl_matrix *cov, gsl_vector *interp,
		  const int final);
#endif
