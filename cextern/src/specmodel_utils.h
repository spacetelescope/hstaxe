/**
 * File: specmodel_utils.h
 * Header file for specmodel_utils.c
 *
 */
#ifndef _SPECMODEL_UTILS_H
#define _SPECMODEL_UTILS_H
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "fitsio.h"
#include "aXe_grism.h"


// interpolation type used for the
// spectral models
#define SMODEL_INTERP_TYPE gsl_interp_linear

/*
 * Struct: energy_distrib
 */
typedef struct
{
  int npoints;
  double *wavelength;
  double *flux;

  gsl_interp *interp;
  gsl_interp_accel *accel;

}
energy_distrib;

/*
 * Struct: spectral_models
 */
typedef struct
{
  int n_models;

  energy_distrib **SEDlist;
}
spectral_models;

typedef struct
{
  int dim_x;               // the number of pixels in x-direction
  int dim_y;               // the number of pixels in y-direction
  double xmean;           // mean x-coordinate of the object
  double ymean;           // mean y-coordinate of the object

  gsl_matrix *modimage;   // matrix with the direct image
}
  dirim_emission;

typedef struct
{
  int n_models;              // the number of direct image object

  dirim_emission **obj_list; // the list of direct image objects
}
object_models;

extern object_models *
load_object_models(const char object_models_file[]);

extern dirim_emission **
load_diremission_list(const char object_models_file[], const int n_models);

extern dirim_emission *
load_diremission(const char object_models_file[], const int n_extension);

extern void
free_object_models(object_models *objmodels);

extern void
free_dirim_emission(dirim_emission *diremission);

extern void
print_object_models(const object_models *objmodels);

extern void
print_dirim_emission(const dirim_emission *diremission);

extern dirim_emission *
get_dirim_emission(const object_models *objmodels, const int objspec);

extern int
has_aperture_dirim(const object_models *objmodels, const object *actobject);

extern double
get_diremission_value(const dirim_emission *diremission,
		      const double xpos, const double ypos);

extern double
bilin_interp_matrix(const gsl_matrix *modimage, const double x, const double y);

extern spectral_models *
load_spectral_models(const char spectral_models_file[]);

extern void
print_spectral_models(const spectral_models *smodels);

extern energy_distrib **
load_SED_list(const char spectral_models_file[], const int n_models);

extern energy_distrib *
load_SED_from_fitsext(const char spectral_models_file[], fitsfile *s_models);

extern int
get_num_extensions(const char spectral_models_file[]);

extern void
free_spectral_models(spectral_models *smodels);

extern void
free_specmodels_SEDs(spectral_models *smodels);

extern energy_distrib *
get_model_sed(const spectral_models *spec_mod, const int modspec);

#endif
