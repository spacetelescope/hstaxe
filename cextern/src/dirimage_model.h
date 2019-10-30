/**
 */
#ifndef _DIRIMAGE_MODEL_H
#define _DIRIMAGE_MODEL_H

// interpolation type used for the
// total passband curves
#define TPASS_INTERP_TYPE gsl_interp_linear

extern int 
compute_dirimage_model(char dirim_file_path[], char conf_file_path[], char tpass_file_path[],
		       char specmod_file_path[], char objmod_file_path[], char aper_file_path[],
		       const double model_scale, const double tel_area, const double lambda_psf,
		       observation *obs, char map_file_path[]);

extern interpolator *
get_filter_sensitivity(const char tpass_file[], const double tel_area);

extern gsl_matrix *
make_dirimage(object **oblist, dirobject **dirlist, const px_point npixels,
	      const double lambda_psf, interpolator *tpass);


extern double
get_cps_for_dirobject(interpolator *tpass, dirobject *actdir);

extern double
integrate_interpolator(interpolator *combine);

extern interpolator *
combine_tpass_SED(interpolator *tpass, dirobject *actdir);

extern interpolator *
get_combined_tpass_SED(gsl_vector *indep_data, interpolator *tpass, dirobject *actdir);
#endif
