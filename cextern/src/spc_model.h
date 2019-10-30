#ifndef _SPC_MODEL_H
#define _SPC_MODEL_H

// this file needs to define beamspec
#include "model_utils.h"
// observation is defined here
#include "aXe_grism.h"
// calib_function defined here
#include "spc_wl_calib.h"
// flux_cube is defined here
#include "spc_fluxcube.h"


extern int
compute_gauss_cont(char grism_file[], char OAF_file[], char CONF_file[],
		   const char specmod_file[], const double model_scale,
		   const int inter_type, const double lambda_psf, observation *obs,
		   const char PET_file[], char map_file[], const int store);

extern int
compute_gaussdirim_cont(char grism_file[], char OAF_file[], char CONF_file[],
			const char specmod_file[], const char objmod_file[],
			const double model_scale, const int inter_type,
			const double lambda_psf, observation *obs,
			const char PET_file[], char map_file[], const int store);

extern int
compute_gauss_dirim(char grism_file[], char OAF_file[], char CONF_file[],
		   const double model_scale,  const int inter_type,
		   const double lambda_psf, observation *obs,
		   const char PET_file[], char map_file[], const int store);

extern int
compute_fcube_cont(char grism_file[], char OAF_file[], char fcube_file[],
		   char CONF_file[], const double model_scale, const int inter_type,
		   observation *obs, const char PET_file[], char map_file[],
		   const int store);

extern int
compute_geometr_cont(char OAF_file[], observation *obs,
		     const char PET_file[], char map_file[],
		     const int store);

extern int
fill_contam_info(const char PET_file[], beamspec  **speclist,
		 const gsl_matrix *all_models, char model_name[]);

extern gsl_matrix *
make_model_image(const px_point npixels, observation *obs, beamspec **speclist);

extern beamspec **
make_fcube_spectra(object **oblist, dirobject **dirlist,
		   const px_point npixels, char CONF_file[],
		   const flux_cube *fcube, const int inter_type);

extern beamspec **
alloc_beamlist_from_dirlist(object **oblist, dirobject **dirlist,
			    const px_point npixels, aperture_conf *conf);

extern gsl_matrix *
make_gauss_dirim(object **oblist, dirobject **dirlist, const double lambda_psf,
		  const px_point npixels, char CONF_file[], observation *obs);

extern beamspec **
make_gauss_spectra(object **oblist, dirobject **dirlist, const double lambda_psf,
		   const px_point npixels, char CONF_file[]);

extern beamspec **
make_gauss_spectra2(object **oblist, dirobject **dirlist, const double lambda_psf,
		    const px_point npixels, char CONF_file[]);

extern int
get_index_for_tracepoint(const tracedata *acttrace, const double dx);

extern int
fill_pixel_in_speed(const dirobject *actdir, const tracedata *acttrace,
		    const d_point dpixel, const spectrum *resp,
		    beamspec *actspec, const calib_function  *wl_calibration);

extern int
no_diffuse_spectrum(int ix, int iy, double cps, beamspec *actspec);

extern int
diffuse_spectrum(double ddx, double ddy, double cps, beamspec *actspec);

extern int
diffuse_spectrumII(double ddx, double ddy, double cps, beamspec *actspec);

beamspec *
dimension_beamspec(dirobject *actdir, object *actobject,
		   const px_point npixels, const aperture_conf * conf, int j);
#endif
