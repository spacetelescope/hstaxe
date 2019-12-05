/**
 * File: model_utils.h
 * Header file for model_utils.c
 */
#ifndef _MODEL_UTILS_H
#define _MODEL_UTILS_H
#include "specmodel_utils.h"
#include "aXe_grism.h"
//calib_function structure defined here
#include "spc_wl_calib.h"
//aper_conf structure is defined here
#include "aper_conf.h"


#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"
#define LIGHTVEL 2.99792458  // light velocity/1.0E+08

#define NSUB 1.0    // half the number of subsampling steps
                    // done in the gaussian (take only integers!!)
#define NDIFF 3.0   // half the number of subsampling steps when
                    // diffusing the spectrum

#define MINPSF 0.1


/*
 * Struct: dirobject
 */
typedef struct dirobject
{
  int ID;		// the ID of the object

  d_point refpoint;     // the refpoint of the object

  d_point xy_off[MAX_BEAMS]; // possible offsets appy to the beams of the direct object

  d_point drzscale;     // the relative scale along the major and minor axis respectively

  int ix_min;           // smallest x-value where the direct object is defined
  int ix_max;           // largest  x-value where the direct object is defined
  int iy_min;           // smallest y-value where the direct object is defined
  int iy_max;           // largest  y-value where the direct object is defined

  int bb_sed;

  energy_distrib *SED;
  dirim_emission *dirim;
  //  double (*dpsf) (const double, const double);
}
dirobject;

/*
 * Struct: beamspec
 */
typedef struct
{
  int objectID;         // the object ID of the modelled beam
  int beamID;           // the beam ID of the modelled beam

  d_point model_ref;    // the start coo's of the beam model in the image

  gsl_matrix * model;   // the model data
}
beamspec;



/*
 * Struct: tracedata
 */
typedef struct
{
  int npoints;
  double dx_start;

  gsl_vector *dx;
  gsl_vector *dy;
  gsl_vector *xi;
  gsl_vector *lambda;
  gsl_vector *dlambda;
  gsl_vector *flux;
  gsl_vector *gvalue;
}
tracedata;

extern dirim_emission *
model_gauss_dirim(dirobject *actdir, beam actbeam, aperture_conf *conf, double psf_offset);

extern calib_function *
get_calib_function(beamspec *actspec, dirobject *actdir, char CONF_file[], const aperture_conf * conf);

extern spectrum *
get_throughput_spec(beamspec *actspec, char CONF_file[]);

extern tracedata *
compute_tracedata(const beam actbeam, const dirobject *actdir,
		  const calib_function *wl_calibration, const beamspec *actspec);

extern tracedata *
compute_short_tracedata(const aperture_conf *conf, const beam actbeam,
		const dirobject *actdir, const calib_function *wl_calibration,
		const beamspec *actspec);

extern double
get_valid_tracedata(tracedata *acttrace, const calib_function *wl_calibration);

extern void
select_tracedata(tracedata *acttrace, const calib_function *wl_calibration, const int nentries);

extern void
fill_fluxfrom_SED(const dirobject *actdir, tracedata *acttrace);

extern energy_distrib *
make_sed_from_beam(const beam onebeam, const int int_type, const int objID);

extern void
print_tracedata(tracedata *acttrace);

extern dirobject **
oblist_to_dirlist(char grism_file[], char CONF_file[], const  px_point npixels,
		  object  **oblist, spectral_models *spec_mod,
		  const double model_scale, const int int_type);

extern dirobject **
oblist_to_dirlist2(char grism_file[], char CONF_file[], const  px_point npixels,
		   object  **oblist, spectral_models *spec_mod, object_models *obj_mod,
		   const double model_scale, const int int_type);

extern dirobject *
fill_dirobject(const object *actobject, const  px_point npixels,
	       gsl_matrix *drzcoeffs, const double model_scale, const int max_offset);

extern dirobject *
fill_dirobj_fromdirim(const object *actobject, object_models *objmodels);

extern gsl_vector *
get_refpoint_ranges(const object *actobject);

extern void
fill_spectrum(const object *actobject, dirobject *actdir,
	      spectral_models *spec_mod, const int int_type);

extern dirobject *
get_dirobject_from_list(dirobject ** dirlist, const int ID);

extern d_point
get_dirobject_meanpos(dirobject *actdir);

extern beam
get_beam_for_beamspec(object **oblist, const int nobjects, const beamspec *actspec);

extern beamspec *
get_beamspec_from_list(beamspec **speclist, const int aperID, const int beamID);

extern void
print_dirobject(const dirobject * actdir);

extern int
check_interp_type(const int inter_type, const int n_flux, const int ID);

extern void
free_dirlist (dirobject ** dirlist);

extern void
free_dirlist_specmodels(dirobject ** dirlist, spectral_models *spec_mod);

extern void
free_enerdist (energy_distrib *sed);

extern void
free_speclist(beamspec **speclist);

extern void
free_tracedata(tracedata *acttrace);

//extern double
//get_flux_from_dirobject(const dirobject *actdir, double in_wave);

extern double
get_flux_from_SED(const energy_distrib *sed, double in_wave);

extern double
get_flambda_from_magab(double mag, double lambda);

extern double
get_aveflux_from_SED(const energy_distrib *sed, double in_wave, double int_wave);

extern void
fill_gaussvalues(const d_point dpixel, const beam actbeam,
		 const dirobject *actdir, const double lambda_ref,
		 const aperture_conf * conf, const double psf_offset,
		 tracedata *acttrace);

extern beam
get_newbeam(const beam actbeam, const double dpsf);

extern double
get_dpsf(const double lambda_ref, const double lambda,
	 const aperture_conf *conf, const beam actbeam);

extern double
get_sub_emodel_value(const d_point dpixel, const beam actbeam,
		     const d_point drzscale);
extern double
get_emodel_value(const d_point dpixel, const beam actbeam,
		 const d_point drzscale);

extern double
get_dpsf_WFC(const double lambda_ref, const double lambda);

extern double
get_dpsf_HRC(const double lambda_ref, const double lambda);

extern double
get_dpsf_SBC(const double lambda_ref, const double lambda);

extern double
get_polyN_gsl (const double x, const gsl_vector *params);
#endif
