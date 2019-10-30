/**
 * The interface to the binning and weight routines
 */
#ifndef _SPC_OPTIMUM_H
#define _SPC_OPTIMUM_H 1

#define NSPARE_PIX 5

extern gsl_matrix  *
create_weightimage(ap_pixel *ap_p, const beam actbeam,
		   const aperture_conf *conf, const double exptime,
		   const double sky_cps);

extern drzstamp *
compute_modvar(ap_pixel *ap_p, const beam actbeam,
	       const drzstamp_dim dimension);

void
prepare_inv_variance(ap_pixel *ap_p, ap_pixel *bg_p, const int dobck,
		     const aperture_conf *conf, const double exptime,
		     const double sky_cps, const double xi_shift);

extern void
compute_object_ivar(ap_pixel *ap_p, const double rdnoise,
		   const double exptime, const double sky_cps);

void
compute_total_ivar(ap_pixel *ap_p, const ap_pixel *bg_p, const double rdnoise,
		   const double exptime, const double sky_cps);

extern drzstamp *
alloc_drzstamp(const drzstamp_dim dimension);

drzstamp_dim
get_all_dims(const ap_pixel *ap_p, const ap_pixel *bg_p, const beam actbeam,
	     const int dobck);

extern drzstamp_dim
get_resample_dims(const ap_pixel *ap_p, const beam actbeam);

extern drzstamp_dim
get_default_dim(void);

extern drzstamp *
get_default_modvar(void);

extern gsl_matrix *
get_default_weight(void);

extern gsl_matrix *
comp_allweight(drzstamp *modvar);

#endif
