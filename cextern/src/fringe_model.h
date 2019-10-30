/**
 */
#ifndef _FRINGE_MODEL_H
#define _FRINGE_MODEL_H


/*
 * Struct: optical_property
 */
typedef struct
{
  double trans_upper;         // the transmission coefficient
                              // to the next upper layer
  double reflect_upper;       // the reflection coefficient
                              // to the next upper layer
  double sign_upper;          // sign for the reflection coefficient
                              // to the next upper layer
  double thickness;           // thickness of the layer
  double double_attenuation;  // the amplitude attenuation when light of a
                              // certain wavelength traverses the layer twice
  gsl_complex double_phshift; // the phase shift when light of a
                              // certain wavelength traverses the layer twice
  double trans_lower;         // the transmission coefficient
                              // to the next lower layer
  double reflect_lower;       // the reflection coefficient
                              // to the next lower layer
  double sign_lower;           // sign for the reflection coefficient
                              // to the next lower layer
}
  optical_property;


extern gsl_matrix  *
compute_fringe_amplitude(fringe_conf *fconf);

extern gsl_vector **
evaluate_wavelength_steps(fringe_conf *fconf);

extern interpolator *
redefine_filter_throughput(const int lower, const int upper,
			   fringe_conf *fconf);

extern gsl_matrix  *
alloc_fringe_image(const ccd_layers *opt_layers);

extern gsl_vector **
get_calibration_data(void);

gsl_vector **
get_PET_calibration_data(void);

extern double
eval_linear_interp(linear_interp *lin_int, const double xval);

extern double
get_layer_thickness(const ccd_layer *opt_layer, const int ii, const int jj);

extern gsl_complex
get_complex_refindex(const ccd_layer *opt_layer, const double lambda);

extern double
compute_reflection(const gsl_complex refract_l1, const gsl_complex refract_l2);

extern double
compute_transmission(const gsl_complex refract_l1,
		     const gsl_complex refract_l2);

extern double
compute_attenuation(const gsl_complex refract, const double thickness,
		    const double phase_number);

extern gsl_complex
compute_pshift(const gsl_complex refract, const double thickness,
	       const double phase_number);

extern optical_property *
alloc_optprops_list(const fringe_conf *fconf);

extern void
print_optprops_list(const optical_property *optprops, const int num_entries);

extern void
free_optprops_list(optical_property *optprops);

extern void
init_optprops_list(const fringe_conf *fconf, const double lambda_mean,
		   optical_property *optprops);

extern void
fill_optprops_thickness(const ccd_layers *opt_layers, const int ii,
			const int jj, optical_property *optprops);

extern void
fill_optprops_all(const ccd_layers *opt_layers, const double lambda,
		  optical_property *optprops);

extern double
fringe_contrib_single(const optical_property *optprops,
		      const fringe_conf *fconf);
#endif
