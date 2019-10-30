/**
 */
#ifndef _IPIXCORR_UTILS_H

#define _IPIXCORR_UTILS_H

// interpolation type used for the
// complex refractive indices
#define IPCORR_INTERP_TYPE gsl_interp_linear

// interpolation type used for the
// complex refractive indices
#define NLINCORR_INTERP_TYPE gsl_interp_linear

// maximum beamID to be corrected
#define CORRMAX    0

// maximum number of iterations
// for the zero point search
#define MAX_ITER 1000

// extension for the interval in the zero
// point search to include the zero
#define INTERVEXT    10

// stepsize in x to compute the
// correction factor and the wavelength
// on the object trace
#define XSTEPSIZE    0.5

// the ipc fit to the
// data must be lower than this number
// to initialize a correction of the spectrum
#define MAXQUALITY 1.5



typedef struct
{
  int       aperID;     // the aperture  ID
  int       beamID;     // the beam ID
  int       nelems;     // the number of spectral elements
  spectrum  *fgr_spec;  // the spectrum structure for the foreground
  spectrum  *bck_spec;  // the spectrum structure for the background
  spectrum  *obj_spec;  // the spectrum structure for the object
}
  full_spectr;



/**
 * Structure: tlength_pars
 */
typedef struct
{
  beam    *actbeam;
  double   tlength;
}
tlength_pars;


extern double
get_intpix_corr(beam *actbeam, gsl_vector *cdisp, interpolator *ipcorr,
		float lambda);

extern double 
get_yfract_for_xvalue(const beam *actbeam, double xvalue);

extern double 
get_xvalue_from_tlength(beam *actbeam, double tlength);

extern d_point 
get_xinterv_from_beam(const beam *actbeam);

extern double
zero_tlength (double x, void *params);

extern double 
get_tlength_from_dispersion(gsl_vector *cdisp, float lambda);

extern gsl_vector *
condense_dispersion(gsl_vector *disp_input);

extern int 
fitting_ipc_corr(beam act_beam, char conf_file_path[], interpolator *ipcorr,
		 full_spectr *SPC, observation *obs, gsl_matrix *data_matrix,
		 char ipc_file_path[], beam *beam_ptr);

extern void
intpix_corr_beam(beam actbeam, char conf_file_path[], interpolator *ipcorr,
		 full_spectr *SPC);

extern gsl_vector *
get_ipc_lambdas(const beam actbeam, char conf_file_path[], gsl_vector *xvalues);

extern gsl_vector *
get_ipc_cvalues(const beam actbeam, interpolator *ipcorr, gsl_vector *xvalues);

extern gsl_vector *
get_ipc_xvalues(const beam actbeam);

extern interpolator *
get_ipclambda(beam actbeam, char conf_file_path[], interpolator *ipcorr);

extern void
apply_corr_pet(interpolator *ipclambda, ap_pixel *PET);

extern void
intpix_corr_pet(beam actbeam, char conf_file_path[],
		interpolator *ipcorr, ap_pixel *PET);

extern int 
is_pointlike(beam actbeam, int spec_OAF, double max_ext);

extern fitsfile *
get_SPC_opened(char SPCname[], int mode);

extern full_spectr *
get_ALL_from_next_in_SPC(fitsfile *SPC_ptr, int *aperID, int *beamID);

extern void
free_full_spectr(full_spectr *act_spectr);

extern spectrum *
get_spectrum_from_SPC(fitsfile *SPC_ptr, char count_col[], char error_col[],
		      int nelems);

extern double *
get_dcolumn_from_SPC_opened(fitsfile *SPC_ptr, char colname[], int nelems);

extern long *
get_lcolumn_from_SPC_opened(fitsfile *SPC_ptr, char colname[], int nelems);

extern interpolator *
create_nlincor(void);

extern void
nlin_corr_beam(interpolator *nlincorr, double adcgain, full_spectr *SPC);

extern double
get_nlin_corr(interpolator *nlincorr, const double lambda, const double cps);
#endif
