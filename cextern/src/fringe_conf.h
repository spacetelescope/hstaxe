/**
 */
#ifndef _FRINGE_CONF_H
#define _FRINGE_CONF_H

// switch to let the debugging
// information print to the screen
//#define DEBUGFCONF 1

// the minimum number of layers
#define  NLAYERS_MIN 1

// the default fringe step in AA
#define DEFAULT_FRINGE_STEP 1.0

// the default number of steps
// to resolve the pixel throughput
#define DEFAULT_NUM_STEPS 100

// the default maximal resolution
#define DEFAULT_MAX_DISPERSION 1.0e+4

// the default maximal resolution
#define DEFAULT_MIN_LAMBDA 0.0

// the default maximal resolution
#define DEFAULT_MAX_LAMBDA 2.0e+4

// threshold for filter values
// to be considered relevant
#define FILTER_THRESHOLD 1.0e-05

// interpolation type used for the filter
// throughput values
#define FILTER_INTERP_TYPE gsl_interp_linear

// interpolation type used for the
// complex refractive indices
#define REFRAC_INTERP_TYPE gsl_interp_linear

/*
 * Struct: interpolator
 */
typedef struct
{
  double           xmin;     // minimum in the independent variable
  double           xmax;     // maximum in the independent variable
  int              nvals;    // the naumber of data values
  double           *xvals;   // array for the independent values
  double           *yvals;   // array for the dependent values
  gsl_interp_accel *acc;     // spline accelerator
  gsl_interp       *interp;  // interpolation structure
}
  interpolator;

/*
 * Struct: linear_interp
 */
typedef struct
{
  double     xmin;      // minimum in the independent variable
  double     xmax;      // maximum in the independent variable
  int        num_elem ; // the number of elements in the arrays
  int        act_index; // the actual position in the interpolator
  gsl_vector *xvals;    // the array of independent values
  gsl_vector *yvals;    // the array of dependent value
}
  linear_interp;


/*
 * Struct: ccd_layer
 */
typedef struct
{
  double       thickness;      // thicknes of the CCD in $[\mu m]$
  gsl_matrix   *thickness2D;   // the layer thickness at every pixel in [mum]
  interpolator *re_refraction; // real part of refractive index
  interpolator *im_refraction; // imaginary part of refractive index
}
  ccd_layer;


/*
 * Struct: ccd_layers
 */
typedef struct
{
  int          num_layers;  // the number of CCD layers
  ccd_layer    **opt_layer; // array with optical layers in a CCD
  interpolator *substrate;  // to represent the CCD substrate
}
  ccd_layers;


/*
 * Struct: fringe_conf
 */
typedef struct
{
  double       fringe_amp;      // the fringe amplitude
  double       fringe_phase;    // phase at the first layer
  double       fringe_step;     // step interval to examine
  int          num_steps;       // number of steps for filter function
  gsl_vector   *fringe_range;   // minimum and maximum wavelength for fringing
  double       max_dispersion;  // maximal dispersion considered in fringing
  interpolator *filter_through; // structure for the filter throughput
  ccd_layers   *opt_layers;     // the different optical layers in a CCD
}
  fringe_conf;

extern ccd_layer *
load_CCD_layer(const char refr_table[], const char thickness[]);

extern void
free_CCD_layer(ccd_layer *opt_layer);

extern ccd_layers *
load_CCD_layers(char fring_conf_path[]);

extern void
free_CCD_layers(ccd_layers *opt_layers);

extern fringe_conf *
load_fringe_conf(char fring_conf_path[]);

extern void
check_fringe_conf(fringe_conf *fconf);

extern void
free_fringe_conf(fringe_conf *fconf);

extern interpolator *
create_interp_ftable(const char table_name[], const int fdunum,
		     char xcol[], char ycol[],
		     const gsl_interp_type *T);

extern interpolator *
create_interp(const int nvals, const gsl_interp_type *interp_type,
	      double *xvals, double *yvals);

extern void 
print_interp(interpolator *interp);

extern double
eval_interp(interpolator *interp, const double xval);

extern void 
free_interp(interpolator *interp);

extern linear_interp *
create_linint_ftable(const char table_name[], const int fdunum,
		     char xcol[], char ycol[]);

extern void
free_linint(linear_interp *lin_int);

#endif
