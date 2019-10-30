#include "spc_flatfield.h"

#define MAX_FLATS 10 /* maximal number of flatfields for xxx_multi functions */
#define SQR(x) ((x)*(x))

/**
  Creates a simulation of a flat field using random numbers.  
  Use only for debugging.
*/
gsl_matrix *
simulate_flatfield (flatfield_d * const flat, double stdev)
{
  int i, j;
  gsl_matrix *flat_pixels;
  
  gsl_rng *r = gsl_rng_alloc (gsl_rng_ran0);
  
  flat_pixels = gsl_matrix_alloc (flat->w, flat->h);
  for (i = 0; i < (int)flat_pixels->size1; i++)
    {
      for (j = 0; j < (int)flat_pixels->size2; j++)
	{
	  gsl_matrix_set (flat_pixels, i, j,
			  1. + gsl_ran_gaussian (r, stdev));
	}
    }
  gsl_rng_free (r);
  
  return flat_pixels;
}

/** 
    Load and return a gsl cube whose planes (z-direction) contain the 
    coefficients of a polynomial of the form ff(i,j,lambda) = a0+a1*lambda+...

    @param fname the name of a flat cube to read and load

    @returns a pointer to a NULL terminated array of gsl_matrix
*/
poly_cube_flatfield * load_flat_poly_cube(char *fname)
{
  int i;
  int order = 0;
  poly_cube_flatfield *poly_cube;
    
  /* Getting the number of extension sin the file */
  order = FITSextnum(fname);
    
  /* Allocating the memory */
  poly_cube = (poly_cube_flatfield *) malloc(sizeof(poly_cube_flatfield));
  poly_cube->poly_order = order;
  poly_cube->coeffs = (gsl_matrix **) malloc(poly_cube->poly_order * sizeof(gsl_matrix *));
  
  /* Load the FF coefficients, on eplane at a time */
  for (i=0;i<poly_cube->poly_order;i++)  
    {
      poly_cube->coeffs[i] = FITSimage_to_gsl (fname, i+1, 1);
    }
  /* Load the minimum and maximum wavelengths this cube applies to */
  poly_cube->wmin = get_float_from_keyword(fname, 1, "WMIN");
  poly_cube->wmax = get_float_from_keyword(fname, 1, "WMAX");
  
  return poly_cube;
}



/**
    Evaluates and return the flatfielding value of pixel coordinates x,y
    and at wavelength lambda, using the field dependent polynomial
    description contained in a poly_cube_flatfield structure. Uses abs(lambda)
    for cosmetic reasons.
    If the sampled wavelength falls outside of the wmin<w<wmax range, then a warning
    is issued and no flat-field coefficient is computed, and 1.0 is returned
	[new version now returns the FF value at the closest known wavelength]
    @param lambda the wavelength 

*/
double
poly_cube_flatfield_lambda (const double lambda, const int x, const int y,
		poly_cube_flatfield *poly_cube)
{
  int i;
  double ff=0.0;
  double w; /* nomalized wavelength */

    
  /* if ( (lambda>poly_cube->wmax) || (lambda<poly_cube->wmin) )
     {
     aXe_message (aXe_M_WARN4, __FILE__, __LINE__, "Sampled wavelengh (%f)"
     " is outside of FF cube (%f < ww < %f). Returning 1.0\n",lambda,poly_cube->wmin,poly_cube->wmax);
     w = 1.0;
     for (i=0;i<poly_cube->poly_order;i++) 
     {
     ff = ff + gsl_matrix_get(poly_cube->coeffs[i],x,y)*pow(w,i);
     }
     return ff;
     }*/
    
  w = (fabs(lambda) - poly_cube->wmin)/(poly_cube->wmax-poly_cube->wmin);

  if (lambda<poly_cube->wmin)
    {
      w = .0;
    }        
  if (lambda>poly_cube->wmax)
    {
      w = 1.0;
    }    
  
  for (i=0;i<poly_cube->poly_order;i++) 
    {
      ff = ff + gsl_matrix_get(poly_cube->coeffs[i],x,y)*pow(w,i);
    }
  return ff;
}

/**
    Completely frees up the space allocated to a gsl_matrix
    @param poly_cube a pointer to a poly_cube_flatfield structure

*/
void free_flat_poly_cube(poly_cube_flatfield *poly_cube)
{
  int i=0;
  
  if (poly_cube!=NULL) {
    for (i=0;i<poly_cube->poly_order;i++)
      gsl_matrix_free(poly_cube->coeffs[i]);
    free(poly_cube->coeffs);
    if (poly_cube) free(poly_cube);
    poly_cube=NULL;
  }
}


/**
  computes the relative efficiency of a pixel from a base flat
  field and a polynomial describing the wavelength dependence of
  the efficiency.

  @param lambda the wave length to assume
  @param x x coordinate of the pixel in question relative to the origin of
    the entire image
  @param y ditto for y
  @param flat pointer to the flatfielding structure containing information
    for the ture
  @returns relative efficiency for (x,y)
*/
static void
flat_poly_func (const double lambda, const int x, const int y,
		PIXEL_T * const val, PIXEL_T * const err,
		const flatfield_d * const flat)
{
  int i;
  double *coeffs = flat->data.poly->poly_coeffs;
  double res = 0;
  
  for (i = flat->data.poly->poly_order; i > 0; i--)
    {
      res += coeffs[i];
      res *= lambda;
    }
  res += coeffs[0];
  //fprintf(stderr,"poly: %f %f\n",res,gsl_matrix_get(flat->data.poly->flatfield,
  //  x-flat->ll_x, y-flat->ll_y));
  
  *val =
    res * gsl_matrix_get (flat->data.poly->flatfield, x - flat->ll_x,
				y - flat->ll_y);
  if (flat->data.poly->errors)
    {
      *err =
	gsl_matrix_get (flat->data.poly->errors, x - flat->ll_x,
			y - flat->ll_y);
    }
  else
    {
      *err = 0;
    }
}


/**
  loads (a part) of a flatfield file

  @param flat a pointer to the flatfield descriptor
  @param flat_name the name of the file containing the flat field
  @return a pointer to a gsl_matrix containing the relative efficiencies
    for each pixel in the bbox defined in flat
  @todo do real fits io
*/
gsl_matrix *
load_flatfield (flatfield_d * const flat, const char *const flat_name)
{
  gsl_matrix *flat_pixels;
  
  if (!flat_name)
    {
      flat_pixels = simulate_flatfield (flat, 0.);
    }
  else
    {
      flat_pixels = NULL;
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "FITS loading for flatfields not yet implemented");
    }
  return flat_pixels;
}


/**
  loads (a part) of the error extension from a flatfield file

  @param flat a pointer to the flatfield function containing the
    bounding box
  @param flat_name the name of the file containing the flat field
  @return a pointer to a gsl_matrix containing the relative efficiencies
    for each pixel in the bbox defined in flat
  @todo do real fits io
*/
gsl_matrix *
load_flatfield_errors (flatfield_d * const ff, const char *const flat_name)
{
  gsl_matrix *err_pixels;
  
  if (!flat_name)
    {
      int i, j;
      
      err_pixels =
	gsl_matrix_alloc (ff->data.poly->flatfield->size1,
			  ff->data.poly->flatfield->size2);
      for (i = 0; i < (int)err_pixels->size1; i++)
	{
	  for (j = 0; j < (int)err_pixels->size2; j++)
	    {
	      //gsl_matrix_set(err_pixels, i, j, 
	      //  1-gsl_matrix_get(ff->data.poly->flatfield, i, j));
	      gsl_matrix_set (err_pixels, i, j, .01);
	    }
	}
    }
  else
    {
      err_pixels = NULL;
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "FITS loading for flatfields not yet implemented");
    }
  return err_pixels;
}


/**
  creates a descriptor for a flatfield with a single flat field and
  a polynomial wave length dependence, such that the efficiency
  at wavelength lambda is flat(x,y)*poly(lambda).

  @param curbeam the beam to create the flatfield for, with the bbox filled
    out (e.g., after spc_extract has run)
  @param flat_name name of a flat field file
  @param order the order of the polynomial
  @param coeffs the coefficients of the polynomial (contains order+1 doubles)
  @return an allocated flatfield descriptor
*/
flatfield_d *
make_poly_flatfield (const beam * const curbeam, const char *const flat_name,
		     const int order, const double coeffs[])
{
  flatfield_d *const ff = malloc (sizeof (flatfield_d));
  int i;
  
  if (!ff)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  
  ff->data.poly = malloc (sizeof (polynom_flatfield));
  if (!ff->data.poly)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  ff->data.poly->poly_order = order;
  ff->data.poly->poly_coeffs = malloc ((order + 1) * sizeof (double));
  if (!ff->data.poly->poly_coeffs)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  
  for (i = 0; i <= order; i++)
    {
      ff->data.poly->poly_coeffs[i] = coeffs[i];
    }
  ff->func = flat_poly_func;
  
  ff->ll_x = curbeam->bbox[0].x;
  ff->ll_y = curbeam->bbox[0].y;
  ff->w = curbeam->bbox[1].x - curbeam->bbox[0].x + 1;
  ff->h = curbeam->bbox[1].y - curbeam->bbox[0].y + 1;
  ff->data.poly->flatfield = load_flatfield (ff, flat_name);
  ff->data.poly->errors = load_flatfield_errors (ff, flat_name);
  
  return ff;
}


/**
  frees a polynom flatfield and its assoicated data structures

  @param flat  the flatfield to free
*/
void
free_poly_flatfield (flatfield_d * flat)
{
  gsl_matrix_free (flat->data.poly->flatfield);
  if (flat->data.poly->errors)
    {
      gsl_matrix_free (flat->data.poly->errors);
    }
  free (flat->data.poly->poly_coeffs);
  flat->data.poly->poly_coeffs = NULL;
  free (flat->data.poly);
  flat->data.poly = NULL;
  free (flat);
  flat = NULL;
}


/**
  computes the relative efficiency of a pixel from a set of flatfields
  at various lambdas.

  @param lambda the wave length to assume
  @param x x coordinate of the pixel in question relative to the origin of
    the entire image
  @param y ditto for y
  @param flat pointer to the flatfielding structure containing information
    for the aperture
  @returns relative efficiency for (x,y)
  @todo check if the interpolation is good at the borders of the lambda range
*/
static void
flat_multi_func (const double lambda, const int x, const int y,
		 PIXEL_T * const val, PIXEL_T * const err,
		 const flatfield_d * const flat)
{
  double xvals[MAX_FLATS], yvals[MAX_FLATS];
  //gsl_interp_factory factory = gsl_interp_factory_cspline_natural;
  //gsl_interp_obj *interpolator;
  //gsl_interp_accel *accelerator;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double d_val;
  
  for (int i = 0; i < flat->data.multi->num_flats; i++)
    {
      xvals[i] = flat->data.multi->lambdas[i];
      yvals[i] =
	    gsl_matrix_get (flat->data.multi->flatfields[i], x - flat->ll_x, y - flat->ll_y);
    }
  //interpolator = factory.create(xvals, yvals, flat->data.multi->num_flats); 
  //accelerator = gsl_interp_accel_new();
  //gsl_interp_eval_impl(interpolator, xvals, yvals, 
  //  lambda, accelerator, &d_val);

  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, flat->data.multi->num_flats);
  gsl_spline_init (spline, xvals, yvals, flat->data.multi->num_flats);
  //d_val = gsl_spline_eval (spline, lambda, acc);
  d_val = spline->interp->type->eval(spline->interp->state, spline->x, spline->y, spline->interp->size, lambda, acc, yvals);
  *val = (PIXEL_T) d_val;
  *err = 0;
}


/**
  creates a descriptor for a flatfield with multiple flatfields at
  various wavelengths

  @param curbeam the beam to create the flatfield for, with the bbox filled
    out (e.g., after spc_extract has run)
  @param num_flats the number of (lambda,flatfield) pairs
  @param lambda (and the next num_flats odd parameters) the wave length of the
    following flat field
  @param flat_name (and the next num_flats even parameters) name of
    the flat field
  @return an allocated flatfield descriptor
*/
flatfield_d *
make_multi_flatfield (const beam * const curbeam, const int num_flats,
		      const double lambda, const char *const flat_name, ...)
{
  flatfield_d *ff = malloc (sizeof (flatfield_d));
  va_list args;
  int i;
  
  if (num_flats < 3)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Number of flats must be" " larger than 3.");
    }
  if (num_flats > MAX_FLATS)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Too many flatfields." " Increase MAX_FLATS in "
		   __FILE__ " and recompile me.");
    }
  
  if (!ff)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  ff->ll_x = curbeam->bbox[0].x;
  ff->ll_y = curbeam->bbox[0].y;
  ff->w = curbeam->bbox[1].x - curbeam->bbox[0].x + 1;
  ff->h = curbeam->bbox[1].y - curbeam->bbox[0].y + 1;
  ff->func = flat_multi_func;
  
  ff->data.multi = malloc (sizeof (multi_flatfield));
  if (!ff->data.multi)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  ff->data.multi->num_flats = num_flats;
  ff->data.multi->errors = malloc (sizeof (gsl_matrix *) * num_flats);
  ff->data.multi->flatfields = malloc (sizeof (gsl_matrix *) * num_flats);
  ff->data.multi->lambdas = malloc (sizeof (double) * num_flats);
  if ((!ff->data.multi->errors) || (!ff->data.multi->flatfields)
      || (!ff->data.multi->lambdas))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");
    }
  
  ff->data.multi->lambdas[0] = lambda;
  ff->data.multi->flatfields[0] = load_flatfield (ff, flat_name);
  ff->data.multi->errors[0] = NULL;
  va_start (args, flat_name);
  for (i = 1; i < num_flats; i++)
    {
      ff->data.multi->lambdas[i] = va_arg (args, double);
      ff->data.multi->flatfields[i] =
	load_flatfield (ff, va_arg (args, char *));
      ff->data.multi->errors[i] = NULL;
    }
  va_end (args);
  
  return ff;
}


/**
  frees a flatfield consisting of multiple images and its 
  assoicated data structures

  @param flat  the flatfield to free
*/
void
free_multi_flatfield (flatfield_d * flat)
{
  int i;

  for (i = 0; i < flat->data.multi->num_flats; i++)
    {
      gsl_matrix_free (flat->data.multi->flatfields[i]);
      if (flat->data.multi->errors[i])
	{
	  gsl_matrix_free (flat->data.multi->errors[i]);
	}
    }
  free (flat->data.multi->flatfields);
  flat->data.multi->flatfields = NULL;
  free (flat->data.multi->errors);
  flat->data.multi->errors = NULL;
  free (flat->data.multi->lambdas);
  flat->data.multi->lambdas = NULL;
  free (flat->data.multi);
  flat->data.multi = NULL;
  free (flat);
  flat = NULL;
}

