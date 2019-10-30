/**
 * Subroutines used in the "aXe_AF2PET.c" for the task aXe_AF2PET
 */

#include <stdarg.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include "spc_wl_calib.h"
#include "aXe_grism.h"
#include "aXe_errors.h"
#include "disp_conf.h"


/**
 * The function performs a wavelength calibration on an allready
 * existing ap_pixel table. It fills in the lambda and dlambda fields.
 *
 * @param ap_p  - the ap_pixel table to work on, with the pathlength field
 *                filled out
 * @param calib - the calibration function
 */
void
wl_calib (ap_pixel * ap_p, const calib_function * const calib)
{
  double l1,l2;

  while (ap_p->p_x != -1)
    {
      ap_p->lambda = calib->func (ap_p->xi, calib->order, calib->coeffs);

      l1 = calib->func (ap_p->xi - .5, calib->order, calib->coeffs);
      l2 = calib->func (ap_p->xi + .5, calib->order, calib->coeffs);
      ap_p->dlambda = fabs(l2-l1);
      ap_p++;
    }
}

/*
 * Function: pwise_wl_calib
 * Make for every PET pixel an individual wavelength
 * calibration. The function was developed for FORS2 MXU spectra
 * and makes especially sense for long spectra and/or fast
 * variations of the 2D calibration solution.
 */
void
pwise_wl_calib(const global_disp *gdisp, const d_point pixel,
	       const beam actbeam, const int for_grism, ap_pixel *ap_p, const calib_function *old_calib)
{

  calib_function *calib;
  d_point actpix;

  double l1,l2;
  double tmp;

  int i=0;

  while (ap_p[i].p_x != -1)
    {
      // get the current xy-position of the pixel
      actpix.x = pixel.x + ap_p[i].dist*cos(actbeam.orient);
      actpix.y = pixel.y + ap_p[i].dist*sin(actbeam.orient);

      // create a calibration for that position
      calib = create_calib_from_gsl_vector(for_grism, get_calvector_from_gdisp(gdisp, actpix));

      // determine the wavelength value
      ap_p[i].lambda = calib->func (ap_p[i].xi, calib->order, calib->coeffs);

      // determine the dlambda value
      l1 = calib->func (ap_p[i].xi - .5, calib->order, calib->coeffs);
      l2 = calib->func (ap_p[i].xi + .5, calib->order, calib->coeffs);
      ap_p[i].dlambda = fabs(l2-l1);

      // determine the trace length
      tmp =  old_calib->func (ap_p[i].xi, old_calib->order, old_calib->coeffs);
      ap_p[i].xi = ap_p[i].xi + (ap_p[i].lambda - tmp) / ap_p[i].dlambda;

      // free the calibration
      // structure
      free_calib(calib);

      // enhance the counter
      i++;
    }
}

/**
 *  The function trims down an ap_pixel table to only allowed values
 *  for the trace path length or trace distance. The input ap_pixel
 *  is only changed. Valid entries are moved to the front side,
 *  and a new endmark is set after all valid entries
 *
 *  @param  ap_pixel       the input ap_pixel table
 *  @param  prange         the valid trace range
 *  @param  wl_calibration the calibration structure, needed to get the
 *                         first dispersion coefficient
 *
 *  @return ret            ap_pixel table output
 */
ap_pixel *
prange_cut(ap_pixel * in_p, const gsl_vector * prange,
           const calib_function * wl_calibration)
{

  //long nelem;
  ap_pixel * ret, *cur_ap, *ap_p;
  double prel;
  double trmin=1.0E+06, trmax=-1.0E6;

  // set all PET list pointers
  // to the beginning of the array
  ret    = in_p;
  cur_ap = in_p;
  ap_p   = in_p;

  while (ap_p->p_x != -1)
    {
      // search for new minimum
      // and amximum in tracelength
      if (ap_p->xi > trmax)
        trmax = ap_p->xi;
      if (ap_p->xi < trmin)
        trmin = ap_p->xi;

      // compute the tracelength with
      // respect to the zeropoint
      prel = ap_p->xi - wl_calibration->coeffs[0];

      // check whether the realitive zeropoint is within the allowed range
      if (prel > gsl_vector_get(prange,0) && prel < gsl_vector_get(prange,1))
        {
          // transfer the properties
          // of the good pixel to the next
          // PET element in the list of good PET's
          cur_ap->p_x    = ap_p->p_x;
          cur_ap->p_y    = ap_p->p_y;
          cur_ap->x      = ap_p->x;
          cur_ap->y      = ap_p->y;
          cur_ap->dist   = ap_p->dist;
          cur_ap->xs     = ap_p->xs;
          cur_ap->ys     = ap_p->ys;
          cur_ap->dxs    = ap_p->dxs;
          cur_ap->xi     = ap_p->xi;
          cur_ap->count  = ap_p->count;
          cur_ap->error  = ap_p->error;
          cur_ap->weight = ap_p->weight;
          cur_ap->contam = ap_p->contam;
          cur_ap->dq     = ap_p->dq;

          // forward the pointer
          // of good pixels
          cur_ap++;
        }

      // always forward the
      // general pointer
      ap_p++;
    }

  // make a new endmark
  cur_ap->p_x = -1;
  cur_ap->p_y = -1;
  cur_ap->count = -1;

  // return the pointer
  return ret;
}


/**
 *  Calibration function for the grism: simple polynomial
 *
 * @param xi     the path length of the sample point
 * @param order  the order of the polynomial, i.e. the number of coefficients
 *               minus one
 * @param coeffs a pointer to order+1 doubles, starting with the zeroth
 *               coefficent
 *
 * @return res   wavelength at xi
 */
static double
grism_calib_func (const double xi, const int order,
		  const double *const coeffs)
{
  int i;
  double res = 0;

  for (i = 0; i <= order; i++)
    {
      res += coeffs[i] * pow(xi,i*1.0);
    }

  return res;
}


/**
 * Calibration function for the prism: reciprocal polynomial in the
 * form: lambda = a1 + a2/(xi-a0) + a3/(xi-a0)**2 + a4/(xi-a0)**3
 *
 * @param  xi     the path length of the sample point
 * @param  order  the order of the polynomial, i.e. the number of coefficients
 *                minus one
 * @param  coeffs a pointer to order+1 doubles, starting with the zeroth
 *                coefficent
 * @return res    wavelength at xi
 */
static double
prism_calib_func (const double xi, const int order,
		  const double *const coeffs)
{
  int i;
  double res = 0;

  if (xi==0.) return GSL_NAN;
  for (i = 0; i <= order-1; i++)
    {
      //	    res += coeffs[i] / pow(xi,i*1.0);
      res += coeffs[i+1] / pow((xi-coeffs[0]),i*1.0);
    }

  return res;
}

/**
 * Create a prism or grism calibration function, i.e. a function of one of
 * the forms "sum(a_i 1/x^i, i=0..order)" or "sum(a_i x^i, i=0..order)"
 * for transforming pathlength along the trace to lambda.
 *
 *@param for_grism if true, a the second form of the calibration function
 *                 is used (normal polynom), otherwise the first one.
 *@param order     order of the polynom
 *@param a         a GSL vector of length n containing the n coefficients
 *                 describing the wavelength calibration polynomial.
 *@return calib    an allocated structure containing an allocated array
 *                 for the coefficients.  Free using calib_free.
 */
calib_function *
create_calib_from_gsl_vector (const int for_grism, const gsl_vector * a)
{
  calib_function *calib;
  int i, order;

  // set the polynomial order
  order = a->size - 1;

  // allocate memory
  calib = malloc(sizeof (calib_function));

  // complain if memory
  // allocation failed
  if (!calib)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");

  // store the order of the
  // calibration
  calib->order = order;

  // transfer the calibration functions
  if (for_grism)
    calib->func = &grism_calib_func;
  else
    calib->func = &prism_calib_func;

  // allocate memory; complain if problems occur
  calib->coeffs = malloc ((order + 1) * sizeof (double));
  if (!calib->coeffs)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory");

  // transfer the coefficients
  calib->coeffs[0] = gsl_vector_get (a, 0);
  for (i = 1; i <= order; i++)
      calib->coeffs[i] = gsl_vector_get (a, i);

  // set the allowed range for
  // the trace lengths to NULL
  calib->pr_range = NULL;

  // return the calibration
  // structure
  return calib;
}

/**
 * Function: free_calib
 * Free a calibration function and its assoicated data.
 *
 * @param  calib the calibration structure to be freed
 *
 */
void
free_calib(calib_function * calib)
{
  free (calib->coeffs);
  calib->coeffs = NULL;
  if (calib->pr_range != NULL)
    {
      gsl_vector_free(calib->pr_range);
      calib->pr_range = NULL;
    }
  free (calib);
  calib = NULL;
}
