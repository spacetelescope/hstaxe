/**
 * Various binning and weighting routines for aperture pixel tables.
 *
*/

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "spc_trace_functions.h"
#include "aXe_utils.h"
#include "aXe_grism.h"
#include "spc_spc.h"
#include "spc_driz.h"
#include "aper_conf.h"
#include "spce_output.h"
#include "spc_optimum.h"

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/*
gsl_matrix *
create_weightimage(ap_pixel *ap_p, const beam actbeam,
		   const aperture_conf *conf, const double exptime,
		   const double sky_cps)
{
  drzstamp_dim  dimension;
  drzstamp     *modvar;
  gsl_matrix   *weights;

  dimension = get_resample_dims(ap_p, actbeam);

  if (!dimension.resolution)
    return get_default_weight();

  compute_model_ivar(ap_p, conf->rdnoise, exptime, sky_cps);

  modvar = compute_modvar(ap_p, actbeam, dimension);

  weights = comp_allweight(modvar);

  return weights;
}
*/

/**
 * Function: compute_modvar
 *
 *
 * Parameters:
 * @param ap_p      - the PET table
 * @param actbeam   - beam to compute the modvar structure for
 * @param dimension - dimenasion of the images
 * @param modvar    - the structure to fill
 */
drzstamp *
compute_modvar(ap_pixel *ap_p, const beam actbeam,
	       const drzstamp_dim dimension)
{
  ap_pixel *cur_p;
  ap_pixel *tmp_p;

  quadrangle quad;

  double jacob;

  int jcen, icen;
  int jupp, iupp;
  int jlow, ilow;

  int ii, jj;
  int stpi, stpj;

  double maxarr, arr;
  //int iim, jjm;
  double value, allweig, weig;
  double stpc;
  double totweight;
  drzstamp *modvar;

  // allocate memory for the result structure
  modvar = alloc_drzstamp(dimension);

  // allocate memory
  tmp_p = (ap_pixel *) malloc(sizeof(ap_pixel));

  // go over each pixel
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      // Skip this pixel if it was not actually used
      if (fabs(cur_p->dist)>actbeam.width+0.5)
        continue;

      // transfer values to the temporary pixel
      tmp_p->lambda  = cur_p->xi;
      tmp_p->dist    = cur_p->dist;
      tmp_p->dxs     = cur_p->dxs;
      tmp_p->dlambda = 1.0;

      // create the quadrangle for the current pixel
      quad = get_quad_from_pixel(tmp_p, actbeam.orient, dimension);


      // get the jacobian (well, easy here)
      // the term "cos(cur_p->dxs)" must be there
      // to correct the enlargement necessary
      // to cover the whole lambda-crossdispersion area!
      // NOT COMPLETELY understood
      jacob = cos(cur_p->dxs);

      // get the central pixel (icen, jcen) of the current PET-pixel
      icen = (int) floor(cur_p->xi   - dimension.xstart+.5);
      jcen = (int) floor(cur_p->dist - dimension.ystart+.5);

      // get the uper and lower extend of the quadrangle in x
      iupp = (int)floor(quad.xmax - (double)icen + 0.5)+1;
      ilow = (int)floor(quad.xmin - (double)icen + 0.5);

      // get the uper and lower extend of the quadrangle in x
      jupp = (int)floor(quad.ymax - (double)jcen + 0.5)+1;
      jlow = (int)floor(quad.ymin - (double)jcen + 0.5);

      maxarr=0.0;
      totweight = 0.0;
      // go over the extend in x
      for (ii=ilow;ii<iupp;ii++)
	{
	  // go over the extend in x
	  for (jj=jlow;jj<jupp;jj++)
	    {

	      // get the coordinates of the current output pixel
	      stpi = icen+ii;
	      stpj = jcen+jj;

	      // check whether the current output pixel is within
	      // the stamp image; continue if not
	      if ((stpi>=dimension.xsize)||(stpi<0)
		  ||(stpj>=dimension.ysize)||(stpj<0))
		continue;

	      // get the area which falls onto the current output pixel
	      arr = boxer(stpi,stpj,quad.x,quad.y);

	      if (arr > 0.0)
		{
		  // get the already existing counts and weights
		  stpc = gsl_matrix_get(modvar->counts,stpi,stpj);
		  weig = gsl_matrix_get(modvar->weight,stpi,stpj);

		  // compute the new, total weight of the current output pixel
		  allweig = weig + arr*cur_p->weight;

		  // do a weighted sum of the count value at the current
		  // output pixel
		  value = (stpc*weig + arr*cur_p->model*cur_p->weight*jacob)
		    / (allweig);

		  // store the new count value and the new weight
		  gsl_matrix_set(modvar->counts,stpi,stpj,value);
		  gsl_matrix_set(modvar->weight,stpi,stpj,allweig);
		}
	    }
	}
    }

  // freep the temporary PET pixel
  free(tmp_p);

  // return the result
  return modvar;
}

void
shift_tracelength(ap_pixel *ap_p, const double xi_shift)
{
  ap_pixel *cur_p;

  // go over each PET pixel
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      // apply a shift in trace distance
      cur_p->xi += xi_shift;
    }
}

/**
 * Function: prepare_inv_variance
 * The function prepares the inverse model variances
 * for a data set containing an object and a background PET.
 * This is done in subroutines which compute the
 * theoretical inverse pixel variance following
 * different methods for the object and the background PET.
 *
 * Parameters:
 * @param ap_p     - the object PET table
 * @param bg_p     - the background PET table
 * @param dobck    - integer indicating background subtraction
 * @param conf     - the configuration structure
 * @param exptime  - the exposure time
 * @param sky_cps  - the sky background
 */
void
prepare_inv_variance(ap_pixel *ap_p, ap_pixel *bg_p, const int dobck,
		     const aperture_conf *conf, const double exptime,
		     const double sky_cps, const double xi_shift)
{

  if (dobck)
    {
      // shift both, the object and the background PET
      shift_tracelength(ap_p, xi_shift);
      shift_tracelength(bg_p, xi_shift);

      // compute the inverse variance for the object + background PET
      compute_total_ivar(bg_p, bg_p, conf->rdnoise, exptime, sky_cps);
    }
  else
    {
      // shift both, the object PET
      shift_tracelength(ap_p, xi_shift);

      // compute the inverse variance for the object PET
      compute_object_ivar(ap_p, conf->rdnoise, exptime, sky_cps);
    }
}

/**
 * Function: compute_model_ivar
 * The function computes for each PET pixel the
 * associated inverse variance value. The input
 * for variance and inverse variance are the
 * model value, the contamination value, the
 * constant sky background value and the readnoise.
 * The computed inverse variance value is stored
 * in the weight entry of the PET.
 * The function also offers the possibility to
 * shift the trace distance by a fixed amount.
 *
 * Parameters:
 * @param ap_p     - the PET table
 * @param rdnoise  - the readnoise value
 * @param exptime  - the exposure time
 * @param sky_cps  - the sky background
 * @param xi_shift - shift in trace distance
 */
void
compute_object_ivar(ap_pixel *ap_p, const double rdnoise,
		    const double exptime, const double sky_cps)
{
  double variance;
  double sqr_rdnoise;
  double sqr_exptime;
  ap_pixel *cur_p;

  // square the readnoise
  sqr_rdnoise = rdnoise*rdnoise;

  // square the exposure time
  sqr_exptime = exptime*exptime;

  // go over each PET pixel
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
	// compute the variance
	variance = ((cur_p->model + cur_p->contam + sky_cps) * exptime
		    + sqr_rdnoise) / sqr_exptime;

      // store the inverse variance
      cur_p->weight = 1.0/variance;
    }
}

/**
 * Function: compute_value_ivar
 * The function computes for each PET pixel the
 * associated inverse variance value. The input
 * to derive the variance is the count value and
 * the readnoise.
 * The computed inverse variance value is stored
 * in the weight entry of the PET.
 * The function also offers the possibility to
 * shift the trace distance by a fixed amount.
 *
 * Parameters:
 * @param ap_p    - the PET table
 * @param rdnoise  - the readnoise value
 * @param exptime - the exposure time
 * @param sky_cps - the sky background
 * @param xi_shift - shift in trace distance
 */
void
compute_total_ivar(ap_pixel *ap_p, const ap_pixel *bg_p, const double rdnoise,
		   const double exptime, const double sky_cps)
{
  double variance;
  double sqr_rdnoise;
  double sqr_exptime;

  ap_pixel *cur_p;
  const ap_pixel *bac_p;

  // square the readnoise
  sqr_rdnoise = rdnoise*rdnoise;

  // square the exposure time
  sqr_exptime = exptime*exptime;

  // go over each PET pixel
  bac_p = bg_p;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      // check that the foreground and background
      // PET element describe the same pixel
      if (bac_p->x != cur_p->x || bac_p->y != cur_p->y)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_PET2SPC: Background PET and Object PET "
		     "have different pixel orders in PET's.\n");
      else
      // compute the variance
      variance = ((cur_p->count + bac_p->count + sky_cps) * exptime + sqr_rdnoise) / sqr_exptime;

      // store the inverse variance
      cur_p->weight = 1.0/variance;

      // increase the background pointer
      bac_p++;
    }
}

/**
 * Function: alloc_drzstamp
 * Allocates memory for a drizzle stamp structure
 * according the the dimension specified in the input.
 *
 * Parameters:
 * @param dimension - the dimesions to allocate
 *
 * Returns:
 * @return ret  - the drizzle stamp image allocated
 */
drzstamp *
alloc_drzstamp(const drzstamp_dim dimension)
{
  drzstamp   *ret;
  gsl_matrix *counts;
  gsl_matrix *weight;

  ret = (drzstamp *) malloc(sizeof(drzstamp));

  /* Allocate the stamp matrix*/
  /* Fill stamp with NaN values */
  counts = gsl_matrix_alloc(dimension.xsize,dimension.ysize);
  //  gsl_matrix_set_all(counts, GSL_NAN);
  gsl_matrix_set_all(counts, 0.0);

  /* Allocate the weight matrix*/
  /* Fill weight with 0.0 values */
  weight = gsl_matrix_alloc(dimension.xsize,dimension.ysize);
  gsl_matrix_set_all(weight, 0.0);

  // fill the output structure
  ret->counts = counts;
  ret->weight = weight;

  // return the allocated structure
  return ret;
}

/**
 * Function: get_all_dims
 * The function derives the dimensional information
 * on the foreground and, if available, on the
 * background PET. The dimensions are compared, and
 * in case of inequality an error is thrown.
 * If the dimensions are equal, one of them is
 * returned.
 *
 * Parameters:
 * @param ap_p    - the table of aperture pixels
 * @param bg_p    - the table of background pixels
 * @param actbeam - the beam to compute the dimensions for
 * @param dobck   - integer indicating background subtraction
 *
 * Returns:
 * @return ret  - the dimension structure
 */
drzstamp_dim
get_all_dims(const ap_pixel *ap_p, const ap_pixel *bg_p, const beam actbeam,
	     const int dobck)
{
  drzstamp_dim ret;
  drzstamp_dim bck;

  // get the aperture PET diemnsion
  ret = get_resample_dims(ap_p, actbeam);

  if (dobck)
    {
      // get the background PET dimension
      bck = get_resample_dims(bg_p, actbeam);

      // check whether the two dimensions
      // are equal, throw an error if not
      if (ret.resolution != bck.resolution
	  || ret.xstart  != bck.xstart
	  || ret.ystart  != bck.ystart
	  || ret.xsize   != bck.xsize
	  || ret.ysize   != bck.ysize)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_PET2SPC: Background PET and Object PET "
		     "do not have the same dimension.\n");
    }

  // return one of the dimensional infos
  return ret;
}


/**
 * Function: get_resample_dims
 * The function parses through a PET table and determines
 * the dimension of the image in the coordinates
 * trace distance - cross-dispersion direction. Those
 * dimensions are stored in a special structure and
 * returned to the calling routine.
 *
 * Parameters:
 * @param ap_p    - the table of aperture pixels
 * @param actbeam - the beam to compute the dimensions for
 *
 * Returns:
 * @return ret  - the dimension structure
 */
drzstamp_dim
get_resample_dims(const ap_pixel *ap_p, const beam actbeam)
{
  drzstamp_dim ret;

  const ap_pixel *cur_p;

  int i_xint, i_dist;
  int imin_xint=0, imin_dist=0;
  int imax_xint=0, imax_dist=0;

  int npixel = 0;

  // immediately return empty PET's
  if (ap_p==NULL)
    return get_default_dim();

  // set the tmp pixel to the table start
  cur_p = ap_p;

  // initialize some relevant numbers
  // based on the first pixel
  i_xint = floor(cur_p->xi   + 0.5);
  i_dist = floor(cur_p->dist + 0.5);
  imin_xint = i_xint;
  imax_xint = i_xint;
  imin_dist = i_dist;
  imax_dist = i_dist;

  // go over all pixels, check
  // for new MIN's ans MAX's
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      // Skip the pixel if it was not actually used
      if (fabs(cur_p->dist)>actbeam.width+.5)
	continue;

      // compute the rounded value in (x,y)
      i_xint = floor(cur_p->xi   + 0.5);
      i_dist = floor(cur_p->dist + 0.5);

      // check for extrema in x
      imin_xint = MIN(imin_xint, i_xint);
      imax_xint = MAX(imax_xint, i_xint);

      // check for extrema in y
      imin_dist = MIN(imin_dist, i_dist);
      imax_dist = MAX(imax_dist, i_dist);

      // enhance the pixel counter
      npixel++;
    }


  // in case there were no valid pixels,
  // return a dummy structure
  if (!npixel)
    return get_default_dim();

  // fill the structure with the correct values
  ret.resolution = 1.0;
  ret.xstart     = imin_xint - NSPARE_PIX;
  ret.ystart     = imin_dist - NSPARE_PIX;
  ret.xsize      = imax_xint - imin_xint + 2 * NSPARE_PIX;
  ret.ysize      = imax_dist - imin_dist + 2 * NSPARE_PIX;

  // return the structure
  return ret;
}

/**
 * Function: get_default_dim
 * The function creates and returns the default
 * dimension structure. The dewfault dimension
 * structure contains a 0.0 in all fields.
 *
 * Returns:
 * @return res  - the dimension structure
 */
drzstamp_dim
get_default_dim()
{
  drzstamp_dim ret;

  // just set everything to zero
  ret.resolution=0.0;
  ret.xstart = 0;
  ret.ystart = 0;
  ret.xsize  = 0;
  ret.ysize  = 0;

  // return the result
  return ret;
}

/**
 * Function: get_default_modvar
 * The function computes the default drizzle stamp
 * structure, which in this case is used for
 * tha calculation of the resampled variance and
 * profile image.
 *
 * Returns:
 * @return res  - the resulting drzstamp structure
 */
drzstamp *
get_default_modvar()
{
  gsl_matrix *counts;
  gsl_matrix *weight;
  drzstamp *res;

  res = (drzstamp *) malloc(sizeof(drzstamp));

  /* Create a dummy stamp image */
  counts = gsl_matrix_alloc(10,10);
  weight = gsl_matrix_alloc(10,10);
  /* Fill stamp with 0.0 values */
  gsl_matrix_set_all(counts, 0.0);
  gsl_matrix_set_all(weight, 0.0);

  // fill the output
  res->counts = counts;
  res->weight = weight;

  // return the output
  return res;
}

gsl_matrix *
get_default_weight()
{
  gsl_matrix *res;

  res = gsl_matrix_alloc(10,10);
  gsl_matrix_set_all(res, 0.0);

  // return the output
  return res;
}

/*
 * Function: comp_opt_weight
 * The function uses the model map and and the variance
 * map and computes optimale weigths according to the
 * Hoorne method.
 * Pixels outside of the extraction
 * area are set to 1000.0
 *
 * Parameters:
 * @param mod_map - the matrix with the model values
 * @param var_map - the matrix with the variance values
 * @param ob      - the object structure
 *
 * Returns:
 * @return weight - the matrix with the exposure time weights
 */
gsl_matrix *
comp_allweight(drzstamp *modvar)
{
  gsl_matrix *weight;

  double mod_sum, weight_sum;
  double contr, norm, allweight;
  double mod_val;
  double act_weight=0;

  //int beamInt = 0;
  int i, j;

  // allocate the weight matrix and set the default
  weight = gsl_matrix_alloc(modvar->counts->size1, modvar->counts->size2);
  gsl_matrix_set_all(weight, 1000.0);

  //* go over all columns
  for (i=0; i < (int)modvar->counts->size1; i++)
    {
      mod_sum = 0.0;
      contr = 0.0;
      allweight = 0.0;
      weight_sum=0.0;

      // determine for each column the total model counts
      // and the number of pixels with non-zero model_counts
      for (j=0; j < (int)modvar->counts->size2; j++)
	{
	  mod_sum = mod_sum + gsl_matrix_get(modvar->counts, i, j);
	  contr = contr + 1.0;
	}

      // check whether the column has model values.
      // normalize the model values and compute
      // optimal weights if yes.
      if (mod_sum > 0.0)
	{
	  //* determine the mean model value
	  //	  norm = mod_sum / contr;
	  norm = mod_sum;

	  // go over each row
	  for (j=0; j < (int)modvar->counts->size2; j++)
	    {
	      // normalize the model counts
	      mod_val = gsl_matrix_get(modvar->counts, i, j)/norm;

	      // store the normalized model counts
	      gsl_matrix_set(modvar->counts, i, j, mod_val);

	      // add up the normalization value for the optimal weights
	      weight_sum = weight_sum
		+ mod_val*mod_val*gsl_matrix_get(modvar->weight, i, j);
	    }

	  // finally compute and write the weights:
	  // go over each pixel
	  for (j=0; j < (int)modvar->counts->size2; j++)
	    {
	      if (gsl_matrix_get(modvar->counts, i, j) > 0.0 )
		{
		  // compute and set the individual pixel weight
		  act_weight = gsl_matrix_get(modvar->counts,i,j)
		    *gsl_matrix_get(modvar->weight,i,j)/weight_sum;
		  //		  gsl_matrix_set(weight,i,j,gsl_matrix_get(modvar->counts,i,j)
		  //		 *gsl_matrix_get(modvar->weight,i,j)/weight_sum);
		  gsl_matrix_set(weight, i, j, act_weight);
		  allweight = allweight + gsl_matrix_get(weight, i, j);
		}
	      else
		{
		  // set the default value if no model value is zero
		  gsl_matrix_set(weight, i, j,0.0);
		}
	    }
	}
      // if the column does not have model values at all:
      else
	{
	  // go over each row
	  for (j=0; j < (int)modvar->counts->size2; j++)
	    {
	      // set the inside default value
	      gsl_matrix_set(weight, i,j,1.0);
	      allweight = allweight + gsl_matrix_get(modvar->counts, i, j);

	    }
	}
    }

  // return the weight matrix
  return weight;
}
