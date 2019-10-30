/**
* Various binning and weighting routines for aperture pixel tables.
*/

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "aXe_grism.h"
#include "aper_conf.h"
#include "spc_driz.h"
#include "spce_output.h"
#include "spc_optimum.h"
#include "spce_binning.h"

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define DEBUG_ME 0x200


#define PIXWEIGHT(x1,y1,x2,y2,pix) ((pix)->weight_function((x1), (y1), (x2),\
  (y2),(pix)))
#define NAIVE_VAL_TO_BININD(x) ((int)floor((x)+1e-6))

/**
 * Function: add_to_spec_table
 * adds some count to an entry in the spectrum table
 *
 * Parameters:
 * @param spec the spectrum table to work on
 * @param bin index of spectrum table entry to add cur_p to
 * @param cur_p the ap_pixel to add
 * @param weight weight of the count
 */
static void
add_to_spec_table (spectrum * const spec, const int bin,
		   const ap_pixel * const cur_p, const int quant_cont,
		   const double weight)
{

  spc_entry *const sp_e = spec->spec + bin;

  // some conditions which should not be violated
  if ((bin < 0) || (bin > spec->spec_len))
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		   "Assignment out of spectrum: %d", bin);
      return;
    }
  if (weight < 0)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Weight cannot be negative " "but is %f", weight);
    }

  // check whether the spectral element is new
  // and without values up to now
  if (isnan (sp_e->lambda_mean))
    {

      // initialize the spectral element
      sp_e->lambda_mean = cur_p->lambda;
      sp_e->dlambda = cur_p->dlambda;
      sp_e->lambda_max = cur_p->lambda;
      sp_e->lambda_min = cur_p->lambda;
      sp_e->weight = weight;
      sp_e->count = cur_p->count * weight;
      sp_e->error = fabs(cur_p->error) * weight;
      sp_e->dq = cur_p->dq;

      // initialize the contamination
      // this is different for quantitative and
      // geometrical contamination
      if (quant_cont)
	{
	  if ((int)(sp_e->contam==-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam = cur_p->contam * weight;
	    }
	}
      else
	{
	  if ((int)(sp_e->contam==-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam = cur_p->contam;
	    }
	}
    }
  else
    {

      // update an existing spectral bin

      // find new maxima and minima
      sp_e->lambda_max = MAX (cur_p->lambda, sp_e->lambda_max);
      sp_e->lambda_min = MIN (cur_p->lambda, sp_e->lambda_min);

      // find new mean lambda via weighted summation
      sp_e->lambda_mean =
	(sp_e->lambda_mean * sp_e->weight +
	 cur_p->lambda * weight) / (weight + sp_e->weight);

      // NEWNEWNEW::
      // find dlambda via weighted summation
      sp_e->dlambda =
	(sp_e->dlambda * sp_e->weight +
	 cur_p->dlambda * weight) / (weight + sp_e->weight);

      // add the weight
      sp_e->weight += weight;

      // add the counts
      sp_e->count += cur_p->count * weight;

      // process the error
      sp_e->error = sqrt (SQR (sp_e->error) + SQR (fabs(cur_p->error) * weight));

      // logically XOR the dq
      sp_e->dq = (sp_e->dq | cur_p->dq);

      // update the contamination,
      // take into account the quantitative
      // contamination
      if (quant_cont)
	{
	  if (((int)sp_e->contam==-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam = cur_p->contam * weight;
	    }
	  if (((int)sp_e->contam!=-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam += cur_p->contam * weight;
	    }
	}
      else
	{
	  if (((int)sp_e->contam==-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam = cur_p->contam;
	    }
	  if (((int)sp_e->contam!=-1)&&((int)cur_p->contam!=-1))
	    {
	      sp_e->contam += cur_p->contam;
	    }
	}
    }
}


/**
 * Function: pixweight_x
 * computes a weight for the bin between coordinates (x1,y1) and (x2,y2)
 * if the pi/4<=angle<=3 pi/4.
 *
 * Parameters:
 * @param x1 - x coordinate for the start point of the bin on the trace
 * @param y1 - y coordinate for the start point of the bin on the trace
 * @param x2 - x coordinate for the end point of the bin on the trace
 * @param y2 - y coordinate for the end point of the bin on the trace
 * @param pix the pixel to compute the weight for in the form of a
 *            w_pixel structure filled out by fill_w_pixel
 *
 * Returns:
 * @return sum - the pixel weight
 */
double
pixweight_x (const double x1, const double y1, const double x2,
	     const double y2, const struct w_pixel_s *const pix)
{
  double a, b;
  double sum = 0;

  a = x1 - pix->cota * (y1 - pix->y0);
  b = x2 - pix->cota * (y2 - pix->y0);
  if ((b >= pix->p0) && (a <= pix->p1))
    {
      sum += pix->slope / 2 * (MIN (b, pix->p1) - MAX (a, pix->p0))
	* (MAX (a, pix->p0) - 2 * pix->p0 + MIN (b, pix->p1))
	* sin (pix->angle);
    }
  if ((b >= pix->p1) && (a <= pix->p2))
    {
      sum += pix->fmax * (MIN (b, pix->p2) - MAX (a, pix->p1))
	* sin (pix->angle);
       }
  if ((b >= pix->p2) && (a <= pix->p3))
    {
      sum += pix->slope / 2 * (MAX (a, pix->p2) - MIN (b, pix->p3))
	* (MAX (a, pix->p2) - 2 * pix->p3 + MIN (b, pix->p3))
	* sin (pix->angle);
    }
  return sum;
}


/**
 * Funtion: pixweight_y
 * Computes a weight for the bin between coordinates (x1,y1) and (x2,y2)
 * if the angle<=pi/4 or 3 pi/4<=angle.
 *
 * @param x1 - x coordinate for the start point of the bin on the trace
 * @param y1 - y coordinate for the start point of the bin on the trace
 * @param x2 - x coordinate for the end point of the bin on the trace
 * @param y2 - y coordinate for the end point of the bin on the trace
 * @param pix - the pixel to compute the weight for in the form of a
 *              w_pixel structure filled out by fill_w_pixel
 *
 * Returns:
 * @return sum - the pixel weight
 */
double
pixweight_y (const double x1, const double y1, const double x2,
	     const double y2, const struct w_pixel_s *const pix)
{
  double a, b;
  double sum = 0;

  a = y1 + pix->tana * (pix->x0 - x1);
  b = y2 + pix->tana * (pix->x0 - x2);
  if (a > b)
    {
      double tmp;
      tmp = a;
      a = b;
      b = tmp;
    }

  if ((b >= pix->p0) && (a <= pix->p1))
    {
      sum += pix->slope / 2 * (MIN (b, pix->p1) - MAX (a, pix->p0))
	* (MAX (a, pix->p0) - 2 * pix->p0 + MIN (b, pix->p1))
	* cos (pix->angle);
    }
  if ((b >= pix->p1) && (a <= pix->p2))
    {
      sum += pix->fmax * (MIN (b, pix->p2) - MAX (a, pix->p1))
	* cos (pix->angle);
    }
  if ((b >= pix->p2) && (a <= pix->p3))
    {
      sum += pix->slope / 2 * (MAX (a, pix->p2) - MIN (b, pix->p3))
	* (MAX (a, pix->p2) - 2 * pix->p3 + MIN (b, pix->p3))
	* cos (pix->angle);
    }
  return sum;
}


/**
 * Function: fill_w_pixel
 * precomputes some properties of a given pixel for purposes
 * of computing the weights it contributes to a given xi bin.
 *
 * Parameters:
 * @param pix a pointer to the w_pix structure to fill
 * @param x0 the x coordinate of the pixel's lower left corner
 * @param y0 the y coordinate of the pixel's lower left corner
 * @param size the size of the pixel
 * @param angle the orientation of the object that has generated the
 *        spectrum
 */
/* since we're only interested in weights here, the size of the square
 *  doesn't matter.  Where I left explicit ones, you could write size
 * and get the true area instead of the fraction of the total area.
 * Dunno why you'd want to do this, though.
 */
static void
fill_w_pixel (w_pixel * const pix, const double x0, const double y0,
	      const double angle)
{
  pix->tana = tan (angle);
  pix->cota = 1 / pix->tana;
  pix->angle = angle;
  pix->x0 = x0;
  pix->y0 = y0;

  if ((angle >= M_PI / 4) && (angle <= 3 * M_PI / 4))
    {
      pix->p0 = MIN (x0, x0 - 1 * pix->cota);
      pix->p1 = MAX (x0, x0 - 1 * pix->cota);
      pix->p2 = MIN (x0 + 1, x0 + 1 - 1 * pix->cota);
      pix->p3 = MAX (x0 + 1, x0 + 1 - 1 * pix->cota);
      pix->fmax = 1 / sin (angle);

      if (fabs (pix->p1 - pix->p0) < 1e-7)
	pix->slope = 0;
      else
	pix->slope = pix->fmax / (pix->p1 - pix->p0);

      pix->weight_function = &pixweight_x;
    }
  else
    {			/* angle between 0 and pi/4 or 3*pi/4 and pi */
      pix->p0 = MIN (y0, y0 - 1 * pix->tana);
      pix->p1 = MAX (y0, y0 - 1 * pix->tana);
      pix->p2 = MIN (y0 + 1, y0 + 1 - 1 * pix->tana);
      pix->p3 = MAX (y0 + 1, y0 + 1 - 1 * pix->tana);
      pix->fmax = 1 / cos (angle);

      if (fabs (pix->p1 - pix->p0) < 1e-7)
	pix->slope = 0;
      else
	pix->slope = pix->fmax / (pix->p1 - pix->p0);

      pix->weight_function = &pixweight_y;
    }
}


/**
 * Function: bin_naive
 * computes a spectrum from a table of aperture pixels generated from
 * spc_extract.  This is adds up the pixel values, distributing them
 * over the trace, taking into account the fracton of the pixel that
 * projects onto the given [xi,xi+1] interval.
 * Return NULL if ap_p is NULL
 *
 * Parameters:
 * @param ap_p the table of aperture pixels
 * @param px_width width of a pixel (=1 if not subsampled)
 * @param ob_width the width of the object that has generated the spectrum
 * @param ob_orientation the orientation of the object that has
 *             generated the spectrum
 * @param flags problems that were accumulated in generating ap_p;
 *   possible flags are defined for the warning member of the spectrum struct
 *
 * Returns:
 * @return spec - the 1D spectrum
 */
spectrum *
bin_naive (const ap_pixel * const ap_p, const double ob_width,
	   const double ob_orient, const int quant_cont)
{
  const ap_pixel *cur_p;
  int upper, lower;

  int bin;
  spectrum *spec, *tspec;
  double phi_trace;
  spc_entry *spec_table;
  double d, frac;

  // immediately return empty PET's
  if (ap_p==NULL)
    return NULL;

  // go through the PET,
  // find the minimum and
  // maximum in trace distance
  cur_p = ap_p;
  upper = NAIVE_VAL_TO_BININD (cur_p->xi);
  lower = NAIVE_VAL_TO_BININD (cur_p->xi);
  while (cur_p->p_x != -1)
    {
      bin = NAIVE_VAL_TO_BININD (cur_p->xi);
      upper = MAX (bin, upper);
      lower = MIN (bin, lower);
      cur_p++;
    }

  // check whether the spectrum
  // will ahve a finite length,
  // exit if not
  if (upper == lower)
    {
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		   "Pixel table empty.\n");
      return NULL;
    }

  // Thresh in some headway so we don't need to worry too much about
  //     accessing invalid elements now and then
  lower -= 10;
  upper += 10;
  spec = allocate_spectrum (upper - lower);

  spec_table = spec->spec;
  cur_p = ap_p;

  while (cur_p->p_x != -1)
    {
      // Compute any fractional pixel that might fall within the
      //desired extraction width
      d = ob_width;
      frac = 1.;

      // continue if the pixel is outside
      // of the extraction region
      if (fabs (cur_p->dist) > d + .5)
	{
	  cur_p++;
	  continue;
	}
      if ((fabs (cur_p->dist) >= d - .5) && (fabs (cur_p->dist) <= d + .5))
	{
	  frac = fabs (d - (fabs (cur_p->dist) - 0.5));
	}

      // store the local trace angle
      phi_trace = cur_p->dxs;

      if (1)
	{
	  double xi;
	  w_pixel pix;
	  double sinp = sin (phi_trace);
	  double cosp = cos (phi_trace);
	  double xc, yc;
	  double w;
	  double sum = 0;

	  fill_w_pixel (&pix, cur_p->x,	cur_p->y, ob_orient);

	  xc = cur_p->xs;
	  yc = cur_p->ys;

	  // at cur_p->xi, there has to be some contribution. We go back
	  // collecting, until w is zero
	  for (xi = cur_p->xi;; xi -= 1)
	    {
	      bin = NAIVE_VAL_TO_BININD (xi);
	      w = PIXWEIGHT (xc + (bin - cur_p->xi) * cosp,
			     yc + (bin - cur_p->xi) * sinp,
			     xc + (bin + 1 - cur_p->xi) * cosp,
			     yc + (bin + 1 - cur_p->xi) * sinp,
			     &pix);
	      if (w < 1e-10)
		break;

	      add_to_spec_table (spec, bin - lower, cur_p, quant_cont,
				 w * frac * cur_p->weight);
	      //	      if (cur_p->weight > 10.0){
	      //		fprintf(stdout,"weight: %f, distance: %f, ewidth: %f\n",
	      //			cur_p->weight, cur_p->dist, d+0.5);
	      //	      }
	      sum += w;
	    }

	  /* Now collect contributions upward of cur_p->xi */
	  for (xi = cur_p->xi + 1;; xi += 1)
	    {
	      bin = NAIVE_VAL_TO_BININD (xi);

	      w = PIXWEIGHT (xc + (bin - cur_p->xi) * cosp,
			     yc + (bin - cur_p->xi) * sinp,
			     xc + (bin + 1 - cur_p->xi) * cosp,
			     yc + (bin + 1 - cur_p->xi) * sinp,
			     &pix);
	      if (w < 1e-10)
		break;

	      add_to_spec_table (spec, bin - lower, cur_p, quant_cont,
				 w * frac * cur_p->weight);
	      sum += w;
	    }

	  if (fabs (sum - 1) > 1e-6)
	    {
	      fprintf(stdout,
		      "Weights added up to only %f for pixel from %4d,%4d\n",
		      sum, cur_p->p_x, cur_p->p_y);
	    }
	}
      cur_p++;
    }

  /* Trimming the INDEF beginning and ending values in spectrum */
  tspec = trim_spectrum (spec);
  free_spectrum (spec);
  spec = NULL;

  // return the spectrum
  return tspec;
}

/**
 * Function: bin_optimal
 * computes a spectrum from a table of aperture pixels generated from
 * spc_extract.  This is adds up the pixel values, distributing them
 * over the trace, taking into account the fracton of the pixel that
 * projects onto the given [xi,xi+1] interval.
 * Return NULL if ap_p is NULL
 *
 * Parameters:
 * @param ap_p the table of aperture pixels
 * @param px_width width of a pixel (=1 if not subsampled)
 * @param ob_width the width of the object that has generated the spectrum
 * @param ob_orientation the orientation of the object that has
 *             generated the spectrum
 * @param flags problems that were accumulated in generating ap_p;
 *   possible flags are defined for the warning member of the spectrum struct
 *
 * Returns:
 * @return spec - the 1D spectrum
 */
spectrum *
bin_optimal (const ap_pixel * const ap_p, const beam curbeam,
	     const int quant_cont, const gsl_matrix *weights,
	     const drzstamp_dim dimension,gsl_matrix *coverage)
{

  const ap_pixel *cur_p;
  ap_pixel *tmp_p;
  spectrum *spec;
  spectrum *tspec;
  spc_entry *spec_table;

  quadrangle quad;
  //  drzstamp_dim dimension;

  double jacob, arr;
  double frac, totweight;

  int jcen, icen;
  int jupp, iupp;
  int jlow, ilow;

  int ii, jj;
  int stpi, stpj;

  // return NULL if
  // empty PET
  if (ap_p==NULL)
    return NULL;

  // allocate memory
  tmp_p = (ap_pixel *) malloc(sizeof(ap_pixel));


  // allocate memory for the spectrum
  spec = allocate_spectrum (weights->size1);
  spec_table = spec->spec;


  // go over each PET pixel
  cur_p = ap_p;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {

      // continue if the pixel is outside
      // of the extraction region
      if (fabs (cur_p->dist) > curbeam.width + .5)
	  continue;

      // determine which fraction
      // of the pixel is inside of the extraction area
      if ((fabs (cur_p->dist) >= curbeam.width - .5) && (fabs (cur_p->dist) <= curbeam.width + .5))
	frac = fabs (curbeam.width - (fabs (cur_p->dist) - 0.5));
      else
	frac = 1.;

      // transfer values to the temporary pixel
      tmp_p->lambda  = cur_p->xi;
      tmp_p->dist    = cur_p->dist;
      tmp_p->dxs     = cur_p->dxs;
      tmp_p->dlambda = 1.0;

      // create the quadrangle for the current pixel
      quad = get_quad_from_pixel(tmp_p, curbeam.orient, dimension);

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
	      if ( (stpi>=dimension.xsize)||(stpi<0)||(stpj>=dimension.ysize)||(stpj<0) )
		continue;

	      // get the area which falls onto the current output pixel
	      arr = boxer(stpi,stpj,quad.x,quad.y);

	      // compute the pixel weight from
	      // the various inputs
	      totweight =  arr*frac*jacob*gsl_matrix_get(weights, stpi, stpj);
	      //totweight =  arr*frac*jacob;//*gsl_matrix_get(weights, stpi, stpj);

	      //gsl_matrix_set(coverage, stpi, stpj, gsl_matrix_get(coverage, stpi, stpj) + arr*frac*jacob);

	      // add the pixel to the spectrum
	      if (totweight > 0.0)
		add_to_spec_table (spec, stpi, cur_p, quant_cont,totweight);
	      //	      gsl_matrix_set(coverage, stpi, stpj, sqrt (SQR (gsl_matrix_get(coverage, stpi, stpj)) + SQR (fabs(cur_p->error) * totweight)));
	    }
	}
    }



  /* Trimming the INDEF beginning and ending values in spectrum */
  tspec = trim_spectrum (spec);
  free_spectrum (spec);
  spec = NULL;

  free(tmp_p);

  // return the spectrum
  return tspec;
}



/**
  does a straight forward summation/binning of an aperture pixel table
  with appropriate weights (cf. Hornes 1986)

  @param ap_p the table of aperture pixels
  @param ob_orient the orientation of the object that has
    generated the spectrum
  @param flags problems that were accumulated in generating ap_p;
    possible flags are defined for the warning member of the spectrum struct
  @param n_sub subsampling factor
  @return an allocated spectrum
  @see spectrum
*/

spectrum *
bin_weighted (const ap_pixel * const ap_p, const double ob_orient,
	      const trace_func * const trace, const int n_sub,
	      const int flags)
{
     gsl_vector *binned_table;
     gsl_vector *wei_table;
     gsl_vector *wei2_table;

     double xi, d, w, wei, wei2;
     int xii, num_bin, bin;
     const ap_pixel *cur_p = ap_p;
     double upper = cur_p->xi, lower = cur_p->xi;
     spectrum *spec;

     while (cur_p->p_x != -1)
       {
	    upper = MAX (cur_p->xi, upper);
	    lower = MIN (cur_p->xi, lower);
	    cur_p++;
       }

     lower -= 10;
     upper += 10;

     lower = floor (lower);
     upper = floor (upper + 1);

     num_bin = floor ((upper - lower) / n_sub);
     spec = allocate_spectrum (num_bin);

     binned_table = gsl_vector_alloc (num_bin);
     wei_table = gsl_vector_alloc (num_bin);
     wei2_table = gsl_vector_alloc (num_bin);

     gsl_vector_set_all (binned_table, 0);
     gsl_vector_set_all (wei_table, 0);
     gsl_vector_set_all (wei2_table, 0);

     cur_p = ap_p;
     while (cur_p->p_x != -1)
       {
	    xi = (cur_p->xi - lower) / n_sub;
	    xii = floor (xi);

	    d = fabs(cur_p->dist);
	    w = exp (-d * d / (2 * 6.66));
	    w = 1.;
	    add_to_spec_table (spec, xii - 1,cur_p, 0, w * (1 - (xi - xii)));
	    add_to_spec_table (spec, xii   , cur_p, 0, w * (xi - xii));

	    gsl_vector_set (wei_table, xii - 1,
			    gsl_vector_get (wei_table, xii - 1) + w);
	    gsl_vector_set (wei_table, xii,
			    gsl_vector_get (wei_table, xii) + w);

	    gsl_vector_set (wei2_table, xii - 1,
			    gsl_vector_get (wei2_table, xii - 1) + w * w);
	    gsl_vector_set (wei2_table, xii,
			    gsl_vector_get (wei2_table, xii) + w * w);

	    //fprintf(stderr,"%f,%d,%d\n",xi,xii,xii+1);
	    cur_p++;
       }
     cur_p = ap_p;
     for (bin = 0; bin < num_bin; bin++)
       {
	    wei = gsl_vector_get (wei_table, bin);
	    wei2 = gsl_vector_get (wei2_table, bin);
	    if (wei2 != 0)
	      {
		   spec->spec[bin].count = spec->spec[bin].count * wei / wei2;
	      }
       }

     spec->warning = 0;
     return spec;
}
