/**
* Output routines for stamped images and overlaid traces.
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>

#include "spce_output.h"
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_pathlength.h"
#include "spc_driz.h"
#include "crossdisp_utils.h"


#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define SQR(x) ((x)*(x))


/* transformation from xi to the index into the trace pixel table */
#define XI_QUANT(xi) ((int)rint(((xi)-minxi)/xi_bin_width))

/**
 * Function: stamp_img
 * Produce a gsl array containing a stamp image.
 * Returns NULL is input object is NULL.
 *
 * Parameters:
 * @param ap_p    - a pointer to a -1 terminated ap_pixel array
 * @param width   - extraction width used
 * @param stp_min - xy-coords of the lower left corner on image
 *
 * Returns:
 * @return res  - a pointer to a newly allocated gsl_matrix containing the
 *                spectrum
 */
gsl_matrix * stamp_img (const ap_pixel * const ap_p, float width, d_point *stp_min)
{
  const ap_pixel *cur_p;
  double min_px = 1e32;
  double max_px = -1e32;
  double min_py = 1e32;
  double max_py = -1e32;
  gsl_matrix *res;
  long xsize,ysize,i,j,n=0;

  if (ap_p==NULL) {
    /* Create a dummy stamp image */
    res = gsl_matrix_alloc(10,10);
    /* Fill stamp with 0.0 values */
    gsl_matrix_set_all(res, 0.0);

    return res;
  }

  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      min_px = MIN (min_px, cur_p->p_x);
      max_px = MAX (max_px, cur_p->p_x);
      min_py = MIN (min_py, cur_p->p_y);
      max_py = MAX (max_py, cur_p->p_y);
      n++;
    }

  xsize = max_px-min_px+2;
  ysize = max_py-min_py+2;

  stp_min->x = min_px;
  stp_min->y = min_py;

  if (n>0)
    {
      res = gsl_matrix_alloc(xsize,ysize);
    } else {
      /* Create a dummy stamp image */
      res = gsl_matrix_alloc(10,10);
      /* Fill stamp with 0.0 values */
      gsl_matrix_set_all(res, 0.0);
      return res;
    }

  /* Fill stamp with NaN values */
  gsl_matrix_set_all(res, 0.0);
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      i = (long) floor(cur_p->p_x - min_px+.5);
      j = (long) floor(cur_p->p_y - min_py+.5);
      if ((i<0) || (i>=xsize)) continue;
      if ((j<0) || (j>=ysize)) continue;

      //fprintf(stderr,"i: %d j: %d\n",i,j);
      gsl_matrix_set(res,i,j,cur_p->count);
    }

  return res;
}

/*
 * Function: get_minxy_from_PET
 * The function determines the minimum in
 * x and y of the pixels in a PET
 *
 * Parameters:
 * @param ap_p - the PET vector
 *
 * Returns:
 * @return res - the minima in x, y
 */
d_point
get_minxy_from_PET(const ap_pixel * const ap_p)
{
  const ap_pixel *cur_p;
  double min_px = 1e32;
  //double max_px = -1e32;
  double min_py = 1e32;
  //double max_py = -1e32;
  d_point res;

  // return -1 in case of an empty PET
  if (ap_p==NULL) {
    res.x = -1;
    res.y = -1;
    return res;
  }

  // go over each pixel, determine the maximum
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      min_px = MIN (min_px, cur_p->p_x);
      min_py = MIN (min_py, cur_p->p_y);
    }

  // fill the return
  res.x = min_px;
  res.y = min_py;

  // return the result
  return res;
}


/*
 * Function: get_maxxy_from_PET
 * The function determines the maximum in
 * x and y of the pixels in a PET
 *
 * Parameters:
 * @param ap_p - the PET vector
 *
 * Returns:
 * @return res - the maximum in x, y
 */
d_point
get_maxxy_from_PET(const ap_pixel * const ap_p)
{
  const ap_pixel *cur_p;
  //double min_px = 1e32;
  double max_px = -1e32;
  //double min_py = 1e32;
  double max_py = -1e32;
  d_point res;

  // return -1 in case of an empty PET
  if (ap_p==NULL) {
    res.x = -1;
    res.y = -1;
    return res;
  }

  // go over each pixel, determine the maximum
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      max_px = MAX (max_px, cur_p->p_x);
      max_py = MAX (max_py, cur_p->p_y);
    }

  // fill the return
  res.x = max_px;
  res.y = max_py;

  // return the result
  return res;
}

/**
 * Function: stamp_img_drzprep
 * Creates all drizzle prepare stamp images from a PET
 * vector. In case of an empty PET, dummy stamps are returned.
 *
 * Parameters:
 * @param opt_extr   - flagg to indicate optimal extraction
 * @param ap_p       - a pointer to a -1 terminated ap_pixel array
 * @param width      - extraction width used
 * @param nullval    - the default pixel value
 * @param quant_cont - flag to indicate quantitative cont.
 * @param dimension  - size information
 * @param drzcoeffs  - the drizzle coefficients
 * @param exptime    - the exposure time
 * @param sky_cps    - the sky background
 *
 * Returns:
 * @return drzprep_stamps - a pointer to the new struct with the
 *                          drzprep stamp images
 */
drzprep *
stamp_img_drzprep (const int opt_extr, const ap_pixel *const ap_p, const ap_pixel *const se_p, float width,
		   float nullval, int quant_cont, drzstamp_dim dimension,
		   gsl_matrix *drzcoeffs, double exptime, double sky_cps,
		   double rdnoise, const int bckmode)
{
  drzprep *drzprep_stamps;
  const ap_pixel *cur_p;
  const ap_pixel *xxx_p;

  int i, j;

  double var;
  double corr;
  double sqr_expt;
  double sqr_rdns;

  // allocate space for the output structure
  drzprep_stamps = (drzprep *)malloc(sizeof(drzprep));

  if (ap_p==NULL) {

    // create a dummy count image and set it to default
    drzprep_stamps->counts = gsl_matrix_alloc(10,10);
    gsl_matrix_set_all(drzprep_stamps->counts, nullval);

    // create a dummy error image and set it to default
    drzprep_stamps->error = gsl_matrix_alloc(10,10);
    gsl_matrix_set_all(drzprep_stamps->error, nullval);

    // create a dummy contamination image and set it to default
    drzprep_stamps->cont = gsl_matrix_alloc(10,10);
    gsl_matrix_set_all(drzprep_stamps->cont, nullval);

    // create a dummy model image and set it to default
    drzprep_stamps->model = gsl_matrix_alloc(10,10);
    gsl_matrix_set_all(drzprep_stamps->model, nullval);

    // create a dummy model image and set it to default
    drzprep_stamps->vari = gsl_matrix_alloc(10,10);
    gsl_matrix_set_all(drzprep_stamps->vari, nullval);

    return drzprep_stamps;
  }

  // allocate the count matrix and fill it with default value
  drzprep_stamps->counts  = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
  gsl_matrix_set_all(drzprep_stamps->counts, nullval);

  // allocate the error matrix and fill it with default value
  drzprep_stamps->error  = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
  gsl_matrix_set_all(drzprep_stamps->error, 0.0);

  // allocate the contamination matrix and fill it with default value
  drzprep_stamps->cont  = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
  gsl_matrix_set_all(drzprep_stamps->cont, 0.0);

  // set the two extentions
  // which are only somethimes
  // used to NULL
  drzprep_stamps->model  = NULL;
  drzprep_stamps->vari  = NULL;

  if (opt_extr)
    {
      // allocate the mode; matrix and fill it with default value
      drzprep_stamps->model  = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
      gsl_matrix_set_all(drzprep_stamps->model, 0.0);

      // allocate the mode; matrix and fill it with default value
      drzprep_stamps->vari  = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
      gsl_matrix_set_all(drzprep_stamps->vari, 0.0);
    }

  // square the exposure time for speed
  sqr_expt = exptime*exptime;
  sqr_rdns = rdnoise*rdnoise;

  // go over each PET pixel
  if (opt_extr && bckmode)
    xxx_p = se_p;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {

      // compute the coordinates in the stamp images
      // if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      i = (int)(cur_p->p_x - dimension.xstart);
      j = (int)(cur_p->p_y - dimension.ystart);


      // check whether the coordinates are inside
      if ((i<0) || (i>=dimension.xsize) || (j<0) || (j>=dimension.ysize))
	{
	  if (opt_extr && bckmode)
	    xxx_p++;
	continue;
	}

      // comute the jacobian of the drizzle coefficients
      corr = fabs(get_det_jacobian(i+1, j+1, drzcoeffs, (int)dimension.xsize,(int)dimension.ysize));

      // compute the count value
      gsl_matrix_set(drzprep_stamps->counts,i,j,(cur_p->count)/corr);

      // MODIFIED in Dec. 2010 in order
      // to get error propagation
      // compute the true error value in electrons
      gsl_matrix_set(drzprep_stamps->error,i,j,(cur_p->error*cur_p->error*sqr_expt)/corr/corr);

      // compute the contamination value
      if (quant_cont)
	gsl_matrix_set(drzprep_stamps->cont,i,j,(cur_p->contam)/corr);
      else
	gsl_matrix_set(drzprep_stamps->cont,i,j,(cur_p->contam));

      if (opt_extr)
	{

	  if (bckmode)
	    {
	      // check that the foreground and background
	      // PET element describe the same pixel
	      if (xxx_p->x != cur_p->x || xxx_p->y != cur_p->y)
		aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			     "aXe_DRZPREP: Background PET and Object PET "
			     "have different pixel orders in PET's.\n");

	      // compute the inverse variance value for the object PET
	      // derived from the noise characteristics
	      var = ((xxx_p->model+xxx_p->contam+sky_cps+cur_p->count)*exptime + sqr_rdns)/sqr_expt;
	      //	      var = cur_p->count*exptime / sqr_expt;

	      // compute the model value
	      gsl_matrix_set(drzprep_stamps->model,i,j,(xxx_p->model/corr));
	    }
	  else
	    {
	      // compute the inverse variance value for the object PET
	      // derived from the noise characteristics
	      var = ((cur_p->model+cur_p->contam+sky_cps)*exptime + sqr_rdns)/sqr_expt;

	      // compute the model value
	      gsl_matrix_set(drzprep_stamps->model,i,j,(cur_p->model/corr));
	    }

	  if (var != 0.0)
	    gsl_matrix_set(drzprep_stamps->vari,i,j,1.0/var);


	}

      if (opt_extr && bckmode)
	xxx_p++;
    }

  // return the result structure
  return drzprep_stamps;
}


/**
 * Function: rectified_stamp_img
 * Produce a gsl array containing a 'rectified aperture'. That is,
 * an image of the spectum showing count on an xi,dist grid instead
 * of the traditional image grid x,y. Returns NULL is input object is NULL.
 *
 * Parameters:
 * @param ap_p  - a pointer to a -1 terminated ap_pixel array
 * @param width - extraction width used
 *
 * Returns:
 * @return res - a pointer to a newly allocated gsl_matrix containing the
 * rectified spectrum
 */
gsl_matrix *
rectified_stamp_img (const ap_pixel * const ap_p, float width, d_point *stp_min)
{
  const ap_pixel *cur_p;
  double min_xi = 1e32;
  double max_xi = -1e32;
  double min_dist = 1e32;
  double max_dist = -1e32;
  gsl_matrix *res;
  long xsize,ysize,i,j,n=0;

  if (ap_p==NULL) {
    /* Create a dummy stamp image */
    res = gsl_matrix_alloc(10,10);
    /* Fill stamp with 0.0 values */
    gsl_matrix_set_all(res, 0.0);
    return res;
  }

  i=0;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      i++;
      min_xi = MIN (min_xi, cur_p->xi);
      max_xi = MAX (max_xi, cur_p->xi);
      min_dist = MIN (min_dist, cur_p->dist);
      max_dist = MAX (max_dist, cur_p->dist);
      n++;
    }

  if (i==0) {
    /* Create a dummy stamp image */
    res = gsl_matrix_alloc(10,10);
    /* Fill stamp with 0.0 values */
    gsl_matrix_set_all(res, 0.0);
    return res;
  }

  min_dist = floor(min_dist);
  max_dist = floor(max_dist);
  min_xi = floor(min_xi);
  max_xi = floor(max_xi);

  // compute the array dimensions
  xsize = max_xi-min_xi+2;
  ysize = max_dist-min_dist+2;

  // store the minima in both coords
  stp_min->x = min_xi;
  stp_min->y = min_dist;

  if (n>0)
    {
      res = gsl_matrix_alloc(xsize,ysize);
    } else {
      /* Create a dummy stamp image */
      res = gsl_matrix_alloc(10,10);
      /* Fill stamp with 0.0 values */
      gsl_matrix_set_all(res, 0.0);
      return res;
    }

  res = gsl_matrix_alloc(xsize,ysize);
  /* Fill stamp with NaN values */
  gsl_matrix_set_all(res, GSL_NAN);


  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      i = (long) floor(cur_p->xi - min_xi+.5);
      j = (long) floor(cur_p->dist - min_dist+.5);
      //fprintf(stderr,"i: %d j: %d\n",i,j);
      gsl_matrix_set(res,i,j,cur_p->count);
    }

  return res;
}

/**
 * Function: drizzled_stamp_img
 * Produce a structure containing a drizzled stamp image plux ther
 * associated weight image. The drizzled stamp image has constant
 * dispersion along the x-axis. No drizzling is applied
 * in cross-dispersion direction. The dimensions of the drizzled
 * stampimage were setr previously, and the information
 * is in an input structure. In case that there are no valid pixels,
 * a dummy structure with 10x10 empty pixels is returned.
 *
 * @param ap_p      - pointer to the list of PET pixels
 * @param width     - extraction width used
 * @param dimension - dimensional information on the drizzled stamp
 *
 * @return res      - structure with counts and weights of the drizzled stamp
*/
drzstamp *
drizzled_stamp_img (const  ap_pixel * const ap_p, double width,
		     double orient, const drzstamp_dim dimension)
{
  const ap_pixel *cur_p;

  quadrangle quad;

  drzstamp *res;
  gsl_matrix *counts;
  gsl_matrix *weight;

  int icen, ilow, iupp;
  int jcen, jlow, jupp;

  double xi, jacob=0;
  double value, allweig, weig;

  double maxarr=0.0;
  int iim, jjm;

  //double totweigth;
  //double cos_phi=0.0;

  int ii,jj;
  double arr;
  int stpi,stpj;
  double stpc;

  // allocate space for the result
  res = (drzstamp *) malloc(sizeof(drzstamp));

  // handle an empty PET
  if (ap_p==NULL || !dimension.resolution) {
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


  /* Allocate the stamp matrix*/
  /* Fill stamp with NaN values */
  counts = gsl_matrix_alloc(dimension.xsize,dimension.ysize);
  gsl_matrix_set_all(counts, GSL_NAN);

  /* Allocate the weight matrix*/
  /* Fill weight with 0.0 values */
  weight = gsl_matrix_alloc(dimension.xsize, dimension.ysize);
  gsl_matrix_set_all(weight, 0.0);


  // go over each pixel
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {

      // Skip this pixel if it was not actually used
      if (fabs(cur_p->dist)>width+1.5)
	continue;

      // create the quadrangle for the current pixel
      //      quad = get_quad_from_pixel2(cur_p, dimension);
      quad = get_quad_from_pixel(cur_p, orient, dimension);

      // get the jacobian (well, easy here)
      // the term "cos(cur_p->dxs)" must be there
      // to correct the enlargement necessary
      // to cover the whole lambda-cross dispersion area!
      // NOT COMPLETELY understood
      jacob = dimension.resolution/cur_p->dlambda*cos(cur_p->dxs);

      // get the central pixel (icen, jcen) of the current PET-pixel
      xi = cur_p->lambda/dimension.resolution;
      icen = (int) floor(xi - dimension.xstart+.5);
      jcen = (int) floor(cur_p->dist - dimension.ystart+.5);

      // get the uper and lower extend of the quadrangle in x
      iupp = (int)floor(quad.xmax - (double)icen + 0.5)+1;
      ilow = (int)floor(quad.xmin - (double)icen + 0.5);

      // get the uper and lower extend of the quadrangle in x
      jupp = (int)floor(quad.ymax - (double)jcen + 0.5)+1;
      jlow = (int)floor(quad.ymin - (double)jcen + 0.5);

      maxarr=0.0;
      //      totweigth = 0.0;
      // go over the extend in x
      for (ii=ilow;ii<iupp;ii++) {
      	// go over the extend in x
      	for (jj=jlow;jj<jupp;jj++) {

      	  // get the coordinates of the current output pixel
      	  stpi = icen+ii;
      	  stpj = jcen+jj;

      	  // check whether the current output pixel is within
      	  // the stamp image; continue if not
      	  if ( (stpi>=dimension.xsize)||(stpi<0)||(stpj>=dimension.ysize)||(stpj<0) )
      	    continue;

      	  // get the area which falls onto the current output pixel
      	  arr = boxer(stpi,stpj,quad.x,quad.y);
      	  if (arr > maxarr)
      	    {
      	      maxarr=arr;
      	      iim = ii;
      	      jjm = jj;
      	    }
      	  // get the already existing counts and weights
      	  stpc = gsl_matrix_get(counts,stpi,stpj);
      	  weig = gsl_matrix_get(weight,stpi,stpj);

      	  // initialize the counts, if necessary
      	  if (isnan(stpc) && (arr!=0.0))
      	    stpc = 0.0;

      	  // compute the new, total weight of the current output pixel
      	  //	  allweig = weig + jacob*arr;
      	  allweig = weig + arr;

      	  // do a weighted sum of the count value at the current output pixel
      	  value = (stpc*weig + arr*cur_p->count*jacob) / (allweig);

      	  // store the new count value and the new weight
      	  gsl_matrix_set(counts,stpi,stpj,value);
      	  gsl_matrix_set(weight,stpi,stpj,allweig);
      	  //	  totweigth = totweigth + arr;
      	}
      }
    }

  // fill the output structure
  res->counts = counts;
  res->weight = weight;

  // return the output
  return res;
}

/**
 * Function: get_drzprep_dim
 * The function computes the dimensional properties
 * of the stamp images created in the drizzle prepare task.
 * Those properties fix the image size
 * and the offsets in both coordinates.
 *
 * Parameter:
 * @param ap_p        - the list of beam pixels to find the diemnsions for
 * @param width       - the extraction width of the beam
 * @param boxwidth    - the boxsize
 * @param boxheight   - the boxwidth
 *
 * Returns:
 * @return dimensions - the filled quadrangle sturcture
 */
drzstamp_dim
get_drzprep_dim(const ap_pixel *const ap_p, float width,
		int boxwidth, int boxheight)
{

  drzstamp_dim dimensions;

  const ap_pixel *cur_p;

  double min_px = 1e32;
  double max_px = -1e32;
  double min_py = 1e32;
  double max_py = -1e32;

  // return a dummy if the PET is NULL
  if (ap_p==NULL)
    {
      dimensions.xstart=0.0;
      dimensions.ystart=0.0;

      dimensions.xsize=0;
      dimensions.ysize=0;

      return dimensions;
    }

  // go over each pixel, find minimum and maximum in x,y
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      // if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      min_px = MIN (min_px, cur_p->p_x);
      max_px = MAX (max_px, cur_p->p_x);
      min_py = MIN (min_py, cur_p->p_y);
      max_py = MAX (max_py, cur_p->p_y);
    }

  // compute the dimension of the stamps
  dimensions.xsize = max_px-min_px+1;
  dimensions.ysize = max_py-min_py+1;

  // correct the dimesions if necessary
  if (boxwidth > dimensions.xsize){
    dimensions.xsize = (long)boxwidth;
  }
  if (boxheight > dimensions.ysize){
    dimensions.ysize = boxheight;
  }

  // fix the start coordinates
  dimensions.xstart= min_px;
  dimensions.ystart= min_py;

  // return the structure
  return dimensions;
}

/**
 * Function: get_stamp_dim
 * The function computes the dimensional properties
 * of the drizzled stamp image associated to a beam.
 * Those properties fix the image size, dispersion
 * and the offsets in both coordinates
 *
 * Parameters:
 * @param ap_p        - the list of beam pixels to find the diemnsions for
 * @param width       - the extraction width of the beam
 * @param conf        - the configuration structure
 * @param beamID      - the the beam ID
 *
 * Returns:
 * @return dimensions - the filled quadrangle sturcture
 */
drzstamp_dim
get_stamp_dim(const ap_pixel *const ap_p, float width,
	      aperture_conf *conf, const int beamID, d_point *stp_min)
{
  drzstamp_dim dimensions;

  const ap_pixel *cur_p;

  double min_lam  = 1e32;
  double max_lam  = -1e32;
  double min_xi   = 1e32;
  double max_xi   = -1e32;
  double min_dist = 1e32;
  double max_dist = -1e32;
  double min_dlam = 1e32;
  double max_dlam = -1e32;

  //double resampwidth = 0.0;

  //long xsize,ysize;
  long npixel=0;
  //long j;

  // check whether there are entries;
  // if not, return a dummy structure
  if (ap_p==NULL)
    {
      dimensions.resolution=0.0;

      dimensions.xstart=0.0;
      dimensions.ystart=0.0;

      dimensions.xsize=0;
      dimensions.ysize=0;

      return dimensions;
    }


  // reset the pixel counter
  npixel=0;

  // go over each pixel
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {

      // Skip the pixel if it was not actually used
      if (fabs(cur_p->dist)>width+.5)
	continue;

      // search for minimum/maximum in wavelength
      min_lam = MIN (min_lam, cur_p->lambda);
      max_lam = MAX (max_lam, cur_p->lambda);

      // search for minimum/maximum in crossdispersion direction
      min_dist = MIN (min_dist, cur_p->dist);
      max_dist = MAX (max_dist, cur_p->dist);

      // search for minimum/maximum in dispersion
      min_dlam = MIN (min_dlam, cur_p->dlambda);
      max_dlam = MAX (max_dlam, cur_p->dlambda);

      // enhance the pixel counter
      npixel++;
    }

  // in case there were no valid pixels,
  // return a dummy structure
  if (!npixel)
    {
      dimensions.resolution=0.0;

      dimensions.xstart=0.0;
      dimensions.ystart=0.0;

      dimensions.xsize=0;
      dimensions.ysize=0;

      return dimensions;
    }

  // fill the resolution value
  // if the drizzle resolution is defined in the
  // configuration file and it is a first order beam
  if (conf->drz_resol && !beamID)
    // take the drizzle value from the config file
    dimensions.resolution = (double)conf->drz_resol;
  // otherwise
  else
    // take the average resoution
    dimensions.resolution = (max_dlam+min_dlam)/2.0;

  // round up and off in crossdispersion direction
  min_dist = floor(min_dist);
  max_dist = ceil(max_dist);

  // round up and off in dispersion direction
  min_xi = floor(min_lam/dimensions.resolution);
  max_xi = ceil(max_lam/dimensions.resolution);

  // derive and store the stamp image dimensions
  dimensions.xsize = max_xi-min_xi+2;
  dimensions.ysize = max_dist-min_dist+2;

  // store the start values
  dimensions.xstart = min_xi;
  dimensions.ystart = min_dist;

  stp_min->x = min_xi * dimensions.resolution;
  stp_min->y = min_dist;

  // return the structure
  return dimensions;
}

/**
 * Function: get_quad_from_pixel
 * The subroutine creates a quadrangle for a PET pixel
 * in the coordinate system of a drizzled stam image.
 * The quadrangle later serves as an input to the boxer
 * routine to resample the PET pixel.
 * Form the math point of view, the generation
 * of the quadrange corresponds to executing an affine
 * transformation of the pixel (x,y) corners into
 * the coordinate system spanned by the wavelength
 * and the crossdispersion direction.
 *
 * Parameters:
 * @param cur_p     - the PET pixel to construct the quadrangle for
 * @param dimension - the dimensions of the stamp image
 *
 * Returns:
 * @return quad - the filled quadrangle sturcture
 */
quadrangle
get_quad_from_pixel(const ap_pixel *cur_p, const double orient,
		    const drzstamp_dim dimension)
{
  quadrangle quad;

  double dxi, ddist;
  double x, y;

  double phi_1, phi_2;
  double cos_phi_1, sin_phi_1, tan_phi_2;
  double term_1, term_2;

  // nomenclature: see contamination.dvi
  phi_1 = orient - 1.5707963267948966;  // crossdispersion direction corresponds to beta; numeric=90deg in rad
  phi_2 = cur_p->dxs - orient + 1.5707963267948966; // traceangle corresponds to alpha

  // nomenclature: see contamination.dvi
  cos_phi_1 = cos(phi_1);
  sin_phi_1 = sin(phi_1);
  tan_phi_2 = tan(phi_2);

  // nomenclature: see contamination.dvi
  term_1 = cos_phi_1*tan_phi_2 + sin_phi_1;
  term_2 = cos_phi_1           - sin_phi_1*tan_phi_2;

  // Bottom left corner (-.5, -.5):
  // get the coos of that corner in the stamp image coo system
  x = -0.5;
  y = -0.5;
  dxi   = +x*cos_phi_1 + y*sin_phi_1;
  ddist = -x*term_1    + y*term_2;
  /*** IMPORTANT:
   * the division by "cos(cur_p->dxs)" in the equations below is
   * is necessary to cover the whole lambda-crossdisp plane.
   * However then the counts are not anymore preserved.
   * This must be taken into account in the calling routine!!
   **/
  quad.x[0] = (cur_p->lambda + cur_p->dlambda*dxi/cos(cur_p->dxs))/dimension.resolution - dimension.xstart;
  quad.y[0] = cur_p->dist + ddist - dimension.ystart;
  //  fprintf(stdout, "%f %f\n", quad.x[0], quad.y[0]);

  // Top left corner (-.5, +.5):
  // get the coos of that corner in the stamp image coo system
  x = -0.5;
  y = +0.5;
  dxi   = +x*cos_phi_1 + y*sin_phi_1;
  ddist = -x*term_1    + y*term_2;
  //  quad.x[1] = cur_p->xi + dxi/cos(cur_p->dxs);
  quad.x[1] = (cur_p->lambda + cur_p->dlambda*dxi/cos(cur_p->dxs))/dimension.resolution - dimension.xstart;
  quad.y[1] = cur_p->dist + ddist - dimension.ystart;

  // Top right corner (+.5, +.5):
  // get the coos of that corner in the stamp image coo system
  x = +0.5;
  y = +0.5;
  dxi   = +x*cos_phi_1 + y*sin_phi_1;
  ddist = -x*term_1    + y*term_2;
  quad.x[2] = (cur_p->lambda + cur_p->dlambda*dxi/cos(cur_p->dxs))/dimension.resolution - dimension.xstart;
  quad.y[2] = cur_p->dist + ddist - dimension.ystart;

  // Bottom right corner (+.5, -.5):
  // get the coos of that corner in the stamp image coo system
  x = +0.5;
  y = -0.5;
  dxi   = +x*cos_phi_1 + y*sin_phi_1;
  ddist = -x*term_1    + y*term_2;
  quad.x[3] = (cur_p->lambda + cur_p->dlambda*dxi/cos(cur_p->dxs))/dimension.resolution - dimension.xstart;
  quad.y[3] = cur_p->dist + ddist - dimension.ystart;

  // get the maximum and minimum of the quadrangle
  // in each dimension
  quad.xmax = MAX(quad.x[3],MAX(quad.x[2],MAX(quad.x[1],quad.x[0])));
  quad.xmin = MIN(quad.x[3],MIN(quad.x[2],MIN(quad.x[1],quad.x[0])));
  quad.ymax = MAX(quad.y[3],MAX(quad.y[2],MAX(quad.y[1],quad.y[0])));
  quad.ymin = MIN(quad.y[3],MIN(quad.y[2],MIN(quad.y[1],quad.y[0])));

  return quad;
}

/**
    Produce a gsl array containing a 'drizzled aperture'. That is,
    an image of the spectum showing count on an xi,dist grid instead
    of the traditional image grid x,y. Returns NULL is input object is NULL.
	[modifies ap_p!!]
    @param ap_p a pointer to a -1 terminated ap_pixel array
    @param width extraction width used
	@param resampwidth A per pixel of output drizzled spectra
    @return a pointer to a newly allocated gsl_matrix containing the
    rectified spectrum
*/
gsl_matrix *
drizzled_stamp_img_orig (const  ap_pixel * const ap_p, float width, aperture_conf *conf)
{
  const ap_pixel *cur_p;

  double min_xi   = 1e32;
  double max_xi   = -1e32;
  double min_dist = 1e32;
  double max_dist = -1e32;
  double min_dlam = 1e32;
  double max_dlam = -1e32;

  double resampwidth = 0.0;

  gsl_matrix *res;

  long xsize,ysize,i,j,n=0;

  double xi;

  if (ap_p==NULL) {
    /* Create a dummy stamp image */
    res = gsl_matrix_alloc(10,10);
    /* Fill stamp with 0.0 values */
    gsl_matrix_set_all(res, 0.0);
    return res;
  }

  i=0;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist)>width+.5) continue; /* Skip this pixel if it was not actually used */
      i++;

      xi = cur_p->lambda;

      min_xi = MIN (min_xi, xi);
      max_xi = MAX (max_xi, xi);

      min_dist = MIN (min_dist, cur_p->dist);
      max_dist = MAX (max_dist, cur_p->dist);

      min_dlam = MIN (min_dlam, cur_p->dlambda);
      max_dlam = MAX (max_dlam, cur_p->dlambda);

      n++;
    }

  if (conf->drz_resol)
    {
      resampwidth = (double)conf->drz_resol;
    }
  else
    {
      resampwidth = (max_dlam+min_dlam)/2.0;
    }

  fprintf(stdout, "resamplewidth: %f\n", resampwidth);

  if (i==0) {
    /* Create a dummy stamp image */
    res = gsl_matrix_alloc(10,10);
    /* Fill stamp with 0.0 values */
    gsl_matrix_set_all(res, 0.0);
    return res;
  }

  min_dist = floor(min_dist);
  max_dist = floor(max_dist);
  min_xi = floor(min_xi/resampwidth);
  max_xi = floor(max_xi/resampwidth);

  xsize = max_xi-min_xi+2;
  ysize = max_dist-min_dist+2;

  if (n>0)
    {
      res = gsl_matrix_alloc(xsize,ysize);
    } else {
      /* Create a dummy stamp image */
      res = gsl_matrix_alloc(10,10);
      /* Fill stamp with 0.0 values */
      gsl_matrix_set_all(res, 0.0);
      return res;
    }

  res = gsl_matrix_alloc(xsize,ysize);
  /* Fill stamp with NaN values */
  gsl_matrix_set_all(res, GSL_NAN);


  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      int ii,jj;
      double x,y,xp,yp,arr;
      double xx[4],yy[4];
      int stpi,stpj;
      double stpc;

      /* Skip this pixel if it was not actually used */
      if (fabs(cur_p->dist)>width+.5)
	continue;


      //cur_p->dxs=0.;
      //cur_p->xi = 100.5 + min_xi;
      //cur_p->dist = 10. + min_dist;
      xi = cur_p->lambda/resampwidth;

      i = (long) floor(xi - min_xi+.5);
      j = (long) floor(cur_p->dist - min_dist+.5);

      /* The center of this PET entry should fall in pixel (i.j) */
      /* Bottom left corner (-.5, -.5) */
      x = -0.5;
      y = -0.5;
      xp =  x*cos(cur_p->dxs) + y*sin(cur_p->dxs);
      yp = -x*sin(cur_p->dxs) + y*cos(cur_p->dxs);
      xx[0] = xp+xi - min_xi;
      yy[0] = yp+cur_p->dist - min_dist;

      /* Top left corner (-.5, +.5) */
      x = -0.5;
      y = +0.5;
      xp =  x*cos(cur_p->dxs) + y*sin(cur_p->dxs);
      yp = -x*sin(cur_p->dxs) + y*cos(cur_p->dxs);
      xx[1] = xp+ xi - min_xi;
      yy[1] = yp+cur_p->dist - min_dist;

      /* Top right corner (+.5, +.5) */
      x = +0.5;
      y = +0.5;
      xp =  x*cos(cur_p->dxs) + y*sin(cur_p->dxs);
      yp = -x*sin(cur_p->dxs) + y*cos(cur_p->dxs);
      xx[2] = xp+ xi - min_xi;
      yy[2] = yp+cur_p->dist - min_dist;

      /* Bottom right corner (+.5, -.5) */
      x = +0.5;
      y = -0.5;
      xp =  x*cos(cur_p->dxs) + y*sin(cur_p->dxs);
      yp = -x*sin(cur_p->dxs) + y*cos(cur_p->dxs);
      // fprintf(stderr,"1 %f %f\n",cur_p->xi - min_xi,cur_p->dist - min_dist);
      xx[3] = xp+ xi - min_xi;
      yy[3] = yp+cur_p->dist - min_dist;

      for (ii=-1;ii<=1;ii++) {
	for (jj=-1;jj<=1;jj++) {
	  //fprintf(stderr,"PET: XI: %f DIST: %f DXS: %f cos:%f\n",cur_p->xi,cur_p->dist,cur_p->dxs,cos(cur_p->dxs));
	  arr = boxer((int) floor(xi - min_xi)+ii,(int) floor(cur_p->dist - min_dist)+jj,xx,yy);

	  //fprintf(stderr,"STP: i:%d j:%d arr: %f\n",(int) floor(cur_p->xi - min_xi)+ii,(int) floor(cur_p->dist - min_dist)+jj,arr);
	  stpi = (int) floor(xi - min_xi)+ii;
	  stpj = (int) floor(cur_p->dist - min_dist)+jj;
	  if ( (stpi>=xsize)||(stpi<0)||(stpj>=ysize)||(stpj<0) ) continue;
	  stpc = gsl_matrix_get(res,stpi,stpj);
	  if (isnan(stpc) && (arr!=0.0)) stpc = 0.0;
	  gsl_matrix_set(res,stpi,stpj,stpc+arr*cur_p->count);
	}
      }


    }

  return res;
}


/**
 * Function: interpolate_over_NaN
 * Function to interpolate over and replace NaN values in a rectified
 * stamp image. This is done by interpolating in the vertical direction.
 * The original gsl_matrix is free'd and is replaced by a new gsl_matrix pointer.
 * Does nothing if input data is NULL.
 *
 * Parameters:
 * @param data - the gsl_matrix where NaN values are to be interpolated over
 */
void interpolate_over_NaN (gsl_matrix *data)
{
  double *x,*y;
  int i,j,n;
  gsl_interp_accel *acc;
  gsl_spline *spline ;

  if (data==NULL) return;

  x = malloc(sizeof(double)*data->size2*data->size1);
  y = malloc(sizeof(double)*data->size2*data->size1);

  for (i=0;i<(int)data->size1;i++) {
    n = 0;
    for (j=0;j<(int)data->size2;j++) {
      if ( !isnan(gsl_matrix_get(data,i,j)) )  {
	       x[n] = j;
	       y[n] = gsl_matrix_get(data,i,j);
	       n++;
      }
    }
    //if (n<data->size2) fprintf(stderr,"%d elements non-NaN found out of %d.\n",n,data->size2);
    if (n<3) continue;

    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, n);

    /* Interpolating */
    gsl_spline_init (spline, x, y, n);

    for (j=0;j<(int)data->size2;j++) {
      double xi = j;
      //double yi = gsl_spline_eval(spline, xi, acc);
      // revert to older extrapolation allowing functionality
      double yi = spline->interp->type->eval(spline->interp->state, spline->x, spline->y, spline->interp->size, x[j], acc, &y[j]);
      //if (n<data->size2)  fprintf(stderr,"%d %g => %d %g\n",j,gsl_matrix_get(data,i,j),j,yi);
      gsl_matrix_set(data, i, j, yi);
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
  }

  free(x);
  free(y);
}

/*
 * Function: free_drzprep
 * The function releases the memory allocated
 * to the various components of a drzprep structure
 *
 * Parameters:
 * @param: drzprep_stamps - the stamp image to free
 */
void
free_drzprep(drzprep *drzprep_stamps)
{
  // release the individual matrixes
  gsl_matrix_free(drzprep_stamps->counts);
  gsl_matrix_free(drzprep_stamps->error);
  gsl_matrix_free(drzprep_stamps->cont);
  if (drzprep_stamps->model)
    gsl_matrix_free(drzprep_stamps->model);
  if (drzprep_stamps->vari)
    gsl_matrix_free(drzprep_stamps->vari);

  // release the struct
  free(drzprep_stamps);
}


/**
 * Function: free_stamp_img
 * Function to free the memeory allocated by a stamp image
 *
 * Parameters:
 * @param rstamp - the stamp image
 */
void free_stamp_img(gsl_matrix *rstamp)
{
  // simply free it and "basta"!
  gsl_matrix_free(rstamp);
}

/**
 * Function: free_drzstamp
 * Function to free the memeory allocated by a
 * drizzled stamp image.
 *
 * Parameters:
 * @param rstamp  - the structure for the drizzled stamp image
 */
void free_drzstamp(drzstamp *rstamp)
{
  // release each component
  gsl_matrix_free(rstamp->counts);
  gsl_matrix_free(rstamp->weight);

  // release the rest
  free(rstamp);
}

/**
  Find the indices of pixels corresponding to the trace (i.e. the
  one with the smallest dist in each column).  If we get
  traces with a slope of significantly more than one, a similar
  function would have to be written operating on rows.

  @param ap_p the aperture pixel table
  @return a gsl_matrix with the indices of the trace pixels in ap_p
*/
gsl_vector_int *
get_trace_inds (const ap_pixel * const ap_p)
{
  const ap_pixel *cur_p;
  int trace_resolution = 10000;	/* Number of xi bins in which to follow the
				   trace */
  gsl_vector_int *trace_inds = gsl_vector_int_alloc (trace_resolution);
  gsl_vector_int *trace_inds2;
  int n, i;
  double distmin;

  gsl_vector_int_set_all (trace_inds, -1);

  /* Find the minimum pixel-to-trace distance in this aperture */
  distmin = 1e32;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      distmin = MIN (distmin, fabs(cur_p->dist));
    }

  n = 0;
  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {
      if (fabs(cur_p->dist) < (distmin + .25))
	{
	  gsl_vector_int_set (trace_inds, n, cur_p - ap_p);
	  n++;
	}
    }
  trace_inds2 = gsl_vector_int_alloc (n);
  n = 0;
  for (i = 0; i < (int)trace_inds->size; i++)
    {
      if (gsl_vector_int_get (trace_inds, i) != -1)
	{
	  gsl_vector_int_set (trace_inds2, n,
			      gsl_vector_int_get (trace_inds, i));
	  n++;
	}
    }

  return trace_inds2;
}
