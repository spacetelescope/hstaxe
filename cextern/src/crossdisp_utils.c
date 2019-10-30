/**
 */
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multiroots.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_trace_functions.h"
#include "disp_conf.h"
#include "trace_conf.h"
#include "crossdisp_utils.h"


/*
 * Function: print_state
 * The function prints a status report
 * on a multi-dimensional root finding
 * solver onto the screen
 *
 * Parameters:
 * @param iter - the iteration number
 * @param s    - the multi-d solver
 */

void print_state (size_t iter, gsl_multiroot_fsolver *s)
{
  printf ("iter = %3zu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
  // print: current x-value, y-value, x-error, y-error
}


/*
 * Function: drizzle_distort
 * This function can be passed to a gsl multi-d function.
 * It defines the function:
 * x'=f'(x,y) = f(x,y)-x_0
 * y'=g'(x,y) = g(x,y)-y_0
 * with f(x,y) and g(x,y) the equations to undistort
 * the coordinates x and y, respectively and the
 * constants x_0 and y_0.
 * Solving for f'(x,y)=g'(x,y)=0.0
 * means to find (x,y) such that
 * f(x,y) = x_0 and g(x,y) = y_0.
 *
 * Parameters:
 * @param x      - the input x,y
 * @param params - the parameters for the function
 * @param f      - the result values in x,y
 *
 * Returns:
 * @return GSL_SUCCESS - always
 */
int
drizzle_distort(const gsl_vector * x, void *params, gsl_vector * f)
{
  // extract the components from the parameter-structure
  gsl_matrix *coeffs = ((drz_pars *) params)->coeffs;
  px_point  npixels = ((drz_pars *) params)->npixels;
  d_point offset = ((drz_pars *) params)->offset;

  d_point xy_in;
  d_point xy_out;

  // extract the input
  xy_in.x = gsl_vector_get(x,0);
  xy_in.y = gsl_vector_get(x,1);

  // undistort the input
  xy_out = get_drz_position(xy_in, coeffs, npixels);

  // subtract the input
  xy_out.x = xy_out.x-offset.x;
  xy_out.y = xy_out.y-offset.y;

  // set the output
  gsl_vector_set(f, 0, xy_out.x);
  gsl_vector_set(f, 1, xy_out.y);

  // return always the same
  return GSL_SUCCESS;
}


/*
 * Function: distort_point
 * The function finds the distorted coordinates
 * for a given undistorted coordinate and the
 * transformations to get the undistorted coordinates.
 * This function makes the inverse transformation
 * to the drizzle transformation.
 *
 * Parameters:
 * @param coeffs   - the drizzle coefficients
 * @param pixmax   - the image dimensions
 * @param xy_image - the undistorted (x,y)
 *
 * Returns:
 * @return xy_ret - the distorted (x,y)
 */
d_point
distort_point(gsl_matrix *coeffs, const px_point pixmax, d_point xy_image)
{
  const gsl_multiroot_fsolver_type *msolve_type;
  gsl_multiroot_fsolver *msolve;

  int status;
  size_t iter = 0;

  //const size_t n = 2;

  gsl_multiroot_function mult_func;

  drz_pars   *drzpars;
  gsl_vector *xy_in = gsl_vector_alloc(2);
  gsl_vector *xy_out = gsl_vector_alloc(2);
  d_point xy_ret;

  // set up the parameters for the
  // multi-d function
  drzpars = (drz_pars *) malloc(sizeof(drz_pars));
  drzpars->coeffs = coeffs;
  drzpars->offset = xy_image;
  drzpars->npixels = pixmax;

  // set up the multi-d function
  mult_func.f = &drizzle_distort;
  mult_func.n = 2;
  mult_func.params = drzpars;

  // set the starting coordinates
  gsl_vector_set(xy_in, 0, xy_image.x);
  gsl_vector_set(xy_in, 1, xy_image.y);

  // allocate and initialize the multi-d solver
  msolve_type = gsl_multiroot_fsolver_dnewton;
  msolve = gsl_multiroot_fsolver_alloc (msolve_type, 2);
  gsl_multiroot_fsolver_set (msolve, &mult_func, xy_in);

  //  print_state (iter, msolve);

  // iterate
  do
  {
    // count the number of iterations
    iter++;

    // do an iteration
    status = gsl_multiroot_fsolver_iterate (msolve);

    //    print_state (iter, msolve);

    // check if solver is stuck
    if (status)
      break;

    // evaluate the iteration
    status = gsl_multiroot_test_residual (msolve->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);
  // chek for the break conditions

  // transfer the result to the return struct
  xy_ret.x = gsl_vector_get(msolve->x,0);
  xy_ret.y = gsl_vector_get(msolve->x,1);

  // deallocate the different structures
  gsl_multiroot_fsolver_free (msolve);
  gsl_vector_free(xy_in);
  gsl_vector_free(xy_out);

  // return the result
  return xy_ret;
}


/*
 * Function: undistort_point
 *
 * Parameters:
 * @param coeffs   - the drizzle coefficients
 * @param pixmax   - the image dimensions
 * @param xy_image - the undistorted (x,y)
 *
 * Returns:
 * @return xy_new - the undistorted (x,y)
 */
d_point
undistort_point(gsl_matrix *coeffs, const px_point pixmax, d_point xy_image)
{
  d_point xy_new;

  xy_new = get_drz_position(xy_image, coeffs, pixmax);

  return xy_new;
}


/*
 * Function: get_crossdisp_matrix
 * The function extracts the drizzle coefficients stored in the
 * keyword headers of an image extension and creates a matrix
 * to store those coefficients. The matrix is returned.
 * Currently the drizzle coefficients are ALWAYS stored in the
 * primary header of the grism images. For this reason
 * all keywords are read from extension '1' (hardcoded).
 *
 * Parameters:
 * @param  filename   - the image filename
 * @param  sci_numext - the extension number
 *
 * Returns:
 * @return ret - the matrix with the drizzle coefficients
 */
gsl_matrix  *
get_crossdisp_matrix(char * filename, int sci_numext)
{
  int ncoeffs;
  int i;
  gsl_matrix * ret;
  //px_point    pixmax;
  char templt[FLEN_CARD];
  float acoeff, tmp;

  // get the number of coefficients
  tmp = get_float_from_keyword(filename, 1, "DRZCNUM");
  //  tmp = get_float_from_keyword(filename, sci_numext, "DRZCNUM");
  if (isnan(tmp)){
    // in case that the keyword does not exist,
    // find a nice way out of the door
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "Could not find drizzle keywords in: %s\n", filename);
    ret = gsl_matrix_alloc(1,1);
    gsl_matrix_set(ret, 0,0,GSL_NAN);
    return ret;
  }
  else {
    // cast the value to int
    ncoeffs = (int)tmp;
  }

  // allocate the matrix, set it to the default
  ret = gsl_matrix_alloc(2,ncoeffs);
  gsl_matrix_set_all (ret, 0.0);

  // go over the number of coefficients
  for (i=0; i<ncoeffs; i++){

    // set the name for the x-coefficient, get it and store it
    sprintf (templt, "DRZ%1iX%02i", sci_numext, i+1);
    //    acoeff = get_float_from_keyword(filename, sci_numext, templt);
    acoeff = get_float_from_keyword(filename, 1, templt);
    gsl_matrix_set(ret, 0,i,acoeff);

    // set the name for the y-coefficient, get it and store it
    sprintf (templt, "DRZ%1iY%02i", sci_numext, i+1);
    //    acoeff = get_float_from_keyword(filename, sci_numext, templt);
    acoeff = get_float_from_keyword(filename, 1, templt);
    gsl_matrix_set(ret, 1,i,acoeff);

  }

  // return the matrix
  return ret;
}


/*
 * Function: get_axis_scales
 * The function computes the distance scale along the major and
 * minor half axis of an object. The distance scale is subject
 * to changes due to the image distortion.
 *
 * Parameters:
 * @param actbeam   - the beam to get the drizzle scales for
 * @param drzcoeffs - the drizzle coefficients
 * @param pixmax    - the image dimensions
 *
 * Returns:
 * @return ret - the drizzle scales along major and minor half axis
 */
d_point
get_axis_scales(beam actbeam, gsl_matrix * drzcoeffs, px_point pixmax){

  d_point dispnt;
  d_point drz_ref;
  d_point drz_dis;

  d_point ret;

  // find a dummy point along the major half axis
  dispnt.x = actbeam.refpoint.x + cos(actbeam.aorient);
  dispnt.y = actbeam.refpoint.y + sin(actbeam.aorient);

  // compute the undistorted coos for the refpoint
  drz_ref = get_drz_position(actbeam.refpoint, drzcoeffs, pixmax);

  // compute the undistorted coos for the dummy point
  drz_dis = get_drz_position(dispnt, drzcoeffs, pixmax);

  // compute the scale along the major axis
  ret.x = sqrt((drz_ref.x-drz_dis.x)*(drz_ref.x-drz_dis.x) +
               (drz_ref.y-drz_dis.y)*(drz_ref.y-drz_dis.y));

  // find a dummy point along the minor half axis
  dispnt.x = actbeam.refpoint.x + cos(actbeam.aorient + M_PI/2.0);
  dispnt.y = actbeam.refpoint.y + sin(actbeam.aorient + M_PI/2.0);

  // compute the undistorted coos for the dummy point
  drz_dis = get_drz_position(dispnt, drzcoeffs, pixmax);

  // compute the scale along the minor axis
  ret.y = sqrt((drz_ref.x-drz_dis.x)*(drz_ref.x-drz_dis.x) +
               (drz_ref.y-drz_dis.y)*(drz_ref.y-drz_dis.y));

  // return the result
  return ret;
}


/*
 * Function: get_crossdisp_scale
 * The function computes the drizzle scale in cross-
 * dispersion direction for a given point and the
 * drizzle parameters of an image.
 *
 * Parameters:
 * @param trace     - the object trace
 * @param refpnt    - the reference point
 * @param drzcoeffs - the drizzle coefficients
 * @param pixmax    - the image dimensions
 *
 * Returns:
 * @return crscale - the drizzle scale in crossdispersion
 *                   direction
 */
double
get_crossdisp_scale(trace_func *trace, d_point refpnt, gsl_matrix * drzcoeffs,
                    px_point pixmax)
{
  double crscale;
  double phi_trace;
  //int beamID = 0;
  d_point dispnt;
  d_point drz_ref;
  d_point drz_dis;

  // get the trace angle
  phi_trace = atan2 (trace->deriv (0.0, trace->data), 1.0);

  // compute a dummy point along the cross-dispersion direction
  dispnt.x = refpnt.x + cos(phi_trace + M_PI/2.0);
  dispnt.y = refpnt.y + sin(phi_trace + M_PI/2.0);

  // get the undistorted positions of the refpoint and
  // the dummy point
  drz_ref = get_drz_position(refpnt, drzcoeffs, pixmax);
  drz_dis = get_drz_position(dispnt, drzcoeffs, pixmax);

  // compute the drizzle scale
  crscale = sqrt((drz_ref.x-drz_dis.x)*(drz_ref.x-drz_dis.x) +
                 (drz_ref.y-drz_dis.y)*(drz_ref.y-drz_dis.y));

  // return the drizzle scale
  return crscale;
}


/*
 * Function: evaln
 * The function evaluates the drizzle coefficients to return
 * either the x- or y-value in the undistorted system.
 * This is a translation of a fortran function with
 * the identical name inside the drizzle package.
 *
 * Parameters:
 * @param x         - the x-value in the distorted system
 * @param y         - the y-value in the distorted system
 * @param drzcoeffs - the drizzle coefficients
 * @param row       - which row of the coefficients to evaluate
 *
 * Returns:
 * @return ret - the undistorted coordinate value
 */
double
evaln(double x, double y, gsl_matrix * drzcoeffs, int row){

  double ret=0.0;
  int order;
  int n,m,nc,ncoeffs;
  order = 0;
  ncoeffs = 1;

  // derive the number of 'orders'
  while(ncoeffs+order+2 <= (int)drzcoeffs->size2){
    ncoeffs = ncoeffs+order+2;
    order++;
  }

  // add each individual term
  nc = 0;
  for (n=1; n < order+2; n++){
    for (m=1; m < n+1; m++){
      ret = ret + gsl_matrix_get(drzcoeffs,row,nc) * pow(x,(double)(n-m)) * pow(y,(double)(m-1));
      nc++;
    }
  }

  // give back the result
  return ret;
}


/*
 * Function: devalndx
 * The function numerically differentiates the
 * function 'evaln' with respect to the x-axis.
 *
 * Parameters:
 * @param x         - the x-value in the distorted system
 * @param y         - the y-value in the distorted system
 * @param drzcoeffs - the drizzle coefficients
 * @param row       - which row of the coefficients to evaluate
 *
 * Returns:
 * @return ret - the differential
 */
double
devalndx(double x, double y, gsl_matrix * drzcoeffs, int row){
  double ret=0.0;
  int order;
  int n,m,nc,ncoeffs;

  // derive the number of 'orders'
  order = 0;
  ncoeffs = 1;
  while(ncoeffs+order+2 <= (int)drzcoeffs->size2){
    ncoeffs = ncoeffs+order+2;
    order++;
  }

  // add each individual term
  nc = 0;
  for (n=1; n < order+2; n++){
    for (m=1; m < n; m++){
      ret = ret + gsl_matrix_get(drzcoeffs,row,nc) * (double)(n-m) * pow(x,(double)(n-m-1)) * pow(y,(double)(m-1));
      nc++;
    }
    nc++;
  }

  // return the result
  return ret;
}


/*
 * Function: devalndy
 * The function numerically differentiates the
 * function 'evaln' with respect to the y-axis.
 *
 * Parameters:
 * @param x         - the x-value in the distorted system
 * @param y         - the y-value in the distorted system
 * @param drzcoeffs - the drizzle coefficients
 * @param row       - which row of the coefficients to evaluate
 *
 * Returns:
 * @return ret - the differential
 */
double
devalndy(double x, double y, gsl_matrix * drzcoeffs, int row){
  double ret=0.0;
  int order;
  int n,m,nc,ncoeffs;

  // derive the number of 'orders'
  order = 0;
  ncoeffs = 1;
  while(ncoeffs+order+2 <= (int)drzcoeffs->size2){
    ncoeffs = ncoeffs+order+2;
    order++;
  }

  // add each individual term
  nc = 0;
  for (n=1; n < order+2; n++){
    nc++;
    for (m=2; m < n+1; m++){
      ret = ret + gsl_matrix_get(drzcoeffs,row,nc) * pow(x,(double)(n-m)) * (double)(m-1) * pow(y,(double)(m-2));
      nc++;
    }
  }

  // return the result
  return ret;
}


/*
 * Function: get_jacobian
 *
 * Parameters:
 * @param i         -
 * @param j         -
 * @param drzcoeffs - the drizzle coefficients
 * @param width     -
 * @param height    -
 *
 * Returns:
 * @return ret -
 */
gsl_matrix *
get_jacobian(int i, int j, gsl_matrix * drzcoeffs, int width, int height){

  double xr, yr;
  gsl_matrix * ret;

  ret = gsl_matrix_alloc(2,2);
  gsl_matrix_set_all (ret, 0.0);

  xr = (double)i - ((double)(width/2)+1.0);
  yr = (double)j - ((double)(height/2)+1.0);

  gsl_matrix_set(ret, 0, 0, devalndx( xr, yr, drzcoeffs, 0));
  gsl_matrix_set(ret, 0, 1, devalndy( xr, yr, drzcoeffs, 0));
  gsl_matrix_set(ret, 1, 0, devalndx( xr, yr, drzcoeffs, 1));
  gsl_matrix_set(ret, 1, 1, devalndy( xr, yr, drzcoeffs, 1));

  return ret;
}


/*
 * Function: get_det_jacobian
 *
 * Parameters:
 * @param i         -
 * @param j         -
 * @param drzcoeffs - the drizzle coefficients
 * @param width     -
 * @param height    -
 *
 * Returns:
 * @return ret -
 */
double
get_det_jacobian(int i, int j, gsl_matrix * drzcoeffs, int width, int height){

  gsl_matrix * jacob;
  double ret;

  jacob = get_jacobian(i,j,drzcoeffs,width,height);

  ret = gsl_matrix_get(jacob, 0, 0) * gsl_matrix_get(jacob, 1, 1)
       -gsl_matrix_get(jacob, 0, 1) * gsl_matrix_get(jacob, 1, 0);

  gsl_matrix_free(jacob);

  return ret;
}


/*
 * Function: get_drz_position
 * The function transforms distorted coordinates image to undistorted
 * image coordinates. As an assumption, the undistorted image
 * has the same dimension as the distorted.
 *
 * Parameters:
 * @param in_point  - the coordinates in the distorted frame
 * @param drzcoeffs - the drizzle coefficients
 * @param pixmax    - the image dimension of the distorted frame
 *
 * Returns:
 * @return: ou_point - the coordinates in the undistorted frame
 */
d_point
get_drz_position(d_point in_point, gsl_matrix *drzcoeffs, px_point pixmax)
{
  d_point ou_point;
  double xr, yr;
  double xn, yn;

  // transform the input coos into the drizzle system
  xr = in_point.x - ((double)(pixmax.x/2)+1.0);
  yr = in_point.y - ((double)(pixmax.y/2)+1.0);

  //  ou_point.x = evaln(xr,yr,drzcoeffs,0);
  //  ou_point.y = evaln(xr,yr,drzcoeffs,1);
  // get the undistorted coos in the drizzle system
  xn = evaln(xr,yr,drzcoeffs,0);
  yn = evaln(xr,yr,drzcoeffs,1);

  // transform back into the image system,
  // assuming identical image dimensions
  ou_point.x = xn + ((double)(pixmax.x/2)+1.0);
  ou_point.y = yn + ((double)(pixmax.y/2)+1.0);

  // return the undistorted coos
  return ou_point;
}


/*
 * Function: get_drz_position_free
 * The function transforms distorted coordinates image to undistorted
 * image coordinates. The size of the input and output image are free.
 *
 * Parameters:
 * @param in_point  - the coordinates in the distorted frame
 * @param drzcoeffs - the drizzle coefficients
 * @param pix_in    - the image dimension of the distorted frame
 * @param pix_out   - the image dimension of the distorted frame
 *
 * Returns:
 * @return: ou_point - the coordinates in the undistorted frame
 */
d_point
get_drz_position_free(d_point in_point, gsl_matrix *drzcoeffs,
                      px_point pix_in, px_point pix_out)
{
  d_point ou_point;
  double xr, yr;
  double xn, yn;

  // transform the input coos into the drizzle system
  xr = in_point.x - ((double)(pix_in.x/2)+1.0);
  yr = in_point.y - ((double)(pix_in.y/2)+1.0);

  //  ou_point.x = evaln(xr,yr,drzcoeffs,0);
  //  ou_point.y = evaln(xr,yr,drzcoeffs,1);
  // get the undistorted coos in the drizzle system
  xn = evaln(xr,yr,drzcoeffs,0);
  yn = evaln(xr,yr,drzcoeffs,1);

  // transform back into the image system,
  // assuming identical image dimensions
  ou_point.x = xn + ((double)(pix_out.x/2)+1.0);
  ou_point.y = yn + ((double)(pix_out.y/2)+1.0);

  // return the undistorted coos
  return ou_point;
}


trace_func *
get_tracefunc_at(char *filename, d_point p){

  int beamID = 0;
  trace_func  *trace;
  tracestruct *tstruct;

  tstruct = get_tracestruct_at_pos (filename, beamID, p);
  trace   = vector_to_trace_polyN ( tstruct->pol );

  return trace;
}
