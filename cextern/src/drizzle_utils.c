/**
* Drizzle utilities
*/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <math.h>

#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_trace_functions.h"
#include "crossdisp_utils.h"
#include "aper_conf.h"
#include "trace_conf.h"
#include "disp_conf.h"
#include "spc_wl_calib.h"
#include "drizzle_utils.h"


d_point
get_refwave_position(dispstruct * disp, trace_func * trace, d_point refpix,
                     aperture_conf *conf)
{
  double a=0.0;
  double b, c, c_star;
  double *cf;
  double rot;
  double dx, dy, dtr;
  d_point res;

  // check for a dispersion solution
  // higher that quadratic order
  if (disp->pol->size > 3)
    // this can not be solved; give an error
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "Order of dispersion solution: %i!\n"
                "At most qudratic solutions are allowed!\n",
                disp->pol->size-1);
  else if (disp->pol->size > 2)
    // store the quadratic term, if it exists
    a = gsl_vector_get(disp->pol, 2);  // ddlambda-term
  else
    // set the qudratic term to zero
    a=0.0;

  // store the constant and linear terms
  b = gsl_vector_get(disp->pol, 1);  // dlambda-term
  c = gsl_vector_get(disp->pol, 0);  // lambda0-term

  // store the dydx value
  cf = trace->data;

  // dydx -->  rotation angle
  rot = get_rotation_angle(cf[2]);

  // determine the wavelength difference
  c_star = c - conf->drz_lamb0;

  // transform the wavelength difference
  // into a path length difference
  if (a)
    dtr = (-1.0*b+sqrt(b*b - 4.0*a*c_star))/(2.0*a);
  else
    dtr = c_star/b;

  // path length difference --> dx, dy
  dx = dtr*cos(rot);
  dy = dtr*sin(rot);

  // compute the absolute values
  // for the reference position
  // in the stamp images
  res.x = refpix.x + dx;
  res.y = refpix.y + dy;

  // return the result
  return res;
}

double
get_drizzle_width(object *ob, int beamnum,trace_func * trace)
{
  double drizzle_width, drizzle_orient;
  double orig_width, orig_orient;
  double rotation, angle, oangle, factor;
  double *cf;
  beam *b = ob->beams+beamnum;

  cf = trace->data;
  rotation = get_rotation_angle(cf[2]);

  orig_orient = b->orient;
  orig_width  = b->width;

  drizzle_orient = orig_orient-rotation;

  angle  = drizzle_orient / M_PI * 180.;
  oangle = orig_orient / M_PI * 180.;

  drizzle_width = fabs(orig_width * sin (drizzle_orient));
  factor = drizzle_width/orig_width;

  if (angle > 180.0)
    angle = angle-180.0;
  if (angle < 30.0 || angle > 150.0)
     fprintf (stdout,
              "Angle: %4.0f --> %4.0f, Width: %5.1f --> %5.1f, %5.3f\n",
              oangle, angle, orig_width, drizzle_width, factor);

  return drizzle_width;
}

gsl_matrix  *
get_drizzle_coeffs(dispstruct * disp, trace_func * trace,
                   int boxwidth, int boxheight, int trlength,
                   double relx, double rely, aperture_conf *conf,
                   double orient, dispstruct * outdisp, double cdcorr,
                   double sprefreso, double spmeanreso)
{
  double a, b, c;
  double ao, bo, co, bref, cref;
  double a11, a12, a21, a22;
  double xr, yr;
  //double lambda0;
  //double dlambda;
  double tmp, rotation;
  double *cf;
  double shear_term;

  gsl_matrix * ret = gsl_matrix_alloc(2,11);
  gsl_matrix * rotcoeffs;

  gsl_matrix_set_all (ret, 0.0);


  /* get the dispersion at the objects point */
  if (disp->pol->size > 2)
    a = gsl_vector_get(disp->pol, 2);
    //    fprintf (stdout, "quadratic solution:  %f\n", a);
  else
    a = 0.0;

  b = gsl_vector_get(disp->pol, 1);  // dlambda-term
  c = gsl_vector_get(disp->pol, 0);  // lambda0-term

  /* get the dispersion at the mean reference point */
  if (outdisp->pol->size > 2){
    ao = gsl_vector_get(outdisp->pol, 2);
    //    fprintf (stdout, "quadratic solution:  %f, old: %f\n", a, ao);
  }
  else{
    ao = 0.0;
  }
  bo = gsl_vector_get(outdisp->pol, 1);  // dlambda-term
  co = gsl_vector_get(outdisp->pol, 0);  // lambda0-term
  if (conf->drz_resol == 0.0) {
    conf->drz_resol = bo + ao*((double)trlength)/2.0;
  }
  if (conf->drz_lamb0 == 0.0) {
    conf->drz_lamb0 = co;
  }

  //  fprintf (stdout, "lambda terms:  4785.0 <-> %f, 24.0 <-> %f\n", cref, bref);

  cref = conf->drz_lamb0;
  bref = conf->drz_resol;
  bref = sprefreso;

  // Bugfix on Sept. 15th 2010:
  // Check whether the ratio of the average dispersion on the grism images
  // and the dispersion on the drizzled images is outsidethe range 0.95 < ratio < 1.05.
  // If yes, adjust "trlength" which, from now on, is the length of the
  // drizzled image. This fix is necessary to allow a different sampling
  // in axedrizzle. The range was introduced such that the bug-fixed
  // version delivers identical results when using the default sampling.
  if ((spmeanreso / sprefreso > 1.05) || (spmeanreso / sprefreso < 0.95))
    trlength = (int)(spmeanreso / sprefreso * (double)trlength) + 0.5;

  /* put the below lines to go BACK to
   * the old computation of the length of
   * the aXedrizzled images:
  fprintf(stdout, "New spectral length: %i", trlength);
  trlength = (int)ceil(bo/bref*sqrt(pow((double)boxwidth,2.0)
                                    +pow((double)boxheight,2.0)));
  trlength = (int)ceil(spmeanreso/sprefreso*sqrt(pow((double)boxwidth,2.0)
                                                 +pow((double)boxheight,2.0)));
  fprintf(stdout, "<--> old spectral length: %i\n", trlength);
  */

  //fprintf (stdout,"Drizzled to resolution: %e, mean resolution for object: %e\n", sprefreso, spmeanreso);


  gsl_matrix_set(ret, 0,10,trlength);
  //  gsl_matrix_set(ret, 1,10,bref);

  cf = trace->data;
  rotation = get_rotation_angle(cf[2]);
  rotcoeffs = get_coeffs_back(trace);
  a11 = gsl_matrix_get(rotcoeffs,0,0);
  a12 = gsl_matrix_get(rotcoeffs,0,1);
  a21 = gsl_matrix_get(rotcoeffs,1,0);
  a22 = gsl_matrix_get(rotcoeffs,1,1);

  // the following lines have to be changed
  // in order to go from a integer
  // center definition to floating point center
  // definition.
  xr = (relx+1.0)-((double)(boxwidth/2)+1);  // +1.0 to compensate for "sp_sex.c" around line 537
  yr = (rely+1.0)-((double)(boxheight/2)+1); // +1.0 to compensate for "sp_sex.c" around line 537

  if (b < 0.0)
    shear_term = -tan(orient-rotation);
  else
    shear_term = tan(orient-rotation);

  // Transformations for X, without the lambda-terms, only rotation:
  //  tmp = xoffs-(((double)trlength/2.0)+1)-a11*xr-a12*yr;
  //  gsl_matrix_set(ret, 0,0,tmp);  // constant term

  //  tmp = a11;
  //  gsl_matrix_set(ret, 0,1,tmp);  // x-term

  //  tmp = a12;
  //  gsl_matrix_set(ret, 0,2,tmp);  // y-term
  //----------------------------------------------------------------------------------




  //------------------------------------------------------------------
  // Mathematica based run
  tmp = 1.0/bref*(c - cref - a11*b*xr + a*a11*a11*xr*xr - a12*b*yr + 2*a*a11*a12*xr*yr + a*a12*a12*yr*yr) +
    (a21*xr + a22*yr)/shear_term - ((double)(trlength/2)+1.0) + conf->drz_xstart;
  gsl_matrix_set(ret, 0,0,tmp);  // constant term

  tmp = 1.0/bref*(a11*b - 2*a*a11*a11*xr - 2*a*a11*a12*yr) - a21/shear_term;
  gsl_matrix_set(ret, 0,1,tmp);  // x-term

  tmp = 1.0/bref*(a12*b - 2*a*a12*a12*yr - 2*a*a11*a12*xr) - a22/shear_term;
  gsl_matrix_set(ret, 0,2,tmp);  // y-term

  tmp = a*a11*a11/bref;
  gsl_matrix_set(ret, 0,3,tmp);  // x^2-term

  tmp = 2*a*a11*a12/bref;
  gsl_matrix_set(ret, 0,4,tmp); ; // xy-term

  tmp = a*a12*a12/bref;
  gsl_matrix_set(ret, 0,5,tmp);  // y^2-term

  //------------------------------------------------------------------


  // Transformations for Y:
  tmp = -1.0*(a21*xr+a22*yr) * cdcorr; // -1.5  is some kind of a fudge factor, no idea where it comes from
  gsl_matrix_set(ret, 1,0,tmp);  // constant term

  tmp = a21 * cdcorr;
  gsl_matrix_set(ret, 1,1,tmp);  // x-term

  tmp = a22 * cdcorr;
  gsl_matrix_set(ret, 1,2,tmp);  // y-term

  gsl_matrix_free(rotcoeffs);

  return ret;
}

gsl_matrix *
get_coeffs_back(trace_func * trace)
{
  double *cf;
  double rotation;
  gsl_matrix * ret = gsl_matrix_alloc(2,2);

  //  ret = gsl_matrix_alloc (2,2);
  gsl_matrix_set_all (ret, 0.0);

  cf = trace->data;

  rotation = get_rotation_angle(cf[2]);

  gsl_matrix_set(ret, 0, 0,      cos(rotation));
  gsl_matrix_set(ret, 0, 1,      sin(rotation));
  gsl_matrix_set(ret, 1, 0, -1.0*sin(rotation));
  gsl_matrix_set(ret, 1, 1,      cos(rotation));

  return ret;
}

double
get_rotation_angle(double dxdy){
  double rotation=0.0;

  rotation = atan2(dxdy, 1.0);

  return rotation;
}

// Start functions for an alternative approach to store the object information
objectobs **
malloc_objectobs(){
  objectobs **allobjects;

  allobjects = (objectobs **) malloc (NMAXOBJ*sizeof(objectobs *));

  return allobjects;
}

void
free_objectobs(objectobs **allobjects){
  //  int i, nobjs=0;
  //  while (allobjects[nobjs])
  //    nobjs++;
  //  for (i=0; i<nobjs; i++){
  //      free(allobjects[i]);
  //      //    allobjects[i] = NULL;
  //  }
  free(allobjects);
}

int
add_observation(char * filename, char *  conf_file, objectobs ** allobjects, int nobjects,
                object ** oblist, int list_size, px_point pixmax, int sci_numext)
{

  int i, j, dec;
  int beamID = 0;
  int get_scale=0;
  int tlength;
  gsl_matrix * drzcoeffs;
  double cdscale;
  double drzscale;
  dispstruct * disp;
  calib_function *wl_calib;
  double l1, l2, spreso;

  drzcoeffs = get_crossdisp_matrix(filename, sci_numext);

  if (drzcoeffs->size1 > 1 && drzcoeffs->size2)
    get_scale = 1;

  if (get_scale){
    drzscale =  (double)get_float_from_keyword(filename, 1, "DRZSCALE");
  }
  for (i=0; i < list_size; i++){
    if (oblist[i]->beams[beamID].ignore == 0){
      dec = 0;

      if (get_scale)
        cdscale = drzscale * get_crossdisp_scale(oblist[i]->beams[beamID].spec_trace,
                                                 oblist[i]->beams[beamID].refpoint, drzcoeffs, pixmax);
      else
        cdscale=1.0;

      disp = get_dispstruct_at_pos(conf_file, 1, beamID,oblist[i]->beams[beamID].refpoint);
      wl_calib = create_calib_from_gsl_vector(1, disp->pol);
      l1 =  wl_calib->func (0.0,  wl_calib->order,  wl_calib->coeffs);
      l2 =  wl_calib->func (1.0,  wl_calib->order,  wl_calib->coeffs);
      spreso = fabs(l2-l1);

      // compute the tracelength
      tlength = get_beam_trace_length(oblist[i]->beams[beamID]);

      for (j=0; j < nobjects; j++){
        if (allobjects[j]->OBJID == oblist[i]->ID){
          add_obs_to_allobj(allobjects[j], oblist[i], pixmax, cdscale, spreso, tlength);
          dec = 1;
        }
      }
      if (dec == 0) {
        nobjects = add_obj_to_allobj(allobjects, nobjects, oblist[i], pixmax, cdscale, spreso, tlength);
      }
      free_dispstruct(disp);
      free_calib(wl_calib);
    }
  }
  gsl_matrix_free(drzcoeffs);
  return nobjects;
}

void
add_obs_to_allobj(objectobs *actobject, object * actobs, px_point pixmax, double cdscale, double spreso, int tlength)
{

  int xmin, xmax, ymin, ymax;
  //double m, b;
  gsl_vector_int * xvec;
  gsl_vector_int * yvec;
  px_point bbox;
  px_point mins;

  double *gaga;

  xvec = gsl_vector_int_alloc (4);
  yvec = gsl_vector_int_alloc (4);

  /* Store the refpoint */
  actobject->refpoint[actobject->nobs].x = actobs->beams[0].refpoint.x;
  actobject->refpoint[actobject->nobs].y = actobs->beams[0].refpoint.y;

  //*************************************************
  // patch to correct the reference point in case
  // that the the trace descritpion
  // does have a non negligeable first order term!
  gaga = actobs->beams[0].spec_trace->data;
  actobject->refpoint[actobject->nobs].y = actobject->refpoint[actobject->nobs].y + gaga[1];
  //*************************************************

  gsl_vector_int_set(xvec, 0, actobs->beams[0].corners[0].x);
  gsl_vector_int_set(xvec, 1, actobs->beams[0].corners[1].x);
  gsl_vector_int_set(xvec, 2, actobs->beams[0].corners[2].x);
  gsl_vector_int_set(xvec, 3, actobs->beams[0].corners[3].x);
  xmin = gsl_stats_int_min(xvec->data, 1, 4);
  xmax = gsl_stats_int_max(xvec->data, 1, 4);
  gsl_vector_int_set(yvec, 0, actobs->beams[0].corners[0].y);
  gsl_vector_int_set(yvec, 1, actobs->beams[0].corners[1].y);
  gsl_vector_int_set(yvec, 2, actobs->beams[0].corners[2].y);
  gsl_vector_int_set(yvec, 3, actobs->beams[0].corners[3].y);
  ymin = gsl_stats_int_min(yvec->data, 1, 4);
  ymax = gsl_stats_int_max(yvec->data, 1, 4);

  bbox = recalc_bbox(xmin, xmax, ymin, ymax, pixmax);
  mins = recalc_mins(xmin, xmax, ymin, ymax, pixmax);

  actobject->width[actobject->nobs]    = bbox.x;
  actobject->height[actobject->nobs]   = bbox.y;
  actobject->tlength[actobject->nobs]  = tlength;
  actobject->objwidth[actobject->nobs] = actobs->beams[0].width;
  actobject->orient[actobject->nobs]   = actobs->beams[0].orient;
  actobject->cdscale[actobject->nobs]   = cdscale;
  actobject->spreso[actobject->nobs]   = spreso;

  actobject->relrefpt[actobject->nobs].x = actobs->beams[0].refpoint.x - (double)mins.x;
  actobject->relrefpt[actobject->nobs].y = actobs->beams[0].refpoint.y - (double)mins.y;

  actobject->nobs = actobject->nobs + 1;

  gsl_vector_int_free(xvec);
  gsl_vector_int_free(yvec);
}

int
add_obj_to_allobj(objectobs **allobjects, int nobjects, object * actobs,
                   px_point pixmax, double cdscale, double spreso, int tlength)
{
  int xmin, xmax, ymin, ymax;
  //double m, b;

  gsl_vector_int * xvec;
  gsl_vector_int * yvec;

  px_point bbox;
  px_point mins;

  double *gaga;

  //  object *ob = malloc (sizeof (object));
  objectobs *objobs = malloc (sizeof (objectobs));

  xvec = gsl_vector_int_alloc (4);
  yvec = gsl_vector_int_alloc (4);

  allobjects[nobjects] = objobs;
  /* Store the refpoint */
  allobjects[nobjects]->refpoint[0].x = actobs->beams[0].refpoint.x;
  allobjects[nobjects]->refpoint[0].y = actobs->beams[0].refpoint.y;

  //*************************************************
  // patch to correct the reference point in case
  // that the the trace descritpion
  // does have a non negligeable first order term!
  gaga = actobs->beams[0].spec_trace->data;
  allobjects[nobjects]->refpoint[0].y = allobjects[nobjects]->refpoint[0].y + gaga[1];
  //*************************************************

  allobjects[nobjects]->OBJID = actobs->ID;
  allobjects[nobjects]->nobs = 1;
  allobjects[nobjects]->pointer = 0;


  gsl_vector_int_set(xvec, 0, actobs->beams[0].corners[0].x);
  gsl_vector_int_set(xvec, 1, actobs->beams[0].corners[1].x);
  gsl_vector_int_set(xvec, 2, actobs->beams[0].corners[2].x);
  gsl_vector_int_set(xvec, 3, actobs->beams[0].corners[3].x);
  xmin = gsl_stats_int_min(xvec->data, 1, 4);
  xmax = gsl_stats_int_max(xvec->data, 1, 4);

  gsl_vector_int_set(yvec, 0, actobs->beams[0].corners[0].y);
  gsl_vector_int_set(yvec, 1, actobs->beams[0].corners[1].y);
  gsl_vector_int_set(yvec, 2, actobs->beams[0].corners[2].y);
  gsl_vector_int_set(yvec, 3, actobs->beams[0].corners[3].y);
  ymin = gsl_stats_int_min(yvec->data, 1, 4);
  ymax = gsl_stats_int_max(yvec->data, 1, 4);

  bbox = recalc_bbox(xmin, xmax, ymin, ymax, pixmax);
  mins = recalc_mins(xmin, xmax, ymin, ymax, pixmax);

  allobjects[nobjects]->width[0]   = bbox.x;
  allobjects[nobjects]->height[0]  = bbox.y;
  allobjects[nobjects]->tlength[0] = tlength;
  allobjects[nobjects]->objwidth[0] = actobs->beams[0].width;
  allobjects[nobjects]->orient[0]   = actobs->beams[0].orient;
  allobjects[nobjects]->cdscale[0]  = cdscale;
  allobjects[nobjects]->spreso[0]   = spreso;

  allobjects[nobjects]->relrefpt[0].x = actobs->beams[0].refpoint.x - (double)mins.x;
  allobjects[nobjects]->relrefpt[0].y = actobs->beams[0].refpoint.y - (double)mins.y;

  ++nobjects;

  gsl_vector_int_free(xvec);
  gsl_vector_int_free(yvec);

  return nobjects;
}

void print_objectobs(objectobs **allobjects, int nobjects){
  int i, j;

  for (i = 0; i < nobjects; i++){
   fprintf(stdout, "OBJID: %i, no of observ.: %i\n", allobjects[i]->OBJID, allobjects[i]->nobs);
    for (j = 0; j < allobjects[i]->nobs; j++){
      fprintf(stdout, "No %i: width %f, height: %f, \n", j, allobjects[i]->relrefpt[j].x, allobjects[i]->relrefpt[j].y);
    }
  }
}
void print_objectobs2(objectobs allobjects[], int nobjects){
  int i, j;

  for (i = 0; i < nobjects; i++){
   fprintf(stdout, "OBJID: %i, no of observ.: %i\n", allobjects[i].OBJID, allobjects[i].nobs);
    for (j = 0; j < allobjects[i].nobs; j++){
      fprintf(stdout, "No %i: width %f, height: %f, \n", j, allobjects[i].relrefpt[j].x, allobjects[i].relrefpt[j].y);
    }
  }
}
int make_refpoints( char * conf_file, char * filename, px_point pixmax, objectobs **allobjects, int nobjects)
{
  int i, j;
  int nobs;
  int max_width, max_height, omax;
  double xmean, ymean;
  //double xdata[NMAXOBS];
  //double ydata[NMAXOBS];
  //  objectobs  oneobject;
  objectobs *oneobject; // = malloc (sizeof (objectobs));

  gsl_vector * xvec;
  gsl_vector * yvec;
  aperture_conf *conf;
  gsl_matrix * drzcoeffs;
  //double cdscale;
  double drzscale;
  trace_func  *trace;
  int beamID=0;
  dispstruct * disp;
  calib_function *wl_calib;
  double l1, l2;

  conf = get_aperture_descriptor(conf_file);
  get_extension_numbers(filename, conf,conf->optkey1,conf->optval1);
  drzcoeffs = get_crossdisp_matrix(filename,conf->science_numext);
  drzscale =  (double)get_float_from_keyword(filename, 1, "DRZSCALE");

  for (i=0; i < nobjects; i++){
    oneobject = allobjects[i];
    nobs = oneobject->nobs;

    xvec = gsl_vector_alloc (nobs);
    yvec = gsl_vector_alloc (nobs);
    for (j=0; j < nobs; j++){
      gsl_vector_set(xvec, j, oneobject->refpoint[j].x);
      gsl_vector_set(yvec, j, oneobject->refpoint[j].y);
    }
    xmean = gsl_stats_mean(xvec->data, 1, nobs);
    ymean = gsl_stats_mean(yvec->data, 1, nobs);

    max_width  = gsl_stats_int_max(oneobject->width,  1, nobs);
    max_height = gsl_stats_int_max(oneobject->height, 1, nobs);
    omax       = gsl_stats_max(oneobject->objwidth, 1, nobs);

    allobjects[i]->mean_refpoint.x = xmean;
    allobjects[i]->mean_refpoint.y = ymean;
    allobjects[i]->max_width  = max_width;
    allobjects[i]->max_height = max_height;
    allobjects[i]->owidthmax  = omax;

    allobjects[i]->max_tlength = gsl_stats_int_max(oneobject->tlength, 1, nobs);

    /*
     * Look whether the cross dispersion scale is given in the configuration.
     * If not determine the cross dispersion scale at the mean reference point an store it
     * store also the mean correction in for to correct the width.
     *
     */
    if (conf->drz_scale < 1.0e-16){
      trace = get_tracefunc_at(conf_file, allobjects[i]->mean_refpoint);
      allobjects[i]->cdrefscale = drzscale * get_crossdisp_scale(trace, allobjects[i]->mean_refpoint, drzcoeffs, pixmax);
    }
    else{
      allobjects[i]->cdrefscale = conf->drz_scale;
    }
    for (j=0; j < nobs; j++){
      gsl_vector_set(xvec, j, oneobject->cdscale[j]);
    }
    xmean = gsl_stats_mean(xvec->data, 1, nobs);
    allobjects[i]->cdmeanscale  = xmean;

    /*
     * look whether the wavelength dispersion is given in the conig file
     * if not determine the wavelength dispersion at the mean reference point amd store
     * it. Store also the mean wavelength dispersion to correct the length
     */
    if (conf->drz_resol < 1.0e-16){
      disp = get_dispstruct_at_pos(conf_file, 1, beamID, allobjects[i]->mean_refpoint);
      wl_calib = create_calib_from_gsl_vector(1, disp->pol);
      l1 =  wl_calib->func (0.0,  wl_calib->order,  wl_calib->coeffs);
      l2 =  wl_calib->func (1.0,  wl_calib->order,  wl_calib->coeffs);
      allobjects[i]->sprefreso = fabs(l2-l1);
    }
    else{
       allobjects[i]->sprefreso = conf->drz_resol;
    }
    for (j=0; j < nobs; j++){
      gsl_vector_set(xvec, j, oneobject->spreso[j]);
    }
    xmean = gsl_stats_mean(xvec->data, 1, nobs);
    allobjects[i]->spmeanreso  = xmean;

    gsl_vector_free(xvec);
    gsl_vector_free(yvec);
  }

  gsl_matrix_free(drzcoeffs);
  free(conf);

  return 0;
}

d_point
get_mean_refpoint(objectobs **allobjects, int nobjects, int ID,int * boxwidth,
                  int * boxheight, double * relx, double * rely,
                  double * objwidth, double *orient, double * cdref,
                  double *cdscale, double *cdmeanscale, double *sprefreso, double *spreso,
                  double *spmeanreso, int *tlength)
{
  int i, dec=0;
  d_point mpoint;

  for (i=0; i < nobjects; i++){
    if (allobjects[i]->OBJID == ID){
      mpoint.x = allobjects[i]->mean_refpoint.x;
      mpoint.y = allobjects[i]->mean_refpoint.y;
      *boxwidth  = allobjects[i]->max_width;
      *boxheight = allobjects[i]->max_height;
      *tlength   = allobjects[i]->max_tlength;
      *objwidth  = allobjects[i]->owidthmax;
      *cdref     = allobjects[i]->cdrefscale;
      *cdmeanscale     = allobjects[i]->cdmeanscale;
      *spmeanreso= allobjects[i]->spmeanreso;
      *sprefreso = allobjects[i]->sprefreso;
      *orient    = allobjects[i]->orient[allobjects[i]->pointer];
      *cdscale   = allobjects[i]->cdscale[allobjects[i]->pointer];
      *spreso    = allobjects[i]->spreso[allobjects[i]->pointer];
      *relx      = allobjects[i]->relrefpt[allobjects[i]->pointer].x;
      *rely      = allobjects[i]->relrefpt[allobjects[i]->pointer].y;
      ++allobjects[i]->pointer;
      dec=1;
    }
  }
  return mpoint;
}
px_point
recalc_bbox(int xmin, int xmax, int ymin, int ymax, px_point pixmax)
{
  px_point ret;

  if (xmax > pixmax.x)
    xmax = pixmax.x;
  if (ymax > pixmax.y)
    ymax = pixmax.y;
  if (xmin < 1)
    xmin = 1;
  if (ymin < 1)
    ymin = 1;


  ret.x = xmax - xmin + 2;
  ret.y = ymax - ymin + 2;

  return ret;
}
/*
I am not sure what this routine is doing,
but one should make it better.
*/
px_point
recalc_mins(int xmin, int xmax, int ymin, int ymax, px_point pixmax)
{
  double m, b;
  px_point ret;

  m = ((double)ymin - (double)ymax)/((double)xmax - (double)xmin);
  b = (double)ymin - m * (double)xmax;
  if (xmin < 1){
    ymax = (int)ceil(1.0*m + b);
    xmin=1;
  }
  if (xmax > pixmax.x){
    ymin = (int)ceil((float)pixmax.x*m + b);
    xmax = pixmax.x;
  }
  if (ymin < 1){
    xmax = (int)ceil((1-b)/m);
    ymin=1;
  }
  if (ymax > pixmax.y){
    xmin = (int)ceil(((float)pixmax.y-b)/m);
    ymax=pixmax.y;
  }
  ret.x = xmin;
  ret.y = ymin;

  return ret;
}

int
get_beam_trace_length(const beam actbeam)
{
  int trlength=0;
  double diagonal_1 = 0.0;
  double diagonal_2 = 0.0;
  double xdiff, ydiff;

  // compute the x- and y- differences for one corner pair
  xdiff = (float)actbeam.corners[0].x - (float)actbeam.corners[2].x;
  ydiff = (float)actbeam.corners[0].y - (float)actbeam.corners[2].y;

  // compute the diagonal distance
  diagonal_1 = xdiff*xdiff + ydiff*ydiff;

  // compute the x- and y- differences for the other corner pair
  xdiff = (float)actbeam.corners[1].x - (float)actbeam.corners[3].x;
  ydiff = (float)actbeam.corners[1].y - (float)actbeam.corners[3].y;

  // compute the diagonal distance
  diagonal_2 = xdiff*xdiff + ydiff*ydiff;

  // compute the tracelength
  // from the larger diagonal
  if (diagonal_1 > diagonal_2)
    trlength = (int)ceil(sqrt(diagonal_1)+3.0);
  else
    trlength = (int)ceil(sqrt(diagonal_2)+3.0);

  // return the result
  return trlength;
}
