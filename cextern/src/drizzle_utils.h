/**
 *
 * File: drizzle_utils.h
 * Subroutines to calculate the
 * drizzle coefficients.
 *
 */

#ifndef _DRIZZLE_UTILS_H
#define _DRIZZLE_UTILS_H


#define NMAXOBS 1100
#define NMAXOBJ 20000

#define DRZMAX 0

typedef struct
{
  int OBJID;
  int nobs;
  int max_width;
  int max_height;
  int pointer;
  int max_tlength;
  double owidthmax;
  double cdrefscale;
  double cdmeanscale;
  double sprefreso;
  double spmeanreso;
  d_point mean_refpoint;
  int width[NMAXOBS];
  int height[NMAXOBS];
  int tlength[NMAXOBS];
  double objwidth[NMAXOBS];
  double orient[NMAXOBS];
  double cdscale[NMAXOBS];
  double spreso[NMAXOBS];
  d_point refpoint[NMAXOBS];    /* Pixel coordinates of the reference point */
  d_point relrefpt[NMAXOBS];
}
objectobs;

gsl_matrix *get_drizzle_coeffs(dispstruct * disp, trace_func * trace,
                               int boxwidth, int boxheight, int trlength,
                               double relx, double rely, aperture_conf *conf,
                               double orient, dispstruct * outdisp,
                               double cdcorr, double spref,double spmeanreso);
gsl_matrix *get_coeffs_back(trace_func * trace);
double get_rotation_angle(double dxdy);
int add_observation(char * filename, char * conf_file,
                    objectobs **allobjects, int nobjects, object ** oblist,
                    int list_size, px_point pixmax, int sci_numext);
void add_obs_to_allobj(objectobs * actobject, object * actobs, px_point pixmax,
                       double cdscale, double spreso, int tlength);
int add_obj_to_allobj(objectobs **allobjects, int nobjects, object * actobs,
                      px_point pixmax, double cdscale, double spreso, int tlength);
int make_refpoints(char * conf_file, char * filename, px_point pixmax,
                   objectobs **allobjects, int nobjects);
d_point get_mean_refpoint(objectobs **allobjects, int nobjects, int ID,
                          int * boxwidth, int * boxheight, double * relx,
                          double * rely, double * objwidth, double * orient,
                          double * cdref, double * cdscale,
                          double * cdmeanscale, double *sprefreso,
                          double *spreso, double *spmeanreso, int *tlength);
px_point recalc_bbox(int xmin, int xmax, int ymin, int ymax, px_point pixmax);
px_point recalc_mins(int xmin, int xmax, int ymin, int ymax, px_point pixmax);
d_point get_refwave_position(dispstruct * disp, trace_func * trace,
                             d_point refpix, aperture_conf *conf);
double get_drizzle_width(object *ob, int beamnum,trace_func * trace);
objectobs **malloc_objectobs(void);
void free_objectobs(objectobs **allobjects);
void print_objectobs(objectobs **allobjects, int nobjects);
void print_objectobs2(objectobs allobjects[], int nobjects);

extern int
get_beam_trace_length(const beam actbeam);
#endif
