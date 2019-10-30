/**
 * Header files for the routines in spc_sex.c
 */
#ifndef SPC_DEF_H
#define SPC_DEF_H

#include "spc_CD.h"

#define CATBUFFERSIZE 10240
#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN_DIFFANGLE 3.0

/**
 * Structure: ellipse
 *  A description of an ellipse with a,b, theta.
 */
typedef struct
{
  double a;
  double b;
  double theta;
}
ellipse;



/**
 * Structure: SexObject
 *  A SeXtractor 2.0 object type which contains
 *  various information about the object
 */
typedef struct
{
  int number;                   /* Sextractor ID of the object */
  d_point xy_image;             /* i,j coordinates of the barycenter of the object in image X_IMAGE, Y_IMAGE */
  sky_coord xy_world;           /* Celestial coordinates of the barycenter of the object X_WORLD, Y_WORLD */
  ellipse el_image;             /* Ellipsoidal description of object in image A_IMAGE, B_IMAGE, THETA_IMAGE */
  ellipse el_world;             /* Ellipsoidal description of object in world coordinates A_WORLD */
                                /* B_WORLD, THETA_WORLD */
  double mag_auto;              /* Object magnitude MAG_AUTO */
  d_point backwindow;
  int modspec;
  int modimage;
  gsl_vector * lambdas;
  gsl_vector * magnitudes;
}
SexObject;


extern SexObject *
create_SexObject(const colinfo *actinfo, char *line, const gsl_vector * waves,
                  const gsl_vector *  cnums, const px_point backwin_cols, const
                 px_point modinfo_cols, const int magcol, const int fullinfo);

extern void
SexObject_fprintf (FILE * output, SexObject * o);

extern void
SexMags_to_beamFlux(SexObject * sobj, beam *actbeam);

extern void
SexObject_to_slitgeom(const aperture_conf *conf, const SexObject * sobj,
                      const double trace_angle, beam *actbeam);

extern void
fill_corner_ignore(SexObject * sobj, observation * const obs, aperture_conf *conf,
                   int beamID, float dmag, beam *actbeam);

extern d_point
check_object_size(const aperture_conf *conf, const SexObject *sobj, const int beamID);

extern void
SexObject_to_beam(SexObject * sobj,  observation * const obs, aperture_conf *conf, char conffile[],
                  float mfwhm, float dmag, int auto_reorient, int bck_mode,
                  int beamID, beam *actbeam);

extern object *
SexObject_to_objectII(SexObject * sobj, observation * const obs,
                      aperture_conf *conf, char conffile[], float mfwhm,
                      float dmag, int auto_reorient, int bck_mode);

extern object *
SexObject_to_object (SexObject * sobj, observation * const obs,
                     aperture_conf *conf, char conffile[], float mfwhm,
                     float dmag, int auto_reorient, int bck_mode);

extern object **
SexObjects_to_oblist (SexObject ** sobjs, observation * const obs,
                      aperture_conf *conf, char conffile[], float mfwhm,
                      float dmag, int auto_reorient, int back_mode);

extern object **
SexObjects_to_oblistII (SexObject ** sobjs, observation * const obs,
                        aperture_conf *conf, char conffile[], float mfwhm,
                        float dmag, int auto_reorient, int bck_mode);

extern int
check_conf_for_slitlessgeom(const aperture_conf *conf, const int auto_reorient);
//check_conf_for_slitlessgeom(const aperture_conf *conf, const int slitless_geom);

extern int
fill_object_bbox (observation * const obs, beam * b, const float m_width,
                  const int dxmin, const int dxmax);

extern int
size_of_sextractor_catalog (char filename[]);

extern SexObject **
get_SexObject_from_catalog (char filename[], const double lambda_mark);

extern void
ABC_image_to_el (d_point *a, d_point *b, d_point *c, ellipse *el);

extern void
ABC_image_to_el2 (d_point *a, d_point *b, d_point *c, ellipse *el);

extern void
el_to_ABC_world (ellipse *el, sky_coord *a, sky_coord *b, sky_coord *c);

extern void
el_to_ABC_world2 (ellipse *el, sky_coord *a, sky_coord *b, sky_coord *c);

extern void
fill_missing_image_coordinates (SexObject * o,
                                struct WorldCoor *from_wcs,int overwrite);

extern void
fill_all_missing_image_coordinates (SexObject ** sobjs,
                                    struct WorldCoor *from_wcs,
                                    int overwrite);

extern void
fill_missing_WCS_coordinates (SexObject * o, struct WorldCoor *from_wcs,
                              int overwrite);

extern void
fill_all_missing_WCS_coordinates (SexObject ** sobjs,
                                  struct WorldCoor *from_wcs,
                                  int overwrite);

extern void
compute_new_image_coordinates (SexObject * o, struct WorldCoor *to_wcs);

extern void
compute_new_image_sexobject (SexObject * o, struct WorldCoor *to_wcs, int th_sky);

extern void
compute_all_new_image_coordinates (SexObject ** sobjs,
                                   struct WorldCoor *to_wcs);

extern void
catalog_to_wcs (char grismfile[], int hdunum, char infile[], char outfile[],
                struct WorldCoor *from_wcs, struct WorldCoor *to_wcs,
                int distortion, int overwrite_wcs, int overwrite_img);

extern void
catalog_to_wcs_nodim (char infile[], char outfile[],
                      struct WorldCoor *to_wcs,int overwrite_wcs,
                      int overwrite_img);

extern gsl_vector *
sobs_to_vout(const SexObject *sobs);

extern void
free_SexObjects (SexObject ** sobjs);
#endif
