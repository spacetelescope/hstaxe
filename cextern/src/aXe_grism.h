/*
 */


#ifndef _aXe_GRISM_H
#define _aXe_GRISM_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "spc_trace_functions.h"
//#include "aXe_errors.h"

#define RELEASE "See github"

#define gsl_matrix              gsl_matrix_float
#define gsl_matrix_get          gsl_matrix_float_get
#define gsl_matrix_set          gsl_matrix_float_set
#define gsl_matrix_set_all      gsl_matrix_float_set_all
#define gsl_matrix_alloc        gsl_matrix_float_alloc
#define gsl_matrix_free         gsl_matrix_float_free
#define gsl_matrix_fprintf      gsl_matrix_float_fprintf
#define gsl_matrix_add_constant gsl_matrix_float_add_constant
#define gsl_matrix_scale        gsl_matrix_float_scale
#define gsl_matrix_min          gsl_matrix_float_min
#define gsl_matrix_max          gsl_matrix_float_max

/* WARNING: If you ever change PIXEL_T, make sure you also change
   the FITS loader */
#define PIXEL_T float           /* type of pixel data -- this should be in
                                   sync with the GSL definitions */
#define MAX_BEAMS 27            /* Max. number of spectrum orders we handle */
                                // now equals the alphabet...

#define BEAM(x) (char)((int)x+65)       /* Defines names for 3 different beams */

#define MAXCHAR 255             /* Number of characters in small strings and FITS keywords */

/**
        A set of two relative (pixel) offsets.
*/
typedef struct
{
  int dx0, dx1;
}
px_offset;

/**
        A set of two relative double offsets.
*/
typedef struct
{
  double dx0, dx1;
}
dx_offset;

/**
  A point with integer (pixel) coordinates.
*/
typedef struct
{
  int x, y;
}
px_point;


/**
  A point with double coordinates.
*/
typedef struct
{
  double x, y;
}
d_point;


/**
  An aggregation of the grism exposure, its errors, the background,
  and flatfields
*/
typedef struct
{
  gsl_matrix *grism;            /* The grism image */
  gsl_matrix *pixerrs;  /* (Absolute) errors associated with the grism image */
  gsl_matrix *dq;        /* Data quality associated with grism image */
}
observation;


/**
   A beam, i.e., a description of the image of a spectrum of a given order,
   including a reference point, the corners of a convex quadrangle
   bounding the image, and the parametrization of the spectrum trace
*/
typedef struct
{
  int ID;                       /* An integer ID for this beam
                                   (see enum beamtype */
  d_point  refpoint;            /* Pixel coordinates of the reference point */
  px_point corners[4];          /* Bounding box of spectrum as an arbitrary
                                   quadrangle */
  px_point bbox[2];             /* Bounding box of corners, will be filled
                                   in by spc_extract */
  trace_func *spec_trace;       /* Parametrization of the spectrum trace */
  double width;                 /* largest axis of object in pixels */
  double orient;                /* orientation of largest axis in radian,
                                   0==horizontal */
  int modspec;
  int modimage;

  double awidth;                /* the object width along the extraction direction */
  double bwidth;                /* the object width perpendicular to ewidth
                                   (which is closer to the dispersion direction */
  double aorient;               /* orientation of awidth (CCW) =0 along the x-axis */

  double slitgeom[4];          /* containst the values for the slit length, angle and
                                   so on */

  gsl_vector *flux;              /* the flux of the object at various wavelengths */

  int ignore;           /* This beam should be ignored if this is not set to 0 */
}
beam;


/**
  An object, describing the width and orientation of the object,
  the beams in which the spectrum images lie, and the images
  comprising an observation of the object.
*/
typedef struct
{
     int ID;                    /* An integer ID number for this object */
     beam beams[MAX_BEAMS];     /* table of beams (orders of spectrum) */
     int nbeams;                /* number of beams in beams member */
     observation *grism_obs;
}
object;


/**
  An aperture pixel, i.e., a pixel within a beam.  The structure
  contains the coordinates, the abscissa of the section point between
  the object axis and the spectrum trace, the distance between
  the pixel and the section point along the object axis, the path
  length of the section point along the trace, and the photon count for
  this pixel.  A table of these is produced by the
  make_spc_table function.
*/
typedef struct
{
  int p_x, p_y;         /* the absolute coordinates of the
                                   source pixel. A -1 in p_x signifies
                                   the end of a list of ap_pixels */
  double x, y;          /* the logical coordinates of this (logical)
                                   pixel relative to the beam's reference
                                   point */
  double dist;          /* distance between point and section */
  double xs;                    /* abscissa of section point relative to reference point */
  double ys;                    /* y coord of section point relative to reference point */

  // traditionally (at least aXe-1.3-1.6 "dxs" was called "Projected size of pixel along the trace".
  // however in the code (spc_extract.c, function "handle_one_pixel()" it becomes clear
  // that this quantity indeed is the local trace angle
  double dxs;

  double xi;                    /* path length of the section along the
                                   spectrum trace (computed from xs) */
  double lambda;                /* path length of the section along the
                                   spectrum trace (computed from xi) */
  double dlambda;               /* dpath length of the section along the
                                   spectrum trace */
  double count;         /* intensity for this pixel */
  double weight;     /* extraction weight for this pixel */
  double error;         /* absolute error of count */
  double contam;     /* flag whether this pixel is contaminated by another spectrum */
  double model;     /* the pixel flux computed for one of the emission models */
  long dq;           /* DQ bit information */
}
ap_pixel;

typedef struct
{
  ap_pixel *ap;                 /* A pointer to an ap_pixel array */
  object *obj;                  /* A pointer to the corresponding object
                                   structure */
  beam *b;                      /* A pointer to the beam in the object
                                   that this group of pixel correspond to */
}
PET_entry;


/* See comment to spectrum structure */
typedef enum
{
     SPC_W_CONTAM,
     SPC_W_BORDER
}
spec_warning;


/**
  A point within a spectrum, consisting of the minimum, maximum, and
  (weighted) mean lambda of the points that contributed to the
  bin, as well as an error estimate for the value (that's in there,
  too, of course)
*/
typedef struct
{
  double lambda_min, lambda_mean;
  double lambda_max, dlambda;/* wavelength information */
  double error;         /* error in counts in this bin */
  double count;         /* number of counts in this bin */
  double weight;                /* weight of this bin */
  double flux;          /* flux in physical units in this bin */
  double ferror;                /* error in physical units in this bin */
  double contam;     /* spectral contamination flag */
  //  int    contam;     /* spectral contamination flag */
  long   dq;        /* DQ information */
}
spc_entry;


/**
  A spectrum, consisting of a gsl_array containing the data, a start
  value (center wave length of lowest bin), a bin size, and warn flags,
  which may be <ul><li>SPC_W_CONTAM -- Spectrum is contaminated</li>
  <li>SPC_W_BORDER -- Spectrum lies partially beyond the image borders.
  </li></ul>
*/
typedef struct
{
  spc_entry *spec;      /* the actual spectrum */
  double lambdamin;
  double lambdamax;
  int spec_len;         /* Number of items in spec */
  spec_warning warning; /* problems with this spectrum */
}
spectrum;

#define DEBUG_FILE 0x01
#define DEBUG 0x0000

extern void
print_ap_pixel_table (const ap_pixel * ap_p);

extern ap_pixel *
make_spc_table (object * const ob,
                const int beamorder, int *const flags);

extern ap_pixel *
make_gps_table (object * const ob,
                const int beamorder, int *const flags,
                int xval, int yval);

#endif /* !_aXe_GRISM_H */
