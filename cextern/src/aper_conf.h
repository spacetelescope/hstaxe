/*
 */
#ifndef _APER_CONF_H
#define _APER_CONF_H

#include "aXe_grism.h"
#include "spc_cfg.h"
#include "disp_conf.h"

#define BCK_BEAM_EXT 1


/**
	The basic description of a beam span
	and its numerical label
*/
typedef struct
{
  int ID;
  px_offset offset;
  float mmag_extract;
  float mmag_mark;
  float psf_offset;
}
beam_conf;

/**
	The basic description of the dispersed apertures
	in coordinates relative to the direct object position.
	Contains MAX_BEAMS beam_conf
	@see beam_conf
*/
typedef struct
{
  char instrument[MAXCHAR];  /* Name of the instrument */
  char camera[MAXCHAR];      /* Name of the camera*/
  int  nbeams;               /* Number of beams defined */
  char science_ext[MAXCHAR]; /* FITS extension name of the data */
  char errors_ext[MAXCHAR];  /* FITS extension name of the error */
  char dq_ext[MAXCHAR];      /* FITS extension name of the data quality */
  int  science_numext;       /* Extension number of the data */
  int  errors_numext;        /* Extension number of the error */
  int  dq_numext;            /* Extension number of the data quality */
  int  dqmask;               /* Data Quality Mask */
  char exptimekey[MAXCHAR];  /* Name of EXPTIME FITS keyword */
  char gainkey[MAXCHAR];     /* Name of the GAIN FITS keyword */
  char optkey1[MAXCHAR];     /* optional key */
  char optval1[MAXCHAR];     /* optional key value */
  char FFname[MAXCHAR];      /* The name of the flatfield file */
  char IPIXfunc[MAXCHAR];    /* Name of the intrapixel correction function */
  int  refx;	             /* Optional reference x pixel of the field dependence */
  int  refy;                 /* Optional reference y pixel of the field dependence */
  beam_conf beam[MAX_BEAMS];
  /* aXe-1.4 :*/
  float drz_resol;           /* The resolution in the drizzled images in [A/pixel]*/
  float drz_scale;           /* The scale in crossdispersion direction in "/pixel */
  float drz_lamb0;           /* The starting wavelength for the drizzled spectrum */
  float drz_xstart;          /* The x-position the starting wavelength is drizzled to */
  float drz_pfrac;           /* The pixfrac for the drizzling */
  char  drz_kernel[MAXCHAR]; /* The kernel for the drizzling */
  /* aXe-1.5 :*/
  gsl_vector *psfcoeffs;
  gsl_vector *psfrange;
  float rdnoise;

  double pobjsize;
  double smfactor;

}
aperture_conf;


extern aperture_conf *
get_aperture_descriptor (char *filename);

extern void
aperture_conf_fprintf (FILE * file, aperture_conf * conf);

extern void
get_extension_numbers(char filename[], aperture_conf * conf,
		      char keyword[],char keyval[]);

extern void
free_aperture_conf(aperture_conf * conf);

double
get_psf_offset(aperture_conf * conf, const beam actbeam);

double
get_max_offset(aperture_conf * conf);

extern void
extend_config_beams(aperture_conf *conf);

#endif
