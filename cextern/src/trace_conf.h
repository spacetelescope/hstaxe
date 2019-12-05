/**
 * File: trace_conf.h
 */

#ifndef _TRACE_CONF_H
#define _TRACE_CONF_H

#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "aXe_grism.h"
#include "aXe_errors.h"
#include "spc_cfg.h"
#include "disp_conf.h"


/**
	A structure to contain the coefficient of a 2D
	polynomial that can be used to compute the given
	coefficient of the polynomial dispersion relation
	as at a particular i.j location in the image
*/
typedef struct
{
     gsl_vector *pol;		/* A vector containing the 2D polynomial coefficients */
     d_point offset;		/* X and Y offsets to apply to input coordinates */
     d_point cpoint;		/* The detector location at which this structure was computed */
     int ID;			/* The ID of the beam for which this coefficient is defined */
     char file[MAXCHAR];
}
tracestruct;


extern int
get_beam_trace_norder (char *filename, int beamID);

extern gsl_vector *
get_beam_trace_order (char *filename, int beamID, int order);

extern float
get_trace_coeff_at_pos (char *filename, int beamID, int order, d_point p);

extern gsl_vector *
get_trace_coeffs_at_pos (char *filename, int beamID, d_point p);

extern gsl_vector *
get_beam_trace_xoff (char *filename, int beamID);

extern gsl_vector *
get_beam_trace_yoff (char *filename, int beamID);

extern float
get_trace_xoff_at_pos (char *filename, int beamID, d_point p);

extern float
get_trace_yoff_at_pos (char *filename, int beamID, d_point p);

extern float
eval_trace_off_at_pos (gsl_vector *coeffs, d_point p, int beamID);

extern tracestruct *
get_tracestruct_at_pos (char *filename, int beamID, d_point p);

extern void
tracestruct_fprintf (FILE * file, tracestruct * trace);

extern void
free_tracestruct (tracestruct * trace);
#endif
