/**
 */

#ifndef _DISP_CONF_H
#define _DISP_CONF_H
#include <gsl/gsl_vector.h>
#include "aXe_grism.h"


/**
 * A structure to contain the coefficient of a 2D
 * polynomial that can be used to compute the given
 * coefficient of the polynomial dispersion relation
 * as at a particular i.j location in the image
 */
typedef struct
{
  gsl_vector *pol; /* A vector containing the 2D polynomial coefficients */
  int for_grism;   /* If true, pol is of form sum(a_i x^i, i=0..order)
                      otherwise, it is of form sum(a_i 1/x^i, i=0..order) */
  d_point cpoint;  /* The detector location at which this structure
                      was computed */
  int ID;          /* The ID of the beam for which this coefficient
                      is defined */
  char file[MAXCHAR];
}
dispstruct;

typedef struct
{
  char file[MAXCHAR];

  int ID;
  int for_grism;

  int n_order;
  gsl_vector **all_coeffs;
}
global_disp;

extern float
get_beam_mmag_extract (char *filename, int beamID);

extern float
get_beam_mmag_mark (char *filename, int beamID);

extern int
get_beam_disp_norder (char *filename, int beamID);

extern gsl_vector *
get_beam_disp_order (char *filename, const int for_grism,
                                 int beamID, int order);

extern float
get_disp_coeff_at_pos (char *filename, const int for_grism, int beamID,
                             int order, d_point p);

extern gsl_vector *
get_disp_coeffs_at_pos (char *filename, const int for_grism,
                                    int beamID, d_point p);
extern void
free_dispstruct(dispstruct *p);

extern dispstruct *
get_dispersion_at_pos (char *filename, int beamID, d_point p);

extern void
dispstruct_fprintf (FILE * file, dispstruct * disp);

extern dispstruct *
get_dispstruct_at_pos (char *filename, const int for_grism,
                                   int beamID, d_point p);

extern int
check_for_grism (char *filename, int beamID);

extern int
check_disp_coeff (char *filename,const int for_grism,int beamID,int order);

extern gsl_vector *
get_prange (char *filename, int beamID);

extern global_disp *
get_global_disp(char *filename, const int for_grism, int beamID);

extern gsl_vector *
get_calvector_from_gdisp(const global_disp *gdisp, const d_point p);

extern double
get_2Ddep_value(const gsl_vector *coeffs, const d_point p);

extern void
print_global_disp(const global_disp *gdisp);

extern int
check_disp_order(const gsl_vector *tmp, const int beamID, const int order);

extern void
free_global_disp(global_disp *gdisp);

#endif
