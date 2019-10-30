#ifndef SPC_FITSCARDS_H
#define SPC_FITSCARDS_H

#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "aXe_grism.h"
#include "fitsio.h"
#include "disp_conf.h"
#include "aper_conf.h"

typedef struct
{
  int n;          /* Number of cards in the **cards array */
  char **cards;   /* A pointer to an array of strings     */
} FITScards;

extern FITScards *
allocate_FITScards(int n);

extern void
free_FITScards(FITScards *cards);

extern FITScards *
get_WCS_FITScards(const double rval1, const double delta1, const double rval2);

extern FITScards *
beam_to_FITScards(object *ob, int beamnum);

extern FITScards *
stpmin_to_FITScards(d_point stp_min);

extern FITScards *
drzinfo_to_FITScards(object *ob, int beamnum, d_point outref,
		     aperture_conf *conf, gsl_matrix *drizzcoeffs,
		     int trlength, double relx, double rely,
		     double objwidth, d_point refwave_pos,
		     float sky_cps, double drizzle_width,
		     double cdcorr, double spcorr);

extern FITScards *
nicbck_info_to_FITScards(const double skypix_frac, const double scale_factor,
			 const double offset);

extern FITScards *
dispstruct_to_FITScards(dispstruct *disp);

extern void
update_contam_model(fitsfile *fptr, char model_name[]);

extern int
check_quantitative_contamination(fitsfile *fptr);

extern void
transport_cont_modname(char grism_file_path[], char PET_file_path[]);

#endif
