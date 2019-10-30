/**
 */
#ifndef _INOUT_APER_H

#define _INOUT_APER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "aXe_errors.h"
#include "spc_cfg.h"
#include "spc_trace_functions.h"

#define APER_MAXLINE 14

extern int
nbeams_from_char_array2 (char **apers, int num);

extern gsl_vector_int *
nbeams_from_char_array (char **apers, int num);

extern int 
object_list_to_file (object * const *oblist, char *filename,
                     int leaveout_ignored);

extern int
get_beam_from_aper_file (char *filename, int aperID, int beamID,beam * b);

extern gsl_vector_int *
aper_file_aperlist (char *filename);

extern int 
aper_file_apernum (char *filename);

extern object *
get_aperture_from_aper_file (char *filename, int aperID);

extern object **
file_to_object_list (char filename[], observation * obs);

extern char **
return_next_aperture(FILE *input);

extern object **
file_to_object_list_seq (char filename[], observation * obs);

extern int 
find_object_in_object_list(object **oblist, const int ID);

extern beam
find_beam_in_object_list(object **oblist, const int objID,
                         const int beamID);

extern beam *
find_beamptr_in_object_list(object **oblist, const int objID,
                            const int beamID);

extern void
refurbish_object_list(object **oblist, const int new_default,
                      const int old_value, const int new_value);

extern int
object_list_size(object **oblist);

extern int
get_beamspec_size(object **oblist);
#endif
