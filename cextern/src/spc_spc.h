#ifndef _SPC_SPC_H

#define _SPC_SPC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include "aXe_grism.h"
#include "aXe_utils.h"
#include "fitsio.h"

#define SPCTOL 1.0e-5

// interpolation type used for the
// response function
#define RESP_FUNC_INTERP_TYPE gsl_interp_linear

extern spectrum *
allocate_spectrum (const int numbin);

extern void
free_spectrum (spectrum * const spec);

extern void
fprintf_spectrum (FILE * output, const spectrum * const sp);

extern spectrum *
subtract_spectra (spectrum * a, spectrum * b);

extern int
get_ID_index_to_SPC (char filename[], char ID[]);

extern void
add_ID_index_to_SPC (char filename[], char ID[], int hdunum);

extern int
add_ID_to_SPC (char filename[], int N, char ID[]);

extern int
find_ID_in_SPC (char filename[], char ID[]);

extern void
create_SPC (char filename[], int overwrite);

extern int
get_SPC_colnum (fitsfile * input, char colname[]);

extern void
add_spectra_to_SPC (char filename[], spectrum * obj_spec,
		    spectrum * bck_spec, spectrum * sobj_spec,
		    int aperID, int beamID);

extern void
add_data_to_SPC (spectrum * spec, char countcolname[],
		 char errorcolname[], char weightcolname[], char ID[],
		 char filename[], int hdunum, long N);

extern spectrum *
trim_spectrum (spectrum * spc);

extern spectrum *
empty_counts_spectrum_copy (spectrum * a);

extern void
add_ID_to_SPC_opened (fitsfile *output, int N, char ID[]);

extern void
add_spectra_to_SPC_opened (fitsfile *input, spectrum * obj_spec,
			   spectrum * bck_spec, spectrum * sobj_spec,
			   int aperID, int beamID);

extern void
add_data_to_SPC_opened (spectrum * spec, char countcolname[],
			char errorcolname[], char weightcolname[],
			char ID[], fitsfile *input, long N);

extern fitsfile *
create_SPC_opened (char filename[], int overwrite);

#endif
