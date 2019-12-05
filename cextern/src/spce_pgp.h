/**
 * File: spce_pgp.h
 * Interface to PGplot routines.
 *
 */

#ifndef _SPCE_PGP_H

#define _SPCE_PGP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_pathlength.h"
#include "spce_output.h"

extern void
pgp_stamp_image (const ap_pixel * const ap_p, const int n_sub,
		 const object * ob, int beamnum, char filename[],
		 const int negative);
#endif /* !_SPCE_PGP_H */
