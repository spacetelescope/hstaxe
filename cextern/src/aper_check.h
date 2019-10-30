/**
 * APER_CHECK - create a FITS image showing the locations of
 *               apertures.
*/

#ifndef _APER_CHECK_H
#define _APER_CHECK_H
#include <gsl/gsl_matrix.h>
#include "aXe_grism.h"


typedef struct
{
  gsl_matrix *img;
}
aXe_mask;

extern aXe_mask *
aXe_mask_init (observation * ob);

extern void
add_ap_p_to_aXe_mask (ap_pixel * ap_p, aXe_mask * mask);

extern void
mark_trace_in_aXe_mask(ap_pixel * ap_p, aXe_mask *mask);
#endif
