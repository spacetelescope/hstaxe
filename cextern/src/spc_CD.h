#ifndef SPC_CD_H
#define SPC_CD_H

#include <string.h>

#include "fitsio.h"
#include "aXe_grism.h"
#include "wcs.h"
#include "aXe_errors.h"


/**
    A point in world coordinates.

*/
typedef struct
{
  // RA
  double ra;
  // dec
  double dec;
  // equinox
  char *equinox;
}
sky_coord;

extern sky_coord *
ij_to_radec (struct WorldCoor *wcs, d_point ij);

extern d_point *
radec_to_ij (struct WorldCoor *wcs, sky_coord radec);

extern char *
get_fits_header (char filenamep[], int hdunum);

extern struct WorldCoor *
get_wcs (char filename[], int hdunum);
#endif
