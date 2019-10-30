#include "spc_FITScards.h"

/**
    Routine to return a pointer to a FITScard structure with
    room allocated for n FITS header cards strings.

    @param n the number of card to allocated memory for

    @return a pointer to an allocated FITScards structure
*/
FITScards *allocate_FITScards(int n)
{
    int i;
    FITScards *res;

    res = malloc(sizeof(FITScards));

    // Allocate enough room (n+1 for NULL termination) for all the cards
    res->cards = (char **) malloc(sizeof(char *)*n);
    if ( res->cards==NULL) {
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "allocate_FITScards: Out of memory.");
    }
    for(i=0;i<n;i++) {
        res->cards[i] = (char *)malloc(sizeof(char)*FLEN_CARD);
        if (res->cards[i]==NULL) {
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "allocate_FITScards: Out of memory.");
        }
    }

    res->n = n;

    return res;
}


/**
    Routine to free a FITScards structure.

    @param cards a pointer to an allocated FITScards structure whose
    allocated memory is to be free'd.
*/
void free_FITScards(FITScards *cards)
{
    int i;
    for (i=0;i<cards->n;i++)
      free(cards->cards[i]);

    free(cards->cards);
    free(cards);
}

/*
 * The function fills the WCS header keywords for
 * the drizzled stamp images. Parameters not given
 * in the input have fixed values in echa stamp image.
 *
 * @param rval1  - the CRPIX1 value
 * @param delta1 - the CDELT1 value
 * @param rval2  - the CRPIX2 value
 *
 * @return cards - a list of filled keyword cards
 */
FITScards *
get_WCS_FITScards(const double rval1, const double delta1, const double rval2)
{
  char templt[FLEN_CARD];
  int i=0,keytype, f_status=0;

  FITScards *cards;

  cards = allocate_FITScards(10);

  i=0;

  // fill the cards for the first (dispersion) axis
  sprintf(templt,"CTYPE1 = 'WAVE' / Grating dispersion function");
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CUNIT1 = 'Angstrom' / Dispersion units");
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CRPIX1 = %f / [pixel] Reference pixel",1.0);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CRVAL1 = %f / [Angstrom] Reference value",rval1);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CDELT1 = %f / [Angstrom/pixel] dispersion",delta1);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

  // fill the cards for the second (cross-dispersion) axis
  sprintf(templt,"CTYPE2 = 'CRDIST' / Cross-dispersion distance");
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CUNIT2 = 'Pixel' / Cross-dispersion units");
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CRPIX2 = %f / [pixel] Reference pixel",1.0);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CRVAL2 = %f / [pixel] Reference value",rval2);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CDELT2 = %f / [pixel/pixel] increment",1.0);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

  return cards;
}


/**
    Routine to generate FITS header cards containing the information
    contained in a beam structure

    The following information gets written about the beam:

    OBJECTID, Object ID, ob->ID
    NBEAM,Number of beam in object, ob->nbeam
    BEAMID, This beam ID, b->ID
    REFPNTX, x coord. of reference point, b->refpoint.x
    REFPNTY, y coord. of reference point, b->refpoint.y
    C0X, x coord. of 1st bounding quadrangle corner, b->corners[0].x
    C0Y, y coord. of 1st bounding quadrangle corner, b->corners[0].y
    C1X, x coord. of 2nd bounding quadrangle corner, b->corners[1].x
    C1Y, y coord. of 2nd bounding quadrangle corner, b->corners[1].y
    C2X, x coord. of 3rd bounding quadrangle corner, b->corners[2].x
    C2Y, y coord. of 3rd bounding quadrangle corner, b->corners[2].y
    C3X, x coord. of 4th bounding quadrangle corner, b->corners[3].x
    C3Y, y coord. of 4th bounding quadrangle corner, b->corners[3].y
    BB0X, x coord. of 1st bounding box corner, b->bbox[0].x
    BB0Y, y coord. of 1st bounding box corner, b->bbox[0].y
    BB1X, x coord. of 2nd bounding box corner, b->bbox[1].x
    BB1Y, y coord. of 2nd bounding box corner, b->bbox[1].y
    WIDTH, Width of object in pixel, b->width
    ORIENT, Orientation of object in deg. wrt y-axis, CCW, b->orient
    IGNORE, Ignore flag for this beam, b->ignore

    @param ob a pointer to an existing object structure containing the beam we are interested in
    @param b a pointer to an existing beam in the passed object ob

    @return a pointer to an allocated FITScards structure

*/
FITScards *beam_to_FITScards(object *ob, int beamnum)
{
    char templt[FLEN_CARD];
    int i=0,keytype, f_status=0;
    beam *b = ob->beams+beamnum;

    FITScards *cards;

    cards = allocate_FITScards(20);

    i=0;
    sprintf(templt,"OBJECTID = %d / Object ID",ob->ID);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"NBEAM = %d / Number of beam in object",ob->nbeams);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BEAMID = %d / This beam ID",b->ID);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"REFPNTX = %f / x coord. of reference point",b->refpoint.x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"REFPNTY = %f / y coord. of reference point",b->refpoint.y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C0X = %d / x coord. of 1st bounding quadrangle corner",b->corners[0].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C0Y = %d / y coord. of 1st bounding quadrangle corner",b->corners[0].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C1X = %d / x coord. of 2nd bounding quadrangle corner",b->corners[1].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C1Y = %d / y coord. of 2nd bounding quadrangle corner",b->corners[1].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C2X = %d / x coord. of 3rd bounding quadrangle corner",b->corners[2].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C2Y = %d / y coord. of 3rd bounding quadrangle corner",b->corners[2].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C3X = %d / x coord. of 4th bounding quadrangle corner",b->corners[3].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C3Y = %d / y coord. of 4th bounding quadrangle corner",b->corners[3].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB0X = %d / x coord. of 1st bounding box corner",b->bbox[0].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB0Y = %d / y coord. of 1st bounding box corner",b->bbox[0].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB1X = %d / x coord. of 2nd bounding box corner",b->bbox[1].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB1Y = %d / y coord. of 2nd bounding box corner",b->bbox[1].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"WIDTH = %f / Width of object in pixel",b->width);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    /* Angle reference frame is converted from aXe's to SeXtractor */
    sprintf(templt,"ORIENT = %f / Orientation of object in deg. wrt y-axis, CCW",((b->orient) / M_PI * 180. - 180.));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"IGNORE = %d / Ignore flag for this beam",b->ignore);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

    return cards;
}

/**
*/
FITScards *stpmin_to_FITScards(d_point stp_min)
{
    char templt[FLEN_CARD];
    int i=0,keytype, f_status=0;

    FITScards *cards;

    cards = allocate_FITScards(2);

    i=0;
    sprintf(templt,"XINISTP = %f / starting x-coordinate of stamp image",stp_min.x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"YINISTP = %f / starting y-coordinate of stamp image",stp_min.y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

    return cards;
}


/**
    Routine to generate FITS header cards containing the information
    for drizzling and beam information contained in a beam structure

    The following information gets written about the beam:

    OBJECTID, Object ID, ob->ID
    NBEAM,Number of beam in object, ob->nbeam
    BEAMID, This beam ID, b->ID
    REFPNTX, x coord. of reference point, b->refpoint.x
    REFPNTY, y coord. of reference point, b->refpoint.y
    C0X, x coord. of 1st bounding quadrangle corner, b->corners[0].x
    C0Y, y coord. of 1st bounding quadrangle corner, b->corners[0].y
    C1X, x coord. of 2nd bounding quadrangle corner, b->corners[1].x
    C1Y, y coord. of 2nd bounding quadrangle corner, b->corners[1].y
    C2X, x coord. of 3rd bounding quadrangle corner, b->corners[2].x
    C2Y, y coord. of 3rd bounding quadrangle corner, b->corners[2].y
    C3X, x coord. of 4th bounding quadrangle corner, b->corners[3].x
    C3Y, y coord. of 4th bounding quadrangle corner, b->corners[3].y
    BB0X, x coord. of 1st bounding box corner, b->bbox[0].x
    BB0Y, y coord. of 1st bounding box corner, b->bbox[0].y
    BB1X, x coord. of 2nd bounding box corner, b->bbox[1].x
    BB1Y, y coord. of 2nd bounding box corner, b->bbox[1].y
    WIDTH, Width of object in pixel, b->width
    ORIENT, Orientation of object in deg. wrt y-axis, CCW, b->orient
    IGNORE, Ignore flag for this beam, b->ignore

    @param ob a pointer to an existing object structure containing the beam we are interested in
    @param b a pointer to an existing beam in the passed object ob

    @return a pointer to an allocated FITScards structure

*/
FITScards *
drzinfo_to_FITScards(object *ob, int beamnum, d_point outref,
		     aperture_conf *conf, gsl_matrix *drizzcoeffs,
		     int trlength, double relx, double rely,
		     double objwidth, d_point refwave_pos, float sky_cps,
		     double drizzle_width, double cdref, double spcorr)
{
    //double dx, dy;
    char templt[FLEN_CARD];
    int i=0,keytype, f_status=0;
    beam *b = ob->beams+beamnum;

    FITScards *cards;

    cards = allocate_FITScards(56);

    i=0;
    sprintf(templt,"OBJECTID = %d / Object ID",ob->ID);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"NBEAM = %d / Number of beam in object",ob->nbeams);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BEAMID = %d / This beam ID",b->ID);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"REFPNTX = %f / x coord. of reference point",b->refpoint.x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"REFPNTY = %f / y coord. of reference point",b->refpoint.y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C0X = %d / x coord. of 1st bounding quadrangle corner",b->corners[0].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C0Y = %d / y coord. of 1st bounding quadrangle corner",b->corners[0].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C1X = %d / x coord. of 2nd bounding quadrangle corner",b->corners[1].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C1Y = %d / y coord. of 2nd bounding quadrangle corner",b->corners[1].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C2X = %d / x coord. of 3rd bounding quadrangle corner",b->corners[2].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C2Y = %d / y coord. of 3rd bounding quadrangle corner",b->corners[2].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C3X = %d / x coord. of 4th bounding quadrangle corner",b->corners[3].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"C3Y = %d / y coord. of 4th bounding quadrangle corner",b->corners[3].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB0X = %d / x coord. of 1st bounding box corner",b->bbox[0].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB0Y = %d / y coord. of 1st bounding box corner",b->bbox[0].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB1X = %d / x coord. of 2nd bounding box corner",b->bbox[1].x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"BB1Y = %d / y coord. of 2nd bounding box corner",b->bbox[1].y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"WIDTH = %f / Width of object in pixel",b->width);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZWIDTH = %f / Width of object in pixel",drizzle_width);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"LENGTH= %i / Length of object in pixel",trlength);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"OWIDTH= %f / Maximum width of object in pixel",objwidth);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    /* Angle reference frame is converted from aXe's to SeXtractor */
    sprintf(templt,"ORIENT = %f / Orientation of object in deg. wrt y-axis, CCW",((b->orient) / M_PI * 180. - 180.));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"IGNORE = %d / Ignore flag for this beam",b->ignore);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"OREFPNTX = %f / x coord. of reference point",outref.x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"OREFPNTY = %f / y coord. of reference point",outref.y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"RREFPNTX = %f / x coord. of refpoint in stamp image",relx);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"RREFPNTY = %f / y coord. of refpoint in stamp image",rely);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"WAVPNTX = %f / x coord. of reference point", refwave_pos.x);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"WAVPNTY = %f / y coord. of reference point", refwave_pos.y);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"XOFFS = %f / x-offset after drizzling",conf->drz_xstart);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"LAMBDA0 = %f / wavelength drizzled to",conf->drz_lamb0);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DLAMBDA = %f / wavelength increment drizzled to",conf->drz_resol);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"SKY_CPS = %f / electrons per second in the background",sky_cps);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"CDSCALE = %f / pixel scale in cross dispersion",cdref);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"SPCORR = %f / correction factor resolution",spcorr);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"SLITWIDT = %f / slit width [pix]",b->slitgeom[2]);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ00 = %e / drizzle coefficient (0, 0)",gsl_matrix_get(drizzcoeffs, 0, 0));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ01 = %e / drizzle coefficient (0, 1)",gsl_matrix_get(drizzcoeffs, 0, 1));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ02 = %e / drizzle coefficient (0, 2)",gsl_matrix_get(drizzcoeffs, 0, 2));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ03 = %e / drizzle coefficient (0, 3)",gsl_matrix_get(drizzcoeffs, 0, 3));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ04 = %e / drizzle coefficient (0, 4)",gsl_matrix_get(drizzcoeffs, 0, 4));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ05 = %e / drizzle coefficient (0, 5)",gsl_matrix_get(drizzcoeffs, 0, 5));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ06 = %e / drizzle coefficient (0, 6)",gsl_matrix_get(drizzcoeffs, 0, 6));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ07 = %e / drizzle coefficient (0, 7)",gsl_matrix_get(drizzcoeffs, 0, 7));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ08 = %e / drizzle coefficient (0, 8)",gsl_matrix_get(drizzcoeffs, 0, 8));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ09 = %e / drizzle coefficient (0, 9)",gsl_matrix_get(drizzcoeffs, 0, 9));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ10 = %e / drizzle coefficient (1, 0)",gsl_matrix_get(drizzcoeffs, 1, 0));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ11 = %e / drizzle coefficient (1, 1)",gsl_matrix_get(drizzcoeffs, 1, 1));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ12 = %e / drizzle coefficient (1, 2)",gsl_matrix_get(drizzcoeffs, 1, 2));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ13 = %e / drizzle coefficient (1, 3)",gsl_matrix_get(drizzcoeffs, 1, 3));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ14 = %e / drizzle coefficient (1, 4)",gsl_matrix_get(drizzcoeffs, 1, 4));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ15 = %e / drizzle coefficient (1, 5)",gsl_matrix_get(drizzcoeffs, 1, 5));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ16 = %e / drizzle coefficient (1, 6)",gsl_matrix_get(drizzcoeffs, 1, 6));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ17 = %e / drizzle coefficient (1, 7)",gsl_matrix_get(drizzcoeffs, 1, 7));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ18 = %e / drizzle coefficient (1, 8)",gsl_matrix_get(drizzcoeffs, 1, 8));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"DRZ19 = %e / drizzle coefficient (1, 9)",gsl_matrix_get(drizzcoeffs, 1, 9));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

    return cards;
}

FITScards *
nicbck_info_to_FITScards(const double skypix_frac, const double scale_factor,
			 const double offset)
{
    char templt[FLEN_CARD];
    int i=0,keytype, f_status=0;

    FITScards *cards;

    cards = allocate_FITScards(3);

    i=0;
    sprintf(templt,"SKY_SCAL = %e / Scale factor for background", scale_factor);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"SKY_OFFS = %e / Offset value for background", offset);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
    sprintf(templt,"F_SKYPIX = %e / Fraction of sky pixels", skypix_frac);
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);

    return cards;
}

/**
	A helper function which generates FITS keywords containing
	the content of a disp structure.

	@param disp A pointer to a dispstruct structure

    @return a pointer to an allocated FITScards structure
*/
FITScards *dispstruct_to_FITScards(dispstruct *disp)
{
  char templt[FLEN_CARD];
  int i=0,j,keytype, f_status=0;
  int numcards=0;
  FITScards *cards;

  numcards += disp->pol->size;
  numcards +=1 ; /* for for_grism */
  //numcards +=2 ; /* for offset */
  numcards +=2 ; /* for cpoint */
  cards = allocate_FITScards(numcards);

  i=0;
  for (j=0; j<(int)disp->pol->size;j++) {
    sprintf(templt,"DLDX%d = %g / Wavelength solution parameter",j,gsl_vector_get(disp->pol,j));
    fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  }

  sprintf(templt,"FORGRISM = %d / If 1 then grism polynomial used. Prism otherwise.",disp->for_grism);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CPOINTX = %f / X position where wavelength disp. was computed",disp->cpoint.x);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  sprintf(templt,"CPOINTY = %f / X position where wavelength disp. was computed",disp->cpoint.y);
  fits_parse_template (templt, cards->cards[i++], &keytype, &f_status);
  return cards;
}

/**
 * The function updates the keyword "CONTAM" in a fits file.
 * Update means that an existing entry is replace.
 * It is used to update this keyword in a PET files
 * within the command aXe_PETCONT.
 *
 * @param fptr       - pointer to the fits file
 * @param model_name - the name of the model
 *
 */
void
update_contam_model(fitsfile *fptr, char model_name[])
{

 //char templt[FLEN_CARD];
 char keyw[FLEN_KEYWORD];
 char comment[FLEN_COMMENT];

 int  f_status=0;

 // fill the keyword valiable
 sprintf(keyw,"CONTAM");
 sprintf(comment,"Name of the contamination model");

 // use the fitsio routine to update the keyword
 // react in case of an error
 fits_update_key(fptr, TSTRING, keyw, model_name,comment, &f_status);
 if (f_status)
   {
     ffrprt (stderr, f_status);
     aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		  "update_contam_model: Could not update the contamination model name.\n");
   }
}

/**
 * The function reads the keyword "CONTAM" in an open FITS file
 * and ananlyzes its content. In case that the content indicates
 * a quantitative contamination model, 1 is returned. If not
 * the value 0 is given back.
 *
 * @param fptr - pointer to the FITS-file
 */
int
check_quantitative_contamination(fitsfile *fptr)
{

  char keyval[FLEN_CARD];
  char keyw[FLEN_KEYWORD];
  char comment[FLEN_COMMENT];

  int  f_status=0;
  int ret;

  // fill the keyword valiable
  sprintf(keyw,"CONTAM");

  // use the fitsio routine to read the keyword
  // react in case of an error
  fits_read_key(fptr, TSTRING, keyw, keyval, comment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		   "check_quantitative_contamination:\nCould not read keyword %s from fits-file.\n", keyw);
    }

  // check whether the keyword was there
  if (strlen(keyval) >0)
    {
      // check the keyword content agaist the
      // names of the non quantitative contamination models
      if (!strcmp(keyval,"GEOM"))
	ret=0;
      // set the quantitative contamination flag
      else
	ret=1;
    }
  // in case that the keyword does not exist,
  // the flag is not set
  else
    {
      ret=0;
    }

  // return the flag
  return ret;
}

/**
 * The function transports the keywird entry for
 * the keyword "CONTAM" from one FITS file to another.
 * As the parameter names suggest, this works particularly
 * form drizzle grism images to PET's
 *
 * @param grism_file_path - the name of the input FITS file
 * @param PET_file_path   - the name of the output FITS file
 *
 */
void
transport_cont_modname(char grism_file_path[], char PET_file_path[])
{

  fitsfile *in_fptr;
  fitsfile *ou_fptr;

  char keyval[FLEN_CARD];
  char keyw[FLEN_KEYWORD];
  char comment[FLEN_COMMENT];

  int  f_status=0;

  // open the input FITS file
  fits_open_file (&in_fptr, grism_file_path, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "transport_cont_modname: Could not open file %s .\n", grism_file_path);
    }

  // name the keyword variable
  sprintf(keyw,"CONTAM");

  // read the contamination keyword from the input file
  fits_read_key(in_fptr, TSTRING, keyw, keyval, comment, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "transport_cont_modname: Could not read keyword %s from %s.\n", keyw, grism_file_path);
    }

  // close the input FITSfile
  fits_close_file (in_fptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "transport_cont_modname: Could not close file %s .\n", grism_file_path);
    }

  // uptdate that keyword in the output FITS-file
  // in case that it is not empty
  if (strlen(keyval) > 0)
    {

      // open the output FITS file
      fits_open_file (&ou_fptr, PET_file_path, READWRITE, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "transport_cont_modname: Could not open file %s .\n", PET_file_path);
	}

      // use the function to update that keyword
      update_contam_model(ou_fptr, keyval);


      // close the output FITS file
      fits_close_file (ou_fptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "transport_cont_modname: Could not close file %s .\n", PET_file_path);
	}
    }
}
