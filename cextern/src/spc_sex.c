/**
 * Subroutines used mostly for working on catalogs
 * in the tasks SEX2GOL and GOL2AF.
 *
 * Howard Bushouse, STScI, 04-Mar-2011, version 1.4
 * Changed FATAL error to WARN3 when get_valid_entries returns no valid
 * magnitudes for an object.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

#include "aXe_utils.h"
#include "aXe_grism.h"
#include "spc_cfg.h"
#include "disp_conf.h"
#include "aper_conf.h"
#include "trace_conf.h"
#include "spc_CD.h"
#include "spc_wl_calib.h"
#include "aper_conf.h"
#include "specmodel_utils.h"
#include "model_utils.h"
#include "spc_fluxcube.h"
#include "spc_model.h"
#include "crossdisp_utils.h"
#include "spc_utils.h"
#include "spc_sex.h"

// define SQUARE
#define SQR(x) ((x)*(x))

/**
 *
 *  Function: create_SexObject
 *  Allocates and creates a SExtractor object from
 *  a line of a SExtractor catalog.
 *
 *  Parameters:
 *  @param  actinfo  - the header structure with the column names
 *  @param  line     - the line to be transformed into a SexObject
 *  @param  waves    - vector with the walelengths of the magnitude columns
 *  @param  cnums    - vector with the column numbers of the magnitude columns
 *  @param  magcol   - column number of the magnitude column selected
 *                     for mag_auto
 *  @param  fullinfo - deecides whether a column is optional or not
 *
 *  Returns:
 *  @return o        - the new created SexObject
 */
SexObject *
create_SexObject(const colinfo *actinfo, char *line, const gsl_vector * waves,
                 const gsl_vector *  cnums, const px_point backwin_cols,
                 const px_point modinfo_cols, const int magcol, const int fullinfo)
{
  gsl_vector *v;
  gsl_vector *mags;
  gsl_vector *wavs;
  SexObject  *o;
  int         i=0;
  // transform the data in the line to a vector
  lv1ws (line);
  v = string_to_gsl_array (line);

  // allocate space for the SexObject
  o = (SexObject *) malloc (sizeof (SexObject));
  if (o == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "create_SexObject: Could not allocate memory.");

  // succesively fill the data into the SexObject.
  // for 'fullinfo=1' missing data results in an error.
  // otherwise NaN is stored into the SexObject
  o->number = (int) get_col_value (actinfo, "NUMBER", v, 1);

  o->xy_image.x     = get_col_value (actinfo, "X_IMAGE", v,fullinfo);
  o->xy_image.y     = get_col_value (actinfo, "Y_IMAGE", v,fullinfo);

  o->xy_world.ra    = get_col_value (actinfo, "X_WORLD", v,1);
  o->xy_world.dec   = get_col_value (actinfo, "Y_WORLD", v,1);

  o->el_image.a     = get_col_value (actinfo, "A_IMAGE", v,fullinfo);
  o->el_image.b     = get_col_value (actinfo, "B_IMAGE", v,fullinfo);
  o->el_image.theta = get_col_value (actinfo, "THETA_IMAGE", v,fullinfo);

  o->el_world.a     = get_col_value (actinfo, "A_WORLD", v,1);
  o->el_world.b     = get_col_value (actinfo, "B_WORLD", v,1);

  // Some clarification is needed for the next two lines.
  // to have get_col_value and get_col_value2 is not
  // satisfying!! PROBLEM
  o->el_world.theta = get_col_value2 (actinfo, "THETA_WORLD", v,0);
  if (o->el_world.theta <-900000){
    //  o->el_world.theta = get_col_value (actinfo, "THETA_WORLD", v,0);
    //  if (isnan(o->el_world.theta)){
    o->el_world.theta = get_col_value (actinfo, "THETA_SKY", v,1);
  }

  if (backwin_cols.x !=-1)
    {
      o->backwindow.x  = get_col_value (actinfo, "BACKWINDOW", v,1);
      o->backwindow.y  = -1;
      // old code for FORS2-MXU
      // o->backwindow.y     = get_col_value (actinfo, "BACKWIN_LOW", v,1);

    }
  else
    {
      o->backwindow.x     = -1;
      o->backwindow.y     = -1;
    }

  // transfer the index for the model spectrum
  if (modinfo_cols.x > -1)
    o->modspec = (int)get_col_value(actinfo, "MODSPEC", v,1);
  else
    o->modspec = -1;

  // transfer the index value for the object shape
  if (modinfo_cols.y > -1)
    o->modimage = (int)get_col_value(actinfo, "MODIMAGE", v,1);
  else
    o->modimage = -1;

  // check whether the MAG_AUTO column exists
  if ((int)waves->size == 1 && (int)gsl_vector_get (waves, 0) == 0)
    {


      // fill the MAG_AUTO value
      // put the magnitude and wavelength vectors to NULL
      o->magnitudes = NULL;
      o->lambdas    = NULL;
      o->mag_auto   = get_col_value (actinfo, "MAG_AUTO", v, 1);
    }
  else
    {
      // create the vectors for the magnitude values
      mags = gsl_vector_alloc ((int)waves->size);
      wavs = gsl_vector_alloc ((int)waves->size);

      // fill the magnitude values into the vectors and store them in
      // the SexObject.
      for (i=0; i < (int)waves->size; i++){
        gsl_vector_set(mags, i, gsl_vector_get (v, (int)gsl_vector_get (cnums, i)-1));
        gsl_vector_set(wavs, i, gsl_vector_get (waves, i));
      }

      // fill the MAG_AUTO value
      o->magnitudes = mags;
      o->lambdas    = wavs;
      o->mag_auto   = gsl_vector_get (v, magcol-1);
    }

  return o;
}


/**
 *  Function: SexObject_fprintf
 *  Function to print the attributes of a SexObbject
 *  to an output stream.
 *
 *  Paramters:
 *  @param output - an ouput file or stream
 *  @param o      - a pointer to a SexObject
 */
void
SexObject_fprintf (FILE * output, SexObject * o)
{
     fprintf (output, "NUMBER: %d\n", o->number);

     fprintf (output, "X_IMAGE: %e\n", o->xy_image.x);
     fprintf (output, "Y_IMAGE: %e\n", o->xy_image.y);

     fprintf (output, "X_WORLD: %e\n", o->xy_world.ra);
     fprintf (output, "Y_WORLD: %e\n", o->xy_world.dec);

     fprintf (output, "A_IMAGE: %e\n", o->el_image.a);
     fprintf (output, "B_IMAGE: %e\n", o->el_image.b);
     fprintf (output, "THETA_IMAGE: %e\n", o->el_image.theta);

     fprintf (output, "A_WORLD: %e\n", o->el_world.a);
     fprintf (output, "B_WORLD: %e\n", o->el_world.b);
     fprintf (output, "THETA_WORLD: %e\n", o->el_world.theta);

     fprintf (output, "MAG_AUTO: %e\n", o->mag_auto);

}


/**
 *  Function: sobs_to_vout
 *  The function transforms a SexObject into
 *  a vector.
 *
 *  Parameters:
 *  @param  sobs - the SexObject to be transformed
 *
 *  Returns:
 *  @return vout - the vector with the Seobject data
 */
gsl_vector *
sobs_to_vout(const SexObject *sobs)
{

  int nentries, i;
  int count=0;

  // determine the size of the vector
  if (sobs->magnitudes){
    nentries = 11+sobs->magnitudes->size;
  }
  else{
    nentries = 12;
  }

  if(sobs->backwindow.x != -1)
    nentries = nentries+1;
    // old FORS2-MXU code:
    //nentries = nentries+2;

  if(sobs->modspec != -1)
    nentries++;
  if(sobs->modimage != -1)
    nentries++;

  // allocate space for the vector
  gsl_vector *vout = gsl_vector_alloc (nentries);

  // fill the vector
  gsl_vector_set (vout, count++, sobs->number);
  gsl_vector_set (vout, count++, sobs->xy_world.ra);
  gsl_vector_set (vout, count++, sobs->xy_world.dec);
  gsl_vector_set (vout, count++, sobs->el_world.a);
  gsl_vector_set (vout, count++, sobs->el_world.b);
  gsl_vector_set (vout, count++, sobs->el_world.theta);
  gsl_vector_set (vout, count++, sobs->xy_image.x);
  gsl_vector_set (vout, count++, sobs->xy_image.y);
  gsl_vector_set (vout, count++, sobs->el_image.a);
  gsl_vector_set (vout, count++, sobs->el_image.b);
  gsl_vector_set (vout, count++, sobs->el_image.theta);

  if(sobs->backwindow.x != -1)
    {
      gsl_vector_set (vout, count++, sobs->backwindow.x);
      // old code for FORS2-MXU
      //gsl_vector_set (vout, count++, sobs->backwindow.y);
    }

  if(sobs->modspec != -1)
      gsl_vector_set (vout, count++, sobs->modspec);

  if(sobs->modimage != -1)
      gsl_vector_set (vout, count++, sobs->modimage);

  // see whether there is MAG_AUTO
  if (sobs->magnitudes){
    // fill in the magnitude values
    for (i=0; i < (int)sobs->magnitudes->size; i++)
      gsl_vector_set (vout, count++, gsl_vector_get(sobs->magnitudes, i));
  }
  else{
    // fill in the MAG_AUTO value
    gsl_vector_set (vout, count++, sobs->mag_auto);
  }

  // return the vector
  return vout;
}

/**
 *  Function: catalog_to_wcs
 *  This function reads a sextractor catalog and using an input and an output
 *  CD matrix, converts the world ra and dec coordinates to the output CD
 *  matrix's image pixel coordinates. Non existing (i.e. NaN) world RA and
 *  DEC values are generated using
 *  the pixel coordinates and the input CD matrix
 *
 *  Parameters:
 *  @param - infile a pointer to a char array containing the name of the
 *         input catalog
 *  @param - outfile a pointer to a char array containing the name of the
 *         output catalog
 *  @param - from_wcs a pointer to an exisiting WoldCoor structure containing
 *         the input catalog CD matrix
 *  @param - to_wcs a pointer to an exisiting WoldCoor structure containing
 *         the output catalog CD matrix
 *  @param - overwrite_wcs if set to 1 then input wold coordinate are
 *         recomputed first using the from_wcs and the x,y image values
 *  @param - overwrite_img if set to 1 then the x,y coordinates of the input
 *         image are first computed using the from_wcs and the
 *          world ra and dec info
 */
void
catalog_to_wcs (char grismfile[], int hdunum, char infile[], char outfile[],
                struct WorldCoor *from_wcs, struct WorldCoor *to_wcs,
                int distortion, int overwrite_wcs, int overwrite_img)
{
  FILE *fin, *fout;
  gsl_vector *v;
  gsl_vector *v_out;
  gsl_vector *waves;
  gsl_vector *cnums;
  gsl_matrix *coeffs=NULL;
  char Buffer[CATBUFFERSIZE];
  char line[CATBUFFERSIZE], str[CATBUFFERSIZE];
  int i, hasmags=0, magcencol = 0;
  SexObject *o;
  colinfo * actcatinfo;
  px_point    pixmax;
  px_point    backwin_cols;
  px_point    modinfo_cols;


  actcatinfo = get_sex_col_descr (infile);
  hasmags = has_magnitudes(actcatinfo);
  if (!hasmags)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "No magnitudes in file %s", infile);
    }
  else{
    waves = gsl_vector_alloc (hasmags);
    cnums = gsl_vector_alloc (hasmags);
    hasmags = get_magcols(actcatinfo, waves, cnums);
    magcencol = get_magauto_col(waves, cnums, 800.0);
  }

  backwin_cols = has_backwindow(actcatinfo);
  modinfo_cols = has_modelinfo(actcatinfo);


  if (!(fin = fopen (infile, "r")))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Cannot open catalog file %s", infile);
    }
  if (!(fout = fopen (outfile, "w")))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Cannot open catalog file %s", outfile);
    }

  if (check_worldcoo_input(actcatinfo, 0))
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "Catalogue %s does not have all WCO columns\n", infile);

  if (check_imagecoo_input(actcatinfo))
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "Catalogue %s does not have all image coo columns\n", infile);

  make_GOL_header(fout, actcatinfo, waves, cnums, backwin_cols, modinfo_cols);

  if (distortion)
    coeffs = get_crossdisp_matrix(grismfile, hdunum);

  pixmax = get_npixels (grismfile, hdunum);

  while (fgets (line, CATBUFFERSIZE, fin))
    {
      if (!(line_is_valid (actcatinfo, line)))
        continue;
       if (line[0] == ';')
         continue;

      /* Create a vector containing the information */
      lv1ws (line);
      v = string_to_gsl_array (line);

      //      o = create_SexObject (sex_col_desc, line);
      o = create_SexObject (actcatinfo, line, waves, cnums, backwin_cols, modinfo_cols, magcencol,  1);

      if ( (from_wcs!=NULL)&&(to_wcs!=NULL) )
        {
          //      if (distortion){
          //        fprintf(stdout, "(%f,%f) &&-->&& ", o->xy_image.x, o->xy_image.y);
          //        o->xy_image = undistort_point(coeffs, pixmax, o->xy_image);
          //      }
          fill_missing_WCS_coordinates (o, from_wcs, overwrite_wcs);
          fill_missing_image_coordinates (o, from_wcs, overwrite_img);
          compute_new_image_coordinates (o, to_wcs);
          if (distortion){
            o->xy_image = distort_point(coeffs, pixmax, o->xy_image);
            //      fprintf(stdout, "(%f,%f)\n\n", o->xy_image.x, o->xy_image.y);
          }
        }

      v_out = sobs_to_vout(o);

      sprintf (Buffer, "%8.5g ", gsl_vector_get (v_out, 0));
      for (i = 1; i < (int)v_out->size; i++)
        {
          sprintf (str, " %8.5g ", gsl_vector_get (v_out, i));
          strcat (Buffer, str);
        }
      strcat (Buffer, "\n");

      fputs (Buffer, fout);

      gsl_vector_free (v_out);
      gsl_vector_free (v);
  }

  if (coeffs)
    gsl_matrix_free(coeffs);

  fclose (fin);
  fclose (fout);
}

/**
 * Function:catalog_to_wcs_nodim
 * This function reads a sextractor catalog and using an input and an output
 * CD matrix, converts the world ra and dec coordinates to the output
 * CD matrix's image pixel coordinates. Non existing (i.e. NaN) world RA
 * and DEC values are generated using the pixel coordinates and the input
 * CD matrix
 *
 * Parameters:
 * @param infile    - a pointer to a char array containing the name of the
 *                    input catalog
 * @param outfile   - a pointer to a char array containing the name of the
 *                    output catalog
 * @param from_wcs  - a pointer to an exisiting WoldCoor structure containing
 *                    the input catalog CD matrix
 * @param to_wcs    - a pointer to an exisiting WoldCoor structure
 *                    containing the output catalog CD matrix
 * @param overwrite_wcs - if set to 1 then input wold coordinate are
 *                        recomputed first using the from_wcs and the
 *                        x,y image values
 * @param overwrite_img - if set to 1 then the x,y coordinates of the
 *                        input image are first computed using the
 *                        from_wcs and the world ra and dec info
 */
void
catalog_to_wcs_nodim (char infile[], char outfile[],
                      struct WorldCoor *grism_wcs, int overwrite_wcs,
                      int overwrite_img)
{
  FILE *fin, *fout;
  gsl_vector *v;
  gsl_vector *v_out;
  gsl_vector *waves;
  gsl_vector *cnums;
  char Buffer[CATBUFFERSIZE];
  char line[CATBUFFERSIZE], str[CATBUFFERSIZE];
  int i;
  SexObject *o;
  int compute_imcoos=0;
  int th_sky=0;
  int checksum, hasmags=0, magcencol=0;
  colinfo * actcatinfo;
  px_point    backwin_cols;
  px_point    modinfo_cols;

  actcatinfo = get_sex_col_descr (infile);
  hasmags = has_magnitudes(actcatinfo);
  if (!hasmags)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "No magnitudes in file %s", infile);
    }
  else{
    waves = gsl_vector_alloc (hasmags);
    cnums = gsl_vector_alloc (hasmags);
    hasmags = get_magcols(actcatinfo, waves, cnums);
    magcencol = get_magauto_col(waves, cnums, 800.0);
  }

  backwin_cols = has_backwindow(actcatinfo);
  modinfo_cols = has_modelinfo(actcatinfo);

  if (!(fin = fopen (infile, "r")))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Cannot open catalog file %s", infile);
    }
  if (!(fout = fopen (outfile, "w")))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Cannot open catalog file %s", outfile);
    }

  //  checksum = check_worldcoo_input(sex_col_desc);
  th_sky=1;
  checksum = check_worldcoo_input(actcatinfo, 1);
  if (checksum)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "Catalogue %s does not have all necessary columns\n", infile);

  //  if (check_imagecoo_input(sex_col_desc))
  if (check_imagecoo_input(actcatinfo))
    compute_imcoos=1;

  make_GOL_header(fout, actcatinfo, waves, cnums, backwin_cols, modinfo_cols);

  while (fgets (line, CATBUFFERSIZE, fin))
    {

      /* If line is not a valid catalog entry, just continue */
      if (!(line_is_valid (actcatinfo, line)))
        continue;

      /* if line starts with ";", just continue */
      if (line[0] == ';')
        continue;

      /* Create a vector containing the information */
      lv1ws (line);
      v = string_to_gsl_array (line);

      o = create_SexObject (actcatinfo, line, waves, cnums, backwin_cols, modinfo_cols, magcencol,  0);


      if (compute_imcoos)
        compute_new_image_sexobject (o, grism_wcs, th_sky);
      v_out = sobs_to_vout(o);


      sprintf (Buffer, "%8.5g ", gsl_vector_get (v_out, 0));
      for (i = 1; i < (int)v_out->size; i++)
        {
          sprintf (str, " %8.5g ", gsl_vector_get (v_out, i));
          strcat (Buffer, str);
        }
      strcat (Buffer, "\n");

      fputs (Buffer, fout);

      gsl_vector_free (v_out);
      gsl_vector_free (v);
    }

  fclose (fin);
  fclose (fout);

  free (actcatinfo);

}

/**
 * Function: SexMags_to_beamFlux
 * Fill an beam structure with the flux and shape information information
 * contained in a SexObject structure.
 *
 * Parameters:
 * @param sobj     - pointer to a  SexObject
 * @param actbeam  - the beam structure to be filled
 *
 * Returns:
 * @return -
 */
void
SexMags_to_beamFlux(SexObject * sobj, beam *actbeam)
  {
    int nvalid=0;
    int jj=0;
    int j=0;

    double fval=0.0;

    gsl_vector *flux;

    // set 'default' shapes
    actbeam->awidth  = -1.0;
    actbeam->bwidth  = -1.0;
    actbeam->aorient = -1.0;

    // set a 'default' flux
    actbeam->flux = NULL;

    /* if flux/wavelength information exists */
    if (sobj->magnitudes)
      {

        // get the number of valid entries in 'sobj->magnitudes'
        nvalid = get_valid_entries(sobj->magnitudes);
        if (!nvalid)
          aXe_message (aXe_M_WARN3, __FILE__, __LINE__,
                       "No valid magnitude for object %i !\n", sobj->number);

        // allocate the flux vector
        flux = gsl_vector_alloc (2*(nvalid));

        // initialize a counter
        jj=0;

        //  go over all magnidute values
        for (j=0; j < (int)sobj->magnitudes->size; j++)
          {
            // check whether the current entry is valid
            if (is_valid_entry(gsl_vector_get(sobj->magnitudes, j)))
              {
                // get the current entry, already converted to flux
                fval = get_flambda_from_magab(gsl_vector_get(sobj->magnitudes, j), gsl_vector_get(sobj->lambdas, j));

                // set the wavelength value and the flux value
                gsl_vector_set(flux, 2*jj, gsl_vector_get(sobj->lambdas, j));
                gsl_vector_set(flux, 2*jj+1, fval);

                // enhance the counter
                jj++;
              }
          }

        // transfer the flux values
        // to the beam
        actbeam->flux = flux;

        // transfer the shape values
        // to the beam
        actbeam->awidth  = sobj->el_image.a;
        actbeam->bwidth  = sobj->el_image.b;
        actbeam->aorient = sobj->el_image.theta;
      }
  }

/**
 * Function: SexObject_to_slitgeom
 * Fill a beam structure with the optimized geometry information.
 * The quantities are computed from the object shape and
 * the trace angle.
 *
 * Parameters:
 * @param sobj        - pointer to a  SexObject
 * @param trace_angle - pointer to a  SexObject
 * @param actbeam     - the beam structure to be filled
 *
 * Returns:
 * @return -
 */
void
SexObject_to_slitgeom(const aperture_conf *conf, const SexObject * sobj,
                      const double trace_angle, beam *actbeam)
  {
    double theta=0.0;
    double orient=0.0;
    double tmp_angle, cos_tmp_angle;
    double A11, A12, A22;

    d_point obj_size;

    // convert the SExtractor angle to rad in the right quadrant
    orient = (180.0 + sobj->el_image.theta) / 180.0 * M_PI;
    while (orient > M_PI)
      orient = orient - M_PI;

    // compute angle between object
    // orientation and trace
    theta = orient - trace_angle;

    // check the object size, possibly setting
    // it to the size of point-like object
    obj_size = check_object_size(conf, sobj, actbeam->ID);

    /*
    // determine the three matrix elements
    A11 = SQR(cos(theta) / sobj->el_image.a) + SQR(sin(theta) / sobj->el_image.b);
    A12 = cos(theta) * sin(theta) * (1.0/(SQR(sobj->el_image.a)) - 1.0/SQR(sobj->el_image.b));
    A22 = SQR(sin(theta) / sobj->el_image.a) + SQR(cos(theta) / sobj->el_image.b);

    // compute a temporary angle
    // and its cosine
    tmp_angle     = atan(A12 / A11);
    cos_tmp_angle = cos(tmp_angle);

    // compute the slit length
    if (cos_tmp_angle > 0.01)
      actbeam->slitgeom[0] = sqrt(A11)*sobj->el_image.a*sobj->el_image.b / cos_tmp_angle;
    else
      actbeam->slitgeom[0] = sqrt(A11)*sobj->el_image.a*sobj->el_image.b / 0.01;

    // compute the slit orientation
    actbeam->slitgeom[1] = tmp_angle + trace_angle + M_PI_2;

    // compute the slit length
    actbeam->slitgeom[2] = 1.0 / sqrt(A11);

    // compute a modified B_IMAGE value, defined to keep the object area constant
    actbeam->slitgeom[3] = sobj->el_image.a*sobj->el_image.b / actbeam->slitgeom[0];
    */

    // determine the three matrix elements
    A11 = SQR(cos(theta) / obj_size.x) + SQR(sin(theta) / obj_size.y);
    A12 = cos(theta) * sin(theta) * (1.0/(SQR(obj_size.x)) - 1.0/SQR(obj_size.y));
    A22 = SQR(sin(theta) / obj_size.x) + SQR(cos(theta) / obj_size.y);

    // compute a temporary angle
    // and its cosine
    tmp_angle     = atan(A12 / A11);
    cos_tmp_angle = cos(tmp_angle);

    // compute the slit length
    if (cos_tmp_angle > 0.01)
      actbeam->slitgeom[0] = sqrt(A11)*obj_size.x*obj_size.y / cos_tmp_angle;
    else
      actbeam->slitgeom[0] = sqrt(A11)*obj_size.x*obj_size.y / 0.01;

    // compute the slit orientation
    actbeam->slitgeom[1] = tmp_angle + trace_angle + M_PI_2;

    // compute the slit width
    actbeam->slitgeom[2] = 1.0 / sqrt(A11);

    // compute a modified B_IMAGE value, defined to keep the object area constant
    actbeam->slitgeom[3] = obj_size.x*obj_size.y / actbeam->slitgeom[0];
  }

/**
 * Function: fill_corner_ignore
 * The method computes and fills the beam corner information
 * and the the ignore flag into a beam structure.
 *
 * Parameters:
 * @param sobj     - pointer to a  SexObject
 * @param obs      - a pointer to the data array containing the image
 * @param conf     - pointer to the configuration structure
 * @param beamID   - the beam ID
 * @param dmag     - number of magnitudes to add to the magnitudes cutoffs
 * @param actbeam  - the beam structure to be filled
 *
 * Returns:
 * @return -
 */
void
fill_corner_ignore(SexObject * sobj, observation * const obs, aperture_conf *conf,
                   int beamID, float dmag, beam *actbeam)
  {
    double mmag_extract;
    double mmag_mark;

    // get the extraction and the mark margnitudes
    mmag_extract = conf->beam[beamID].mmag_extract + dmag;
    mmag_mark    = conf->beam[beamID].mmag_mark    + dmag;

    // make a default for the
    // ignore flag
    actbeam->ignore = 0;

    /* magnitude cut */
    // PROBLEM: there is a logical flaw inside
    //          whhat happend when (mag > mag_mark) && (mag <  mmag_extract) ????
    //          of course this has only relevance when mag_mark < mmag_extract
    //      if ((sobj->mag_auto <= mmag_mark)&&(sobj->mag_auto <= mmag_extract))
    // That's now an easy fix, but makes sense..
    if (sobj->mag_auto <= mmag_extract)
        actbeam->ignore = 0; // this object will be extracted

    if ((sobj->mag_auto <= mmag_mark)&&(sobj->mag_auto > mmag_extract))
        actbeam->ignore = 2; // this object will be not be extracted

    if ((sobj->mag_auto > mmag_mark)&&(sobj->mag_auto > mmag_extract))
        actbeam->ignore = 1; // this object will be ignored

    // fill the boundary box values;
    // returns if the beams is completely out of the image
    if (!fill_object_bbox (obs, actbeam, 2, conf->beam[beamID].offset.dx0,
                           conf->beam[beamID].offset.dx1))
      actbeam->ignore = 1;
  }

/**
 * Function: set_extraction_parameters
 * The function computes and fills the extraction width and the
 * extraction orientation into a beam structure. Depending on the
 * input parameters and the object shape, different methods
 * are deployed.
 *
 * Parameters:
 * @param sobj     - pointer to a  SexObject
 * @param bck_mode - pointer for background mode
 * @param mfwhm    - the fwhm multiplicator constant to apply to
 *                   determine the width of the aperture box for the object.
 * @param dmag     - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set then this task tries to optimize
 *                        the orientation of the
 *                        extraction slit - Should in general be left at 1 so that strange
 *                        geomety is avoided. If set to 2 the extraction is
 *                        forced to be vertical (90 deg.)
 * @param trace_angle - the trace angle
 * @param actbeam     - the beam structure to be filled
 *
 * Returns:
 * @return -
 */
void
set_extraction_parameters(SexObject * sobj, int bck_mode, float mfwhm,
                          int auto_reorient, double trace_angle, beam *actbeam)
  {
    double orient;
    double theta;
    double theta_deg;
    double trace_dist;
    double dya;
    double dyb;

    //int iturn=0;

    // convert the orientation to rad;
    // turn into the right quadrant
    orient = (180.0 + sobj->el_image.theta) / 180.0 * M_PI;
    while (orient > M_PI)
      orient = orient - M_PI;

    // compute the relative between the trace angle
    // and object orientation
    theta = orient - trace_angle;

    // convert theta to deg;
    // convert it into the range:
    // 0.0 < theta_deg < 180.0
    theta_deg = theta / M_PI * 180.;
    while (theta_deg < 0.0)
      theta_deg += 180.0;
    while (theta_deg > 180.0)
      theta_deg -= 180.0;

    // that's the 'normal', slanted extaction
    if (auto_reorient == 0)
      {
        // take the major half axis value * mfwhm
        // as width and its orientation
        // as extraction direction
        actbeam->width  = sobj->el_image.a * mfwhm;
        actbeam->orient = orient;

        // give a warning for small angles
        // between the orientation direction and
        // the trace angle.
        if (fabs(theta_deg )< MIN_DIFFANGLE || fabs(theta_deg-180.0) < MIN_DIFFANGLE)
          aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
                       "aXe_GOL2AF: Object ID: %i : The angle between the extraction orientation and "
                       "the trace is less than %5.2f degrees. You may get severe problems"
                       " down to core dumps in the 1D extraction later on!\n", sobj->number, MIN_DIFFANGLE);
      }
    // thats the slanted, adjusted extraction
    else if (auto_reorient == 1)
      {
        // use the optimized parameters
        actbeam->width  = actbeam->slitgeom[0] * mfwhm;
        actbeam->orient = actbeam->slitgeom[1];
      }
    // thats the perpendicular extraction
    else if(auto_reorient == 2)
      {
        // extraction angle is prependicular to trace
        actbeam->orient = trace_angle + M_PI_2;
        if (mfwhm < 0.0)
          {
            // take the mfwhm value as width;
            actbeam->width  = -1.0 * mfwhm;
          }
        else
          {
            // compute the projections of the major and minor
            // half axis to the extraction direction
            dya = fabs (sobj->el_image.a * sin (theta));
            dyb = fabs (sobj->el_image.b * sin (M_PI_2 - theta));

            // use the larger for extraction width
            actbeam->width = mfwhm*MAX(dya,dyb);
          }
      }

    // this piece of code was specifically added
    // for NICMOS HLA, inorder to nelarge the background
    // area around bright, point-like objects.
    if (sobj->backwindow.x != -1 && bck_mode)
      {
        // compute the projected distance to the trace
        trace_dist = sin(actbeam->orient - trace_angle) * actbeam->width;

        // check whether the distance is too small
        if (trace_dist < sobj->backwindow.x)
          // elongate the width such that the projected distance
          // equals the minimum value given in the SExobject
          actbeam->width = sobj->backwindow.x / sin(actbeam->orient - trace_angle);
      }
  }

d_point
check_object_size(const aperture_conf *conf, const SexObject *sobj, const int beamID)
  {
    double obj_size = 0.0;

    d_point ret;

    // set the input values
    // as default
    ret.x = sobj->el_image.a;
    ret.y = sobj->el_image.b;

    // check whether a minimum
    // size is defined
    if (conf->pobjsize < 0.0)
      // return the value
      return ret;

    // compute the object size
    obj_size = sobj->el_image.a * sobj->el_image.b;

    // check whether the object is point-like
    if (obj_size < SQR(conf->pobjsize))
      {
        // transfer the point-like
        // sizes to the return
        ret.x = conf->pobjsize;
        ret.y = conf->pobjsize;
      }

    // return the resulting
    // object size
    return ret;
  }


/**
 * Function: SexObject_to_beam
 * Fill an beam structure with the information contained in a SexObject
 * structure plus further information derived via the confiugration file.
 *
 * Parameters:
 * @param sobj     - pointer to a  SexObject
 * @param obs      - a pointer to the data array containing the image
 * @param conf     - pointer to the configuration structure
 * @param conffile - the name of the aperture configuration file
 * @param mfwhm    - the fwhm multiplicator constant to apply to
 *                   determine the width of the aperture box for the object.
 * @param dmag     - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set then this task tries to optimize
 *                        the orientation of the
 *                        extraction slit - Should in general be left at 1 so that strange
 *                        geomety is avoided. If set to 2 the extraction is
 *                        forced to be vertical (90 deg.)
 * @param bck_mode - pointer for background mode
 * @param beamID   - the beam ID
 * @param actbeam  - the beam structure to be filled
 *
 * Returns:
 * @return -
 */
void
SexObject_to_beam(SexObject * sobj, observation * const obs, aperture_conf *conf,
                  char conffile[], float mfwhm, float dmag, int auto_reorient,
                  int bck_mode, int beamID, beam *actbeam)
  {
    double trace_angle;

    d_point pixel;

    tracestruct *trace;

    // set the beam ID
    actbeam->ID = conf->beam[beamID].ID;

    // set the model template ID's
    actbeam->modspec  = sobj->modspec;
    actbeam->modimage = sobj->modimage;

    // Adjust for posibly non (0,0) ref point of the 2D field dependence
    // get the geometrical description of the trace at position "pixel"
    pixel.x = sobj->xy_image.x - 1.0 - conf->refx;
    pixel.y = sobj->xy_image.y - 1.0 - conf->refy;
    trace = get_tracestruct_at_pos (conffile, conf->beam[beamID].ID, pixel);

    // transfer trace information to the beam
    actbeam->spec_trace = vector_to_trace_polyN(trace->pol);

    // compute the local trace angle
    trace_angle = atan2(actbeam->spec_trace->deriv (0, actbeam->spec_trace->data),1.0);

    // set the reference point
    actbeam->refpoint.x = sobj->xy_image.x - 1.0 + trace->offset.x;
    actbeam->refpoint.y = sobj->xy_image.y - 1.0 + trace->offset.y;

    // fill the flux and shape information
    SexMags_to_beamFlux(sobj, actbeam);

    // fill in the slit geometry
    SexObject_to_slitgeom(conf, sobj, trace_angle, actbeam);

    // actually set the extraction parameters
    set_extraction_parameters(sobj,bck_mode,mfwhm, auto_reorient, trace_angle, actbeam);

    // set the beam corners and the ignore flag
    fill_corner_ignore(sobj, obs, conf, beamID, dmag, actbeam);
  }

/**
 * Function: SexObject_to_objectII
 * Fill an object structure with the information contained in a SexObject
 * structure This function uses a configuration file to compute the full
 * list of beams for each object.
 *
 * Parameters:
 * @param sobj     - pointer to a  SexObject
 * @param obs      - a pointer to the data array containing the image
 * @param conf     - pointer to the configuration structure
 * @param conffile - the name of the aperture configuration file
 * @param mfwhm    - the fwhm multiplicator constant to apply to
 *                   determine the width of the aperture box for the object.
 * @param dmag     - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set then this task tries to optimize
 *                        the orientation of the
 *                        extraction slit - Should in general be left at 1 so that strange
 *                        geomety is avoided. If set to 2 the extraction is
 *                        forced to be vertical (90 deg.)
 * @param bck_mode - pointer for background mode
 *
 * Returns:
 * @return a pointer to a newly allocated object structure
 */
object * SexObject_to_objectII(SexObject * sobj, observation * const obs,
                               aperture_conf *conf, char conffile[], float mfwhm, float dmag,
                               int auto_reorient, int bck_mode)
  {
    int i=0;

    object *ob;

    // allocate an object
    ob = (object *)malloc(sizeof(object));

    // store the object specific
    // information
    ob->ID = sobj->number;
    ob->nbeams = conf->nbeams;
    ob->grism_obs = obs;

    // go over all beams
    for (i = 0; i < conf->nbeams; i++)
        // fill the current beam
        SexObject_to_beam(sobj, obs, conf, conffile, mfwhm, dmag, auto_reorient,
                          bck_mode, i, &(ob->beams[i]));

    // return the object
    return ob;
  }

/**
 * Function: SexObject_to_object
 * Fill an object structure with the information contained in a SexObject
 * structure This function uses a configuration file to compute the full
 * list of beams for each object.
 *
 * Parameters:
 * @param sobj - pointer to a  SexObject
 * @param obs  - a pointer to the data array containing the image
 * @param conffile - the name of the aperture configuration file
 * @param maxmag  - upper magniture bound. If object has a magnitude
 *                  greater than this has its ignore flag set to 1,
 *                  zero otherwise.
 * @param mfwhm - the fwhm multiplicator constant to apply to
 *                determine the width of the aperture box for the object.
 * @param dmag - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set then this task tries to optimize
 *                        the orientation of the
 *  extraction slit - Should in general be left at 1 so that strange
 *                    geomety is avoided. If set to 2 the extraction is
 *                    forced to be vertical (90 deg.)
 *
 * Returns:
 * @return a pointer to a newly allocated object structure
 *
 */
// object *
// SexObject_to_object (SexObject * sobj, observation * const obs,
//                      aperture_conf *conf, char conffile[], float mfwhm,
//                      float dmag,
//                      //              char conffile[], float mfwhm, float dmag,
//                      int auto_reorient, int bck_mode)
// {
//   int i, dx0, dx1;
//   beam *b;
//   object *ob = malloc (sizeof (object));
//   d_point pixel;
//   tracestruct *trace;
//   //  aperture_conf *conf;
//   float mmag_extract, mmag_mark;
//   //  conf = get_aperture_descriptor (conffile);
//   int j = 0, jj=0, nvalid=0;
//   double fval=0.0;
//   gsl_vector *flux;
//   float aposang, paposang;
//   float dya, dyb;
//   int iturn;

//   ob->ID = sobj->number;


//   for (i = 0; i < conf->nbeams; i++)
//     {
//       dx0 = conf->beam[i].offset.dx0;
//       dx1 = conf->beam[i].offset.dx1;
//       mmag_extract = conf->beam[i].mmag_extract+dmag;
//       mmag_mark = conf->beam[i].mmag_mark+dmag;
//       b = &(ob->beams[i]);
//       b->ID = conf->beam[i].ID;

//       //---------------------------------------
//       // some code for FORS2 MXU
//       // b->backwindow.x = sobj->backwindow.x;
//       // b->backwindow.y = sobj->backwindow.y;

//       b->modspec  = sobj->modspec;
//       b->modimage = sobj->modimage;

//       /* Adjust for posibly non (0,0) ref point of the 2D field dependence */
//       pixel.x = sobj->xy_image.x-1 - conf->refx;
//       pixel.y = sobj->xy_image.y-1 - conf->refy;
//       /* get the geometrical description of the trace at position "pixel" */
//       trace =
//         get_tracestruct_at_pos (conffile, conf->beam[i].ID, pixel);

//       b->spec_trace = vector_to_trace_polyN ( trace->pol );

//       b->width   = sobj->el_image.a * mfwhm;


//       /* if flux/wavelength information exists */
//       if (sobj->magnitudes)
//         {

//           // get the number of valid entries in 'sobj->magnitudes'
//           nvalid = get_valid_entries(sobj->magnitudes);
//           if (!nvalid)
//             aXe_message (aXe_M_WARN3, __FILE__, __LINE__,
//                          "No valid magnitude for object %i !\n", sobj->number);

//           // allocate the vector for the flux
//           flux = gsl_vector_alloc (2*(nvalid));

//           // fill the flux vector with values
//           jj=0;
//           for (j=0; j < (int)sobj->magnitudes->size; j++)
//             {
//               if (is_valid_entry(gsl_vector_get(sobj->magnitudes, j)))
//                 {

//                   fval = get_flambda_from_magab(gsl_vector_get(sobj->magnitudes, j), gsl_vector_get(sobj->lambdas, j));

//                   gsl_vector_set(flux, 2*jj, gsl_vector_get(sobj->lambdas, j));
//                   //              gsl_vector_set(flux, 2*jj+1, gsl_vector_get(sobj->magnitudes, j));
//                   gsl_vector_set(flux, 2*jj+1, fval);
//                   jj++;
//                 }
//             }
//           b->flux = flux;

//           /* fill the structual parameters of the object */
//           b->awidth  = sobj->el_image.a;
//           b->bwidth  = sobj->el_image.b;
//           b->aorient = sobj->el_image.theta;
//         }
//       // if flux/wavelength information does not exists
//       else
//         {
//           // set everything to dummy values
//           b->flux = NULL;

//           b->awidth  = -1.0;
//           b->bwidth  = -1.0;
//           b->aorient = -1.0;
//         }

//       /* Convert from SeXtractor angle reference point to aXe's */
//       b->orient = (180 + sobj->el_image.theta) / 180. * M_PI;
//       while (b->orient > M_PI)
//         b->orient = b->orient - M_PI;

//       // calculate the angle between
//       // the extraction direction and the trace
//       aposang =
//         b->orient - atan2(b->spec_trace->deriv (0, b->spec_trace->data),1.0);

//       // transform the angle to degrees
//       paposang = aposang / M_PI * 180.;

//       // give a warning for small angles
//       // between the orientation direction and
//       // the trace angle.
//       // BUGFIX comparison should be to paposang not the logical comparison MLS
//       if ( ( (fabs(paposang) < MIN_DIFFANGLE) || (fabs(paposang-180.0) < MIN_DIFFANGLE) ) && (auto_reorient == 0) )
//         aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
//                      "aXe_GOL2AF: Object ID: %i: The angle between the extraction orientation and "
//                      "the trace is less than %5.2f degrees. You may get severe problems"
//                      " down to core dumps in the 1D extraction later on!\n", ob->ID, MIN_DIFFANGLE);


//        Case when the extraction angle is modified to that the extraction 
//       /* proceeds along the semi-axis which projects farther away from the
//        * trace i.e. we avoid trying to extract spectra in a direction nearly
//        * parallel the trace */
//       if (auto_reorient==1)
//         {

//           aposang =
//             b->orient - atan2(b->spec_trace->deriv (0,
//                                                     b->spec_trace->data),1.0);

//           /* February 2004 introduced to make a hardstop in the range for
//              allowed angles. */
//           if (aposang > 0.0){
//             paposang = aposang / M_PI * 180.;
//           }
//           else{
//             paposang = aposang / M_PI * 180. + 360.0;
//           }
//           if (paposang > 180.0)
//             paposang = paposang-180.0;
//           /* the hard stop is here */
//           if (paposang > 30.0 && paposang < 150.0){
//             iturn = 0;}
//           else{
//             iturn = 1;}
//           /*end introduction */


//           /* Compute how far each axes extends away from the trace */
//           dya = fabs (sobj->el_image.a * sin (aposang));
//           dyb = fabs (sobj->el_image.b * sin (M_PI / 2. - aposang));

//           /* Select the broader of the two axes */
//           if (dya > dyb && iturn == 0)
//             {
//               b->width = sobj->el_image.a * mfwhm;
//             }
//           else
//             {
//               b->width = sobj->el_image.b * mfwhm;
//               b->orient = b->orient + M_PI / 2.;
//             }
//         }

//         /* Case when we force the extraction to be vertical,
//            the extraction width is recomputed */
//       if (auto_reorient==2) {
//         // new system for fixed extraction:
//         //    - extraction direction perpendicular
//         //      to the trace
//         //    - the mfwhm is used as a fixed extraction widt in pixels
//         b->orient =
//           atan2(b->spec_trace->deriv (0,b->spec_trace->data),1.0)+ M_PI / 2.0;

//         if (mfwhm < 0.0)
//           {
//             b->orient =
//               atan2(b->spec_trace->deriv (0.0,b->spec_trace->data),1.0)+ M_PI / 2.0;
//             b->width  = -mfwhm;

//           }
//         else
//           {
//             aposang = b->orient
//               - atan2(b->spec_trace->deriv (0,b->spec_trace->data),1.0);

//             /* Compute how far each axes extends away from the trace */
//             dya = fabs (sobj->el_image.a * sin (aposang));
//             dyb = fabs (sobj->el_image.b * sin (M_PI / 2. - aposang));
//             b->width = mfwhm*MAX(dya,dyb);
//             b->orient =
//               atan2(b->spec_trace->deriv (0,b->spec_trace->data),1.0)+ M_PI / 2.0;
//           }

//       }

//       if (sobj->backwindow.x != -1 && bck_mode)
//         b->width = MAX(b->width, sobj->backwindow.x);

//       /* On coordinate systems:
//          the sexrtractor image positions are given in the iraf system,
//          which means the value of the lower left pixel is associated
//          with the coordinate (1.0,1.0).
//          aXe works in a kind of 'matrix system', where the value of
//          the lower left pixel is stored in the matrix indices (0,0)
//          or (0.0,0.0) seen as a coordinate system.
//          To transform into this system 1.0 is subracted from
//          both sextractor coordinates.
//        */
//       b->refpoint.x = sobj->xy_image.x-1.0 + trace->offset.x;
//       b->refpoint.y = sobj->xy_image.y-1.0 + trace->offset.y;

//       /* magnitude cut */
//       // PROBLEM: there is a logical flaw inside
//       //          whhat happend when (mag > mag_mark) && (mag <  mmag_extract) ????
//       //          of course this has only relevance when mag_mark < mmag_extract
//       //      if ((sobj->mag_auto <= mmag_mark)&&(sobj->mag_auto <= mmag_extract))
//       // That's now an easy fix, but makes sense..
//       if (sobj->mag_auto <= mmag_extract)
//         {
//           b->ignore = 0; /* This object will be extracted */
//         }
//       if ((sobj->mag_auto <= mmag_mark)&&(sobj->mag_auto > mmag_extract))
//         {
//           b->ignore = 2; /* This object will be not be extracted */
//         }
//       if ((sobj->mag_auto > mmag_mark)&&(sobj->mag_auto > mmag_extract))
//         {
//           b->ignore = 1; /* This object will be ignored */
//         }

//       if (!fill_object_bbox (obs, b, 2, dx0, dx1))
//         b->ignore = 1;

//     b->slitgeom[0] = -1.0;
//     b->slitgeom[1] = -1.0;
//     b->slitgeom[2] = -1.0;
//     b->slitgeom[3] = -1.0;
//     }

//   ob->nbeams = conf->nbeams;
//   ob->grism_obs = obs;

//   return ob;
// }

/**
 * Function: SexObjects_to_oblist
 * Produces an object list from the data contained in an array of SexObjects
 *
 * Parameters:
 * @param sobjs - an NULL terminated array containing pointers to SexObjects
 * @param obs - a pointer to the data array containing the image
 * @param conffile - the name of the aperture configuration file
 * @param mmag_extract - upper magniture bound. Any object with a magnitude greater than
 *  this has its ignore flag set to 1.
 * @param mmag_mark - upper magniture bound. Any object with a magnitude greater than
 *  this has its ignore flag set to 2.
 * @param mfwhm - the fwhm multiplicator constant to apply to determine the width of the
 *  aperture box for the object.
 * @param dmag - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set to 1 then this task tries to optimize the orientation of the
 *  extraction slit. Should in general be left at 1 so that strange geomety is avoided. If set to 2
 *  the extraction is forced to be vertical (90 deg.)
 *
 *  Return:
 *  @return a pointer to a NULL terminated object array.
 */
object **
SexObjects_to_oblist (SexObject ** sobjs, observation * const obs,
                      aperture_conf *conf, char conffile[], float mfwhm,
                      float dmag,
                      //                      char conffile[], float mfwhm, float dmag,
                      int auto_reorient, int bck_mode)
{
     int i, nobjs = 0;
     object **oblist;

     /* Find the number of SexObjects in sobjs */
     while (sobjs[nobjs])
          nobjs++;
     /* Allocate enough room for a new object list */
     oblist = (object **) malloc ((nobjs + 1) * sizeof (object *));

     for (i = 0; i < nobjs; i++)
       {
         //fprintf(stdout, "Using the old routine...\n");
         oblist[i] =
         SexObject_to_objectII(sobjs[i], obs, conf, conffile, mfwhm, dmag, auto_reorient, bck_mode);
        }
     oblist[nobjs] = NULL;

     return oblist;
}

/**
 * Function: SexObjects_to_oblistII
 * Produces an object list from the data contained in an array of SexObjects
 *
 * Parameters:
 * @param sobjs - an NULL terminated array containing pointers to SexObjects
 * @param obs - a pointer to the data array containing the image
 * @param conffile - the name of the aperture configuration file
 * @param mmag_extract - upper magniture bound. Any object with a magnitude greater than
 *  this has its ignore flag set to 1.
 * @param mmag_mark - upper magniture bound. Any object with a magnitude greater than
 *  this has its ignore flag set to 2.
 * @param mfwhm - the fwhm multiplicator constant to apply to determine the width of the
 *  aperture box for the object.
 * @param dmag - number of magnitudes to add to the magnitudes cutoffs
 * @param auto_reorient - if set to 1 then this task tries to optimize the orientation of the
 *  extraction slit. Should in general be left at 1 so that strange geomety is avoided. If set to 2
 *  the extraction is forced to be vertical (90 deg.)
 *
 *  Return:
 *  @return a pointer to a NULL terminated object array.
 */
object **
SexObjects_to_oblistII (SexObject ** sobjs, observation * const obs,
                      aperture_conf *conf, char conffile[], float mfwhm,
                      float dmag,
                      int auto_reorient, int bck_mode)
{
     int i, nobjs = 0;
     object **oblist;
     fflush(stdout);
     /* Find the number of SexObjects in sobjs */
     while (sobjs[nobjs])
          nobjs++;
     /* Allocate enough room for a new object list */
     oblist = (object **) malloc ((nobjs + 1) * sizeof (object *));

     for (i = 0; i < nobjs; i++)
       {
         //fprintf(stdout, "Using the new routine...\n");
         oblist[i] =
         SexObject_to_objectII(sobjs[i], obs, conf, conffile, mfwhm, dmag, auto_reorient, bck_mode);
        }
     oblist[nobjs] = NULL;
     return oblist;
}

/**
 * Function: check_conf_for_slitlessgeom
 * The functions checks whether a smoothed flux conversion
 * is possible or not. In case that keywords in the configuration
 * files are missing, it is NOT possible, and 0 is returned.
 *
 * Parameters:
 * @param conf          - the configuarion file structure
 * @param auto_reorient - integer indicating the extraction method
 *
 * Returns
 * @return is_possible - the paramters used in the smoothing
 */
int
check_conf_for_slitlessgeom(const aperture_conf *conf, const int auto_reorient)
//check_conf_for_slitlessgeom(const aperture_conf *conf, const int slitless_geom)
  {
    int is_possible=1;

    if (auto_reorient == 1 && conf->pobjsize < 0.0)
      // change the switch
      is_possible = 0;

    // return the pointer
    return is_possible;
  }


/**
 * Function: fill_object_bbox
 * Fill up the aperture part of the beam of an object structure by looking at the trace polynomial
 * description associated with this object.
 *
 * Parameters:
 * @param obs - a pointer to the data array containing the image (used for bound checking)
 * @param b - A pointer to a beam whose aperture must be filled.
 * @param m_width - The width in pixel of the extraction box
 * @param dxmin - How far (in pixel) to follow the trace on the left side of the spectrum
 * @param dxmax - How far (in pixel) to follow the trace on the right hand side of the spectrum
 *
 * Returns:
 * @return 1 if bounding box appears to be valied (i.e. non-empty)
 */
int
fill_object_bbox (observation * const obs, beam * b, const float m_width,
                  const int dxmin, const int dxmax)
{
     d_point pmin, pmax;
     d_point pminl, pminh, pmaxl, pmaxh;
     float w = b->width/2. * m_width + 2.0;
     float wcos, wsin;
     float area;

     int xmin, xmax, xact;
     double yact;

     pmin.x = b->refpoint.x + dxmin;
     pmin.y =
          b->refpoint.y + b->spec_trace->func (dxmin, b->spec_trace->data);

     pmax.x = b->refpoint.x + dxmax;
     pmax.y =
          b->refpoint.y + b->spec_trace->func (dxmax, b->spec_trace->data);


     if (b->spec_trace->type > 1)
       {
         if (dxmin < dxmax)
           {
             xmin = dxmin;
             xmax = dxmax;
           }
         else
           {
             xmin = dxmax;
             xmax = dxmin;
           }

         for (xact = xmin; xact < xmax; xact++)
           {
             yact = b->refpoint.y + b->spec_trace->func (xact, b->spec_trace->data);

             if (yact < pmin.y)
               pmin.y = yact;
             if (yact > pmax.y)
               pmax.y = yact;
           }
         wcos = w * cos (b->orient);
         wsin = w * sin (b->orient);

         if (wcos > 0)
           {
             pmin.x = pmin.x - wcos;
             pmax.x = pmax.x + wcos;
           }
         else
           {
             pmin.x = pmin.x + wcos;
             pmax.x = pmax.x - wcos;
           }
         if (wsin > 0)
           {
             pmin.y = pmin.y - wsin;
             pmax.y = pmax.y + wsin;
           }
         else
           {
             pmin.y = pmin.y + wsin;
             pmax.y = pmax.y - wsin;
           }
         b->corners[0].x = pmin.x;
         b->corners[0].y = pmin.y;

         b->corners[1].x = pmax.x;
         b->corners[1].y = pmin.y;

         b->corners[2].x = pmax.x;
         b->corners[2].y = pmax.y;

         b->corners[3].x = pmin.x;
         b->corners[3].y = pmax.y;
       }
     else
       {

         wcos = w * cos (b->orient);
         wsin = w * sin (b->orient);

         pminl.x = pmin.x - wcos;
         pminl.y = pmin.y - wsin;

         pminh.x = pmin.x + wcos;
         pminh.y = pmin.y + wsin;

         pmaxl.x = pmax.x - wcos;
         pmaxl.y = pmax.y - wsin;

         pmaxh.x = pmax.x + wcos;
         pmaxh.y = pmax.y + wsin;

         b->corners[0].x = pminl.x;
         b->corners[0].y = pminl.y;

         b->corners[1].x = pmaxl.x;
         b->corners[1].y = pmaxl.y;

         b->corners[2].x = pmaxh.x;
         b->corners[2].y = pmaxh.y;

         b->corners[3].x = pminh.x;
         b->corners[3].y = pminh.y;
       }

     /* Compute the area of this aperture */
     area =
          0.5 * (abs (b->corners[0].y - b->corners[3].y) *
                 abs (b->corners[0].x - b->corners[1].x) +
                 abs (b->corners[2].y -
                       b->corners[1].y) * abs (b->corners[3].x -
                                               b->corners[2].x));
     if (area == 0)
       {
            fprintf (stderr, "aper debug: area is zero!\n");
            return 0;
       }
     return 1;
}

/**
 * Function: size_of_sextractor_catalog
 * A utility function which parses a Sextractor catalog file,
 * and returns the number of valid catalog entries found in the file
 *
 * Parameters:
 * @param filename - a pointer pointing to a char array containing the
 *                   list of a sextractor object output catalog.
 *                   Ignores rows starting with a ";"
 */
int
size_of_sextractor_catalog (char filename[])
{
  FILE *input;
  char Buffer[CATBUFFERSIZE];
  gsl_vector *v;
  int catsize;
  int num = 0;
  colinfo * actcatinfo;
  actcatinfo = get_sex_col_descr (filename);
  //  catalog_header = get_sex_col_descr (filename);
  catsize = actcatinfo->numcols;
  if (!(input = fopen (filename, "r")))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Could not open Sextractor catalog" "file %s,\n",
                   filename);
    }


  while (fgets (Buffer, CATBUFFERSIZE, input))
    {
      if (Buffer[0] == ';')
        continue;
      lv1ws (Buffer);
      v = string_to_gsl_array (Buffer);
      if (v==NULL) continue;
      if ((int)v->size == catsize)
        num++;
    }
  return num;
}


/**
 * Function: get_SexObject_from_catalog
 * Parses a Sextractor 2.0 catalog file, and outputs a NULL terminated
 * array of SexObjects pointers. Ignores rows starting with a ;
 *
 * Parameters:
 * @param filename a pointer pointing to a char array containing the
 *                 list of a sextractor object output catalog
 *
 * Returns:
 * @return a pointer to an array of SeXObject pointer. NULL terminated
 */
SexObject **
get_SexObject_from_catalog (char filename[], const double lambda_mark)
{
  FILE *input;
  char Buffer[CATBUFFERSIZE];
  gsl_vector *v;
  gsl_vector *waves;
  gsl_vector *cnums;
  size_t hasmags=0;
  size_t nobjs;
  size_t catsize;
  int i;
  size_t magcencol=0;
  SexObject **sobjs, *sobj;
  colinfo * actcatinfo;
  px_point  backwin_cols;
  px_point  modinfo_cols;

  actcatinfo = get_sex_col_descr (filename);
  hasmags = has_magnitudes(actcatinfo);
  waves = gsl_vector_alloc (hasmags);
  cnums = gsl_vector_alloc (hasmags);
  hasmags = get_magcols(actcatinfo, waves, cnums);
  magcencol = get_magauto_col(waves, cnums, lambda_mark);
  backwin_cols = has_backwindow(actcatinfo);
  modinfo_cols = has_modelinfo(actcatinfo);


  //  catsize = count_keys (catalog_header);
  catsize = actcatinfo->numcols;
  nobjs = size_of_sextractor_catalog (filename);
  if (!(input = fopen (filename, "r")))
       {
         aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                      "Could not open Sextractor catalog" "file %s,\n",
                      filename);
       }

  /* Allocate enough room for nobjs+1 SexObject  pointers */
  sobjs = (SexObject **) malloc ((nobjs + 1) * sizeof (SexObject *));
  if (!sobjs)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory. Couldn't allocate Sextractor Object");
    }

  i = 0;
  while (fgets (Buffer, CATBUFFERSIZE, input))
    {
      if (Buffer[0] == ';')
        continue;
      v = string_to_gsl_array (Buffer);
      if (v==NULL) continue;
      if (v->size == catsize)
        {
          // BIG BUG!!!
          //               sobj = create_SexObject (catalog_header, Buffer);
          sobj = create_SexObject (actcatinfo, Buffer, waves, cnums, backwin_cols, modinfo_cols, magcencol, 0);
          sobjs[i++] = sobj;
        }
    }
  sobjs[i] = NULL;
  return sobjs;
}

/**
 * Function: el_to_ABC_world
 * Convert a wold coordinate ellipse into a set of 3 points on the
 * tangent plane point a is at pos of obj, semi major-axis end is point b,
 * semi minor axis is point c theta is measured from x-axis to semi-major
 * axis, counter clock wise.
 */
void
el_to_ABC_world (ellipse *el, sky_coord *a, sky_coord *b, sky_coord *c)
{

     a->ra = 0.0;
     a->dec = 0.0;

     b->ra = el->a * cos (el->theta / 180. * M_PI);
     b->dec = el->a * sin (el->theta / 180. * M_PI);

     c->ra = -1.0 * el->b * sin (el->theta / 180. * M_PI);
     c->dec = el->b * cos (el->theta / 180. * M_PI);
}

/*
 * Function: el_to_ABC_world2
 *
 */
void
el_to_ABC_world2 (ellipse *el, sky_coord *a, sky_coord *b, sky_coord *c)
{
  b->ra  = a->ra  - (el->a * cos (el->theta / 180. * M_PI))/cos(a->dec / 180. * M_PI);
  b->dec = a->dec + (el->a * sin (el->theta / 180. * M_PI));

  c->ra  = a->ra  - (el->b * cos ((el->theta-90.0)/180.*M_PI))/cos(a->dec / 180. * M_PI);
  c->dec = a->dec + (el->b * sin ((el->theta-90.0)/180.*M_PI));
}

/**
 * Function: ABC_image_to_el
 * Convert a set of three image coordinates into an elipse structure
 *
 */
void
ABC_image_to_el (d_point *a, d_point *b, d_point *c, ellipse *el)
{
  if (b->x != 0.0)
    {
      el->theta = atan (b->y / b->x) / M_PI * 180.0;
    }
  else
    {
      el->theta = 90.0;
    }

  el->a = sqrt ((b->x - a->x) * (b->x - a->x) + (b->y - a->y) * (b->y - a->y));
  el->b = sqrt ((c->x - a->x) * (c->x - a->x) + (c->y - a->y) * (c->y - a->y));
}

/**
 * Function: ABC_image_to_el2
 * Convert a set of three image coordinates into an elipse structure
 *
 */
void
ABC_image_to_el2 (d_point *a, d_point *b, d_point *c, ellipse *el)
{
  el->theta = atan ((b->y - a->y )/ (b->x - a->x)) / M_PI * 180.0;

  el->a = sqrt ((b->x - a->x) * (b->x - a->x) + (b->y - a->y) * (b->y - a->y));
  el->b = sqrt ((c->x - a->x) * (c->x - a->x) + (c->y - a->y) * (c->y - a->y));
}

/**
 *  Function: fill_missing_WCS_coordinates
 *  Update a SexObject structure and replace all missing WCS coordinates
 *  by re-computing them using the world one and the from_wcs wcs.
 *  Re-computes:
 *      peak_world_x
 *     peak_world.y
 *       xy_world.x
 *      xy_world.y
 *      el_world.a
 *      el_world.b
 *      el_world.theta
 *
 *  @param o - a pointer to an exsting SexObject. Should have been checked for validity already
 *  @param from_wcs - a pointer to an existing  WoldCoord wcs structure
 */
void
fill_missing_WCS_coordinates (SexObject * o, struct WorldCoor *from_wcs, int overwrite)
{
  sky_coord as, bs, cs;
  d_point a, b, c;
  //int offscl;

  /***************************************************************/
  /* If the WCS coordinates are missing, generate them using     */
  /* from_wcs and sextractor image coordinates                   */
  /***************************************************************/

  if ((overwrite) || (isnan (o->xy_world.ra)) || (isnan (o->xy_world.dec)))
    {
      /* recomputes xy_image  from the world values */
      /*pix2wcs (from_wcs, o->xy_image.x, o->xy_image.y,
               &(o->xy_world.ra),&(o->xy_world.dec), &offscl);
      */
      /* recomputes xy_image  from the world values */
      // new version requested since wcstools-3.6
      pix2wcs (from_wcs, o->xy_image.x, o->xy_image.y,
               &(o->xy_world.ra),&(o->xy_world.dec));
    }

  if ((isnan (o->el_world.a)) || (isnan (o->el_world.b))
      || (isnan (o->el_world.theta)) || (overwrite))
    {

      /* recomputes el_image them from the world values */

      /* Note the following is an unconventional use of the et_to*()
         routines since we store image coords in as, bs, and cs */
      el_to_ABC_world (&(o->el_image), &as, &bs, &cs);  /* Convert ellipse into a set of 3 coordinates */

      /* Compute the end point of the axes in the image coord system */
      /*      pix2wcs (from_wcs, as.ra, as.dec, &(a.x), &(a.y), &offscl);
      pix2wcs (from_wcs, bs.ra, bs.dec, &(b.x), &(b.y), &offscl);
      pix2wcs (from_wcs, cs.ra, cs.dec, &(c.x), &(c.y), &offscl);
      */
      /* Compute the end point of the axes in the image coord system */
      // new version requested since wcstools-3.6
      pix2wcs (from_wcs, as.ra, as.dec, &(a.x), &(a.y));
      pix2wcs (from_wcs, bs.ra, bs.dec, &(b.x), &(b.y));
      pix2wcs (from_wcs, cs.ra, cs.dec, &(c.x), &(c.y));

      ABC_image_to_el (&a, &b, &c, &(o->el_world));     /* convert 3 coordinates into ellipse */
    }
}

/**
 * Function: fill_missing_image_coordinates
 *  Update a SexObject structure and replace all missing image coordinates
 *  by re-computing them using the world one and the from_wcs wcs.
 *  Re-computes:
 *      peak_image_x
 *      peak_image.y
 *      xy_image.x
 *      xy_image.y
 *      el_image.a
 *      el_image.b
 *      el_image.theta
 *
 * Parameters:
 * @param o - a pointer to an exsting SexObject. Should have been
 *            checked for validity already
 * @param from_wcs - a pointer to an existing  WoldCoord wcs structure
 */
void
fill_missing_image_coordinates (SexObject *o, struct WorldCoor *from_wcs, int overwrite)
{
  sky_coord as, bs, cs;
  d_point a, b, c;
  int offscl;


  /************************/
  /* If the image coordinates are missing, generate them using from_wcs and sextractor world
     coordinates */
  /************************/
  if ((overwrite) || (isnan (o->xy_image.x)) || (isnan (o->xy_image.y)))
    {
      /* recomputes xy_image  from the world values */
      wcs2pix (from_wcs, o->xy_world.ra, o->xy_world.dec,
               &(o->xy_image.x), &(o->xy_image.y), &offscl);

    }


  if ((isnan (o->el_image.a)) || (isnan (o->el_image.b))
      || (isnan (o->el_image.theta)) || (overwrite))
    {

      /* recomputes el_image them from the world values */

      el_to_ABC_world (&(o->el_world), &as, &bs, &cs);  /* Convert ellipse into a set of 3 coordinates */

      /* Compute the end point of the axes in the image coord system */
      wcs2pix (from_wcs, as.ra, as.dec, &(a.x), &(a.y), &offscl);
      wcs2pix (from_wcs, bs.ra, bs.dec, &(b.x), &(b.y), &offscl);
      wcs2pix (from_wcs, cs.ra, cs.dec, &(c.x), &(c.y), &offscl);

      /* DISABLED !!*/
      //            ABC_image_to_el (&a, &b, &c, &(o->el_image));       /* convert 3 coordinates into ellipse */
    }
}

/**
 * Function: fill_all_missing_image_coordinates
 * Replaces all NaN values in a set of SexObjects with ones computed
 * from the available information and the associated WCS.
 *
 * Parameters:
 * @param sobjs    - an NULL terminated array containing pointers to SexObjects
 * @param from_wcs - a pointer to an existing  WoldCoord wcs structure
 * @param overwrite - forces the re-computation of all image coordinates
 *
 */
void
fill_all_missing_image_coordinates (SexObject ** sobjs,
                                    struct WorldCoor *from_wcs, int overwrite)
{
  int i, nobjs = 0;

  /* Find the number of SexObjects in sobjs */
  while (sobjs[nobjs])
    nobjs++;

  for (i = 0; i < nobjs; i++)
    {
      fill_missing_image_coordinates (sobjs[i], from_wcs, overwrite);
    }

}

/**
 * Function: fill_all_missing_WCS_coordinates
 *  Replaces all NaN values in a set of SexObjects with ones computed from the
 *  available information and the associated WCS.
 *
 * Parameters:
 *  @param sobjs - an NULL terminated array containing pointers to SexObjects
 *  @param from_wcs - a pointer to an existing  WoldCoord wcs structure
 *  @param overwrite - forces the re-computation of all image coordinates
 *
 */
void
fill_all_missing_WCS_coordinates (SexObject ** sobjs,
                                    struct WorldCoor *from_wcs, int overwrite)
{
  int i, nobjs = 0;

  /* Find the number of SexObjects in sobjs */
  while (sobjs[nobjs])
    nobjs++;

  for (i = 0; i < nobjs; i++)
    {
      fill_missing_WCS_coordinates (sobjs[i], from_wcs, overwrite);
    }

}

/**
 * Function: compute_new_image_coordinates
 *  Function that uses the exisiting world coordinates of a SexObject to
 *  compute new pixel coordinates using a (new) WCS.
 *
 * Parameters:
 *  @param o - a pointer to an exsting SexObject. Should have been checked
 *   for validity already
 *  @param to_wcs - a pointer to an existing  WoldCoord wcs structure
 *
 */
void
compute_new_image_coordinates (SexObject * o, struct WorldCoor *to_wcs)
{
  sky_coord as, bs, cs;
  d_point a, b, c;
  int offscl;


  /******************************************************************/
  /* uses the existing world coordinates of a SexObject to compute  */
  /* new pixel coordinates using a (new) WCS.                       */
  /*******************************************************************/

  /* recomputes xy_image  from the world values */
  wcs2pix (to_wcs, o->xy_world.ra, o->xy_world.dec, &(o->xy_image.x),
           &(o->xy_image.y), &offscl);


  /* recomputes el_image them from the world values */
  el_to_ABC_world (&(o->el_world), &as, &bs, &cs);      /* Convert ellipse into a set of 3 coordinates */

  /* Compute the end point of the axes in the image coord system */
  wcs2pix (to_wcs, as.ra, as.dec, &(a.x), &(a.y), &offscl);
  wcs2pix (to_wcs, bs.ra, bs.dec, &(b.x), &(b.y), &offscl);
  wcs2pix (to_wcs, cs.ra, cs.dec, &(c.x), &(c.y), &offscl);


  /* DISABLED */
  //     ABC_image_to_el (&a, &b, &c, &(o->el_image));  /* convert 3 coordinates into ellipse */
}

/**
 * Function: compute_new_image_sexobject
 * Function that uses the exisiting world coordinates of a SexObject to
 * compute new pixel coordinates using a (new) WCS.
 *
 * Parameters:
 * @param o - a pointer to an exsting SexObject. Should have been checked
 *  for validity already
 *  @param to_wcs - a pointer to an existing  WoldCoord wcs structure
 */
void
compute_new_image_sexobject (SexObject * o, struct WorldCoor *to_wcs, int th_sky)
{
  sky_coord  bs, cs;
  d_point a, b, c;
  int offscl;

  /******************************************************************/
  /* uses the existing world coordinates of a SexObject to compute  */
  /* new pixel coordinates using a (new) WCS.                       */
  /******************************************************************/

  /* recomputes xy_image  from the world values */
  wcs2pix (to_wcs, o->xy_world.ra, o->xy_world.dec, &(o->xy_image.x),
           &(o->xy_image.y), &offscl);

  /* Convert ellipse into a set of 3 coordinates */
  el_to_ABC_world2 (&(o->el_world), &(o->xy_world), &bs, &cs);

  /* Compute the end point of the axes in the image coord system */

  wcs2pix (to_wcs, o->xy_world.ra, o->xy_world.dec, &(a.x), &(a.y), &offscl);
  wcs2pix (to_wcs, bs.ra, bs.dec, &(b.x), &(b.y), &offscl);
  wcs2pix (to_wcs, cs.ra, cs.dec, &(c.x), &(c.y), &offscl);


  /* convert 3 coordinates into ellipse */
  ABC_image_to_el2 (&a, &b, &c, &(o->el_image));
  if (th_sky)
    o->el_image.theta = o->el_image.theta - 90.0;
}

/**
 * Function: compute_all_new_image_coordinates
 * Replaces all SexObjects image coordinates with new ones computed using  the
 * the passed WCS.
 *
 * Parameters:
 * @param sobjs - a NULL terminated array containing pointers to SexObjects
 * @param to_wcs - a pointer to an existing  WoldCoord wcs structure
 */
void
compute_all_new_image_coordinates (SexObject ** sobjs,
                                   struct WorldCoor *to_wcs)
{
     int i, nobjs = 0;

     /* Find the number of SexObjects in sobjs */
     while (sobjs[nobjs])
          nobjs++;

     for (i = 0; i < nobjs; i++)
       {
            compute_new_image_coordinates (sobjs[i], to_wcs);
       }

}


/**
 * Function: free_SexObjects
 *    Free a NULL terminated array of SexObject
 *
 * Parameters:
 *  @param sobjs - a NULL terminated array containing pointers to SexObjects
 *
 */
void
free_SexObjects (SexObject ** sobjs)
{
  int i, nobjs = 0;

  /* Find the number of SexObjects in sobjs */
  while (sobjs[nobjs])
    nobjs++;

  for (i = 0; i < nobjs; i++)
    {
      if (sobjs[i]->magnitudes){
        gsl_vector_free (sobjs[i]->lambdas);
        gsl_vector_free (sobjs[i]->magnitudes);
      }
      free (sobjs[i]);
      sobjs[i] = NULL;
    }
  free (sobjs);
  sobjs = NULL;
}
