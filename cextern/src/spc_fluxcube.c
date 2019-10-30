/**
 * Subroutines to calculate the
 * various contamination models
 *
 */

#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

#include "fitsio.h"
#include "inout_aper.h"
#include "aXe_grism.h"
#include "spce_PET.h"
#include "spc_wl_calib.h"
#include "aXe_errors.h"
#include "fringe_conf.h"
#include "spc_resp.h"
#include "spce_pathlength.h"
#include "aper_conf.h"
#include "trace_conf.h"
#include "specmodel_utils.h"
#include "model_utils.h"
#include "spc_fluxcube.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define SQR(x) ((x)*(x))

/**
 * The function creates and loads a new flux_cube structure
 * The data is taken from the file specified in the input.
 *
 * @param  fcube_file - the file name of the fluxcube
 *
 * @return ret        - the new flux_cube structure
 */
flux_cube *
load_fluxcube(const char fcube_file[])
{
  flux_cube *fcube;

  //gsl_matrix *test;

  int nflux, n_ext=0;
  int i;


  // get the number of extensions
  n_ext = FITSextnum(fcube_file);
  if (n_ext <3)
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "aXe_PETCONT: " "the fluxcube file %s has only %i extensions!\n",
                fcube_file, n_ext);

  // determine the number of flux images and allocate
  // the space for the fluxcube structure
  nflux = n_ext-2;
  fcube = alloc_fluxcube(nflux);

  // load XOFFS and YOFFS
  load_offsets(fcube_file, fcube);

  // load the segmentation image into the structure
  fcube->segmentation = load_segmentation(fcube_file);
  //  gsl_to_FITSimage (fcube->segmentation, "gogo.fits", 1, "MOO");

  // load the fluximages
  for (i=0; i < nflux; i++)
    {
      fcube->fluxims[i] = load_fluximage(fcube_file, i+3);
    }

  // fill the number of fluximages
  fcube->n_fimage = nflux;

  // order the fluximages and store the
  // vector with the ordered indices
  fcube->fimage_order = order_fluxims(fcube);

  return fcube;
}

/**
 * The function allocates space for a new
 * fluxcube structure.
 *
 * @param  nflux - the number of fluximages
 *
 * @return ret   - the fluxcube structure
 */
flux_cube *
alloc_fluxcube(const int nflux)
{
   flux_cube *fcube;

   // allocate space for the flux_cube structure
   fcube = (flux_cube *)malloc(sizeof(flux_cube));

   // allocate space for the array of fluximages
   fcube->fluxims      = (flux_image **)malloc(nflux * sizeof(flux_image *));

   // return the result
   return fcube;
}

/**
 * The function fills the keywords XOFFS and YOFFS
 * from the primary header of the fitscube file
 * into the according items of the fitscube structure.
 * In case that those keywords are not present,
 * they are set to 0.
 *
 * @param  fcube_file - the file name of the fluxcube
 * @param  fcube      - the fluxcube structure
 *
 * @return 1          - always true
 */
int
load_offsets(const char fcube_file[], flux_cube *fcube)
{
  char tmp[MAXCHAR];

  double dtmp=0.0;

  // avoid problems with the 'const' declaration
  strcpy(tmp, fcube_file);

  // read XOFFS from the primary header
  dtmp = get_float_from_keyword(tmp, 1, "XOFFS");

  // store the keyword value or 0 (if not present) in the fcube structure
  if (isnan(dtmp))
    fcube->xoffs = 0;
  else
    fcube->xoffs = (int)dtmp;

  // read XOFFS from the primary header
  dtmp = get_float_from_keyword(tmp, 1, "YOFFS");

  // store the keyword value or 0 (if not present) in the fcube structure
  if (isnan(dtmp))
    fcube->yoffs = 0;
  else
    fcube->yoffs = (int)dtmp;

  return 1;
}

/**
 * The function creates and loads the segmentation image
 * of a flux_cube file into a gsl_matrix_int structure.
 *
 * @param  fcube_file - the file name of the fluxcube
 *
 * @return ret        - pointer to the gsl_matrix_int structure
 */
gsl_matrix_int *
load_segmentation(const char fcube_file[])
{

  gsl_matrix_int *segmentation;
  gsl_matrix     *dvals;

  int i,j;

  double diff;

  // load the segmentation image into a gsl_matrix_(double)
  dvals = FITSimage_to_gsl(fcube_file, 2, 1);

  // add 0.5 to get proper rounding with (int)
  gsl_matrix_add_constant (dvals, 0.5);

  // allocate space for the integer matrix
  segmentation = gsl_matrix_int_alloc(dvals->size1,dvals->size2);

  // fill the integer matrix with values from the double matrix
  for (i=0; i < (int)dvals->size1; i++)
    {
      for (j=0; j < (int)dvals->size2; j++)
        {
          gsl_matrix_int_set(segmentation, i, j, (int)gsl_matrix_get(dvals, i, j));

          diff = gsl_matrix_get(dvals, i, j) - (double)gsl_matrix_int_get(segmentation, i, j);
          if (diff > 0.6 || diff < 0.4)
            fprintf(stdout, "diff %f; %f --> %i, (%i, %i)\n", diff, gsl_matrix_get(dvals, i, j), gsl_matrix_int_get(segmentation, i, j), i, j);
        }
    }

  // free the space for the double matrix
  gsl_matrix_free(dvals);

  return segmentation;
}

/**
 * The function creates and loads a fluximage structure
 * from an extention of a fluxcube file.
 *
 * @param  fcube_file - the file name of the fluxcube
 * @param  hdunum     - the extension number to load in a fluximage
 *
 * @return fimage     - pointer to the fluximage structure loaded
 */
flux_image *
load_fluximage(const char fcube_file[], int hdunum)
{
  flux_image *fimage;

  char tmp[MAXCHAR];

  double wavelength;

  strcpy(tmp, fcube_file);

  // allocate space for the fluximage structure
  fimage = (flux_image *)malloc(sizeof(flux_image));

  // load the array from the specified extension number
  fimage->flux = FITSimage_to_gsl(fcube_file, hdunum, 1);

  // get the keyword "WAVELENG' from the specified extension number
  wavelength = (double)get_float_from_keyword(tmp, hdunum, "WAVELENG");

  // report an error in case that the extension does not have that keyword.
  // Otherwise store its content in the fluximage structure
  if (isnan(wavelength))
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "aXe_PETCONT: " "Missing keyword WAVELENGTH in fluxcube %s, extension %i!\n",
                fcube_file, hdunum);
  else
    fimage->wavelength = wavelength;

  return fimage;
}


/**
 * The function analyzes the fluximages in a
 * fluxcube structure and sorts their indices in
 * ascending wavelength. The vector with the
 * ordered indices is given back.
 *
 * @param  fcube - the pointer to the fluxcube structure
 *
 */
gsl_vector_int *
order_fluxims(flux_cube *fcube)
{
  gsl_vector     *fimage_wavelength;

  gsl_vector_int *fimage_order;

  int i, j, k;

  // allocate space for the order vector
  fimage_order = gsl_vector_int_alloc(fcube->n_fimage);

  // put a default order in the vector
  for (i=0; i<fcube->n_fimage; i++)
    gsl_vector_int_set(fimage_order, i, i);

  // order the wavelength
  if (fcube->n_fimage > 1)
    {

      // create a vector for the wavelengths
      fimage_wavelength = gsl_vector_alloc(fcube->n_fimage);

      // set the first entry for the wavelength
      gsl_vector_set(fimage_wavelength, 0, fcube->fluxims[0]->wavelength);

      // go over all fluximages
      for (i=1; i<fcube->n_fimage; i++)
        {

          // find the correct position of the fluximage
          j=0;
          while (j < i && gsl_vector_get(fimage_wavelength, j) <= fcube->fluxims[i]->wavelength)
            j++;

          // check whether the correct poition of the fluximage is at the end
          if (j == i)
            {

              // append the fluximage at the end
              gsl_vector_set(fimage_wavelength, j, fcube->fluxims[i]->wavelength);
              gsl_vector_int_set(fimage_order, j, i);
            }
          else
            {

              // move all fluximages one index up
              for (k=i; k > j; k--)
                {
                  gsl_vector_set(fimage_wavelength, k, gsl_vector_get(fimage_wavelength, k-1));
                  gsl_vector_int_set(fimage_order, k, gsl_vector_int_get(fimage_order, k-1));
                }

              // place the flux image at the correct position
              gsl_vector_set(fimage_wavelength, j, fcube->fluxims[i]->wavelength);
              gsl_vector_int_set(fimage_order, j, i);
            }
        }

      // release the space for the wavelength vector
      gsl_vector_free(fimage_wavelength);
    }

  // return the order vector
  return fimage_order;
}

/**
 * The function releases the space allocated
 * for a fluxcube structure.
 *
 * @param  fcube - the pointer to the fluxcube structure
 *
 */
void
free_fluxcube(flux_cube *fcube)
{
  int i;

  // free the segmentation matrix
  gsl_matrix_int_free(fcube->segmentation);

  // free the wavelength order vector
  gsl_vector_int_free(fcube->fimage_order);

  // go over the fluximages, and free them
  for (i=0; i < fcube->n_fimage; i++)
    free_fluximage(fcube->fluxims[i]);

  // free the fluximage vector
  free(fcube->fluxims);

  // free the fluximage order vector
  //  gsl_vector_int_free(fcube->fimage_order);

  // free the fluxcube
  free(fcube);

  // set the fluxcube to NULL
  fcube=NULL;
}

/**
 * The function releases the space allocated
 * for a fluximage structure.
 *
 * @param  fimage - the pointer to the fluximage structure
 *
 */
void
free_fluximage(flux_image *fimage)
{
  // free the flux matrix
  gsl_matrix_free(fimage->flux);

  // free the fluximage
  free(fimage);

  // set the fluximage to NULL
  fimage=NULL;
}

/**
 * The function transforms a pixel coordinate point in a
 * flt image into an pixel coordinate point in the associated
 * fluxcube image
 *
 * @param  fcube       - the fluxcube structure
 * @param  flt_point   - pixel coordinate point in the flt-image
 *
 * @return fcube_point - the pixel coordinate in the fluxcube
 */
d_point
flt_to_fcube_trans(const flux_cube *fcube, const d_point flt_point)
{
  d_point fcube_point;

  fcube_point.x = flt_point.x - (double)fcube->xoffs;
  fcube_point.y = flt_point.y - (double)fcube->yoffs;

  return fcube_point;
}

/**
 * The function transforms an image coordinate point in a
 * fluxcube into an image coo point in the associated
 * flt image
 *
 * @param  fcube       - the fluxcube structure
 * @param  fcube_point - image-coo in the fcube
 *
 * @return flt_point   - the image coo in the flt
 */
d_point
fcube_to_flt_trans(const flux_cube *fcube, const d_point fcube_point)
{
  d_point flt_point;

  flt_point.x = fcube_point.x + (double)fcube->xoffs;
  flt_point.y = fcube_point.y + (double)fcube->yoffs;

  return flt_point;
}


/**
 * The function converts the information stored in the segmenation
 * matrix of a fluxcube and the object list to a list of direct objects
 *
 * @param  fcube   - the fluxcube structure
 * @param  oblist  - the object list
 *
 * @return dirlist - the dirobject list created
 */
dirobject **
fluxcube_to_dirlist(const flux_cube *fcube, object  **oblist)
{
  dirobject **dirlist;
  dirobject  *actdir;

  //px_point fcube_point;
  px_point flt_point;

  d_point tmp1, tmp2;

  int iact=0;
  int nobjects=0;
  int ndir=0;
  int i, j;
  int objindex;

  //int itmp=0;

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // allocate space for the dirobject list
  dirlist = (dirobject **) malloc((nobjects+1) * sizeof(dirobject *));
  if (dirlist == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "fluxcube_to_dirlist:" " Could not allocate"
                 " memory for pointers to %i dirobject objects", nobjects+1);

  // initialize the number of direct object
  // and set the first item to NULL
  ndir=0;
  dirlist[ndir]=NULL;

  //  if (nobjects < 1)
  //    return dirlist;

  // go over each column in the segmentation image
  for (i=0; i < (int)fcube->segmentation->size1; i++)
    {
      // go over each row in the segmentation image
      for (j=0; j < (int)fcube->segmentation->size2; j++)
        {

          // get the pixel value in the segmentation image
          iact = gsl_matrix_int_get(fcube->segmentation, i, j);

          // do something if the pixel is part of an object
          if (iact > 0)
            {

              // transform the cooridnates from the fcube-system
              // into the flt system
              tmp1.x = (double)i;
              tmp1.y = (double)j;
              tmp2 = fcube_to_flt_trans(fcube, tmp1);
              flt_point.x = (int)tmp2.x;
              flt_point.y = (int)tmp2.y;

              // try to get a dirobject for that object (=pixel value)
              actdir = get_dirobject_from_list(dirlist, iact);

              // check whether a dirobject just exists
              if (actdir != NULL)
                {
                  // udate an existing object
                  update_dirobject(actdir, flt_point);
                }
              else
                {

                  // find the the object index for a (hypotetical) new dirobject
                  objindex = find_object_in_object_list(oblist, iact);
                  // check whether object exists in oblist
                  if (objindex >-1)
                    {
                      // create a new dirobject and enhance the counter
                      dirlist[ndir++] = dirobject_from_segpoint(flt_point, oblist[objindex]);

                      // set the actually last dirobject to NULL
                      // to get a properly closed list at all time
                      dirlist[ndir]=NULL;
                    }
                }
            }
        }
    }


  // return the direct object list
  return dirlist;
}

/**
 * The function creates a a dirobject on the basis
 * of an object and the pixel coos of a point
 * which is part of that object.
 *
 * @param  flt_point - the fluxcube structure
 * @param  actobject - the object list
 *
 * @return dirobject - the dirobject created
 */
dirobject *
dirobject_from_segpoint(const px_point flt_point, const object  *actobject)
{

  dirobject  *actdir;

  // allocate space for the dirobject
  actdir = (dirobject *) malloc (sizeof (dirobject));
  if (actdir == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "fill_dirobject:" " Could not allocate"
                 " memory for a dirobject object");

  // transfer refpoint and ID
  actdir->ID = actobject->ID;
  actdir->refpoint.x = actobject->beams[0].refpoint.x;
  actdir->refpoint.y = actobject->beams[0].refpoint.y;

  // derive and store min/max in x/y for the corners
  actdir->ix_min = flt_point.x;
  actdir->ix_max = flt_point.x;
  actdir->iy_min = flt_point.y;
  actdir->iy_max = flt_point.y;

  // derive the distortion scales along the major an minor axis
  actdir->drzscale.x = 1.0;
  actdir->drzscale.y = 1.0;

  // set the variable
  actdir->bb_sed = 0;

  // make it a dummy
  actdir->SED = NULL;

  // return the new dirobject
  return actdir;
}

/**
 * The function updates a dirobject with a point
 * belonging to the same object. The new point
 * may extend the object size.
 *
 * @param actdir    - the direct object
 * @param flt_point - the new coordinate point
 *
 */
void
update_dirobject(dirobject *actdir, const px_point flt_point)
{
  // check whether the new points extends the
  // old min/max values in x or y.
  // correct the values if necessary.
  actdir->ix_min = MIN(actdir->ix_min, flt_point.x);
  actdir->ix_max = MAX(actdir->ix_max, flt_point.x);
  actdir->iy_min = MIN(actdir->iy_min, flt_point.y);
  actdir->iy_max = MAX(actdir->iy_max, flt_point.y);
}


/**
 * The function computes and stores at the according position
 * the xy offsets of a direct image. This is only important
 * in the fluxcube emission model, since the image coordinates
 * in fluxcubes are 'filter' coordinates and do not take into
 * account the offsets introduce by the grism.
 * For each object and beam those offsets are evaluated
 * and stored in the direct object structure.
 *
 * @param dirlist   - the direct object
 * @param CONF_file - the new coordinate point
 *
 */
void fill_xy_offsets(dirobject **dirlist, char CONF_file[])
{
  aperture_conf   *conf;

  gsl_vector *x_coeffs;
  gsl_vector *y_coeffs;

  d_point m_point;

  //double xoffs, yoffs;

  int beamID=0;
  int i;

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  // go over each beam defined in the configuration
  // file
  for (beamID=0; beamID < conf->nbeams; beamID++)
    {

      // determine the coeeficients for the offsets
      // in x and in y
      x_coeffs = get_beam_trace_xoff (CONF_file, beamID);
      y_coeffs = get_beam_trace_yoff (CONF_file, beamID);

      // go over each direct object
      i=0;
      while (dirlist[i] !=NULL)
        {

          // get the mean direct object position
          m_point = get_dirobject_meanpos(dirlist[i]);

          // evaluate the coefficients for the 2D variable offsets
          // at the mean position
          dirlist[i]->xy_off[beamID].x = eval_trace_off_at_pos (x_coeffs, m_point, beamID);
          dirlist[i]->xy_off[beamID].y = eval_trace_off_at_pos (y_coeffs, m_point, beamID);

          // iterate the counter
          i++;
        }

      // free the vectors for the coefficients
      gsl_vector_free(x_coeffs);
      gsl_vector_free(y_coeffs);
    }
        // release memory
        free_aperture_conf(conf);
}

/**
 * The function fills the flux information from a point in
 * a fluxcube structure into the flux part of a direct object.
 *
 * @param fcube  - the fluxcube structure
 * @param point  - the position to take the flux from
 * @param actdir - the direct object to update
 *
 */
void fill_fluxvalues(const flux_cube *fcube, const px_point point,
                     dirobject *actdir, const int inter_type)
{


  energy_distrib *sed;
  double *sed_wavs;
  double *sed_flux;

  int i;
  int iact;


  sed = (energy_distrib*) malloc(sizeof(energy_distrib));
  sed_wavs = (double*) malloc(fcube->n_fimage*sizeof(double));
  sed_flux = (double*) malloc(fcube->n_fimage*sizeof(double));

  if (actdir->SED)
    free_enerdist (actdir->SED);

  // go over each fluximage in the fluxcube structure
  for (i=0; i < fcube->n_fimage; i++)
    {

      // determine the index of the next fluximage
      iact = gsl_vector_int_get(fcube->fimage_order, i);

      // transfer the wavelength and flux from the fluxcube
      sed_wavs[i] = fcube->fluxims[iact]->wavelength;
      sed_flux[i] = gsl_matrix_get(fcube->fluxims[iact]->flux, point.x, point.y);
    }

  sed->npoints    = fcube->n_fimage;
  sed->wavelength = sed_wavs ;
  sed->flux       = sed_flux;

  if (fcube->n_fimage > 1)
    {
      if (inter_type == 2)
        sed->interp     = gsl_interp_alloc (gsl_interp_polynomial, (size_t)fcube->n_fimage);
      else if (inter_type == 3)
        sed->interp     = gsl_interp_alloc (gsl_interp_cspline, (size_t)fcube->n_fimage);
      else
        sed->interp     = gsl_interp_alloc (gsl_interp_linear, (size_t)fcube->n_fimage);

      sed->accel      = gsl_interp_accel_alloc ();
      gsl_interp_init (sed->interp, sed->wavelength, sed->flux, (size_t)sed->npoints);
    }
  else
    {
      sed->interp = NULL;
      sed->accel= NULL;
    }

  // transfer the SED;
  // announce the SED
  actdir->SED    = sed;
  actdir->bb_sed = 1;
}

/**
 * The function prints the main properties of a
 * fluxcube structure onto the screen.
 * This is mostly for debuggin/development
 * purposes, and not for production
 *
 * @param fcube      - the fluxcube structure
 *
 */
void
print_fluxcube(const flux_cube *fcube)
{
  int i;

  fprintf(stdout, "XOFFS,YOFFS: %i,%i\n", fcube->xoffs, fcube->yoffs);
  fprintf(stdout, "Number of flux-extension: %i\n", fcube->n_fimage);

  for (i=0; i<fcube->n_fimage; i++)
    fprintf(stdout, "Wavelength for flux extension %i: %f\n", i,
            fcube->fluxims[gsl_vector_int_get(fcube->fimage_order,i)]->wavelength);
}
