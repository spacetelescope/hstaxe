/**
 *  Subroutines to calculate the
 *  various contamination models
*/
#include <time.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

#include "inout_aper.h"
#include "spc_model.h"
#include "spce_sect.h"
#include "spc_back.h"
#include "spce_PET.h"
#include "aXe_errors.h"
#include "fringe_conf.h"
#include "spc_resp.h"
#include "spce_pathlength.h"
#include "aper_conf.h"
#include "specmodel_utils.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define SQR(x) ((x)*(x))


/*-------------------
 * Section 1: Gauss contamination and directly
 *            dependent subroutines
 */
/**
 * Function: compute_gauss_cont
 * The subroutine computes ans stores the quantitative contamination
 * using the Gaussian emission model. The emitting object with Gaussian
 * shape and SED as derived from the AB filter magnitudes are composed
 * The individual nbeams are modelled, then the complete contamination
 * image is composed from the beam models.
 * Finally, the contaminating flux is computed by subtracting the
 * beam emission from the contamination image for all PET pixels.
 *
 * Parameters:
 * @param grism_file  - the full name of the grism file
 * @param OAF_file    - the name of the aperture file
 * @param CONF_file   - the full name of configuration file
 * @param specmod_file- fits table with the model spectra
 * @param model_scale - the scale for extension of the direct object area
 * @param inter_type  - the interpolation method for the flux values
 * @param lambda_psf  - wavelength the Gaussian parameters were determined at
 * @param obs         - the observation
 * @param PET_file    - the name of the PET which is modified
 * @param map_file    - the name of the contamination map
 * @param store       - flag whether the contamination image is stored or not
 *
 * Returns:
 * @return status     - returns success or failure
 */
int
compute_gauss_cont(char grism_file[], char OAF_file[], char CONF_file[],
                   const char specmod_file[], const double model_scale,
                   const int inter_type, const double lambda_psf,
                   observation *obs, const char PET_file[], char map_file[],
                   const int store)
{

  object    **oblist;
  dirobject **dirlist;
  beamspec  **speclist;
  spectral_models *spec_mod;

  //dirobject  *actdir;
  gsl_matrix *all_models;
  char model_name[60];

  px_point npixels;

  //time_t timer;

  // load the object list
  fprintf (stdout, "aXe_PETCONT: Loading object aperture list...");
  oblist = file_to_object_list_seq (OAF_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // check whether highres models
  // are given
  if (strlen(specmod_file) > 0)
    // load the spectral models
    spec_mod = load_spectral_models(specmod_file);
  else
    // or set the struct to NULL
    spec_mod = NULL;

  // get the image dimensions
  npixels = get_npixel(obs);

  // create the list of direct objects to be simulated
  dirlist = oblist_to_dirlist(grism_file, CONF_file, npixels, oblist, spec_mod, model_scale, inter_type);

  //  timer=time(NULL);
  //  printf("The current time is %s.\n",asctime(localtime(&timer)));
  speclist = make_gauss_spectra(oblist, dirlist, lambda_psf, npixels, CONF_file);
  //  timer=time(NULL);
  //  printf("The current time is %s.\n",asctime(localtime(&timer)));

  // compose the contamination image from the modelled beams
  all_models = make_model_image(npixels, obs, speclist);

  // check whether the contamination image
  // should be stored
  if (store)
    // store the contamination image
    gsl_to_FITSimage (all_models, map_file, 1, NULL);

  // store the name of the
  // contamination model
  sprintf (model_name, "GAUSS");

  // put the contamination info into the PET
  fill_contam_info(PET_file, speclist, all_models, model_name);

  // release allocated memory
  // in the various structures
  gsl_matrix_free(all_models);
  if (spec_mod != NULL)
    free_spectral_models(spec_mod);
  free_dirlist(dirlist);
  free_speclist(speclist);
  if (oblist !=NULL)
    free_oblist (oblist);

  // return always '1'
  return 1;
}

/**
 * Function: compute_gaussdirim_cont
 * The subroutine computes ans stores the quantitative contamination
 * using the Gaussian emission model. The object shape can be defined
 * in a fits file with object models, and the SED may be defined in a
 * fits table with high resolution spectra. If no detailed object shape
 * or SED is given for an individual object, Gaussian shape and SED as
 * derived from the AB filter magnitudes are composed.
 * The individual nbeams are modelled, then the complete contamination
 * image is composed from the beam models.
 * Finally, the contaminating flux is computed by subtracting the
 * beam emission from the contamination image for all PET pixels.
 *
 * Parameters:
 * @param grism_file  - the full name of the grism file
 * @param OAF_file    - the name of the aperture file
 * @param CONF_file   - the full name of configuration file
 * @param specmod_file- fits table with the model spectra
 * @param objmod_file - fits image with object models
 * @param model_scale - the scale for extension of the direct object area
 * @param inter_type  - the interpolation method for the flux values
 * @param lambda_psf  - wavelength the Gaussian parameters were determined at
 * @param obs         - the observation
 * @param PET_file    - the name of the PET which is modified
 * @param map_file    - the name of the contamination map
 * @param store       - flag whether the contamination image is stored or not
 *
 * Returns:
 * @return status     - returns success or failure
 */
int
compute_gaussdirim_cont(char grism_file[], char OAF_file[], char CONF_file[],
                        const char specmod_file[],  const char objmod_file[],
                        const double model_scale, const int inter_type,
                        const double lambda_psf, observation *obs,
                        const char PET_file[], char map_file[], const int store)
{

  object        **oblist;
  dirobject     **dirlist;
  beamspec      **speclist;
  spectral_models *spec_mod;
  object_models   *obj_mod;

  //dirobject  *actdir;
  gsl_matrix *all_models;
  char model_name[60];

  px_point npixels;

  //time_t timer;

  // load the object list
  fprintf (stdout, "aXe_PETCONT: Loading object aperture list...");
  oblist = file_to_object_list_seq (OAF_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // check whether highres models
  // are given
  if (strlen(specmod_file) > 0)
    // load the spectral models
    spec_mod = load_spectral_models(specmod_file);
  else
    // or set the struct to NULL
    spec_mod = NULL;

  // check whether direct emission models
  // are given
  if (strlen(objmod_file) > 0)
    obj_mod = load_object_models(objmod_file);
  else
    obj_mod = NULL;

  // get the image dimensions
  npixels = get_npixel(obs);

  // create the list of direct objects to be simulated
  dirlist = oblist_to_dirlist2(grism_file, CONF_file, npixels, oblist,
                               spec_mod, obj_mod, model_scale, inter_type);

  //timer=time(NULL);
  //printf("The current time is %s.\n",asctime(localtime(&timer)));
  speclist = make_gauss_spectra2(oblist, dirlist, lambda_psf, npixels, CONF_file);
  //timer=time(NULL);
  //printf("The current time is %s.\n",asctime(localtime(&timer)));

  // compose the contamination image from the modelled beams
  all_models = make_model_image(npixels, obs, speclist);

  // check whether the contamination image
  // should be stored
  if (store)
    // store the contamination image
    gsl_to_FITSimage (all_models, map_file, 1, NULL);

  // store the name of the
  // contamination model
  if (obj_mod != NULL)
    sprintf (model_name, "DIRIM");
  else
    sprintf (model_name, "GAUSS");

  // check whether a PET exists
  if (strlen(PET_file) > 0)
    // put the contamination info into the PET
    fill_contam_info(PET_file, speclist, all_models, model_name);

  // release allocated memory
  // in the various structures
  gsl_matrix_free(all_models);
  if (spec_mod != NULL)
    free_spectral_models(spec_mod);
  free_dirlist(dirlist);
  if (obj_mod != NULL)
    free_object_models(obj_mod);
  free_speclist(speclist);
  if (oblist !=NULL)
    free_oblist (oblist);

  // return always '1'
  return 1;
}

/**
 * Function: make_gauss_spectra2
 * The function creates a spectral model for the beams
 * that are considered in the contamination. The modelled
 * beams are returned as a list.
 *
 * Parameters:
 * @param  oblist     - the object list as input to select beams
 * @param  dirlist    - the direct object list to dimension the models
 * @param  lambda_psf - the wavelength the object psf was determined at
 * @param  npixels    - the dimensions of the model for the whole image
 * @param  CONF_file  - the name of the configuration file
 *
 * Returns:
 * @return speclist  - the list of modelled beams
 */
beamspec **
make_gauss_spectra2(object **oblist, dirobject **dirlist,
                   const double lambda_psf, const px_point npixels,
                   char CONF_file[])
{
  beamspec       **speclist;
  dirobject       *actdir;
  aperture_conf   *conf;
  calib_function  *wl_calibration;
  spectrum        *resp;
  beam             actbeam;
  tracedata       *acttrace;

  //double eval=0.0;
  double psf_offset=0;

  //int nspecs;
  //int i=0;
  //int j=0;
  //int jj=0;
  int ii=0;
  int nobjects;

  //int kk, ll;
  double sval;
  double frac_prev, frac;

  int nx, ny;
  d_point dpixel;

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  // allocate ther list of spectral beams
  speclist = alloc_beamlist_from_dirlist(oblist, dirlist, npixels, conf);

  // go over each beam model
  ii = 0;
  while (speclist[ii] != NULL)
    {
      // get the direct object for the actual model spectrum
      actdir = get_dirobject_from_list(dirlist, speclist[ii]->objectID);

      // get the beam for the actual model spectrum
      actbeam = get_beam_for_beamspec(oblist,nobjects,speclist[ii]);

      // get the psf offset values
      psf_offset = get_psf_offset(conf, actbeam);

      // get the wavelength calibration for the actual model spectrum
      wl_calibration = get_calib_function(speclist[ii], actdir, CONF_file, conf);

      // get the sensitivity data for the actual model spectrum
      resp = get_throughput_spec(speclist[ii], CONF_file);

      // fill the tracedata structure for the model spectrum
      //acttrace = compute_tracedata(actbeam,actdir, wl_calibration,speclist[ii]);
      acttrace = compute_short_tracedata(conf, actbeam,actdir, wl_calibration,speclist[ii]);

      if (acttrace->npoints < 1)
        {
          // release the space for the various structures
          free_calib(wl_calibration);
          free_spectrum(resp);
          free_tracedata(acttrace);

          // give feedback to the screen
          fprintf(stderr, "aXe_PETCONT: skipping object %i beam %c ...", speclist[ii]->objectID, BEAM(speclist[ii]->beamID));

          // enhance the spectral beam counter
          // and go to the next
          ii++;
          continue;
        }

      // fill the flux information int the tracedata
      fill_fluxfrom_SED(actdir, acttrace);

      // give feedback to the screen
      fprintf(stdout, "aXe_PETCONT: modelling object %i beam %c ...", speclist[ii]->objectID, BEAM(speclist[ii]->beamID));

      frac_prev=10.0;
      // go over each pixel in the direct object area
      for (nx=actdir->ix_min; nx<=actdir->ix_max; nx++)
        {
    	  //----------------------------------------------
    	  // compute the percentage that has been computed
    	  // report progress in 10% increments
    	  frac = fmod(100.0*(nx - actdir->ix_min)/(actdir->ix_max - actdir->ix_min),10);
    	  if (frac < frac_prev)
    	  {
            fprintf(stdout, " %i ", (int)100.0*(nx - actdir->ix_min)/(actdir->ix_max - actdir->ix_min));
            fflush(stdout);
    	  }
    	  frac_prev=frac;
    	  //----------------------------------------------

    	  for (ny=actdir->iy_min; ny<=actdir->iy_max; ny++)
            {
              // fill the dpixel structure
              dpixel.x = (double)nx;
              dpixel.y = (double)ny;

              if (actdir->dirim)
                {
                  sval = get_diremission_value(actdir->dirim, dpixel.x - actbeam.refpoint.x, dpixel.y - actbeam.refpoint.y);
                  gsl_vector_set_all (acttrace->gvalue, sval);
                }
              else
                {
                  // check whether a wavelength-dependent
                  // emission profile is given
                  if ((conf->psfcoeffs && conf->psfrange) || psf_offset)
                    {
                      // fill in the wavelength dependend
                      // emission values
                      fill_gaussvalues(dpixel, actbeam, actdir, lambda_psf, conf, psf_offset, acttrace);
                    }
                  else
                    {
                      // do a subsampling over the pixel
                      // to get a more appropriate value for the
                      // emission val
                      sval = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);
                      gsl_vector_set_all (acttrace->gvalue, sval);
                    }
                }

              // insert the spectrum of this direct object pixel in the beam spectrum
              fill_pixel_in_speed(actdir, acttrace, dpixel, resp, speclist[ii], wl_calibration);
            }
        }

      // release the space for the various structures
      fprintf(stdout, " Done\n");
      free_calib(wl_calibration);
      free_spectrum(resp);
      free_tracedata(acttrace);

      // enhance the counter
      ii++;
    }

  // free the memory in the conf structure
  free_aperture_conf(conf);

  // return the list of modelled beams
  return speclist;
}

/**
 * Function: make_gauss_spectra
 * The function creates a spectral model for the beams
 * that are considered in the contamination. The modelled
 * beams are returned as a list.
 *
 * Parameters:
 * @param  oblist     - the object list as input to select beams
 * @param  dirlist    - the direct object list to dimension the models
 * @param  lambda_psf - the wavelength the object psf was determined at
 * @param  npixels    - the dimensions of the model for the whole image
 * @param  CONF_file  - the name of the configuration file
 *
 * Returns:
 * @return speclist  - the list of modelled beams
 */
beamspec **
make_gauss_spectra(object **oblist, dirobject **dirlist,
                   const double lambda_psf, const px_point npixels,
                   char CONF_file[])
{
  beamspec       **speclist;
  dirobject       *actdir;
  aperture_conf   *conf;
  calib_function  *wl_calibration;
  spectrum        *resp;
  beam             actbeam;
  tracedata       *acttrace;

  //double eval=0.0;
  double psf_offset=0;

  int nspecs;
  int i=0;
  int j=0;
  int jj=0,ii=0;
  int nobjects;

  //int kk, ll;
  double sval;

  int nx, ny;
  d_point dpixel;

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  // get the number of beams included in the contamination
  // (mag < mag_mark(BEAM)
  nspecs = get_beamspec_size(oblist);
  speclist = (beamspec  **) malloc((nspecs+1) * sizeof(beamspec  *));
  if (speclist == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "make_model_spectra:" " Could not allocate"
                 " memory for pointers to %i beamspec objects", nspecs+1);

  fprintf(stdout, "aXe_PETCONT: %i beams are included in the contamination.\n",nspecs);

  // go over all objects
  for (i = 0; i < nobjects; i++)
    {

      // go over each beam in an object
      for (j=0; j < oblist[i]->nbeams; j++)
        {

          // check whether the beam is included in the contamination
          if (oblist[i]->beams[j].ignore != 1)
            {

              // get the direct object for the beam
              actdir = get_dirobject_from_list(dirlist, oblist[i]->ID);

              // check whether the beam is inside the image,
              // if yes, allocate space for the beam model
              speclist[jj] = dimension_beamspec(actdir, oblist[i],  npixels, conf, j);

              // increment the counter if space was allocated
              if (speclist[jj] != NULL)
                jj++;

            }
        }
    }

  fprintf(stdout, "aXe_PETCONT: %i beams are modelled.\n",jj);

  // terminate the beamspec list with NULL's
  for (ii=jj; ii < nspecs+1; ii++)
    speclist[ii] = NULL;

  // go over each beam model
  ii = 0;
  while (speclist[ii] != NULL)
    {

      //      fprintf(stdout, "modelling 1: %i\n", ii);

      // get the direct object for the actual model spectrum
      actdir = get_dirobject_from_list(dirlist, speclist[ii]->objectID);

      // get the beam for the actual model spectrum
      actbeam = get_beam_for_beamspec(oblist,nobjects,speclist[ii]);

      //      fprintf(stdout, "modelling 2: %i\n", ii);
      psf_offset = get_psf_offset(conf, actbeam);

      //      fprintf(stdout, "modelling 3: %i\n", ii);

      // get the wavelength calibration for the actual model spectrum
      wl_calibration = get_calib_function(speclist[ii], actdir, CONF_file, conf);

      // get the sensitivity data for the actual model spectrum
      resp = get_throughput_spec(speclist[ii], CONF_file);

      // fill the tracedata structure for the model spectrum
      acttrace = compute_tracedata(actbeam,actdir, wl_calibration,speclist[ii]);
      //      fprintf(stdout, "modelling 4: %i\n", ii);
      if (acttrace->npoints < 1)
        {
          // release the space for the various structures
          //      fprintf(stdout, "modelling 50: %i\n", ii);
          free_calib(wl_calibration);
          free_spectrum(resp);
          free_tracedata(acttrace);
      fprintf(stderr, "aXe_PETCONT: skipping object %i beam %c ...", speclist[ii]->objectID, BEAM(speclist[ii]->beamID));

          ii++;
          continue;
        }

      //      fprintf(stdout, "modelling 51: %i\n", ii);
      // fill the flux information int the tracedata
      fill_fluxfrom_SED(actdir, acttrace);
      //      fprintf(stderr, "modelling 6: %i\n", ii);

      fprintf(stdout, "aXe_PETCONT: modelling object %i beam %c ...", speclist[ii]->objectID, BEAM(speclist[ii]->beamID));

      // go over each pixel in the direct object area
      for (nx=actdir->ix_min; nx<=actdir->ix_max; nx++)
        {
          for (ny=actdir->iy_min; ny<=actdir->iy_max; ny++)
            {

              //fprintf(stderr, "point %i %i %i <-> %i %i <-> %i \n", nx, ny, actdir->ix_min, actdir->ix_max, actdir->iy_min, actdir->iy_max);
              // fill the dpixel structure
              dpixel.x = (double)nx;
              dpixel.y = (double)ny;


              // check whether a wavelength-dependent
              // emission profile is given
              if ((conf->psfcoeffs && conf->psfrange) || psf_offset)
                {
                  // fill in the wavelength dependend
                  // emission values
                  fill_gaussvalues(dpixel, actbeam, actdir, lambda_psf, conf, psf_offset, acttrace);
                }
              else
                {
                  // do a subsampling over the pixel
                  // to get a more appropriate value for the
                  // emission val
                  sval = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);
                  gsl_vector_set_all (acttrace->gvalue, sval);
                }

              //
              fill_pixel_in_speed(actdir, acttrace, dpixel, resp, speclist[ii], wl_calibration);
            }
        }

      // release the space for the various structures
      fprintf(stdout, " Done\n");
      free_calib(wl_calibration);
      free_spectrum(resp);
      free_tracedata(acttrace);

      // enhance the counter
      ii++;
    }

  // fre the memory in the conf structure
  free_aperture_conf(conf);

  // return the list of modelled beams
  return speclist;
}


/*-------------------
 * Section 2: Fluxcube contamination and directly
 *            dependent subroutines
 */
/*
 * Function: compute_fcube_cont
 * The subroutine computes ans stores the quantitative contamination
 * using the Fluxcube emission model. The emitting objects are defined
 * in fluxcube images, which were derived from MultiDrizzled direct images.
 * The individual nbeams are modelled, then the complete contamination
 * image is composed from the beam models.
 * Finally, the contaminating flux is computed by subtracting the
 * beam emission from the contamination image for all PET pixels.
 *
 *
 * Parameters:
 * @param grism_file  - the full name of the grism file
 * @param OAF_file    - the name of the aperture file
 * @param fcube_file  - the name of the fluxcube file
 * @param CONF_file   - the full name of configuration file
 * @param model_scale - the scale for extension of the direct object area
 * @param inter_type  - the interpolation method for the flux values
 * @param obs         - the observation
 * @param PET_file    - the name of the PET which is modified
 * @param map_file    - the name of the contamination map
 * @param store       - flag whether the contamination image is stored or not
 *
 * Returns:
 * @return status     - returns success or failure
 */
int
compute_fcube_cont(char grism_file[], char OAF_file[], char fcube_file[],
                   char CONF_file[], const double model_scale, const int inter_type,
                   observation *obs, const char PET_file[], char map_file[],
                   const int store)
{
  object    **oblist;
  dirobject **dirlist;
  beamspec  **speclist;
  flux_cube  *fcube;

  //dirobject  *actdir;
  gsl_matrix *all_models;

  char model_name[60];

  px_point npixels;

  int i_type;

  fprintf (stdout, "aXe_PETCONT: Loading fluxcube image %s ...", fcube_file);
  fcube = load_fluxcube(fcube_file);
  fprintf (stdout, " Done\n");

  //  load the object list
  fprintf (stdout, "aXe_PETCONT: Loading object aperture list...");
  oblist = file_to_object_list_seq (OAF_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // get the dimension of the grism images
  npixels = get_npixel(obs);

  // generate the list of emitters from the
  // fluxcube image
  dirlist = fluxcube_to_dirlist(fcube, oblist);

  // make the offsets??
  fill_xy_offsets(dirlist, CONF_file);

  // check whether there is enough information
  // for the desired interpolation type
  i_type = check_interp_type(inter_type, fcube->n_fimage, 0);

  // model the beams
  speclist = make_fcube_spectra(oblist, dirlist, npixels, CONF_file, fcube, i_type);

  // compute the contamination image from the
  // modelled beams
  all_models = make_model_image(npixels, obs, speclist);

  // check whether the contamination
  // image should be stored
  if (store)
    // store the contamination image
    gsl_to_FITSimage (all_models, map_file, 1, NULL);

  // store the name of the contamination model
  sprintf (model_name, "FLUXCUBE");

  // check whether a PET exists
  if (strlen(PET_file) > 0)
    // compute and transfer the contamination
    // information fot the PET pixels
    fill_contam_info(PET_file, speclist, all_models, model_name);

  // free the memory allocated
  // in the various structures
  gsl_matrix_free(all_models);
  free_speclist(speclist);
  free_dirlist (dirlist);
  free_fluxcube(fcube);
  if (oblist != NULL)
    free_oblist (oblist);

  // return '1' as a dummy
  return 1;
}


/**
 * Function: make_fcube_spectra
 * The function creates a spectral model for the beams
 * that are considered in the contamination. It uses
 * the fluxcube emission model to calculate the emission
 * in the different beams.
 *
 * Parameters:
 * @param  oblist    - the object list as input to select beams
 * @param  dirlist   - the direct object list to dimension the models
 * @param  npixels   - the dimensions of the model for the whole image
 * @param  CONF_file - the name of the configuration file
 * @param  fcube     - the fluxcube to get the flux data from
 *
 * Returns:
 * @return speclist  - the list of modelled beams
 */
beamspec **
make_fcube_spectra(object **oblist, dirobject **dirlist,
                   const px_point npixels, char CONF_file[],
                   const flux_cube *fcube, const int inter_type)
{
  beamspec       **speclist;
  dirobject       *actdir;
  aperture_conf   *conf;
  calib_function  *wl_calibration;
  spectrum        *resp;
  beam             actbeam;
  tracedata       *acttrace;

  int i;
  int nobjects;

  int nx, ny;
  //d_point dpixel;
  d_point dflt_point;
  px_point fcube_point;

  d_point tmp2;


  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // load the configuration file
  conf = get_aperture_descriptor (CONF_file);

  speclist = alloc_beamlist_from_dirlist(oblist, dirlist, npixels, conf);

  // go over each beam model
  i = 0;
  while (speclist[i] != NULL)
    {

      // get the direct object for the actual model spectrum
      actdir = get_dirobject_from_list(dirlist, speclist[i]->objectID);

      // get the beam for the actual model spectrum
      actbeam = get_beam_for_beamspec(oblist,nobjects,speclist[i]);

      // get the wavelength calibration for the actual model spectrum
      wl_calibration = get_calib_function(speclist[i], actdir, CONF_file, conf);

      // get the sensitivity data for the actual model spectrum
      resp = get_throughput_spec(speclist[i], CONF_file);

      // fill the tracedata structure for the model spectrum
      acttrace = compute_tracedata(actbeam,actdir, wl_calibration,speclist[i]);
      if (acttrace->npoints < 1)
        {
          // release the space for the various structures
          free_calib(wl_calibration);
          free_spectrum(resp);
          free_tracedata(acttrace);
     fprintf(stdout, "aXe_PETCONT: skipping object %i beam %c ...", speclist[i]->objectID, BEAM(speclist[i]->beamID));

          continue;
        }

      fprintf(stdout, "aXe_PETCONT: modelling object %i beam %c ...", speclist[i]->objectID, BEAM(speclist[i]->beamID));

      // go over each pixel in the direct object area
      for (nx=actdir->ix_min; nx<=actdir->ix_max; nx++)
        {
          for (ny=actdir->iy_min; ny<=actdir->iy_max; ny++)
            {

              // transform the flt-coordinates
              // to fcube coordinates
              dflt_point.x = (double)nx;
              dflt_point.y = (double)ny;
              tmp2 = flt_to_fcube_trans(fcube, dflt_point);
              fcube_point.x = (int)tmp2.x;
              fcube_point.y = (int)tmp2.y;

              // check whether the coordinate point actually
              // does belong to the spectral beam
              if ( gsl_matrix_int_get(fcube->segmentation, fcube_point.x,fcube_point.y) == actdir->ID)
                {

                  // transfer the flux information from the current pixel
                  // to the SED of the direct object
                  fill_fluxvalues(fcube, fcube_point, actdir, inter_type);

                  // fill theflux-vector of the tracedata
                  fill_fluxfrom_SED(actdir, acttrace);

                  // compute the contribution of the current pixel to the
                  // current beam object
                  fill_pixel_in_speed(actdir, acttrace, dflt_point, resp, speclist[i], wl_calibration);
                }
            }
        }

      // release the space fpr the various structures
      fprintf(stdout, " Done\n");
      free_calib(wl_calibration);
      free_spectrum(resp);
      free_tracedata(acttrace);

      i++;
    }
  // release memory
  free_aperture_conf(conf);

  // return the spectrum list
  return speclist;
}


/**
 * Function: alloc_beamlist_from_dirlist
 * The function determines which beams of a direct object
 * are included in the contamination model. Then the size for each
 * of the spectral beams is estimated, and the space is allocated.
 *
 * Parameters:
 * @param  oblist    - the object list as input to select beams
 * @param  dirlist   - the direct object list to dimension the models
 * @param  npixels   - the dimensions of the model for the whole image
 * @param  conf      - configuration structure
 *
 * Returns:
 * @return speclist  - the list of modelled beams
 */
beamspec **
alloc_beamlist_from_dirlist(object **oblist, dirobject **dirlist,
                            const px_point npixels, aperture_conf *conf)
{
  beamspec  **speclist;
  //dirobject  *actdir;
  //object     *actobj;

  int nspecs;
  int i, j, jj=0;
  int objindex;

  // get the number of beams included in the contamination
  // (mag < mag_mark(BEAM)
  nspecs = get_beamspec_size(oblist);

  // allocate space for the vector of spectral beams
  speclist = (beamspec  **) malloc((nspecs+1) * sizeof(beamspec  *));
  if (speclist == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "make_model_spectra:" " Could not allocate"
                 " memory for pointers to %i beamspec objects", nspecs+1);

  fprintf(stdout, "aXe_PETCONT: %i beams are included in the contamination.\n",nspecs);

  // go over all objects
  i=0;
  while (dirlist[i] != NULL)
    {

      // find the object structure to the directo object
      objindex = find_object_in_object_list(oblist, dirlist[i]->ID);
      if (objindex < 0)
        {
          i++;
          continue;
        }

      // go over each beam in an object
      for (j=0; j < oblist[objindex]->nbeams; j++)
        {

          // check whether the beam is included in the contamination
          if (oblist[objindex]->beams[j].ignore != 1)
            {
              // check whether the beam is inside the image,
              // if yes, allocate space for the beam model
              //              fprintf(stdout, "Trying %i, beam %c\n", dirlist[i]->ID, BEAM(j));
              speclist[jj] = dimension_beamspec(dirlist[i], oblist[objindex],  npixels, conf, j);

              // increment the counter if space was allocated
              if (speclist[jj] != NULL)
                jj++;

            }
        }
      i++;
    }

  fprintf(stdout, "aXe_PETCONT: %i beams are modelled.\n",jj);

  // terminate the beamspec list with NULL's
  for (i=jj; i < nspecs+1; i++)
    speclist[i] = NULL;

  // return the list of spectral beams
  return speclist;
}


/*-------------------
 * Section 3: Subroutines directly dependent
 *            from Gaussian AND Fluxcube contamination
 */
/**
 * Function: dimension_beamspec
 * The function checks the boundaries of a candidate
 * model spectrum of a beam. In case that the
 * model spectrum is completely outside the frame,
 * a NULL model spectrum is returned.
 * Otherwise space is allocated for the model spectrum and
 * the matrix therein, and the data are filled in.
 * Only the matrix is not defined and filled.
 *
 * Parameters:
 * @param  dirlist   - the list with all dirobjects
 * @param  actobject - one beam of this object is examined
 * @param  npixels   - the size of the CCD
 * @param  conf      - the configutration structure
 * @param  j         - the beam number to be examined
 *
 * Returns:
 * @return actspec   - the beamspec object created
 */
beamspec *
dimension_beamspec(dirobject *actdir, object *actobject,
                   const px_point npixels, const aperture_conf * conf, int j)
{
  beamspec   *actspec;
  trace_func *tracefun;
  gsl_matrix *stamp;
  double     dx0, dy0, dx1, dy1;
  double     xmin, xmax, ymin, ymax;


  // allocate space for the beamspoec
  actspec = (beamspec *)malloc(sizeof(beamspec));
  if (actspec == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "make_model_spectra:" " Could not allocate"
                 " memory for a beamspec object");

  // determine the range of x-values covered by that beam
  // (relative to a refpoint)
  tracefun =  actobject->beams[j].spec_trace;
  dx0 = (double)conf->beam[actobject->beams[j].ID].offset.dx0;
  dx1 = (double)conf->beam[actobject->beams[j].ID].offset.dx1;

  // determine the range of y-values covered by that beam
  // (relative to a refpoint)
  dy0 = tracefun->func (dx0, tracefun->data);
  dy1 = tracefun->func (dx1, tracefun->data);

  // apply the positional corrections;
  // in case of the fluxcube model
  // they are needed to transform
  // the coordinates from the direct
  // image to the grism image system
  dx0 = dx0 + actdir->xy_off[j].x;
  dx1 = dx1 + actdir->xy_off[j].x;
  dy0 = dy0 + actdir->xy_off[j].y;
  dy1 = dy1 + actdir->xy_off[j].y;

  // using the corners of the direct image object as refpoints,
  // translate the x- and y-ranges to the real area subtended
  // bye that model beam on the CCD
  xmin = floor(MIN((double)actdir->ix_min + dx0, (double)actdir->ix_min + dx1)+0.5);
  xmax = floor(MAX((double)actdir->ix_max + dx0, (double)actdir->ix_max + dx1)+0.5);
  ymin = floor(MIN((double)actdir->iy_min + dy0, (double)actdir->iy_min + dy1)+0.5);
  ymax = floor(MAX((double)actdir->iy_max + dy0, (double)actdir->iy_max + dy1)+0.5);

  // check whether the area is completely outside the CCD
  if ((xmax < 0.0) || (ymax < 0.0))
    {
      //      fprintf(stdout, "Beam's too low: %i %i\n",actobject->ID, actobject->beams[j].ID);
      // if yes, release the memory and return NULL
      free(actspec);
      actspec = NULL;
      return actspec;
    }

  // check whether the area is completely outside the CCD
  if ((xmin > (double)(npixels.x-1)) || (ymin > (double)(npixels.y-1)))
    {
      //      fprintf(stdout, "Beam's too high: %i %i\n",actobject->ID, actobject->beams[j].ID);
      // if yes, release the memory and return NULL
      free(actspec);
      actspec = NULL;
      return actspec;
    }

  // check whether part of the are is outside;
  // correct if necessary
  xmin = MAX(xmin-1, 0.0);
  ymin = MAX(ymin-1, 0.0);

  // check whether part of the are is outside;
  // correct if necessary
  xmax = MIN(xmax+1, (double)(npixels.x-1.0));
  ymax = MIN(ymax+1, (double)(npixels.y-1.0));

  // transfer the ID's to the model beam
  actspec->objectID = actobject->ID;
  actspec->beamID   = actobject->beams[j].ID;

  // transfer the coo's of the starting point;
  // leave one pixel for the 'drizzling'
  actspec->model_ref.x = (int)xmin;
  actspec->model_ref.y = (int)ymin;

  // allocate space for the matrix;
  // give two pixels more on each side for 'drizzling';
  // set the matrix to 0.0 and give it to the model beam
  stamp = gsl_matrix_alloc((int)(xmax-xmin)+1, (int)(ymax-ymin)+1);
  if (stamp == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "dimension_beamspec:" " Could not allocate"
                 " memory for %i x %i GSL matrix", (int)(xmax-xmin)+1, (int)(ymax-ymin)+1);
  gsl_matrix_set_all(stamp, 0.0);
  actspec->model = stamp;

  // return the beam created
  return actspec;
}


/**
 * Function: fill_pixel_in_speed
 * The function determines and coadds the contribution of a single pixel in
 * the direct image area to the model spectrum.
 *
 * Parameters:
 * @param  actdir   - the direct object of the model spectrum
 * @param  acttrace - the tracedata of the model spectrum
 * @param  dpixel   - the coordinates of the modelled pixel
 * @param  eval     - the emission value of the source at the modelled pixel
 * @param  resp     - the sensitivity data
 * @param  actspec  - the structure for the model spectrum
 *
 * Returns:
 * @return 1        -
 */
int
fill_pixel_in_speed(const dirobject *actdir, const tracedata *acttrace,
                    const d_point dpixel, const spectrum *resp,
                    beamspec *actspec, const calib_function  *wl_calibration)
{

  double dx;
  double sens;
  double fval;
  double tmp1;
  //double tmp2, tmp3;

  double ddx, ddy;

  int ix, iy;

  int xstart, xend;
  int xact;
  int ipos;
  int nguess;

  // define the dx-range for which a pixel is modelled
  xstart = actspec->model_ref.x;
  xend   = actspec->model_ref.x + actspec->model->size1;

  // nguess is the approximate possition
  // to find a wavelength in the sensitivity table.
  // should speed upt things
  nguess=0;

  // go over the dx-range
  for (xact = xstart; xact <  xend; xact++)
    {

      // compute the actual dx-value
      dx = xact-dpixel.x;

      // get the index to find the data for the actual
      ipos = get_index_for_tracepoint(acttrace, dx);

      // in case that the dx-value is not in the tracedata,
      // something is wrong. exit
      if (ipos <0)
        {
    	  //if (wl_calibration->pr_range == NULL)
          //  aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
          //              "aXe_PETCONT: " "tracepoint is outside the prepared range!\n");
          //else
          continue;
        }

      // compute the y-position in the matrix of the beamspec
      iy = (int)floor(gsl_vector_get(acttrace->dy,ipos) + actdir->xy_off[actspec->beamID].y + 0.5) + (int)dpixel.y - actspec->model_ref.y;

      // the same quantity as a double
      ddy = gsl_vector_get(acttrace->dy,ipos) + actdir->xy_off[actspec->beamID].y + dpixel.y - actspec->model_ref.y;

      // if the y-position is outside, go to the next dx value
      if (iy < 0 || iy > (int)actspec->model->size2-1)
        continue;

      // in case that the actual wavelength is outside the range
      // of the sensitivity data, continue with the next dx-value
      if (gsl_vector_get(acttrace->lambda, ipos) < resp->lambdamin || gsl_vector_get(acttrace->lambda, ipos) > resp->lambdamax)
        continue;

      // compute the x-position in the matrix of the beamspec
      ix = (int)(gsl_vector_get(acttrace->dx,ipos) + actdir->xy_off[actspec->beamID].x) + (int)dpixel.x - actspec->model_ref.x;

      // the same quantity as double
      ddx = gsl_vector_get(acttrace->dx,ipos) + actdir->xy_off[actspec->beamID].x + dpixel.x - actspec->model_ref.x;

      // get the total flux value of the source at the wavelength of the actual dx-value
      fval = gsl_vector_get(acttrace->gvalue, ipos) * gsl_vector_get(acttrace->flux, ipos);


      // get the sensitivity at the wavelength of the actual dx-value
      sens = get_response_value_plus(resp, gsl_vector_get(acttrace->lambda, ipos), &nguess);

      // compute the contribution of the actual dx-value to the model spectrum
      tmp1 = fval * sens * gsl_vector_get(acttrace->dlambda,ipos);

      // double check whether we are inside the image
      if (ix < 0 || iy < 0 || ix > (int)(actspec->model->size1-1) || iy > (int)(actspec->model->size2-1))
        {
          fprintf(stdout, "xval: %i, yval: %i, size1: %zi, size2: %zi\n",ix, iy, actspec->model->size1, actspec->model->size2);
        }
      else
        {
    	  // add the contribution of the actual dx value to
    	  // the model spectrum, using NO diffusion
    	  // no_diffuse_spectrum(ix, iy, tmp1, actspec);

    	  // add the contribution of the actual dx value to
    	  // the model spectrum, using GRID diffusion
    	  // diffuse_spectrum(ddx, ddy, tmp1, actspec);

    	  // add the contribution of the actual dx value to
    	  // the model spectrum, using PROPORTIONAL diffusion
    	  diffuse_spectrumII(ddx, ddy, tmp1, actspec);
        }
    }

  // return a dummy
  return 1;
}


/**
 * Function: get_index_for_tracepoint
 * For a given dx-value, the function derives the index
 * under which data for this dx-value is stored in a tracedata structure.
 * The index is simply returned. If the dx-value is not inside
 * the tracedata, '-1' is returned.
 * In contrast to the first attempt (now 'get_index_for_tracepoint2')
 * this routine uses the minimum dx value stored in the tracedata
 * structure. It becomes much much faster by that!
 *
 * Parameters:
 * @param  acttrace - the tracedata to search an index
 * @param  dx       - the dx-value to identify in the tracedata
 *
 * Returns:
 * @return ipos     - the index where the dx-value is stored
 */
int
get_index_for_tracepoint(const tracedata *acttrace, const double dx)
{
  int ipos;

  // compute the index using the minimum dx value
  // stored in the structure
  ipos = (int)(dx - acttrace->dx_start);

  // check whether the index is within the acceptable range
  // make it -1 if not
  if (ipos < 0 || ipos > (acttrace->npoints -1))
    ipos=-1;

  // return the index
  return ipos;
}

/**
 * Function no_diffuse_spectrum
 * Adds the simulated flux for one pixel exactly into
 * one pixel, using the integer value of the correct
 * position.
 *
 * Parameters:
 * @param  ix      - the integer x-position in the beam
 * @param  iy      - the integer y-position in the beam
 * @param  cps     - the cps-value for one entire pixel
 * @param  actspec - the beam to model
 *
 * Returns:
 * @return 1.0     - return a dummy value
 */
int
no_diffuse_spectrum(int ix, int iy, double cps, beamspec *actspec)
{
  double tmp  = 0.0;

  // add the contribution of the actual dx value to the model spectrum
  tmp = gsl_matrix_get(actspec->model, ix, iy);

  // set the new value
  gsl_matrix_set(actspec->model, ix, iy, cps + tmp);

  // return a dummy
  return 1;
}

/**
 * Function diffuse_spectrum
 * Distribute the flux simulated for the area of one pixel
 * square in the emission model over a pixel square in
 * the beam model.
 * This method uses sub-steps starting from the exact
 * center of the light
 * This smoothes the modelled beam by avoiding integer
 * rounding effects.
 *
 * Parameters:
 * @param  ddx     - the exact x-position in the beam
 * @param  ddy     - the exact y-position in the beam
 * @param  cps     - the cps-value for one entire pixel
 * @param  actspec - the beam to model
 *
 * Returns:
 * @return 1.0     - return a dummy value
 */
int
diffuse_spectrum(double ddx, double ddy, double cps, beamspec *actspec)
{
  double step;
  double offset;
  double ncps;
  double oldvalue;

  double sum=0.0;

  double dx_act;
  double dy_act;

  int irange;
  int kk, ll;

  int ix, iy;

  // convert the number of steps to a local integer
  irange = (int)NDIFF;

  // compute the step size
  step = 1.0/(2.0*(double)NDIFF);

  // compute the initial offset
  offset = step/2.0;

  // compute the incremental flux
  ncps = cps * step * step;

  for (kk=-irange; kk < irange; kk++)
    {
      for (ll=-irange; ll < irange; ll++)
        {
          // determine the actual grid position
          dx_act = ddx + (double)kk * step + offset;
          dy_act = ddy + (double)ll * step + offset;

          // get the current pixel grid value
          ix = floor(dx_act + 0.5);
          iy = floor(dy_act + 0.5);

          // double check whether we are inside the image
          if (ix < 0 || iy < 0 || ix > (int)(actspec->model->size1-1) || iy > (int)(actspec->model->size2-1))
            continue;

          // add the contribution of the actual dx value to the model spectrum
          oldvalue = gsl_matrix_get(actspec->model, ix, iy);
          gsl_matrix_set(actspec->model, ix, iy, oldvalue + ncps);
          sum += ncps;
        }
    }

  // return a dummy
  return 1.0;
}

/**
 * Function diffuse_spectrumII
 * Distribute the flux simulated for the area of one pixel
 * square in the emission model over a pixel square in
 * the beam model.
 * Starting from the exact center of the light, the exact fraction
 * falling on the four affected pixels is comuted and then added
 * to the spectral beam.
 * This smoothes the modelled beam by avoiding integer
 * rounding effects.
 *
 * Parameters:
 * @param  ddx     - the exact x-position in the beam
 * @param  ddy     - the exact y-position in the beam
 * @param  cps     - the cps-value for one entire pixel
 * @param  actspec - the beam to model
 *
 * Returns:
 * @return 1.0     - return a dummy value
 */
int
diffuse_spectrumII(double ddx, double ddy, double cps, beamspec *actspec)
{
  int ix;
  int iy;
  int ix_rem;
  int iy_rem;

  double p=0.0;
  double q=0.0;
  double p_rem=0.0;
  double q_rem=0.0;

  double d_incr = 0.0;
  double oldvalue;

  double sum = 0.0;

  // compute the indices
  // of the pixels involved
  ix     = (int)floor(ddx);
  iy     = (int)floor(ddy);
  ix_rem = ix + 1;
  iy_rem = iy + 1;


  // compute the basic
  // quantities for
  // the increments;
  // it seems to be the wrong way,
  // however we are talking about
  // pixels, not coordinates...
  p_rem = ddx - floor(ddx);
  q_rem = ddy - floor(ddy);
  p     = 1.0 - p_rem;
  q     = 1.0 - q_rem;


  // increment the first quarter
  // double check whether we are inside the image
  if ( ! (ix < 0 || iy < 0 || ix > (int)(actspec->model->size1-1) || iy > (int)(actspec->model->size2-1)))
    {
      // compute the area
      d_incr = p * q;

      // add the fractional contribution to the model spectrum
      oldvalue = gsl_matrix_get(actspec->model, ix, iy);
      gsl_matrix_set(actspec->model, ix, iy, oldvalue + d_incr*cps);
      sum += d_incr;
    }

  // increment the second quarter
  // double check whether we are inside the image
  if ( ! (ix < 0 || iy_rem < 0 || ix > (int)(actspec->model->size1-1) || iy_rem > (int)(actspec->model->size2-1)))
    {
      // compute the area
      d_incr = p * q_rem;

      // add the fractional contribution to the model spectrum
      oldvalue = gsl_matrix_get(actspec->model, ix, iy_rem);
      gsl_matrix_set(actspec->model, ix, iy_rem, oldvalue + d_incr*cps);
      sum += d_incr;
    }

  // increment the third quarter
  // double check whether we are inside the image
  if ( ! (ix_rem < 0 || iy < 0 || ix_rem > (int)(actspec->model->size1-1) || iy > (int)(actspec->model->size2-1)))
    {
      // compute the area
      d_incr = p_rem * q;

      // add the fractional contribution to the model spectrum
      oldvalue = gsl_matrix_get(actspec->model, ix_rem, iy);
      gsl_matrix_set(actspec->model, ix_rem, iy, oldvalue + d_incr*cps);
      sum += d_incr;
    }

  // increment the fourth quarter
  // double check whether we are inside the image
  if ( ! (ix_rem < 0 || iy_rem < 0 || ix_rem > (int)(actspec->model->size1-1) || iy_rem > (int)(actspec->model->size2-1)))
    {
      // compute the area
      d_incr = p_rem * q_rem;

      // add the fractional contribution to the model spectrum
      oldvalue = gsl_matrix_get(actspec->model, ix_rem, iy_rem);
      gsl_matrix_set(actspec->model, ix_rem, iy_rem, oldvalue + d_incr*cps);
      sum += d_incr;
    }


  // return a dummy
  return 1.0;
}


/**
 *
 * Function: make_model_image
 * Sums up the modeled spectral beams to create a
 * spectral model for the whole image.
 *
 * Parameters:
 * @param  npixels    - the dimensions of the input grism image
 * @param  speclist   - the list of modelled beams
 *
 * Returns:
 * @return all_models - the model for the whole image
 */
gsl_matrix *
make_model_image(const px_point npixels, observation *obs, beamspec **speclist)
{

  gsl_matrix *all_models;

  double oldval, addval;

  int i=0;
  int xact, yact;
  int ix, iy;

  // allocate space for the result,
  // set the matrix to 0.0
  //  all_models = gsl_matrix_alloc(npixels.x, npixels.y);
  all_models = obs->grism;
  gsl_matrix_set_all(all_models,0.0);

  fprintf (stdout, "\naXe_PETCONT: Summing up the model beam spectra\n");

  // go over each beam in the list
   while (speclist[i] != NULL)
     {

       fprintf(stdout, "aXe_PETCONT: summing up object %i beam %c ...", speclist[i]->objectID, BEAM(speclist[i]->beamID));

       // go over each pixel in the array of the beam model
       for (xact=0; xact < (int)speclist[i]->model->size1; xact++)
         {
           for (yact=0; yact < (int)speclist[i]->model->size2; yact++)
             {

               // find the coordinates of the pixels in the whole image model
               ix = speclist[i]->model_ref.x + xact;
               iy = speclist[i]->model_ref.y + yact;

               // check for safety reasons whether the coordinates are inside
               if (ix < 0 || iy < 0 || ix > (npixels.x-1) || ix > (npixels.x-1)){
                 fprintf(stdout, "This should not happen!\n");
               }
               else{

                 // summ up the pixel
                 addval = gsl_matrix_get(speclist[i]->model, xact, yact);
                 oldval = gsl_matrix_get(all_models, ix, iy);
                 gsl_matrix_set(all_models, ix, iy, oldval+addval);
               }
             }
         }

       fprintf(stdout, " Done\n");

       // increment the counter
       i++;
     }

   // return the model
  return all_models;
}

/**
 * Function: fill_contam_info
 * The function fills the contamination information into the
 * PET. To do that it transfers the information from the
 * aperture mask matrix into each PET. Self contamination is
 * taken into account by subtracting the value from the modelled
 * spectrum before storing the contamination.
 *
 * Parameters:
 * @param PET_file   - the PET file to add contamination
 * @param speclist   - the list of modelled spectra
 * @param all_models - the contamination image
 *
 * Returns:
 * @return 1         - returns always 1
 */
int
fill_contam_info(const char PET_file[], beamspec **speclist,
                 const gsl_matrix *all_models, char model_name[])
{
  fitsfile *OPET_ptr;
  ap_pixel *PET;
  beamspec *actspec;
  //FITScards *cards;

  char ID[60];

  double c;
  double m;

  int aperID, beamID;
  int f_status=0;
  int status;
  int ix, iy;
  int j;

  // Open the OPET file for reading/writing
  ffopen (&OPET_ptr, PET_file, READWRITE, &f_status);
  if (f_status){
    ffrprt(stdout, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "aXe_PETCONT: Could not open file: %s\n",
                 PET_file);
  }

  // store the contamination model name in the PET
  update_contam_model(OPET_ptr, model_name);

  // report the action
  fprintf (stdout, "\naXe_PETCONT: Writing the contamination into the PET.\n");

  while (1)
    {

      // get the PET at the actual position of the PET file
      PET = get_ALL_from_next_in_PET(OPET_ptr, &aperID, &beamID);

      // leave the while loop if there is the end
      if ((aperID==-1) && (beamID==-1))
        break;

      // report progress
      fprintf(stdout, "aXe_PETCONT: writing contamination for object %i beam %c ... ", aperID, BEAM(beamID));

      // skip the PET if empty
      if (PET==NULL)
        {
          fprintf (stdout, ".Done\n");
          continue;
        }

      // get the model spectrum for the actual PET
      actspec = get_beamspec_from_list(speclist, aperID, beamID);

      // go over all PET entries
      j=0;
      while  (PET[j].p_x != -1)
        {

          // get the value from the contamination matrix
          c = gsl_matrix_get(all_models,PET[j].p_x,PET[j].p_y);
          m = 0.0;

          // correct for self-contamination if the beams was modelled
          if (actspec != NULL){

            // get the coordinates in the model array
            ix = PET[j].p_x - actspec->model_ref.x;
            iy = PET[j].p_y - actspec->model_ref.y;

            // check whether the PET entry is outside of the model array.
            // correct for self contamination if the PET-etry is inside

            if (ix > -1 && iy >-1 && ix < (int)actspec->model->size1 && iy < (int)actspec->model->size2)
              {
                m = gsl_matrix_get(actspec->model, ix, iy);
                c = c - m;
              }
          }

          // well, it should never be below zero
          if (c < 0.0)
            c = 0.0;

          // write the contamination into the PET-entry
          PET[j].contam = c;
          PET[j].model = m;

          // enhance the counter
          j++;
        }

      // write the updated PET into the PET file
      sprintf (ID, "%d%c", aperID, BEAM (beamID));
      add_ALL_to_PET (PET, ID, OPET_ptr,1);

      // free PET memory
      if (PET!=NULL)
        {
          free(PET);
          PET=NULL;
        }
      fprintf (stdout, ".Done\n");
    }

  // close the PET file
  fits_close_file (OPET_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "aXe_PETCONT: " "Error closing PET: %s \n",
                   PET_file);
    }

  // set the status to 'success' and return it
  status=1;
  return status;
}


/*-------------------
 * Section 4: Geometrical contamination
 */
/**
 * Function: compute_geometr_cont
 * The subroutine computes the original, geometrical contamination
 * model where every beams contaminates the area he occupies.
 * Source brightness and source shape are not taken into account.
 * For every beam the sum of all contaminations by othe beams is
 * written into the PET.
 *
 * Parameters:
 * @param OAF_file - the name of the aperture file
 * @param obs      - the observation
 * @param PET_file - the name of the PET which is modified
 * @param map_file - the name of the contamination map
 * @param store    - flagg whether the contamination image is stored or not
 *
 * Returns:
 * @return status  - returns success or failure
 */
int
compute_geometr_cont(char OAF_file[], observation *obs,
                     const char PET_file[], char map_file[],
                     const int store){

  object **oblist;
  ap_pixel *PET;
  int f_status=0;
  fitsfile *OPET_ptr;
  int aperID, beamID, objindex;
  gsl_matrix *aper_mask;
  int status=0;
  char model_name[60];

  //  load the object list
  fprintf (stdout, "aXe_PETCONT: Loading object aperture list...");
  oblist = file_to_object_list_seq (OAF_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // Open the OPET file for reading/writing
  ffopen (&OPET_ptr, PET_file, READWRITE, &f_status);
  if (f_status)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "aXe_PETCONT: Could not open file: %s\n",
                 PET_file);

  // store the contamination model name in the PET
  sprintf (model_name, "GEOM");
  update_contam_model(OPET_ptr, model_name);

  // do something in case that there are objects
  if (oblist!=NULL)
    {

      // make the aperture mask
      aper_mask = aperture_mask(obs,oblist);

      // store the aperture mask into a fits file if requested
      if (store)
        gsl_to_FITSimage (aper_mask, map_file, 1, NULL);

      // loop over all beams
      while (1)
        {

          // get the PET at the actual position of the PET file
          PET = get_ALL_from_next_in_PET(OPET_ptr, &aperID, &beamID);

          // leave the while loop if there is the end
          if ((aperID==-1) && (beamID==-1))
            break;

          // search for the particular beam in 'oblist', continue
          // if PET is empty
          fprintf (stdout, "aXe_PETCONT: BEAM %d%c", aperID, BEAM(beamID));
          objindex =  find_object_in_object_list(oblist,aperID);
          if (PET==NULL)
            {
              fprintf (stdout, ".Done\n");
              continue;
            }

          // transfer the information from the aperture mask
          // into the beam
          {
            int j = 0;
            double c;

            // go over all PET entries
            while  (PET[j].p_x != -1)
              {
                // subtract self contamination, then store the contamination
                c = gsl_matrix_get(aper_mask,PET[j].p_x,PET[j].p_y)-1.0;
                if (c < 0.0)
                  c = 0.0;
                PET[j].contam = c;
                j++;
              }
          }

          // write the updated PET into the PET file
          {
            char ID[60];
            sprintf (ID, "%d%c", oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));
            add_ALL_to_PET (PET, ID, OPET_ptr,1);
          }

          // free PET memory
          if (PET!=NULL)
            {
              free(PET);
              PET=NULL;
            }
          fprintf (stdout, ".Done\n");
        }
    }

  // free obl;ist memory
  if (oblist!=NULL)
    free_oblist (oblist);

  // free observation memory
  free_observation(obs);

  // close the PET file
  fits_close_file (OPET_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "aXe_PETCONT: " "Error closing PET: %s \n",
                   PET_file);
    }

  // set the status to 'success' and return it
  status=1;
  return status;
}


/*-------------------
 * Section 5: Simualting dirdct images with Gaussians
 */
/**
 * Function: compute_gauss_dirim
 *
 * @param grism_file  - the full name of the grism file
 * @param OAF_file    - the name of the aperture file
 * @param obs         - the observation
 * @param PET_file    - the name of the PET which is modified
 * @param CONF_file   - the full name of configuration file
 * @param map_file    - the name of the contamination map
 * @param model_scale - the scale for extension of the direct object area
 * @param store       - flagg whether the contamination image is stored or not
 *
 * @return status     - returns success or failure
 */
int
compute_gauss_dirim(char grism_file[], char OAF_file[], char CONF_file[],
                    const double model_scale,  const int inter_type,
                    const double lambda_psf, observation *obs,
                    const char PET_file[], char map_file[], const int store){

  object    **oblist;
  dirobject **dirlist;

  //dirobject  *actdir;
  gsl_matrix *all_models;

  //char model_name[60];

  px_point npixels;

  //  load the object list
  fprintf (stdout, "aXe_PETCONT: Loading object aperture list...");
  oblist = file_to_object_list_seq (OAF_file, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  npixels = get_npixel(obs);

  dirlist = oblist_to_dirlist(grism_file, CONF_file, npixels, oblist, NULL, model_scale, inter_type);

  all_models = make_gauss_dirim(oblist, dirlist, lambda_psf, npixels, CONF_file, obs);


  if (store)
    gsl_to_FITSimage (all_models, map_file, 1, NULL);


  gsl_matrix_free(all_models);
  free_dirlist (dirlist);
  if (oblist != NULL)
    free_oblist (oblist);

  return 1;
}
/**
 * Function: make_gauss_dirim
 * The function creates a spectral model for the beams
 * that are considered in the contamination. The modelled
 * beams are returned as a list.
 *
 * Parameters:
 * @param  oblist    - the object list as input to select beams
 * @param  dirlist   - the direct object list to dimension the models
 * @param  npixels   - the dimensions of the model for the whole image
 * @param  CONF_file - the name of the configuration file
 *
 * Returns:
 * @return speclist  - the list of modelled beams
 */
gsl_matrix *
make_gauss_dirim(object **oblist, dirobject **dirlist,
                 const double lambda_psf, const px_point npixels,
                 char CONF_file[], observation *obs)
{
  gsl_matrix *all_models;
  dirobject       *actdir;
  beam             actbeam;

  //int nspecs;
  //int i=0;
  //int j=0;
  //int jj=0;
  int ii=0;
  int nobjects;

  //int kk;
  //int ll;
  double sval=0.0;
  double flux=0.0;
  double value=0.0;

  int nx, ny;
  d_point dpixel;

  // determine the number of objects in the object list
  nobjects = object_list_size(oblist);

  // allocate space for the result,
  // set the matrix to 0.0
  //  all_models = gsl_matrix_alloc(npixels.x, npixels.y);
  all_models = obs->grism;
  gsl_matrix_set_all(all_models,0.0);

  // go over each beam model
  ii = 0;
  while (oblist[ii] != NULL)
    {

      // get the direct object for the actual model spectrum
      actdir = get_dirobject_from_list(dirlist, oblist[ii]->ID);

      // get the beam for the actual model spectrum
      actbeam = oblist[ii]->beams[0];

      fprintf(stdout, "aXe_PETCONT: modelling object %i ...", oblist[ii]->ID);

      // go over each pixel in the direct object area
      for (nx=actdir->ix_min; nx<=actdir->ix_max; nx++)
        {
          for (ny=actdir->iy_min; ny<=actdir->iy_max; ny++)
            {

              // fill the dpixel structure
              dpixel.x = (double)nx;
              dpixel.y = (double)ny;

              sval = get_sub_emodel_value(dpixel, actbeam, actdir->drzscale);
              flux = get_flux_from_SED(actdir->SED, 890.0);

              if (nx > -1 && ny > -1 && nx < npixels.x && ny < npixels.y)
                {
                  value = gsl_matrix_get(all_models, nx, ny) + sval*flux;
                  gsl_matrix_set(all_models, nx, ny, value);
                }
            }
        }

      // release the space for the various structures
      fprintf(stdout, " Done\n");

      ii++;
    }

  return all_models;
}
