/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "aper_conf.h"
#include "spce_output.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"
#define AXE_SIMDATA_PATH "AXE_SIMDATA_PATH"

extern int
compute_dirimage_model(char *, char *, char *,
                       char *, char *, char *,
                       const double , const double , const double ,
                       observation *, char *);

int
main (int argc, char *argv[])
{
  char dirim_file[MAXCHAR];
  char dirim_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char tpass_file[MAXCHAR];
  char tpass_file_path[MAXCHAR];

  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char specmod_file[MAXCHAR];
  char specmod_file_path[MAXCHAR];

  char objmod_file[MAXCHAR];
  char objmod_file_path[MAXCHAR];

  char map_file[MAXCHAR];
  char map_file_path[MAXCHAR];

  char *opt;
  aperture_conf *conf;
  observation *obs;

  int index;
  double model_scale=0.0;
  double lambda_psf=0.0;
  double tel_area=0.0;

  FITScards      *cards;

  if ( ((argc < 3) || (opt = get_online_option ("help", argc, argv))) )
    {
      fprintf (stdout,
	       "ST-ECF European Coordinating Facility\n"
	       "aXe_DIRIMAGE Version %s:\n"
	       "\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_DIRIMAGE: Starting...\n");

  // initialize an index for the command arguments
  index = 0;

  // save the first argument as direct image name
  strcpy (dirim_file, argv[++index]);
  build_path (AXE_IMAGE_PATH, dirim_file, dirim_file_path);

  // save the second argument as aXe config file name
  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  // save the third argument as total throughput file
  strcpy (tpass_file, argv[++index]);
  build_path (AXE_SIMDATA_PATH, tpass_file, tpass_file_path);

  // load the aXe configuration file
  conf = get_aperture_descriptor (conf_file_path);

  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(dirim_file_path, conf,conf->optkey1,conf->optval1);

  /* Get or set up the name of the output Aperture File */
  if ((opt = get_online_option ("in_AF", argc, argv)))
    {
      /* get it */
      strcpy (aper_file, opt);
      strcpy (aper_file_path, opt);
    }
  else {
    /* Build aperture file name */
    replace_file_extension (dirim_file, aper_file, ".fits",
			    ".OAF", conf->science_numext);
    build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
  }

  // check whether a name for the spectral
  // models file is given
  if ( (opt = get_online_option ("model_spectra", argc, argv)) )
    {
      // get and set up the filename
      strcpy (specmod_file, opt);
      build_path (AXE_IMAGE_PATH, specmod_file, specmod_file_path);
    }
  else
    {
      // set the filenames to NULL,
      // indicating that they are NOT  used
      strcpy (specmod_file, "");
      strcpy (specmod_file_path, "");
    }

  if ( (opt = get_online_option ("model_images", argc, argv)) )
    {
      // get and set up the name for the image templates
      strcpy (objmod_file, opt);
      build_path (AXE_IMAGE_PATH, objmod_file, objmod_file_path);
    }
  else
    {
      // set the filenames to NULL,
      // indicating that they are NOT  used
      strcpy (objmod_file, "");
      strcpy (objmod_file_path, "");
    }


  /* Build object CONT file name */
  replace_file_extension (dirim_file, map_file, ".fits",
			  ".CONT.fits", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, map_file, map_file_path);

  // determine the extend of the gaussian emission model
  if ( (opt = get_online_option ("model_scale", argc, argv)) )
    model_scale = atof(opt);
  else
    model_scale = 3.0;

  // determine the telescope area
  if ( (opt = get_online_option ("tel_area", argc, argv)) )
    tel_area = atof(opt);
  else
    // set the default HST collecting area
    tel_area = 45238.93;

  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf (stdout, "aXe_DIRIMAGE: Input Aperture file name:            %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_DIRIMAGE: Output CONT file name:               %s\n",
	   map_file_path);
  fprintf (stdout, "aXe_DIRIMAGE: Computing gaussian contamination\n");
  fprintf (stdout, "aXe_DIRIMAGE: Scale factor for the direct image model %f\n", model_scale);
  if ( (strlen(specmod_file_path) > 0) )
    fprintf (stdout, "aXe_DIRIMAGE: Using spectral models in table: %s\n", specmod_file_path);
  if ( (strlen(objmod_file_path) > 0) )
    fprintf (stdout, "aXe_DIRIMAGE: Using direct emission objects in image:: %s\n", specmod_file_path);


  fprintf (stdout, "aXe_DIRIMAGE: ");
  obs = load_sci_image (dirim_file_path, conf->science_numext);


  compute_dirimage_model(dirim_file_path, conf_file_path, tpass_file_path, specmod_file_path,
			 objmod_file_path, aper_file_path, model_scale, tel_area,
			 lambda_psf, obs, map_file_path);

  // copy the header from the grism image
  // to the contamination image
  cards = get_FITS_cards(dirim_file_path, conf->science_numext);
  put_FITS_cards(map_file_path, 2, cards);
  free_FITScards(cards);

  // de-allocate memory
  //free_observation(obs);
  free_aperture_conf(conf);

  // free the observation structure
  free_observation(obs);

  // report on the program end
  fprintf (stdout, "aXe_DIRIMAGE: Done...\n");

  exit (0);
}
