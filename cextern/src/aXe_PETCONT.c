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
#include "spc_model.h"
#include "model_utils.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

int
main (int argc, char *argv[])
{
  char *opt;
  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char grism_file[MAXCHAR];
  char grism_file_path[MAXCHAR];

  char specmod_file[MAXCHAR];
  char specmod_file_path[MAXCHAR];

  char objmod_file[MAXCHAR];
  char objmod_file_path[MAXCHAR];

  char fcube_file[MAXCHAR];
  char fcube_file_path[MAXCHAR];

  char map_file[MAXCHAR];
  char map_file_path[MAXCHAR];

  char PET_file[MAXCHAR];
  char PET_file_path[MAXCHAR];

  aperture_conf *conf;

  observation *obs;

  int index;
  int store_map=0;
  int cont_model=0;
  double model_scale=0.0;
  double lambda_psf=0.0;
  int inter_type=1;

  FITScards      *cards;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "aXe_PETCONT\n"
	       "           Task that populates the CONTAMINATION part of a\n"
	       "           Pixel Extraction Table (PET). Each pixel is examined\n"
	       "           and the number of aperture minus 1 in which they\n"
	       "           appear is written into the PET\n"
	       "\n"
	       "Usage:\n"
	       "      aXe_PETCONT g/prism_filename configuration_filename [option]\n"
	       "\n"
	       "Options:\n"
	       "      -in_AF=[string] - overwrite the automatically generated name\n"
	       "                        of the input aperture file\n"
	       "      -cont_map       - write the contamination map into a FITS file"
	       "\n");
      exit (1);
    }

  fprintf (stdout, "aXe_PETCONT: Starting...\n");

  index = 0;
  strcpy (grism_file, argv[++index]);
  build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);


  conf = get_aperture_descriptor (conf_file_path);

  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);

  /* Get or set up the name of the output Aperture File */
  if ((opt = get_online_option ("in_AF", argc, argv)))
    {
      /* get it */
      strcpy (aper_file, opt);
      strcpy (aper_file_path, opt);
    }
  else {
    /* Build aperture file name */
    replace_file_extension (grism_file, aper_file, ".fits",
			    ".OAF", conf->science_numext);
    build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
  }

  // check whether a name for the spectral
  // models file is given
  if ((opt = get_online_option ("model_spectra", argc, argv)))
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

  // determine the contamination model
  if ((opt = get_online_option ("cont_model", argc, argv)))
    cont_model = atoi(opt);
  else
    cont_model = 1;

  // check whether the direct emission model was chosen and
  // a name for the object models file is given
  if (cont_model == 2)
    {
      if ( (opt = get_online_option ("model_images", argc, argv)) )
	{
	  // get and set up the name for the image templates
	  strcpy (objmod_file, opt);
	  build_path (AXE_IMAGE_PATH, objmod_file, objmod_file_path);
	}
      else
	{
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "For the direct image contamination a file with\nobject templates must be given!\n");
	}
    }
  else
    {
      // set the filenames to NULL,
      // indicating that they are NOT  used
      strcpy (objmod_file, "");
      strcpy (objmod_file_path, "");
    }


  // check whether a PET is there
  if ((get_online_option ("noPET", argc, argv)))
    {
      // set the PET filenames to NULL,
      // indicating that they are NOT  used
      strcpy (PET_file, "");
      strcpy (PET_file_path, "");
    }
      else
    {
      /* Build object PET file name */
      replace_file_extension (grism_file, PET_file, ".fits",
			      ".PET.fits", conf->science_numext);
      build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
    }

  /* Build object CONT file name */
  replace_file_extension (grism_file, map_file, ".fits",
			  ".CONT.fits", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, map_file, map_file_path);

  // determine the extend of the gaussian emission model
  if ((opt = get_online_option ("model_scale", argc, argv)))
    model_scale = atof(opt);
  else
    model_scale = 3.0;

  // determine the wavelength
  // the object extend was determined at
  if ((opt = get_online_option ("lambda_psf", argc, argv)))
    lambda_psf = atof(opt);
  else
    lambda_psf = 800.0;

  // determine the interpolation type
  if ((opt = get_online_option ("inter_type", argc, argv)))
    inter_type = atoi(opt);
  else
    inter_type = 1;

  // store the contamination map?
  if ((get_online_option ("cont_map", argc, argv)))
    store_map=1;
  else
    store_map=0;


  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf (stdout, "aXe_PETCONT: Input Aperture file name:            %s\n",
	   aper_file_path);
  if (strlen(PET_file_path) > 0)
    fprintf (stdout, "aXe_PETCONT: Input PET file name:                 %s\n",
	   PET_file_path);
  if (store_map)
      fprintf (stdout, "aXe_PETCONT: Output CONT file name:               %s\n",
	       map_file_path);
  if (cont_model ==1)
    {
      fprintf (stdout, "aXe_PETCONT: Computing gaussian contamination\n");
      fprintf (stdout, "aXe_PETCONT: Scale factor for the direct image model %f\n", model_scale);
      fprintf (stdout, "aXe_PETCONT: Object parameters determined at %fnm\n", lambda_psf);
      if (strlen(specmod_file_path) > 0)
	fprintf (stdout, "aXe_PETCONT: Using spectral models in table: %s\n", specmod_file_path);
    }
  else if (cont_model ==2)
    {
      fprintf (stdout, "aXe_PETCONT: Computing direct emission object contamination.\n");
      fprintf (stdout, "aXe_PETCONT: Using direct emission objects in image: %s\n", objmod_file_path);
      if (strlen(specmod_file_path) > 0)
	fprintf (stdout, "aXe_PETCONT: Using spectral models in table: %s\n", specmod_file_path);
    }
  else if (cont_model ==3)
    {
      fprintf (stdout, "aXe_PETCONT: Computing fluxcube contamination.\n");
    }
  else if (cont_model ==4)
    {
      fprintf (stdout, "aXe_PETCONT: Computing geometrical contamination.\n");
    }
  if (cont_model == 1 || cont_model == 3)
    {
      if (inter_type == 1)
	fprintf (stdout, "aXe_PETCONT: Using linear flux interpolation.\n");
      else if (inter_type == 2)
	fprintf (stdout, "aXe_PETCONT: Using polynomial flux interpolation.\n");
      else if (inter_type == 3)
	fprintf (stdout, "aXe_PETCONT: Using spline flux interpolation.\n");
    }

  fprintf (stdout, "aXe_PETCONT: ");
  obs = load_sci_image (grism_file_path, conf->science_numext);


  if (cont_model == 1 || cont_model == 2)
    {
      compute_gaussdirim_cont(grism_file_path, aper_file_path, conf_file_path,
			      specmod_file_path, objmod_file_path, model_scale, inter_type,
			      lambda_psf, obs, PET_file_path, map_file_path, store_map);
      //      compute_gauss_cont(grism_file_path, aper_file_path, conf_file_path,
      //			 specmod_file_path, model_scale, inter_type,
      //			 lambda_psf, obs, PET_file_path, map_file_path, store_map);
      //      compute_gauss_dirim(grism_file_path, aper_file_path, conf_file_path,
      //			 model_scale, inter_type, lambda_psf, obs, PET_file_path, map_file_path,
      //			 store_map);
    }
  else if (cont_model == 3)
    {
      replace_file_extension (grism_file, fcube_file, ".fits",
			      ".FLX.fits", conf->science_numext);
      build_path (AXE_IMAGE_PATH, fcube_file, fcube_file_path);

      compute_fcube_cont(grism_file_path, aper_file_path, fcube_file_path,
			 conf_file_path, model_scale, inter_type, obs,
			 PET_file_path, map_file_path, store_map);

    }
  else if (cont_model == 4)
    {
      compute_geometr_cont(aper_file_path, obs, PET_file_path,
			   map_file_path, store_map);
    }

  // copy the header from the grism image
  // to the contamination image
  cards = get_FITS_cards(grism_file_path, conf->science_numext);
  put_FITS_cards(map_file_path, 2, cards);
  free_FITScards(cards);

  //free_observation(obs);
  free_aperture_conf(conf);
  fprintf (stdout, "aXe_PETCONT: Done...\n");
  exit (0);
}
