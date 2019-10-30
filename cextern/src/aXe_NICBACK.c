/*
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_sect.h"
#include "inout_aper.h"
#include "aper_conf.h"
#include "spc_back.h"
#include "nicback_utils.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

int
main(int argc, char *argv[])
{
  char *opt;
  char  grism_image[MAXCHAR];
  char  grism_image_path[MAXCHAR];

  char  conf_file[MAXCHAR];
  char  conf_file_path[MAXCHAR];

  char  master_bck[MAXCHAR];
  char  master_bck_path[MAXCHAR];

  char  corr_bck[MAXCHAR];
  char  corr_bck_path[MAXCHAR];

  char  msk_name[MAXCHAR];
  char  msk_name_path[MAXCHAR];

  char  plist_name[MAXCHAR];
  char  plist_name_path[MAXCHAR];

  char  bck_name[MAXCHAR];
  char  bck_name_path[MAXCHAR];

  //  int   np;
  //  int   interp;

  //  int   makemask=0;
  //  int   nor_flag=0;

  //  object        **oblist;
  observation    *obs;
  aperture_conf  *conf;
  //  background     *backg;

  //  int niter_med;
  //  int niter_fit;

  //  double kappa;
  double exptime;

  if ((argc <= 3) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
	    "aXe_NICBACK Version %s:\n"
	    "\n"
	    "Example:   aXe_BE slim_grismb.fits -np=10 -interp=3 \n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_NICKBACK: Starting...\n");

  // Getting the grism/prism image name
  strcpy(grism_image, argv[1]);
  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);

  // Get the name of the configuration file
  strcpy(conf_file, argv[2]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);

  // Getting the master background name
  strcpy(master_bck, argv[3]);
  build_path(AXE_CONFIG_PATH, master_bck, master_bck_path);

  // check whether a pedestal file was given
  if (argc > 4)
    {
      // transfer file name
      strcpy(corr_bck, argv[4]);
      build_path(AXE_CONFIG_PATH, corr_bck, corr_bck_path);
    }
  else
    {
      // give a default
      strcpy(corr_bck, "");
      strcpy(corr_bck_path, "");
    }

  // read the configuration file
  conf = get_aperture_descriptor(conf_file_path);
  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);

  // fix the mask file name
  replace_file_extension(grism_image, msk_name, ".fits",
			 ".MSK.fits", conf->science_numext);
  build_path(AXE_OUTPUT_PATH, msk_name, msk_name_path);

  // fix the name of the pixel list
  replace_file_extension(grism_image, plist_name, ".fits",
			 ".plis", conf->science_numext);
  build_path(AXE_OUTPUT_PATH, plist_name, plist_name_path);

  // fix the mask file name
  replace_file_extension(grism_image, bck_name, ".fits",
			 ".NBCK.fits", conf->science_numext);
  build_path(AXE_OUTPUT_PATH, bck_name, bck_name_path);


  fprintf(stdout,
	  "aXe_NICBACK: Input data file name:                 %s\n",
	  grism_image_path);
  fprintf(stdout,
	  "aXe_NICBACK: Main configuration file name:         %s\n",
	  conf_file_path);
  fprintf(stdout,
	  "aXe_NICBACK: Master background file name:          %s\n",
	  master_bck_path);
  if (strlen(corr_bck_path) > 0)
    fprintf(stdout,
        "aXe_NICBACK: Correction background file name:      %s\n",
        corr_bck_path);
  fprintf(stdout,
	  "aXe_NICBACK: Grism mask file name:                 %s\n",
	  msk_name_path);
  fprintf(stdout,
	  "aXe_NICBACK: Background image file name:           %s\n",
	  bck_name_path);
  fprintf(stdout, "\n\n");

  fprintf(stdout, "aXe_NICBACK: ");

  // try to get the descriptor 'exptime' from the 'sci'-extension
  exptime = (double)get_float_from_keyword(grism_image_path,
					   conf->science_numext, conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_image_path, 1, conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;

  // load the grism image in an obs structure
  obs = load_image_t(grism_image_path, conf->science_numext,
		   conf->errors_numext, conf->dq_numext,
		   conf->dqmask, exptime, conf->rdnoise);

  // run the subroutine which does everything
  make_nicmos_back(obs, msk_name_path, master_bck_path, bck_name_path,
		   plist_name_path, corr_bck_path);

  // free the memory in the ob structure
  free_observation(obs);
  fprintf(stdout, "aXe_NICBACK: Done...\n");
  exit(0);
}
