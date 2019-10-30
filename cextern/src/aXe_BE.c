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

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

int
main(int argc, char *argv[])
{
  char *opt;
  char  grism_image[MAXCHAR];
  char  grism_image_path[MAXCHAR];

  char  aper_file[MAXCHAR];
  char  aper_file_path[MAXCHAR];

  char  back_file[MAXCHAR];
  char  back_file_path[MAXCHAR];

  char  conf_file[MAXCHAR];
  char  conf_file_path[MAXCHAR];

  int   np;
  int   interp;

  int   makemask=0;
  int   nor_flag=0;

  object        **oblist;
  observation    *obs;
  aperture_conf  *conf;
  background     *backg;

  int niter_med;
  int niter_fit;

  double kappa;
  double exptime;

  int sm_length=0;
  double fwhm=0.0;

  if ((argc <= 2) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
            "ST-ECF European Coordinating Facility\n"
            "aXe_BE Version %s:\n"
            "           aXe task to compute an estimate of the p/grism image background\n"
            "           or to create a mask file for background subtraction (option -msk).\n"
            "           This tasks uses the existing p/grism image and an existing\n"
            "           aperture file defining the position, orientation, and width\n"
            "           of all the beams (orders) that are believed to be in the image.\n"
            "           The values in the regions within each of these beams (orders) are\n"
            "           replaced by the median, average, linear, or n^th order\n"
            "           interpolation of pixels which are immediately above and below a\n"
            "           beam (but not within any other beam). The number of pixels to use\n"
            "           for this is by default se to 10 both below and above each aperture.\n"
            "           The -np option can be used to change this default value.\n"
            "           If the number of points is set to a value which is 0 or less, then\n"
            "           the entire column of an image will be used, ignoring any pixel\n"
            "           which are within any known beam. This option allows for a global\n"
            "           background estimate to be created instead of a local background\n"
            "           estimate.\n\n"
            "           The type of interpolation is controlled by the -interp option:\n"
            "                         -interp= -1    ; Median\n"
            "                         -interp= 0     ; Average\n"
            "                         -interp= 1     ; Linear fit\n"
            "                         -interp= (n>1) ; n^th order polynomial fit\n"
            "\n"
            "           The output of this task is a FITS image containing two extensions:\n"
            "           A SCI extension containing the actual image of the background\n"
            "           A ERR extension containing an estimate of the error in the fit\n"
            "\n"
            "             Input FITS mages are looked for in $AXE_IMAGE_PATH\n"
            "             aXe config file is looked for in $AXE_CONFIG_PATH\n"
            "             All outputs are writen to $AXE_OUTPUT_PATH\n"
            "\n"
            "Usage:\n"
            "      aXe_BE [g/prism image filename] [aXe config filename] [options]\n"
            "\n"
            "Options:\n"
            "           -msk                - create a mask image with the aperture area.\n"
            "                                 marked. Those masks are an input for the task\n"
            "                                 aXe_PREPARE.\n"
            "           -SCI_hdu=[integer]  - overwrite the default from the aXe config file\n"
            "           -ERR_hdu=[integer]  - overwrite the default from the aXe config file\n"
            "           -DQ_hdu=[integer]   - overwrite the default from the aXe config file\n"
            "           -np=[integer]       - The number of points to attempt to use to \n"
            "                                 compute the medianed/averaged/fitted\n"
            "                                 background\n"
            "           -interp=[integer]   - The type of interpolation to perform\n"
            "           -out_BCK=[string]   - Overwrite the default output background\n"
            "                                 filename\n"
            "           -in_AF=[string]     - Overwrites the default input aperture\n"
            "                                 filename\n"
            "\n"
            "Example:   aXe_BE slim_grismb.fits -np=10 -interp=3 \n"
            "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_BE: Starting...\n");

  /* Getting the grism/prism image name */
  strcpy(grism_image, argv[1]);
  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);

  /* Get the name of the configuration file */
  strcpy(conf_file, argv[2]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);

  /* Read the configuration file */
  conf = get_aperture_descriptor(conf_file_path);
  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);

  /* Online parameter overwrite all */
  if ((opt = get_online_option("SCI_hdu", argc, argv))) {
    conf->science_numext = atoi(opt);
  }
  if ((opt = get_online_option("ERR_hdu", argc, argv))) {
    conf->errors_numext = atoi(opt);
  }
  if ((opt = get_online_option("DQ_hdu", argc, argv))) {
    conf->dq_numext = atoi(opt);
  }
  if ((opt = get_online_option("msk", argc, argv)))
    makemask=1;
  if ((opt = get_online_option("nor_flag", argc, argv)))
    nor_flag=1;

  /* Get or set up the name of the background output file */
  if ((opt = get_online_option("out_BCK", argc, argv))) {
    strcpy(back_file, opt);
    strcpy(back_file_path, opt);
  } else {
    if (makemask){
      replace_file_extension(grism_image, back_file, ".fits",
                             ".MSK.fits", conf->science_numext);
      build_path(AXE_OUTPUT_PATH, back_file, back_file_path);
    } else {
      replace_file_extension(grism_image, back_file, ".fits",
                             ".BCK.fits", conf->science_numext);
      build_path(AXE_OUTPUT_PATH, back_file, back_file_path);
    }
  }

  /* Get the number of points to use to compute the background */
  if ((opt = get_online_option("np", argc, argv))) {
    {
      char            str[MAXCHAR];
      strcpy(str, opt);
      sscanf(str, "%d", &np);
    }
  } else {
    np = 0;
  }

  /* Get the type of interpolation to perform */
  if ((opt = get_online_option("interp", argc, argv)))
    {
      {
        char            str[MAXCHAR];
        strcpy(str, opt);
        sscanf(str, "%d", &interp);
      }
    }
  else
    {
      interp = -1;
      //Default to Median
    }

  /* Get or set up the name of the input background aperture file */
  if ((opt = get_online_option("in_AF", argc, argv)))
    {
      strcpy(aper_file, opt);
      strcpy(aper_file_path, opt);
    }
  else
    {
      if (makemask)
        {
          replace_file_extension(grism_image, aper_file, ".fits", ".OAF",
                                 conf->science_numext);
          build_path(AXE_OUTPUT_PATH, aper_file, aper_file_path);
        }
      else
        {
          replace_file_extension(grism_image, aper_file, ".fits", ".BAF",
                                 conf->science_numext);
          build_path(AXE_OUTPUT_PATH, aper_file, aper_file_path);
        }
    }

  if ((opt = get_online_option("niter_med", argc, argv)))
    {
      niter_med = atoi(opt);
    }
  else
    {
      niter_med = 0;
    }

  if ((opt = get_online_option("niter_fit", argc, argv)))
    {
      niter_fit = atoi(opt);
    }
  else
    {
      niter_fit=0;
    }

  // check for the parameter "kappa"
  if ((opt = get_online_option("kappa", argc, argv)))
    kappa = atof(opt);
  else
    kappa=0.0;

  // check for the parameter "smooth_length"
  if ((opt = get_online_option("smooth_length", argc, argv)))
    sm_length = atoi(opt);
  else
    sm_length=0;

  // check for the parameter "fwhm"
  if ((opt = get_online_option("fwhm", argc, argv)))
    fwhm = atof(opt);
  else
    fwhm=0.0;


  fprintf(stdout,
          "aXe_BE: Main configuration file name:         %s\n",
          conf_file_path);
  fprintf(stdout,
          "aXe_BE: Input data file name:                 %s\n",
          grism_image_path);
  fprintf(stdout,
          "aXe_BE: SCI extension number:                 %d\n",
          conf->science_numext);
  fprintf(stdout,
          "aXe_BE: ERR extension number:                 %d\n",
          conf->errors_numext);
  fprintf(stdout,
          "aXe_BE: DQ extension number:                  %d\n",
          conf->dq_numext);
  fprintf(stdout,
          "aXe_BE: Input aperture file name:             %s\n",
          aper_file_path);
  fprintf(stdout,
          "aXe_BE: Output background file name:          %s\n",
          back_file_path);

  fprintf(stdout,
          "aXe_BE: Number of points for computation :    %d\n", np);
  fprintf(stdout,
          "aXe_BE: Interpolation order \n"
          "        (-1=median, 0= cst, 1= linear etc..): %d\n", interp);
  if (niter_med)
    fprintf(stdout,
          "aXe_BE: Number of klippings on median:        %d\n", niter_med);
  if (niter_fit)
    fprintf(stdout,
          "aXe_BE: Number of klippings on fit:           %d\n", niter_fit);
  if (niter_med || niter_fit)
    fprintf(stdout,
          "aXe_BE: Kappa for kappa-sigma klippings:      %f\n", kappa);
  if (sm_length && fwhm)
    {
     fprintf(stdout,
          "aXe_BE: Final smoothing length:               %i\n", sm_length);
     fprintf(stdout,
          "aXe_BE: Gaussian FWHM in smoothing:           %f\n", fwhm);
    }
  if (makemask)
    fprintf(stdout, "aXe_BE: Producing mask file.\n");
  if (nor_flag)
    fprintf(stdout, "aXe_BE: Producing old background image.\n");

  fprintf(stdout, "\n\n");

  fprintf(stdout, "aXe_BE: ");

  //
  // try to get the descriptor 'exptime' from the 'sci'-extension
  //
  exptime = (double)get_float_from_keyword(grism_image_path, conf->science_numext, conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_image_path, 1, conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;

  obs = load_image_t(grism_image_path, conf->science_numext,
                   conf->errors_numext, conf->dq_numext,
                   conf->dqmask, exptime, conf->rdnoise);

  fprintf(stdout, "aXe_BE: Loading the object list...");
  oblist = file_to_object_list_seq(aper_file_path, obs);

  fprintf(stdout, "%d objects loaded.\n", object_list_size(oblist));

  fprintf(stdout, "aXe_BE: Computing the background image...\n");

  if (makemask)
    {
      backg = compute_backsub_mask (obs, oblist);
    }
  else
    {
      if (np>0)
        {
          //backg = compute_fullimg_background2(obs, oblist, np, interp);
          backg = compute_fullimg_background(obs, oblist, np, interp, niter_med,
                                             niter_fit, kappa, nor_flag, sm_length, fwhm);
        }
      else
        {
          backg = compute_fullimg_global_background(obs, oblist, interp, sm_length, fwhm);
        }
    }

  background_to_FITSimage(back_file_path, backg, obs);

  if (oblist != NULL)
    free_oblist(oblist);

  free_observation(obs);
  {
    FITScards      *cards;
    cards = get_FITS_cards(grism_image_path, 1);
    put_FITS_cards(back_file_path, 1, cards);
    free_FITScards(cards);
  }
  fprintf(stdout, "aXe_BE: Done...\n");
  exit(0);
}
