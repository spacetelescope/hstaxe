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
#include "scaleback_utils.h"

#define AXE_IMAGE_PATH  "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

int
main(int argc, char *argv[])
{
  char *opt;

  char  scale_image[MAXCHAR];
  char  scale_image_path[MAXCHAR];
  
  char  conf_file[MAXCHAR];
  char  conf_file_path[MAXCHAR];

  char  grism_image[MAXCHAR];
  char  grism_image_path[MAXCHAR];

  char  grism_mask[MAXCHAR];
  char  grism_mask_path[MAXCHAR];

  char  bck_image[MAXCHAR];
  char  bck_image_path[MAXCHAR];

  char  plist_name[MAXCHAR];
  char  plist_name_path[MAXCHAR];

  int scale_to_master=0;
  int make_plis=0;


  if ((argc <= 3) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
        "aXe_SCALEBCK Version %s:\n"
        "              Task for finding the scaling value between the pixels of\n"
        "              two images. pixels marked in the dq-extension or in a mask\n"
        "              image are excluded. From the extracted list of pixels, the\n"
        "              scale is determined via several kappa-sigma clipping\n"
        "              iterations around the median.\n"
        "              A scaled version of either the input image or the grism\n"
        "              image is written to the disk.\n"
        "\n"
        "Usage:\n"
        "     aXe_SCALEBCK [g/prism image] [g/prism mask] [config file] [master sky] [options]\n"
        "\n"
        "Options:\n"
        "             -toMaster       - scale the g/prism image to the master sky\n"
        "             -make_plis      - create an ASCII pixel list as output\n"
        "\n"
        "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_SCALEBCK: Starting...\n");

  // getting the grism/prism image name
  // get the name mask file for the scale image
  strcpy(grism_image, argv[1]);
  build_path(AXE_IMAGE_PATH,  grism_image, grism_image_path);

  // get the name mask file for the scale image
  strcpy(grism_mask,  argv[2]);
  build_path(AXE_OUTPUT_PATH,  grism_mask, grism_mask_path);

  // get the configuration file name
  strcpy(conf_file,   argv[3]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);

  // get the master sky file name
  strcpy(scale_image, argv[4]);
  build_path(AXE_CONFIG_PATH,  scale_image, scale_image_path);

  // fix the scaled image path
  replace_file_extension(grism_image, bck_image, ".fits",
			 ".SGRISM.fits", 2);
  build_path(AXE_OUTPUT_PATH,  bck_image, bck_image_path);

  // fix the name of the pixel list
  replace_file_extension(grism_image, plist_name, ".fits",
			 ".plis", 2);
  build_path(AXE_OUTPUT_PATH,  plist_name, plist_name_path);

  // check option for scaling target
  if ((opt = get_online_option ("toMaster", argc, argv)))
	  scale_to_master = 1;
  else
      scale_to_master = 0;

  // check option to make pixel list
  if ((opt = get_online_option ("make_plis", argc, argv)))
	  make_plis = 1;
  else
	  make_plis = 0;

  fprintf(stdout,
      "aXe_SCALEBCK: Input grism image:       %s\n",
      grism_image_path);
  fprintf(stdout,
      "aXe_SCALEBCK: Input grism mask:       %s\n",
      grism_mask_path);
  fprintf(stdout,
      "aXe_SCALEBCK: aXe configuration file: %s\n",
      conf_file_path);
  fprintf(stdout,
      "aXe_SCALEBCK: Master sky image:        %s\n",
      scale_image_path);
  fprintf(stdout,
      "aXe_SCALEBCK: Scaled grism image:      %s\n",
      bck_image_path);
  if (make_plis)
    fprintf(stdout,
        "aXe_SCALEBCK: Pixel list name:         %s\n",
        plist_name_path);
  if (scale_to_master)
    fprintf(stdout, "aXe_SCALEBCK: Grism is scaled to master sky!\n");
  else
    fprintf(stdout, "aXe_SCALEBCK: Master sky is scaled to grism!\n");

  // do some line feed
  fprintf(stdout, "\n\n");

  // read the configuration file
  // determine all extensions
  //conf = get_aperture_descriptor(conf_file_path);
  //get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);

  // run the subroutine which does everything
  make_scale_back(grism_image_path, grism_mask_path, conf_file_path,
      scale_image_path, bck_image_path, plist_name_path,
      scale_to_master, make_plis);

  // give feedback and out
  fprintf(stdout, "aXe_SCALEBCK: Done...\n");
  exit(0);
}
