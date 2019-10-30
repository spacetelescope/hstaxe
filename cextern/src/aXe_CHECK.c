/*
 */
#include "aper_check.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aXe_utils.h"
#include "spce_PET.h"
#include "inout_aper.h"
#include "trace_conf.h"
#include "aper_conf.h"
#include "spc_sex.h"
#include "disp_conf.h"
#include "spc_wl_calib.h"
#include "spce_binning.h"
#include "spc_spc.h"
#include "spce_pgp.h"
#include "spce_output.h"
#include "spc_FITScards.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"




int
main(int argc, char *argv[])
{
  char           *opt;

  char            grism_image[MAXCHAR];
  char            grism_image_path[MAXCHAR];

  char            aper_file[MAXCHAR];
  char            aper_file_path[MAXCHAR];

  char            PET_file[MAXCHAR];
  char            PET_file_path[MAXCHAR];

  char            outputroot[MAXCHAR];
  char            outputroot_path[MAXCHAR];

  char            output_path[MAXCHAR];

  //object        **oblist;
  observation    *obs;
  ap_pixel       *PET;
  fitsfile       *PET_ptr;

  int             index;
  //char            label[MAXCHAR];

  int             f_status = 0;

  aXe_mask       *mask;


  int             aperID, beamID;
  //FITScards      *cards;

  if ((argc < 3) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
	    "ST-ECF European Coordinating Facility\n"
	    "Copyright (C) 2002 Nor Pirzkal\n"
	    "aXe_CHECK Version %s:\n"
	    "           Task that generates some Postscript stamp images from the\n"
	    "           content of a Pixel Extraction Table (PET)\n"
	    "           The trace is overplotted and is generated from the group\n"
	    "           of pixels in the PET which have the smallest projected distance\n"
	    "           to the trace in the PET. No analytical a-priori trace function \n"
	    "           is used.\n"
	    "\n"
	    "Usage:\n"
	    "      aXe_CHECK grism_image Pixel_extraction_table config_filename [options]\n"
	    "\n"
	    "Options:\n"
	    "      -outputroot=[string]  - overwrite the automatically generated path\n"
	    "                              and rootname of the output postscript stamp\n"
	    "                              images\n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_CHECK: Starting...\n");

  index = 0;

  strcpy(grism_image, argv[++index]);
  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);

  strcpy(PET_file, argv[++index]);
  build_path(AXE_OUTPUT_PATH, PET_file, PET_file_path);

  strcpy(aper_file, argv[++index]);
  build_path(AXE_OUTPUT_PATH, aper_file, aper_file_path);

  if ((opt = get_online_option("outputroot", argc, argv))) {
    strcpy(outputroot, opt);
  } else {
    replace_file_extension(PET_file, outputroot, ".PET.fits", "",
				       -1);
  }

  build_path(AXE_OUTPUT_PATH, outputroot, outputroot_path);
  sprintf(output_path, "%s.CHECK.fits", outputroot_path);

  fprintf(stdout, "aXe_CHECK: Input file name:            %s\n",
	  grism_image_path);
  fprintf(stdout, "aXe_CHECK: Input Aperture file name:            %s\n",
	  aper_file_path);
  fprintf(stdout, "aXe_CHECK: Input PET file name:                 %s\n",
	  PET_file_path);
  fprintf(stdout, "aXe_CHECK: Name of CHECK file :            %s\n",
	  output_path);

  obs = load_dummy_observation();
  obs = load_image_t(grism_image_path, 1, -1, -1, 0, 0, 0);

  /* Loading the object list */
  //fprintf(stdout, "aXe_CHECK: Loading object aperture list...");
  //oblist = file_to_object_list_seq(aper_file_path, obs);
  //fprintf(stdout, "Done.\n");

  //Open the OPET file for reading
  fits_open_file(&PET_ptr, PET_file_path, READONLY, &f_status);
  if (f_status) {
    ffrprt(stdout, f_status);
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
		"aXe_CHECK: Could not open file: %s",
		PET_file_path);
  }
  mask = aXe_mask_init(obs);

  while (1) {

    /* Get the PET for this object */
    PET = get_ALL_from_next_in_PET(PET_ptr, &aperID, &beamID);
    if ((aperID == -1) && (beamID == -1))
      break;
    fprintf(stdout, "aXe_CHECK: object %d%c", aperID, BEAM(beamID));
    fflush(stdout);
    //objindex = find_object_in_object_list(oblist, aperID);

    {
      //add_ap_p_to_aXe_mask(PET, mask);
      //mark_trace_in_aXe_mask(PET, mask);
    }

    if (PET != NULL)
      free(PET);

    fprintf(stdout, ".Done\n");
    fflush(stdout);
  }
  //free_oblist(oblist);

  fits_close_file(PET_ptr, &f_status);

  gsl_to_FITSimage(mask->img, "tmp.fits", 1, NULL);

  fprintf(stdout, "aXe_CHECK: Done...\n");


	/* Copy the header info from the grism image */
  {
    FITScards      *cards;
    cards = get_FITS_cards(PET_file_path, 1);
    put_FITS_cards("tmp.fits", 1, cards);
    free_FITScards(cards);
  }

  exit(0);
}
