/*
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aXe_grism.h"
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
#include "spc_FITScards.h"
#include "spce_output.h"

#define AXE_IMAGE_PATH   "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH  "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH  "AXE_CONFIG_PATH"
#define AXE_DRIZZLE_PATH "AXE_DRIZZLE_PATH"


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

  char PET_file[MAXCHAR];
  char PET_file_path[MAXCHAR];

  //  char outputroot[MAXCHAR];
  //  char outputroot_path[MAXCHAR];
  char STP_file[MAXCHAR];
  char STP_file_path[MAXCHAR];

  //  char output_path[MAXCHAR];

  aperture_conf *conf;

  object **oblist;
  observation *obs;

  ap_pixel *PET;

  int index, i;
  char label[MAXCHAR];

  fitsfile *PET_ptr, *STP_ptr;
  int f_status = 0;

  int aperID, beamID, objindex;
  FITScards *cards;
  FITScards *xymin_cards;

  gsl_matrix *rstamp;
  drzstamp *drzstmp;
  drzstamp_dim dim;

  int rectified = 0;
  //int drizzled = 0;
  int drizzle=0;
  int drzstamp=0;
  int for_grism=0;

  d_point stp_min;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "ST-ECF European Coordinating Facility\n"
	       "aXe_STAMPS Version %s:\n"
	       "           Task that generates some Postscript stamp images from the\n"
	       "           content of a Pixel Extraction Table (PET)\n"
	       "           The trace is overplotted and is generated from the group\n"
	       "           of pixels in the PET which have the smallest projected distance\n"
	       "           to the trace in the PET. No analytical a-priori trace function \n"
	       "           is used.\n"
	       "\n"
	       "Usage:\n"
	       "      aXe_STAMPS g/prism_filename configuration_filename [options]\n"
	       "\n"
	       "Options:\n"
	       "      -rectified            - produces rectified stamp image following the\n"
	       "                              direction of the extraction process\n"
	       "      -drz                  - use $AXE_DRIZZLE_PATH to locate the grism-, OAF-\n"
	       "                              and PET-files instead of $AXE_IMAGE/OUTPUT_PATH\n"
	       "      -in_AF=[string]       - overwrite the automatically generated name\n"
	       "                              of the input aperture file\n"
	       "      -in_PET=[string]      - overwrite the automatically generated name\n"
	       "                              of the input PET file\n"
	       "      -out_STP=[string]     - overwrite the automatically generated name\n"
	       "                              of the output stamp FITS image\n"
	       "\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_STAMPS: Starting...\n");

  index = 0;
  strcpy (grism_file, argv[++index]);
  if ((opt = get_online_option ("drz", argc, argv))){
    build_path (AXE_DRIZZLE_PATH, grism_file, grism_file_path);
    drizzle=1;
  }
  else{
    build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);
    drizzle=0;
  }

  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);


  conf = get_aperture_descriptor (conf_file_path);

  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);

  /* Build aperture file name */
  if ((opt = get_online_option("in_AF", argc, argv))){
    /* get it */
    strcpy (aper_file, opt);
    strcpy (aper_file_path, opt);
  }
  else {
    replace_file_extension (grism_file, aper_file, ".fits",
			    ".OAF", conf->science_numext);
    if (drizzle)
      build_path (AXE_DRIZZLE_PATH, aper_file, aper_file_path);
    else
      build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
  }

  if ((opt = get_online_option("in_PET", argc, argv))){
    /* get it */
    strcpy (PET_file, opt);
    strcpy (PET_file_path, opt);
  }
  else{
    /* Build object PET file name */
    replace_file_extension (grism_file, PET_file, ".fits",
			    ".PET.fits", conf->science_numext);
    if (drizzle)
      build_path (AXE_DRIZZLE_PATH, PET_file, PET_file_path);
    else
      build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
  }

  if ( (opt = get_online_option ("out_STP", argc, argv)) )
    {
      strcpy (STP_file, opt);
      //      strcpy (STP_file_PATH, opt);
    }
  else
    {
      //      replace_file_extension (PET_file, STP_file, ".PET.fits", ".STP.fits",-1);
      replace_file_extension (grism_file, STP_file, ".fits", ".STP.fits",conf->science_numext);
    }


  if (drizzle)
    build_path (AXE_DRIZZLE_PATH, STP_file, STP_file_path);
  else
    build_path (AXE_OUTPUT_PATH, STP_file, STP_file_path);

  if ((opt=get_online_option("rectified",argc,argv)))
    rectified = 1;
  else if ((opt = get_online_option ("drzstamp", argc, argv)))
    drzstamp=1;


  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameters
  fprintf (stdout, "aXe_STAMPS: Input Image file name:               %s\n",
	   grism_file_path);
  fprintf (stdout, "aXe_STAMPS: Input Aperture file name:            %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_STAMPS: Input PET file name:                 %s\n",
	   PET_file_path);
  fprintf (stdout, "aXe_STAMPS: Name of STP file :                   %s\n",
	   STP_file_path);
  if (rectified)
    fprintf (stdout, "aXe_STAMPS: Producing rectified stamp images.\n");
  else if (drzstamp)
    fprintf (stdout, "aXe_STAMPS: Producing drizzled stamp images.\n");
  else
    fprintf (stdout, "aXe_STAMPS: Producing trace stamp images.\n");


  // Loading the object list
  obs = load_dummy_observation ();
  fprintf (stdout, "\naXe_STAMPS: Loading object aperture list...");
  oblist = file_to_object_list_seq (aper_file_path, obs);
  fprintf (stdout,"%d objects loaded.\n\n",object_list_size(oblist));

  //  Open the OPET file for reading
  fits_open_file (&PET_ptr, PET_file_path, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_STAMPS: Could not open file: %s\n",
		   PET_file_path);
    }

  STP_ptr = create_FITSimage_opened (STP_file_path, 1);

  i = 0;
  if (oblist!=NULL)
    {
      while (1)
	{

	  /* Get the PET for this object */
	  PET = get_ALL_from_next_in_PET(PET_ptr, &aperID, &beamID);
	  if ((aperID==-1) && (beamID==-1)) break;
	  /*fprintf (stdout, "aXe_STAMPS: BEAM %d%c.", aperID, BEAM(beamID));*/
	  objindex =  find_object_in_object_list(oblist,aperID);

	  sprintf (label, "%s.%d%c.ps/CPS", STP_file_path,
		   oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));

	  // determine the minimum in x and y
	  stp_min.x = -1.0;
	  stp_min.y = -1.0;

	  // special treatment for FORS2: dirzzled stamp images
	  if (drzstamp)
	    {

	      for_grism = check_for_grism (conf_file_path, beamID);
	      if (!for_grism)
		aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			     "aXe_STAMPS: Drizzled stamp images are not\n"
			     "supported for prism images. Please choose a different option!\n");

	      // get the dimensions of the drizzled stamp
	      dim = get_stamp_dim(PET, oblist[objindex]->beams[beamID].width, conf, beamID, &stp_min);

	      // fill the drzstamp structure (counts plus weight matrix)
	      drzstmp = drizzled_stamp_img (PET,oblist[objindex]->beams[beamID].width,
					    oblist[objindex]->beams[beamID].orient,  dim);

	      // does this make sense -> no, STP files are only for validation
	      // interpolate_over_NaN (drzstmp->counts);

	      // give the extension name and store the counts
	      sprintf (label, "BEAM_%d%c",oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));
	      gsl_to_FITSimage_opened (drzstmp->counts, STP_ptr ,0,label);

	      // in case that the stamp is no dummy, make and store the WCS header
	      if (dim.resolution)
		{
		  cards = get_WCS_FITScards(dim.xstart*dim.resolution, dim.resolution, dim.ystart);
                  put_FITS_cards_opened(STP_ptr,cards);
                  free_FITScards(cards);
		}

	      // make and store the default header
	      cards = beam_to_FITScards(oblist[objindex],beamID);
              xymin_cards = stpmin_to_FITScards(stp_min);
	      put_FITS_cards_opened(STP_ptr,cards);
              put_FITS_cards_opened(STP_ptr,xymin_cards);
	      free_FITScards(cards);
              free_FITScards(xymin_cards);

	      // give the extension name and store the weight image
	      sprintf (label, "WEIG_%d%c",oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));
	      gsl_to_FITSimage_opened (drzstmp->weight, STP_ptr ,0,label);

	      // in case that the stamp is no dummy, make and store the WCS header
	      if (dim.resolution)
		{
		  cards = get_WCS_FITScards(dim.xstart*dim.resolution, dim.resolution, dim.ystart);
		  put_FITS_cards_opened(STP_ptr,cards);
		  free_FITScards(cards);
		}

	      // make and store the default header
	      cards = beam_to_FITScards(oblist[objindex],beamID);
              xymin_cards = stpmin_to_FITScards(stp_min);
	      put_FITS_cards_opened(STP_ptr,cards);
              put_FITS_cards_opened(STP_ptr,xymin_cards);
	      free_FITScards(cards);
              free_FITScards(xymin_cards);

	      // free the structure with the count and weight matrix
	      free_drzstamp(drzstmp);
	    }
	  else
	    {
	      if (rectified)
		{
		  rstamp = rectified_stamp_img (PET,oblist[objindex]->beams[beamID].width, &stp_min);
		}
	      else
		{
		  rstamp = stamp_img (PET,oblist[objindex]->beams[beamID].width, &stp_min);
		}
	      sprintf (label, "BEAM_%d%c",oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));
	      //interpolate_over_NaN (rstamp);
	      gsl_to_FITSimage_opened (rstamp, STP_ptr ,0,label);
	      cards = beam_to_FITScards(oblist[objindex],beamID);
              xymin_cards = stpmin_to_FITScards(stp_min);
	      put_FITS_cards_opened(STP_ptr,cards);
              put_FITS_cards_opened(STP_ptr,xymin_cards);
	      free_FITScards(cards);
              free_FITScards(xymin_cards);

	      free_stamp_img(rstamp);
	    }

	  if (PET!=NULL)
	    free(PET);

	  /*fprintf (stdout, " Done.\n");*/
	  i++;
	}

    }
  if (oblist!=NULL) free_oblist (oblist);

  fits_close_file (STP_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_STAMPS: Could not" " close STP file: %s\n", STP_file_path);
    }

  fprintf (stdout, "aXe_STAMPS: Done...\n");

  /* Copy the header info from the grism image */
  {
    FITScards *cards;
    cards = get_FITS_cards (PET_file_path, 1);
    put_FITS_cards(STP_file_path,1,cards);
    free_FITScards(cards);
  }
  exit (0);
}
