/**
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_PET.h"
#include "inout_aper.h"
//#include "trace_conf.h"
//#include "aper_conf.h"
//#include "spc_sex.h"
//#include "disp_conf.h"
//#include "spc_wl_calib.h"
//#include "spce_binning.h"
//#include "spc_spc.h"
//#include "spce_pgp.h"
//#include "spce_output.h"
//#include "spc_FITScards.h"
//#include "spc_flatfield.h"
#include "fringe_conf.h"
#include "fringe_utils.h"

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

  char fconf_file[MAXCHAR];
  char fconf_file_path[MAXCHAR];

  char grism_file[MAXCHAR];
  char grism_file_path[MAXCHAR];

  char FF_file[MAXCHAR];
  //char FF_file_path[MAXCHAR];

  char OPET_file[MAXCHAR];
  char OPET_file_path[MAXCHAR];

  char BPET_file[MAXCHAR];
  char BPET_file_path[MAXCHAR];

  aperture_conf *conf=NULL;
  fringe_conf   *fconf=NULL;

  object **oblist;
  observation *obs=NULL;

  ap_pixel *OPET=NULL;
  ap_pixel *BPET=NULL;

  fitsfile *OPET_ptr=NULL;
  fitsfile *BPET_ptr=NULL;

  beam act_beam;

  int index;
  int i;
  int f_status = 0;
  int bckmode = 0;
  int oaperID, obeamID;
  int baperID, bbeamID;

  if ((argc < 4) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "aXe_FRINGECORR Version %s:\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_FRINGECORR: Starting...\n\n");

  // copy the first parameter to the grism filename
  index = 0;
  strcpy (grism_file, argv[++index]);
  build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

  // copy the next parameter to the aXe configuration file
  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  // copy the next parameter to the fringe configuration file
  strcpy (fconf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, fconf_file, fconf_file_path);

  // Determine if we are using the special bck mode
  // In this mode, file names are handled diferently
  if ((opt = get_online_option("bck", argc, argv)))
    bckmode = 1;

  // load the aXe configuration file
  conf = get_aperture_descriptor (conf_file_path);

  // load the fringe configuration file
  fconf = load_fringe_conf(fconf_file_path);

  // Determine where the various extensions are in the FITS file
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);

  // Build aperture file name
  replace_file_extension (grism_file, aper_file, ".fits",
			  ".OAF", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);

  // Build object PET file name
  replace_file_extension (grism_file, OPET_file, ".fits",
			  ".PET.fits", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, OPET_file, OPET_file_path);

  // build the background PET file name
  // if necessary
  if (bckmode)
    {
      replace_file_extension (grism_file, BPET_file, ".fits",
			      ".BCK.PET.fits", conf->science_numext);
      build_path (AXE_OUTPUT_PATH, BPET_file, BPET_file_path);
  }


  // Give a short feedback on the filenames
  // of the various input
  fprintf (stdout, "aXe_FRINGECORR: Grism image file name:          %s\n",
	   grism_file_path);
  fprintf (stdout, "aXe_FRINGECORR: Aperture file name:             %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_FRINGECORR: aXe configuration file name:    %s\n",
	   conf_file_path);
  fprintf (stdout, "aXe_FRINGECORR: Fringe configuration file name: %s\n",
	   fconf_file_path);
  fprintf (stdout, "aXe_FRINGECORR: Object PET file name:           %s\n",
	   OPET_file_path);
  if (bckmode)
    fprintf (stdout, "aXe_FRINGECORR: Background PET file name:       %s\n",
	     BPET_file_path);

  // check whether all necessary information
  // is in place, provide defaults
  check_fringe_conf(fconf);


  fprintf (stdout, "aXe_FRINGECORR: Loading object aperture list...");
  fflush(stdout);
  oblist = file_to_object_list_seq (aper_file_path, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // Open the OPET file for reading/writing
  fits_open_file (&OPET_ptr, OPET_file_path, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_FRINGECORR: Could not open file: %s\n",
		   OPET_file_path);
    }

  // if necessary, open the background PET
  if (bckmode)
    {
      fits_open_file (&BPET_ptr, BPET_file_path, READWRITE, &f_status);
      if (f_status)
	{
	  ffrprt (stdout, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_FRINGECORR: Could not open file: %s\n",
		       BPET_file_path);
	}
    }

  i = 0;
  if ((oblist!=NULL) && (strcmp(FF_file,"None")))
    {
      while (1)
	{
	  // Get the PET for this object
	  OPET = get_ALL_from_next_in_PET(OPET_ptr, &oaperID, &obeamID);
	  if ((oaperID==-1) && (obeamID==-1))
	    break;

	  // Get the PET for this object
	  if (bckmode)
	    {
	      BPET = get_ALL_from_next_in_PET(BPET_ptr, &baperID, &bbeamID);
	      if ((baperID != oaperID) ||  (bbeamID != obeamID))
		{
		  fprintf(stderr, "%i, %i, %i, %i\n", oaperID, baperID, obeamID, bbeamID);
		  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			       "aXe_FRINGECORR: The beam order in the\n"
			       "object PET %s and the background PET\n"
			       "%s does not coincide!\n",
			       OPET_file_path, BPET_file_path);
		}
	    }

	  // check whether there is content in the PET;
	  // skip it if empty
	  if (OPET==NULL)
	    continue;

	  // refer about the actual beam in process
	  fprintf (stdout, "aXe_FRINGECORR: BEAM %d%c",oaperID, BEAM(obeamID));

	  // get the beam information from the OAF
	  act_beam = find_beam_in_object_list(oblist, oaperID, obeamID);
	  fprintf (stdout, " Done. \n");

	  fringe_correct_PET(fconf, act_beam, OPET, BPET);

	  // update PET table with the new FF info
	  {
	    char ID[60];
	    sprintf (ID, "%d%c", oaperID, BEAM (act_beam.ID));
	    add_ALL_to_PET (OPET, ID, OPET_ptr,1);

	    if (bckmode)
	      add_ALL_to_PET (BPET, ID, BPET_ptr,1);
            }

	  // free the memory of the object pixels
	  if (OPET!=NULL)
	    {
	      free(OPET);
	      OPET=NULL;
	    }

	  // free the memory of the background pixels
	  if (BPET!=NULL)
	    {
	      free(BPET);
	      BPET=NULL;
	    }
	}
    }

  // free the object list
  if (oblist!=NULL)
    free_oblist (oblist);

  // close the object PET fits-file
  fits_close_file (OPET_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETFF: " "Error closing PET: %s\n",
		   OPET_file_path);
    }

  // if necessary, close the object PET fits-file
  if (bckmode)
    {
      fits_close_file (BPET_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_PETFF: " "Error closing PET: %s\n",
		       BPET_file_path);
	}
    }

  fprintf (stdout, "aXe_FRINGECORR: Done...\n\n");
  exit (0);
}
