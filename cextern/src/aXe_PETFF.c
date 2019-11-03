/**
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
#include "spce_output.h"
#include "spc_FITScards.h"
#include "spc_flatfield.h"

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

  char FF_file[MAXCHAR];
  char FF_file_path[MAXCHAR];

  char PET_file[MAXCHAR];
  char PET_file_path[MAXCHAR];

  aperture_conf *conf;

  object **oblist;
  observation *obs;

  ap_pixel *PET;

  int index, i;

  fitsfile *OPET_ptr;
  int f_status = 0;
  int bckmode = 0;
  int aperID, beamID, objindex;

  poly_cube_flatfield *FF_poly_cube;

  double exptime;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "ST-ECF European Coordinating Facility\n"
	       "aXe_PETFF Version %s:"
	       "           Task to apply a wavelength dependent flat-fielding to the COUNT and\n"
	       "           ERROR parts of a Pixel Extraction Table (PET).\n"
	       "           The wavelength of a pixel is used in conjunction with a flat-fielding\n"
	       "           data cube containing the coefficents of a polynomial which can be used \n"
	       "           to compute at each pixel (x,y): \n"
	       "                   FF(x,y,lambda) = a_0(x,y) + a_1(x,y)*lambda + .. +a_i * lambda^i\n"
	       "           The coefficients a_0(x,y) are stored in the first data extention of the\n"
	       "           flat-field cube, a_1(x,y) in the second, etc...\n"
	       "           The name of the flat-field cube is read from the aXe configuration file.\n"
	       "\n"
	       "Usage:\n"
	       "      aXe_PETFF g/prism_filename configuration_filename [option]\n"
	       "Options:\n"
	       "             -FFNAME=[string]   - overwrite the default input flat-field cube name\n"
	       "             -bck               - apply FF to the background Pixel Extraction Table (BPET)\n"
	       "\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_PETFF: Starting...\n");

  index = 0;
  strcpy (grism_file, argv[++index]);
  build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  /* Determine if we are using the special bck mode */
  /* In this mode, file names are handled diferently */
  if ((opt = get_online_option("bck", argc, argv)))
    bckmode = 1;

  conf = get_aperture_descriptor (conf_file_path);

  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);

  /* Build aperture file name */
  replace_file_extension (grism_file, aper_file, ".fits",
			  ".OAF", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);

  /* Build object PET file name */
  if (bckmode) {
    replace_file_extension (grism_file, PET_file, ".fits",
			    ".BCK.PET.fits", conf->science_numext);
    build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
  } else {
    replace_file_extension (grism_file, PET_file, ".fits",
			    ".PET.fits", conf->science_numext);
    build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
  }


  /* Build the FF file name */
  sprintf(FF_file,"%s",conf->FFname);
  build_path (AXE_CONFIG_PATH, FF_file, FF_file_path);

  if ((opt = get_online_option ("FFNAME", argc, argv)))
    {
      strcpy(FF_file,opt);
      strcpy(FF_file_path,opt);
    }

  fprintf (stdout, "aXe_PETFF: Input Aperture file name:            %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_PETFF: Input PET file name:                 %s\n",
	   PET_file_path);
  fprintf (stdout, "aXe_PETFF: Input Flat-field cube file name:     %s\n",
	   FF_file_path);

  //
  // try to get the descriptor 'exptime' from the 'sci'-extension
  //
  exptime = (double)get_float_from_keyword(grism_file_path, conf->science_numext, conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_file_path, 1, conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;
  obs = load_image_t(grism_file_path, conf->science_numext,
		   conf->errors_numext, conf->dq_numext,
		   conf->dqmask, exptime, conf->rdnoise);


  /* Loading the object list */
  fprintf (stdout, "aXe_PETFF: Loading object aperture list...");fflush(stdout);
  oblist = file_to_object_list_seq (aper_file_path, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  // Open the OPET file for reading/writing
  fits_open_file (&OPET_ptr, PET_file_path, READWRITE, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETFF: Could not open file: %s\n",
		   PET_file_path);
    }

  i = 0;
  if ((oblist!=NULL) && (strcmp(FF_file,"None")))
    {
      /* Load the FF cube */
      fprintf (stdout, "aXe_PETFF: Loading the FF cube...");fflush(stdout);
      FF_poly_cube = load_flat_poly_cube(FF_file_path);
      fprintf (stdout, "Done.\n");
      while (1)
	{
	  /* Get the PET for this object */
	  PET = get_ALL_from_next_in_PET(OPET_ptr, &aperID, &beamID);
	  if ((aperID==-1) && (beamID==-1)) break;
	  if (PET==NULL) continue; /* PET is empty, skip it */
	  /*fprintf (stdout, "aXe_PETFF: BEAM %d%c", aperID, BEAM(beamID));*/
	  objindex =  find_object_in_object_list(oblist,aperID);

	  /* Compute FF information for each pixel */
	  /* Divide count and error by it */
	  {
	    int j = 0;
	    double ff;
	    while (PET[j].p_x != -1)
	      {
		ff = poly_cube_flatfield_lambda (PET[j].lambda, PET[j].p_x, PET[j].p_y,
						 FF_poly_cube);
		if (ff!=0)
		  {
		    PET[j].count /= ff;
		    PET[j].error /= ff;
		  }
		j++;
	      }
	  }

	  /* update PET table with the new FF info */
	  {
	    char ID[60];
	    sprintf (ID, "%d%c", oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));
	    add_ALL_to_PET (PET, ID, OPET_ptr,1);
            }
	  if (PET!=NULL)
	    {
	      free(PET);
	      PET=NULL;
	    }
	  /*fprintf (stdout, ".Done\n");*/
	  i++;
	}
      /* Free the FF cube */
      free_flat_poly_cube(FF_poly_cube);
    }
  if (oblist!=NULL) free_oblist (oblist);

  fits_close_file (OPET_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETFF: " "Error closing PET: %s\n",
		   PET_file_path);
    }
  free_observation(obs);
  fprintf (stdout, "aXe_PETFF: Done...\n");
  exit (0);
}
