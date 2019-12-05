/**
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "inout_aper.h"
#include "aper_conf.h"
#include "spc_spc.h"
#include "spce_PET.h"
#include "fringe_conf.h"
#include "ipixcorr_utils.h"

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

  char IPC_file[MAXCHAR];
  char IPC_file_path[MAXCHAR];

  char PET_file[MAXCHAR];
  char PET_file_path[MAXCHAR];

  char IPC_func[MAXCHAR];
  char IPC_func_path[MAXCHAR];
  char ID[60];

  aperture_conf *conf;

  object **oblist;
  observation *obs = NULL;

  interpolator *ipcorr;

  int index;
  //int i;

  ap_pixel *PET;

  fitsfile  *PET_ptr;
  fitsfile  *IPC_ptr;
  FITScards *cards;

  int f_status=0;
  int bckmode = 0;

  int aperID=0, beamID=0;
  //int objindex;

  int spec_OAF=0;
  int point_like=0;
  int origname=0;

  beam act_beam;

  double max_ext;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
         "aXe_FRINGECORR Version %s:\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_PETIPC: Starting...\n");

  // reset the parameter index
  index = 0;

  // read in the grism image name
  strcpy(grism_file, argv[++index]);
  build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

  // read in the configuration file name
  strcpy(conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  // load the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // determine the science extension number
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);


  // Determine if we are using the special bck mode
  // In this mode, file names are handled diferently
  if ((opt = get_online_option("bck", argc, argv)))
    bckmode = 1;

  // determine whether the corrected PET should maintain the
  // filename of the original PET
  if ((opt = get_online_option("origname", argc, argv)))
    origname = 1;

  // get or set up the input OAF name
  if ((opt = get_online_option("in_OAF", argc, argv)))
    {
      // copy the parameter value if given
      strcpy(aper_file, opt);

      // mark that there is a special OAF file
      spec_OAF=1;
    }
  else
    {
      // copy the grism name and construct the name
      // if not given as parameter
      strcpy(aper_file, grism_file);
      replace_file_extension (grism_file, aper_file, ".fits",
			      ".OAF", conf->science_numext);
    }
  build_path(AXE_OUTPUT_PATH, aper_file, aper_file_path);


  // get or set up the output PET name
  if ((opt = get_online_option("out_PET", argc, argv)))
    {
      // copy the parameter value if given
      strcpy(IPC_file, opt);
    }
  else
    {
      // copy the grism name and construct the name
      // if not given as parameter
      strcpy(IPC_file, grism_file);

      if (bckmode)
	replace_file_extension(grism_file, IPC_file, ".fits",
			       "_ipc.BCK.PET.fits", conf->science_numext);
      else
	replace_file_extension(grism_file, IPC_file, ".fits",
			       "_ipc.PET.fits", conf->science_numext);
    }
  build_path(AXE_OUTPUT_PATH, IPC_file, IPC_file_path);


  // get or set up the output PET name
  if ((opt = get_online_option("in_PET", argc, argv)))
    {
      // copy the parameter value if given
      strcpy(PET_file, opt);
    }
  else
    {
      // copy the grism name and construct the name
      // if not given as parameter
      strcpy(PET_file, grism_file);
      if (bckmode)
	replace_file_extension(grism_file, PET_file, ".fits",
			       ".BCK.PET.fits", conf->science_numext);
      else
	replace_file_extension(grism_file, PET_file, ".fits",
			       ".PET.fits", conf->science_numext);

    }
  build_path(AXE_OUTPUT_PATH, PET_file, PET_file_path);


  // Build the intrapixel correction file name
  sprintf(IPC_func,"%s",conf->IPIXfunc);

  // overwrite it with the parameter if given
  if ((opt = get_online_option ("IPIXFUNCTION", argc, argv)))
    strcpy(IPC_func,opt);

  if (!strcmp(IPC_file,"None"))
    // check whether the correction function is given;
    // complain if not given!
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_PETIPC: Correction function must be given!!\n");
  else
    // compose the full pathname
    build_path (AXE_CONFIG_PATH, IPC_func, IPC_func_path);


  // get the maximum extension number
  if ((opt = get_online_option ("max_ext", argc, argv)))
    max_ext = atof(opt);
  else
    max_ext = 0.0;

  // check whether either max_ext is defined
  // OR a special OAF file is given
  if (!max_ext && !spec_OAF)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_PETIPC: Either maximum extension\n"
		 "OR a special OAF file must be given!\n");


  fprintf (stdout, "aXe_PETIPC: Grism file name:                     %s\n",
	   grism_file_path);
  fprintf (stdout, "aXe_PETIPC: Configuration file name:             %s\n",
	   conf_file_path);
  fprintf (stdout, "aXe_PETIPC: Correction function:                 %s\n",
	   IPC_func_path);
  fprintf (stdout, "aXe_PETIPC: Input OAF file name:                 %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_PETIPC: Input PET file name:                 %s\n",
	   PET_file_path);
  if (origname)
    fprintf (stdout, "aXe_PETIPC: Output PET file name:                %s\n",
	     PET_file_path);
  else
    fprintf (stdout, "aXe_PETIPC: Output PET file name:                %s\n",
	     IPC_file_path);
  fprintf (stdout, "aXe_PETIPC: Maximal extension for point objects: %f\n",
	   max_ext);

  // load the aperture list
  fprintf (stdout, "aXe_PETIPC: Loading object aperture list...");
  fflush(stdout);
  oblist = file_to_object_list_seq (aper_file_path, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));
  fflush(stdout);


  // Open the PET file for reading
  fits_open_file (&PET_ptr, PET_file_path, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stdout, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETIPC: Could not open file: %s\n",
		   PET_file_path);
    }


  // open the output PET file
  IPC_ptr = create_PET_opened(IPC_file_path,1);

  // copy the keywords from the input
  // to the output PET
  cards = get_FITS_cards_opened(PET_ptr);
  put_FITS_cards_opened(IPC_ptr, cards);
  free_FITScards(cards);

  // create an interpolator reading in the fits file
  ipcorr = create_interp_ftable(IPC_func_path,2,"X","FX",IPCORR_INTERP_TYPE);

  while ((aperID!=-1) || (beamID!=-1))
    {
      // Get the PET for this object
      PET = get_ALL_from_next_in_PET(PET_ptr, &aperID, &beamID);
      if ((aperID==-1) && (beamID==-1))
	break;

      // PET is empty, skip it
      if (PET==NULL)
	{
	  fprintf(stdout, "Empty PET BEAM %d%c...\n", aperID, BEAM(beamID));
	  continue; }

      // select the actual beam from the beam list
      act_beam = find_beam_in_object_list(oblist, aperID, beamID);

      // determine whether the object is pointlike
      point_like = is_pointlike(act_beam, spec_OAF, max_ext);

      if (point_like)
	{
	  fprintf(stdout, "aXe_PETIPC: Correcting BEAM %d%c...\n", aperID, BEAM(beamID));
	  intpix_corr_pet(act_beam, conf_file_path, ipcorr, PET);
	}
      else
	{
	  fprintf(stdout, "aXe_PETIPC: Skipping BEAM %d%c...\n", aperID, BEAM(beamID));
	}

      sprintf(ID, "%d%c", aperID, BEAM (beamID));
      add_ALL_to_PET(PET, ID, IPC_ptr, 0);

      /* Copy header from OPET extension into this SPC extension */
      cards = get_FITS_cards_opened(PET_ptr);
      put_FITS_cards_opened(IPC_ptr,cards);
      free_FITScards(cards);

      if (PET!=NULL)
	free(PET);
    }

  // free the object list
  if (oblist!=NULL)
    free_oblist (oblist);

  // free the interpolator
  free_interp(ipcorr);

  // free the configuration structure
  free_aperture_conf(conf);

  fits_close_file (PET_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETIPC: " "Error closing PET: %s\n",
		   PET_file_path);
    }

  fits_close_file (IPC_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETIPC: " "Error closing PET: %s\n",
		   IPC_file_path);
    }

  // check whether the old
  // PET name should be retained
  if (origname)
    {
      // delete the original
      // PET file
      if (remove(PET_file_path) != 0)
	fprintf (stdout, "aXe_PETIPC: Problems deleting %s!\n", PET_file_path);

      // rename the corrcted PET file
      if (rename(IPC_file_path, PET_file_path) !=0)
	fprintf (stdout, "aXe_PETIPC: Problems renaming %s to %s!\n", IPC_file_path, PET_file_path);


    }

  // say goodbye
  fprintf (stdout, "aXe_PETIPC: Done...\n");
  exit (0);
}
