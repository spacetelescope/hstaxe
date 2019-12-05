/**
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_multifit_nlin.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
//#include "spce_PET.h"
#include "inout_aper.h"
//#include "trace_conf.h"
#include "aper_conf.h"
//#include "spc_sex.h"
//#include "disp_conf.h"
//#include "spc_wl_calib.h"
//#include "spce_binning.h"
#include "spc_spc.h"
//#include "spce_output.h"
//#include "spc_FITScards.h"
#include "fringe_conf.h"
#include "trfit_utils.h"
#include "ipixcorr_utils.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"


int
main (int argc, char *argv[])
{
  char *opt;

  char ipc_aper_file[MAXCHAR];
  char ipc_aper_file_path[MAXCHAR];

  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char grism_file[MAXCHAR];
  char grism_file_path[MAXCHAR];

  char bck_image[MAXCHAR];
  char bck_image_path[MAXCHAR];

  char IPC_file[MAXCHAR];
  char IPC_file_path[MAXCHAR];

  char ipc_file[MAXCHAR];
  char ipc_file_path[MAXCHAR];

  char SPC_file[MAXCHAR];
  char SPC_file_path[MAXCHAR];

  char IPC_func[MAXCHAR];
  char IPC_func_path[MAXCHAR];

  aperture_conf *conf;

  object **oblist;
  observation *obs=NULL, *bck;

  //Initialize obs struct memebers
  obs->grism = NULL;
  obs->pixerrs = NULL;
  obs->dq = NULL;

  interpolator *ipcorr;

  interpolator *nlincorr;

  int index;
  //int i;

  full_spectr *SPC;

  fitsfile  *SPC_ptr;
  fitsfile  *IPC_ptr;
  FITScards *cards;

  int f_status=0;

  int aperID=0, beamID=0;
  //objindex;

  int spec_OAF=0;
  int point_like=0;

  beam act_beam;
  beam *beam_ptr;

  double max_ext=0.0;
  //double adcgain;
  double exptime=NAN;

  gsl_matrix *data_matrix;

  int ipixc=0;
  int nlinc=0;
  int icorr=0;
  int num=0;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
         "aXe_FRINGECORR Version %s:\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_INTPIXCORR: Starting...\n");

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

  // build the name of the modified
  // OAF file
  replace_file_extension (grism_file, ipc_aper_file, ".fits",
			      "_ipc.OAF", conf->science_numext);
  build_path(AXE_OUTPUT_PATH, ipc_aper_file, ipc_aper_file_path);

  // get or set up the output SPC name
  if ((opt = get_online_option("out_SPC", argc, argv)))
    {
      // copy the parameter value if given
      strcpy(IPC_file, opt);
    }
  else
    {
      // copy the grism name and construct the name
      // if not given as parameter
      strcpy(IPC_file, grism_file);
      replace_file_extension(grism_file, IPC_file, ".fits",
			     "_ipc.SPC.fits", conf->science_numext);
    }
  build_path(AXE_OUTPUT_PATH, IPC_file, IPC_file_path);


  // get or set up the output SPC name
  if ((opt = get_online_option("in_SPC", argc, argv)))
    {
      // copy the parameter value if given
      strcpy(SPC_file, opt);
    }
  else
    {
      // copy the grism name and construct the name
      // if not given as parameter
      strcpy(SPC_file, grism_file);
      replace_file_extension(grism_file, SPC_file, ".fits",
			     ".SPC.fits", conf->science_numext);
    }
  build_path(AXE_OUTPUT_PATH, SPC_file, SPC_file_path);


  // form the background image name
  replace_file_extension(grism_file, bck_image, ".fits",
			 ".BCK.fits", conf->science_numext);
  build_path (AXE_OUTPUT_PATH, bck_image, bck_image_path);

  // check whether an intra pixel sensitivity correction
  // shall be applied
  if ((opt = get_online_option("ipixcorr", argc, argv)))
    ipixc = 1;
  // check whether an intra pixel sensitivity correction
  // shall be applied
  if ((opt = get_online_option("nlincorr", argc, argv)))
    nlinc = 1;

  // check whether something is
  // done at all
  if (ipixc || nlinc)
    {
      // Build the intrapixel correction file name
      sprintf(IPC_func,"%s",conf->IPIXfunc);

      // overwrite it with the parameter if given
      if ((opt = get_online_option ("IPIXFUNCTION", argc, argv)))
	strcpy(IPC_func,opt);

      if (!strcmp(IPC_file,"None"))
	// check whether the correction function is given;
	// complain if not given!
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_INTPIXCORR: Correction function must be given!!\n");
      else
	// compose the full pathname
	build_path (AXE_CONFIG_PATH, IPC_func, IPC_func_path);


      // get the maximum extension number
      if ((opt = get_online_option ("max_ext", argc, argv)))
	max_ext = atof(opt);

      // check whether either max_ext is defined
      // OR a special OAF file is given
      if (!max_ext && !spec_OAF)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_INTPIXCORR: Either maximum extension\n"
		     "OR a special OAF file must be given!\n");
    }

  // print on screen what's done
  fprintf (stdout, "aXe_INTPIXCORR: Grism file name:                     %s\n",
	   grism_file_path);
  fprintf (stdout, "aXe_INTPIXCORR: Configuration file name:             %s\n",
	   conf_file_path);
  fprintf (stdout, "aXe_INTPIXCORR: Correction function:                 %s\n",
	   IPC_func_path);
  fprintf (stdout, "aXe_INTPIXCORR: Input OAF file name:                 %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_INTPIXCORR: Input SPC file name:                 %s\n",
	   SPC_file_path);
  fprintf (stdout, "aXe_INTPIXCORR: Output SPC file name:                %s\n",
	   IPC_file_path);
  fprintf (stdout, "aXe_INTPIXCORR: Maximal extension for point objects: %f\n",
	   max_ext);
  if (ipixc)
    fprintf (stdout, "aXe_INTPIXCORR: Intra-pixel sensitivity correction is aplied!\n");
  if (nlinc)
    fprintf (stdout, "aXe_INTPIXCORR: Non-linearity correction is aplied!\n");


  // check whether something is
  // done at all
  if (ipixc || nlinc)
    {
      // check whether the intra pixel
      // sensitivity is applied
      if (ipixc)
	{
	  // determine the exposure time
	  if (isnan(exptime))
	    exptime = (double)get_float_from_keyword(grism_file_path, 1,
						     conf->exptimekey);
	  if (isnan(exptime))
	    exptime = 1.0;

	  // load the image data
	  obs = load_image_t(grism_file_path, conf->science_numext,
			     conf->errors_numext, conf->dq_numext,
			     conf->dqmask, exptime, conf->rdnoise);

	  // load the background data
	  bck = load_image_t(bck_image_path, conf->science_numext,
			     conf->errors_numext, conf->dq_numext,
			     conf->dqmask, exptime, conf->rdnoise);

	  // subtract the background from
	  // the image data
	  data_matrix = bcksub_observation(obs, bck);
	}

      // load the aperture list
      fprintf (stdout, "aXe_INTPIXCORR: Loading object aperture list...");fflush(stdout);
      oblist = file_to_object_list_seq (aper_file_path, obs);
      fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));fflush(stdout);

      // open the input SPC file
      SPC_ptr = get_SPC_opened(SPC_file_path, 0);

      // open the output SPC file
      IPC_ptr = create_SPC_opened (IPC_file_path,1);
      cards = get_FITS_cards_opened(SPC_ptr);
      put_FITS_cards_opened(IPC_ptr, cards);
      free_FITScards(cards);

      if (ipixc)
	// create an interpolator reading in the fits file
	ipcorr = create_interp_ftable(IPC_func_path,2,"X","FX",IPCORR_INTERP_TYPE);
      else
	ipcorr = NULL;

      if (nlinc)
	// for the NICMOS nonlinearity correction
	// create an interpolator for the non-linear correction
	//nlincorr = create_nlincor();
	nlincorr = create_interp_ftable(IPC_func_path,2,"LAMBDA","BPAR",NLINCORR_INTERP_TYPE);
      else
	nlincorr = NULL;

      while ((aperID!=-1) || (beamID!=-1))
	{
	  // load the next extension from the SPC table
	  SPC = get_ALL_from_next_in_SPC(SPC_ptr, &aperID, &beamID);

	  // break the loop if the end was reached
	  if ((aperID==-1) && (beamID==-1))
	    break;

	  // if the full spectrum is NULL,
	  // go to the next beam
	  if (!SPC)
	    continue;

	  // select the actual beam from the beam list
	  act_beam = find_beam_in_object_list(oblist, aperID, beamID);
	  beam_ptr = find_beamptr_in_object_list(oblist, aperID, beamID);

	  // determine whether the object is pointlike
	  point_like = is_pointlike(act_beam, spec_OAF, max_ext);

	  // if the obejct is indeed
	  // a pointlike one
	  if (point_like)
	    {
	      // check whether the intra-pixel sensititvity correction
	      // shall be applied
	      if (ipixc)
		{
		  // form the name of the file which contain
		  // the data
		  replace_file_extension(grism_file, ipc_file, ".fits",
					 ".IPC.dat", aperID);
		  build_path(AXE_OUTPUT_PATH, ipc_file, ipc_file_path);

		  icorr = fitting_ipc_corr(act_beam, conf_file_path, ipcorr,
					   SPC, obs, data_matrix, ipc_file_path, beam_ptr);

		  if (icorr)
		    fprintf(stdout, "aXe_INTPIXCORR: Corrected BEAM %d%c\n", SPC->aperID, BEAM(SPC->beamID));
		  else
		    fprintf(stdout, "aXe_INTPIXCORR: phase accuracy in BEAM %d%c too low\n", SPC->aperID, BEAM(SPC->beamID));

		}

	      // check whether the non-linearity correction
	      // shall be applied
	      if (nlinc)
		{
		  /*
		  if (!strcmp(conf->gainkey,"None"))
		    adcgain=1.0;
		  else
		    adcgain = (double)get_float_from_keyword(grism_file_path, 1, conf->gainkey);
		  */
		  // apply the NICMOS non-linearity correction
		  fprintf(stdout, "aXe_INTPIXCORR: Non-linearity correting BEAM %d%c...\n", SPC->aperID, BEAM(SPC->beamID));
		  // the gain mostly used in NICMOS
		  // is hardcoded here!!!!!
		  nlin_corr_beam(nlincorr, 6.5, SPC);
		}
	    }
	  else
	    {
	      fprintf(stdout, "aXe_INTPIXCORR: Skipping BEAM %d%c...\n", SPC->aperID, BEAM(SPC->beamID));
	    }

	  add_spectra_to_SPC_opened (IPC_ptr, SPC->fgr_spec, SPC->bck_spec,
				     SPC->obj_spec, SPC->aperID, SPC->beamID);

	  /* Copy header from OPET extension into this SPC extension */
	  cards = get_FITS_cards_opened(SPC_ptr);
	  put_FITS_cards_opened(IPC_ptr,cards);
	  free_FITScards(cards);

	  free_full_spectr(SPC);
      	}

      refurbish_object_list(oblist, 1, -100, 0);

      fprintf (stdout, "aXe_INTPIXCORR: Writing aperture file %s ... ", ipc_aper_file);fflush(stdout);
      num = object_list_to_file (oblist, ipc_aper_file_path, 0);
      fprintf (stdout, "%d beams written.\n", num);

      // free the object list
      if (oblist!=NULL)
	free_oblist (oblist);

      if (ipcorr)
	// free the interpolator
	free_interp(ipcorr);

      if (nlincorr)
	// free the interpolator
	free_interp(nlincorr);

      fits_close_file (SPC_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_INTPIXCORR: " "Error closing SPC: %s\n",
		       SPC_file_path);
	}

      fits_close_file (IPC_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_INTPIXCORR: " "Error closing SPC: %s\n",
		       IPC_file_path);
	}
    }
  // say goodbye
  fprintf (stdout, "aXe_INTPIXCORR: Done...\n");
  exit (0);
}
