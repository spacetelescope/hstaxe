#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>

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
#include "spc_optimum.h"
#include "spc_spc.h"
#include "fringe_conf.h"
#include "spc_resp.h"
#include "spc_FITScards.h"

#define AXE_IMAGE_PATH   "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH  "AXE_OUTPUT_PATH"
#define AXE_DRIZZLE_PATH "AXE_DRIZZLE_PATH"
#define AXE_CONFIG_PATH  "AXE_CONFIG_PATH"

int
main (int argc, char *argv[])
{
  char *opt;

  char grism_file[MAXCHAR];
  char grism_file_path[MAXCHAR];

  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char obj_PET_file[MAXCHAR];
  char obj_PET_file_path[MAXCHAR];

  char bck_PET_file[MAXCHAR];
  char bck_PET_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char SPC_file[MAXCHAR];
  char SPC_file_path[MAXCHAR];

  char SPC_opt_file[MAXCHAR];
  char SPC_opt_file_path[MAXCHAR];

  char WHT_file[MAXCHAR];
  char WHT_file_path[MAXCHAR];

  char label[MAXCHAR];

  int i, index, dobck = 0, noflux = 1;

  object **oblist;

  FITScards      *cards;

  ap_pixel *obj_PET = NULL, *bck_PET = NULL;

  observation *obs;

  //tracestruct *trace;
  aperture_conf *conf;

  spectrum *obj_spec = NULL, *bck_spec = NULL, *sobj_spec = NULL;
  spectrum *resp;
  response_function *resp_func;
  calib_function *wl_calibration;

  fitsfile *OPET_ptr, *BPET_ptr;
  int f_status = 0;

  fitsfile *SPC_ptr, *SPC_opt_ptr, *WHT_ptr;

  gsl_matrix   *weights;
  drzstamp     *modvar;

  int obj_aperID, obj_beamID, objindex;
  int bck_aperID, bck_beamID;

  char table[MAXCHAR], table_path[MAXCHAR];
  //char comment[FLEN_COMMENT];
  int empty;
  int drizzle;
  int quant_cont=0;
  int opt_weights=0;
  int for_grism=0;
  int smooth_conv=0;

  d_point smooth_params;

  double exptime;
  double sky_cps;

  drzstamp_dim  dimension;
  gsl_matrix *coverage = NULL;

  if (((argc < 3))
      || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "aXe_PET2SPC Version %s:\n"
	       "             aXe task that produces 1-D, binned spectra using information\n"
	       "             contained in a OPET and a BPET. A 1-D spectrum is generated\n"
	       "             for each beam in each of both the OPET and the BPET files\n"
	       "             The binned background spectra of each beam (order) is then\n"
	       "             subtraced from the corresponding spectra form the OPET. Th\ne"
	       "             background subtraction (and reading a BPET altogether can be\n"
	       "             avoided by using the -noBPET option\n"
	       "             An SPC file, a multi-extension FITS file containing binned\n"
	       "             spectra is produced. Each extension (named after the beam ID,\n"
	       "             e.g. 11B for aperture (object) 11,beam (order) B) contains the\n"
	       "             following columns:\n"
	       "\n"
	       "              LAMBDA      ; the wavelength (in A)\n"
	       "              TCOUNT      ; the total number of counts (in DN)\n"
	       "              TERROR      ; the error in TERRORS (in DN)\n"
	       "              BCOUNT      ; the estimated number of counts from the\n"
	       "                            background (in DN)\n"
	       "              BERROR      ; the error in BCOUNTS (in DN)\n"
	       "              COUNT       ; the estimated number of counts from the\n"
	       "                            object (in DN)\n"
	       "              ERROR       ; the error in COUNTS (in DN)\n"
	       "              FLUX        ; the estimated flux (in erg/s/cm^2/A)\n"
	       "              FERROR      ; the error in FLUX (in erg/s/cm^s/A)\n"
	       "              WEIGHT      ; weight (in pixels)\n"
	       "\n"
	       "             Input FITS mages are looked for in $AXE_IMAGE_PATH\n"
	       "             aXe config file is looked for in $AXE_CONFIG_PATH\n"
	       "             All outputs are writen to $AXE_OUTPUT_PATH\n"
	       "\n"
	       "Usage:\n"
	       "     aXe_PET2SPC [g/prism image filename] [aXe config file name] [options]\n"
	       "\n"
	       "Options:\n"
	       "             -noBPET         - to disable the use of a BPET file\n"
	       "             -noflux         - to disable the flux calibration\n"
	       "             -drz            - use $AXE_DRIZZLE_PATH to locate the grism-, OAF-\n"
	       "                             - and PET-files instead of $AXE_IMAGE/OUTPUT_PATH\n"
	       "             -in_AF=[string] - overwrite the default input Aperture File name\n"
	       "             -OPET=[string]  - overwrite the default input Object PET file name\n"
	       "             -BPET=[string]  - overwrite the default input Background PET\n"
	       "                               file name\n"
	       "             -out_SPC=[string] - overwrite the default output SPC file name\n"
	       "\n"
	       "Example:\n"
	       "       ./aXe_PET2SPC slim_grism.fits SLIM.conf.A.0\n"
	       "\n",RELEASE);
      exit (1);
    }

  // make a general opening statement
  fprintf (stdout, "aXe_PET2SPC: Starting...\n");

  // get the name of the flt-file
  index = 0;
  strcpy (grism_file, argv[++index]);
  if ((opt = get_online_option ("drz", argc, argv)))
    {
      build_path (AXE_DRIZZLE_PATH, grism_file, grism_file_path);
      drizzle=1;
    }
  else
    {
      build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);
      drizzle=0;
    }

  // get the name of the configuration file
  strcpy (conf_file, argv[++index]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  // load the configuration file
  conf = get_aperture_descriptor (conf_file_path);

  // Determine where the various extensions are in the FITS file
  get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);


  // Build aperture file name
  replace_file_extension (grism_file, aper_file, ".fits",
			  ".OAF", conf->science_numext);
  if (drizzle)
    build_path (AXE_DRIZZLE_PATH, aper_file, aper_file_path);
  else
    build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);

  // Build object PET file name
  replace_file_extension (grism_file, obj_PET_file, ".fits",
			  ".PET.fits", conf->science_numext);
  // make the total filename
  if (drizzle)
    build_path (AXE_DRIZZLE_PATH, obj_PET_file, obj_PET_file_path);
  else
    build_path (AXE_OUTPUT_PATH, obj_PET_file, obj_PET_file_path);


  // Build background PET file name
  replace_file_extension (grism_file, bck_PET_file, ".fits",
			  ".BCK.PET.fits", conf->science_numext);

  // make the total filename
  if (drizzle)
    build_path (AXE_DRIZZLE_PATH, bck_PET_file, bck_PET_file_path);
  else
    build_path (AXE_OUTPUT_PATH, bck_PET_file, bck_PET_file_path);

  // make a non-standard AF name if necessary
  if ((opt = get_online_option ("in_AF", argc, argv)))
    {
      strcpy(aper_file,opt);
      strcpy(aper_file_path,opt);
    }

  // make a non-standard PET name if necessary
  if ((opt = get_online_option ("OPET", argc, argv)))
    {
      strcpy(obj_PET_file,opt);
      strcpy(obj_PET_file_path,opt);
    }

  // make a non-standard BPET name if necessary
  if ((opt = get_online_option ("BPET", argc, argv)))
    {
      strcpy(bck_PET_file,opt);
      strcpy(bck_PET_file_path,opt);
    }

  // set the flagg for no background subtraction
  if ((opt = get_online_option ("noBPET",argc,argv)))
    {
      dobck = 0;
      strcpy(bck_PET_file,"None");
      strcpy(bck_PET_file_path, "None");
    }
  else
    {
      dobck = 1;
    }

  // set the flagg for flux
  if ((opt = get_online_option ("noflux",argc,argv)))
    noflux = 1;
  else
    noflux = 0;

  // check for the weights flagg,
  // set the file name if the flagg is set
  if ((opt = get_online_option ("opt_weights",argc,argv)))
    opt_weights = 1;
  else
    opt_weights = 0;

  // read the trigger for smoothing the sensitivity
  if ((opt = get_online_option ("smooth_conv",argc,argv)))
    smooth_conv = 1;
  else
    smooth_conv = 0;

  if ((opt = get_online_option ("out_SPC", argc, argv)))
    {
      strcpy (SPC_file, opt);
      strcpy (SPC_file_path, opt);

      if (opt_weights)
	{
	  replace_file_extension (SPC_file, SPC_opt_file, ".fits",
			      "_opt.SPC.fits", -1);
	  replace_file_extension (SPC_file_path, SPC_opt_file_path, ".fits",
			      "_opt.SPC.fits", -1);
	  replace_file_extension (SPC_opt_file, WHT_file, ".SPC.fits",
				  ".WHT.fits", -1);
	  replace_file_extension (SPC_opt_file_path, WHT_file_path, ".SPC.fits",
				  ".WHT.fits", -1);
	}
    }
  else
    {
      replace_file_extension (obj_PET_file, SPC_file, ".PET.fits",
			      ".SPC.fits", -1);
      if (drizzle)
	build_path (AXE_DRIZZLE_PATH, SPC_file, SPC_file_path);
      else
	build_path (AXE_OUTPUT_PATH, SPC_file, SPC_file_path);

      if (opt_weights)
	{
	  replace_file_extension (obj_PET_file, SPC_opt_file, ".PET.fits",
				  "_opt.SPC.fits", -1);
	  replace_file_extension (SPC_opt_file, WHT_file, ".SPC.fits",
				  ".WHT.fits", -1);
	  if (drizzle)
	    {
	      build_path (AXE_DRIZZLE_PATH, SPC_opt_file, SPC_opt_file_path);
	      build_path (AXE_DRIZZLE_PATH, WHT_file, WHT_file_path);
	    }
	  else
	    {
	      build_path (AXE_OUTPUT_PATH, SPC_opt_file, SPC_opt_file_path);
	      build_path (AXE_OUTPUT_PATH, WHT_file, WHT_file_path);
	    }
	}
    }

  // check the configuration file for the smoothin keywords
  if (!check_conf_for_smoothing(conf, smooth_conv))
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "aXe_PET2SP: Either the configuration file %s does not contain\n"
        "the necessary keywords for the smoothing (POBJSIZE, SMFACTOR),\n"
        "or one of these keywords has an unreasonable value < 0.0!\n",
        conf_file_path);


  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf (stdout, "aXe_PET2SPC: Input configuration file name:   %s\n",
	   conf_file_path);
  fprintf (stdout, "aXe_PET2SPC: Input Object Aperture file name: %s\n",
	   aper_file_path);
  fprintf (stdout, "aXe_PET2SPC: Input Object PET file name:      %s\n",
	   obj_PET_file_path);
  if (dobck)
    fprintf (stdout, "aXe_PET2SPC: Input Background PET file name:  %s\n",
	     bck_PET_file_path);
  fprintf (stdout, "aXe_PET2SPC: Output SPC file name:            %s\n",
	   SPC_file_path);
  if (opt_weights)
    {
      fprintf (stdout, "aXe_PET2SPC: Computing optimal weights.\n");
      fprintf (stdout, "aXe_PET2SPC: Optimized SPC file name:         %s\n",
	       SPC_opt_file_path);
      fprintf (stdout, "aXe_PET2SPC: Output WHT file name:            %s\n",
	       WHT_file_path);
    }
  if (!noflux)
    {
      fprintf (stdout, "aXe_PET2SPC: Performing flux calibration.\n");
      if (smooth_conv)
	fprintf (stdout, "aXe_PET2SPC: Using smoothed sensitivity curves.\n");
    }
  fprintf (stdout, "\n\n");

  //
  // try to get the descriptor 'exptime' from the 'sci'-extension
  //
  exptime = (double)get_float_from_keyword(grism_file_path, conf->science_numext, conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_file_path, 1, conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;

  //
  // try to get the descriptor 'SKY_CPS' from the 'sci'-extension
  //
  sky_cps = (double)get_float_from_keyword(grism_file_path, conf->science_numext, "SKY_CPS");
  if (isnan(sky_cps))
    sky_cps = 0.0;


  //  Open the OPET file for reading
  fits_open_file (&OPET_ptr, obj_PET_file_path, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PET2SPC: Could not open file: %s\n",
		   obj_PET_file_path);
    }

  // check whether the contamination is quantitative
  quant_cont = check_quantitative_contamination(OPET_ptr);

  if (opt_weights && !quant_cont)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_PET2SPC: The optimal extractions needs quantitative contamination! "
		 " Please re-run aXe with Gauss or Fluxcube contamination!");

  obs = load_dummy_observation ();
  /* Loading the object list */
  fprintf (stdout, "aXe_PET2SPC: Loading object aperture list...");
  oblist = file_to_object_list_seq (aper_file_path, obs);
  fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

  if (dobck)
    {
      //  Open the file for reading
      fits_open_file (&BPET_ptr, bck_PET_file_path, READONLY, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_PET2SP: Could not open file: %s\n",
		       bck_PET_file_path);
	}
    }


  /* Copy the header info from the grism image */
  SPC_ptr = create_SPC_opened (SPC_file_path,1);
  cards = get_FITS_cards_opened(OPET_ptr);
  put_FITS_cards_opened(SPC_ptr, cards);
  free_FITScards(cards);

  if (opt_weights)
    {
      cards = get_FITS_cards_opened(OPET_ptr);

      // open the WHT file and add the header keywords
      WHT_ptr = create_FITSimage_opened (WHT_file_path, 1);
      put_FITS_cards_opened(WHT_ptr, cards);

      // open the opt_SPC and add the header keywords
      SPC_opt_ptr = create_SPC_opened (SPC_opt_file_path,1);
      put_FITS_cards_opened(SPC_opt_ptr, cards);

      // delete the header keywords
      free_FITScards(cards);
    }

  // do something only if there
  // exist valid objects
  i = 0;
  if (oblist!=NULL)
    {

      // turn until the end of the PET is reached
      while (1)
	{
	  empty=0;

	  // Get the PET for this object
	  obj_PET = get_ALL_from_next_in_PET(OPET_ptr, &obj_aperID, &obj_beamID);

	  // load the background PET if requested
	  if (dobck)
	    {
	      bck_PET = get_ALL_from_next_in_PET(BPET_ptr, &bck_aperID, &bck_beamID);
	      if ((bck_aperID!=obj_aperID)||(bck_beamID!=obj_beamID))
		aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			     "Background PET and Object PET extensions are not"
			     " in the same order and cannot be combined.\n");
	    }

	  // end of PET reached: the break condition
	  if ((obj_aperID==-1) && (obj_beamID==-1))
	    break;

	  // signal an empty PET
	  if (obj_PET==NULL)
	    empty=1;

	  // give feedback to the screen
	  fprintf (stdout, "aXe_PET2SPC: BEAM %d; %d%c\n", i, obj_aperID, BEAM(obj_beamID));fflush(stdout);

	  // identify the object which matches the PET
	  objindex =  find_object_in_object_list(oblist,obj_aperID);

	  // look whether we are for grisms or prisms
	  for_grism = check_for_grism (conf_file_path, obj_beamID);
	  wl_calibration  = get_calfunc_for_beam(oblist[objindex]->beams[obj_beamID], for_grism, conf_file_path, conf);

	  // compute the object spectrum
	  obj_spec = bin_naive (obj_PET, oblist[objindex]->beams[obj_beamID].width,
				oblist[objindex]->beams[obj_beamID].orient,	quant_cont);

	  // check for the existence of a background PET
	  if (dobck)
	    {
	      // compute the background spectrum
	      bck_spec = bin_naive (bck_PET,
				    oblist[objindex]->beams[bck_beamID].width,
				    oblist[objindex]->beams[bck_beamID].orient, quant_cont);
	    }
	  else
	    {
	      // create a dummy background spectrum
	      bck_spec = empty_counts_spectrum_copy(obj_spec);
	    }


	  // subtract the background spectrum from the
	  // object (or forground) spectrum
	  sobj_spec = subtract_spectra (obj_spec, bck_spec);

	  if(!noflux)
	    {
	      get_troughput_table_name(conf_file_path,
				       oblist[objindex]->beams[obj_beamID].ID,
				       table);
	      if (strcmp(table,"None"))
		{
		  build_path (AXE_CONFIG_PATH, table, table_path);
		  resp=get_response_function_from_FITS(table_path,2);
		  resp_func = create_response_function(table_path);
		  if (resp->spec_len <2)
		    {
		      aXe_message (aXe_M_WARN1, __FILE__, __LINE__,
				   "Throughput table %s contains only %d"
				   " values. No sensitivity curve was applied.\n",
				   table_path,resp->spec_len);
		    }
		  else
		    {
		      fprintf(stdout,"aXe_PET2SPC: Applying sensitivity contained in %s\n",table_path);
		      smooth_params = get_smooth_pars_for_beam(conf, smooth_conv, oblist[objindex]->beams[obj_beamID]);
		      if (smooth_params.x > 0.0)
		        {
		          // apply a smoothed flux conversion
			        apply_smoothed_response(wl_calibration, for_grism, quant_cont, resp_func, smooth_params, sobj_spec);
			       }
		      else
			   {
			       // apply a normal flux conversion
			       apply_response_function(sobj_spec, resp, quant_cont);
			   }
		    }
		  // free the memory of the
		  // response functions
		  free_spectrum(resp);
		  free_response_function(resp_func);
		}
	    }
	  if (empty!=1)
	    {
	      add_spectra_to_SPC_opened (SPC_ptr, obj_spec, bck_spec,
					 sobj_spec, oblist[objindex]->ID, oblist[objindex]->beams[obj_beamID].ID);

	      /* Copy header from OPET extension into this SPC extension */
	      cards = get_FITS_cards_opened(OPET_ptr);
	      put_FITS_cards_opened(SPC_ptr,cards);
	      free_FITScards(cards);

	    }

	  free_spectrum (bck_spec);
	  free_spectrum (sobj_spec);
	  free_spectrum (obj_spec);


	  if (opt_weights)
	    {
	      // get the dimension in trace length
	      // and crossdispersion
	      dimension = get_all_dims(obj_PET, bck_PET,
				       oblist[objindex]->beams[obj_beamID], dobck);

	      // check for empty PET
	      if (!dimension.resolution)
		{
		  // create dummies in case of empty PET's
		  weights = get_default_weight();
		  modvar  = get_default_modvar();
		}
	      else
		{
		  // prepare the PET's by computing the inverse variance.
		  // Also the trace distances are shifted by 0.5
		  // to get a sampling comparable to the unweighted
		  // extraction
		  prepare_inv_variance(obj_PET, bck_PET, dobck, conf, exptime, sky_cps, 0.0);
		  // compute the inverse variance and the profile
		  // image in the trace distance - crossdispersion plane
		  modvar = compute_modvar(obj_PET, oblist[objindex]->beams[obj_beamID], dimension);

		  // compute the optimal weights
		  weights = comp_allweight(modvar);

		}

	      if (dimension.resolution && empty != 1)
		{
		  sprintf (label, "WHT_%d%c", obj_aperID, BEAM (obj_beamID));
		  gsl_to_FITSimage_opened (weights, WHT_ptr ,0,label);

		  // make and store the default header
		  cards = beam_to_FITScards(oblist[objindex],obj_beamID);
		  put_FITS_cards_opened(WHT_ptr,cards);
		  free_FITScards(cards);
		}

	      // create the optimal weighted
	      // foreground spectrum
	      obj_spec = bin_optimal (obj_PET,oblist[objindex]->beams[obj_beamID],
				      quant_cont, weights, dimension, coverage);



	      // check for the presence of a background PET
	      if (dobck)
		{
		  // create the optimal weighted
		  // background spectrum
		  bck_spec = bin_optimal (bck_PET,oblist[objindex]->beams[obj_beamID],
					  quant_cont, weights, dimension, coverage);
		}
	      else
		{
		  // make a dummy background spectrum
		  bck_spec = empty_counts_spectrum_copy(obj_spec);
		}

	      // release memory
	      gsl_matrix_free(weights);
	      free_drzstamp(modvar);


	  // subtract the background spectrum from the
	  // object (or forground) spectrum
	  sobj_spec = subtract_spectra (obj_spec, bck_spec);

	  if(!noflux)
	    {
	      get_troughput_table_name(conf_file_path,
				       oblist[objindex]->beams[obj_beamID].ID,
				       table);
	    if (strcmp(table,"None"))
	      {
		build_path (AXE_CONFIG_PATH, table, table_path);
		resp=get_response_function_from_FITS(table_path,2);
		resp_func = create_response_function(table_path);
		if (resp->spec_len <2)
		  {
		    aXe_message (aXe_M_WARN1, __FILE__, __LINE__,
				 "Throughput table %s contains only %d"
				 " values. No sensitivity curve was applied.\n",
				 table_path,resp->spec_len);
		  }
		else
		  {
		    fprintf(stdout,"aXe_PET2SPC: Applying sensitivity contained in %s\n",table_path);
		    smooth_params = get_smooth_pars_for_beam(conf, smooth_conv, oblist[objindex]->beams[obj_beamID]);
        if (smooth_params.x > 0.0)
          {
            apply_smoothed_response(wl_calibration, for_grism, quant_cont, resp_func, smooth_params, sobj_spec);
          }
		    else
		      {
			apply_response_function(sobj_spec, resp, quant_cont);
		      }
		  }
		// free the memory
		free_spectrum(resp);
		free_response_function(resp_func);
	      }
	    }
	  if (empty!=1)
	    {
	      add_spectra_to_SPC_opened (SPC_opt_ptr, obj_spec, bck_spec,
					 sobj_spec, oblist[objindex]->ID, oblist[objindex]->beams[obj_beamID].ID);

	      /* Copy header from OPET extension into this SPC extension */
	      cards = get_FITS_cards_opened(OPET_ptr);
	      put_FITS_cards_opened(SPC_opt_ptr,cards);
	      free_FITScards(cards);

	    }

	  free_spectrum (bck_spec);
	  free_spectrum (sobj_spec);
	  free_spectrum (obj_spec);
	    }

	  // free memory
	  free_calib (wl_calibration);
	  if (bck_PET!=NULL)
	    free (bck_PET);
	  if (obj_PET!=NULL)
	    free (obj_PET);
	  i++;
        }
    }

  fits_close_file (SPC_ptr, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PET2SPC: " "Error closing SPC: %s\n",
		   SPC_file_path);
    }

  if (opt_weights)
    {
      fits_close_file (WHT_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_PET2SPC: " "Error closing WHT: %s\n",
		       WHT_file_path);
	}
      fits_close_file (SPC_opt_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_PET2SPC: " "Error closing PSC: %s\n",
		       SPC_opt_file_path);
	}
    }

  if (oblist!=NULL)
    free_oblist (oblist);

  fprintf (stdout, "aXe_PET2SPC: Done...\n");
  exit (0);
}
