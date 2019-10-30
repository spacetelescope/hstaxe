/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "spc_FITScards.h"
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_PET.h"
#include "inout_aper.h"
#include "trace_conf.h"
#include "aper_conf.h"
#include "disp_conf.h"
#include "drizzle_utils.h"
#include "crossdisp_utils.h"
#include "spc_cfg.h"
#include "inima_utils.h"
#include "spce_output.h"

#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"


int
main (int argc, char *argv[])
{
  char *opt;
  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char list_file[MAXCHAR];
  char list_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char grism_file[MAXCHAR];
  char grism_file_path[MAXCHAR];

  char PET_file[MAXCHAR];
  char PET_file_path[MAXCHAR];

  char SEC_file[MAXCHAR];
  char SEC_file_path[MAXCHAR];

  char outputroot[MAXCHAR];
  char outputroot_path[MAXCHAR];

  char output_path[MAXCHAR];

  aperture_conf *conf;

  object **oblist;
  observation *obs;

  ap_pixel *pri_PET;
  ap_pixel *sec_PET;

  d_point     pixel;
  d_point     refwave_pos;
  d_point     outref;
  d_point     minxy_PET;
  px_point    pixmax;
  dispstruct  *disp, *outdisp;
  trace_func  *trace;
  gsl_matrix  *drizzcoeffs;
  //gsl_vector  *lambdaref = gsl_vector_alloc(2);

  FILE *in_file;

  objectobs  **allobjects;

  int nobjects, nobjects2;
  int boxwidth, boxheight,trlength;
  double relx, rely, objwidth, orient;
  //double min_px, min_py;
  double drizzle_width;
  double cdref, cdscale, cdcorr, cdmeanscale;
  double sprefreso, spreso, spmeanreso, spcorr;
  double sky_cps;
  double exptime;

  double *gaga;

  int index, i, for_grism, in_index;
  char label[MAXCHAR];

  fitsfile *PET_ptr, *SEC_ptr, *DPP_ptr;
  int f_status = 0;

  int aperID, beamID, objindex;
  FITScards *cards;

  drzstamp_dim dimension;
  drzprep *drzprep_stamps;

  axe_inputs *input_list;

  //int rectified = 0;
  //int drizzled  = 0;
  int bckmode   = 0;
  int backpet   = 1;
  //int usemode   = 0;
  int quant_cont= 0;
  int opt_extr  = 0;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "aXe_DRZPREP Version %s:\n"
	       "           aXe task to produce a set of Drizzle PrePare (DPP) files\n"
	       "           for a set of images given in an image list. A DPP-file is\n"
	       "           a multi extension fits file with an pixel stamp image and a\n"
	       "           contamination stamp image for each aperture in the\n"
	       "           grism image. Uses the PET file to derive the pixel/contamination\n"
	       "           values for the stamp images and the aperture files\n"
	       "           (AF's) to define a common geometry for the individual objects,\n"
	       "           which are usually located on several images/PET's. The task also\n"
	       "           derives and stores keywords such that aXe_DRIZZLE can create\n"
	       "           coadded images for all objects in the set of DPP-files.\n"
	       "\n"
	       "           The image list is looked for in the current directory.\n"
	       "           The images listed in the image list is looked for in\n"
	       "           $AXE_IMAGE_PATH. The OAF/BAF- and PET-files computed from\n"
	       "           those images is looked for in $AXE_OUTPUT_PATH\n"
	       "           The DPP's are written to $AXE_OUTPUT_PATH\n"
	       "\n"
	       "Usage:\n"
	       "      aXe_DRZPREP [image list filename] [aXe_config1,aXe_config2,..] [options]\n"
	       "\n"
	       "Options:\n"
	       "           -bck  - generating a background DPP from a background PET\n"
	       "\n",RELEASE);
      exit (1);
    }

  fprintf (stdout, "aXe_DRZPREP: Starting...\n\n");

  index = 0;
  strcpy (list_file, argv[++index]);
  strcpy (list_file_path, list_file);

  strcpy (conf_file, argv[++index]);

  // build the aXe inputs structure
  input_list = get_axe_inputs(list_file_path, conf_file);
  //print_axe_inputs(input_list);

  /* Determine if we are using the special bck mode */
  /* In this mode, file names are handled differently */
  if ((opt = get_online_option("bck", argc, argv)))
    bckmode = 1;
  if ((opt = get_online_option("backpet", argc, argv)))
    backpet = 1;
  if ((opt = get_online_option("opt_extr", argc, argv)))
    opt_extr = 1;


  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf (stdout, "aXe_DRZPREP: Input Image List:           %s\n",
	   list_file_path);
  fprintf (stdout, "aXe_DRZPREP: Configuration file names:   %s\n",
	   conf_file);
  if (opt_extr)
    fprintf (stdout, "aXe_DRZPREP: Storing optimal weighting extentions.\n");
  fprintf(stdout, "\n\n");


  nobjects = 0;
  nobjects2 = 0;

  allobjects = malloc_objectobs();

  fprintf (stdout, "aXe_DRZPREP: Checking input files ...");
  for (in_index=0; in_index < input_list->nitems; in_index++)
    {
      // copy the config file name and the grism file
      // name from the struct to loacal variables
      strcpy(conf_file, input_list->axe_items[in_index].config_file);
      strcpy(grism_file, input_list->axe_items[in_index].grism_file);

      fprintf(stdout, "(%i): %s, %s\n", in_index, conf_file, grism_file);

      // construct the full path to the current configuration fille
      build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);
      for (i=0; i<= DRZMAX; i++)
	{
	  for_grism = check_for_grism (conf_file_path,i);
	  if (!for_grism)
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "aXe_DRZPREP: Beam %c in configuration file %s has no grism dispersion and can not be drizzled.\n",
			 BEAM(i), conf_file_path);
	}

      conf = get_aperture_descriptor (conf_file_path);

      // compose the full name to the grism image
      build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

      // check whether the grism file exists
      in_file = fopen(grism_file_path, "r");
      if (!in_file)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_DRZPREP: File %s does not exist!\n", grism_file_path);
      fclose (in_file);

      get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);


      // find the full name to the PET file
      if (bckmode)
	{
	  /* Build object PET file name */
	  replace_file_extension (grism_file, PET_file, ".fits",
				  ".BCK.PET.fits", conf->science_numext);
	      build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
	}
      else
	{
	  /* Build object PET file name */
	  replace_file_extension (grism_file, PET_file, ".fits",
				  ".PET.fits", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
	}

      //  Open the OPET file for reading
      fits_open_file (&PET_ptr, PET_file_path, READONLY, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_DRZPREP: Could not open file: %s\n",
		       PET_file_path);
	}

      if (!bckmode)
	// check whether there is quantitative contamination
	quant_cont = check_quantitative_contamination(PET_ptr);
      else
	quant_cont=1;

      // check whether optimal extraction can be done
      if (opt_extr && !quant_cont)
	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		     "aXe_DRZPREP: PET file: %s was assembled without quantitative contamination.\n",
		     "Optimal extraction is not possible.", PET_file_path);

      // close the PET file
      fits_close_file (PET_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_DRZPREP: Could not" " close PET file: %s\n", PET_file_path);
	}

      if (bckmode && opt_extr)
	{
	  /* Build object PET file name */
	  replace_file_extension (grism_file, PET_file, ".fits",
				  ".PET.fits", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);

	  //  Open the OPET file for reading
	  fits_open_file (&PET_ptr, PET_file_path, READONLY, &f_status);
	  if (f_status)
	    {
	      ffrprt (stderr, f_status);
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "aXe_DRZPREP: Could not open file: %s\n",
			   PET_file_path);
	    }

	  // check whether there is quantitative contamination
	  quant_cont = check_quantitative_contamination(PET_ptr);

	  if (!quant_cont)
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "aXe_DRZPREP: PET file: %s was assembled without quantitative contamination.\n",
			 "Optimal extraction is not possible.", PET_file_path);

	  // close the PET file
	  fits_close_file (PET_ptr, &f_status);
	  if (f_status)
	    {
	      ffrprt (stderr, f_status);
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "aXe_DRZPREP: Could not" " close PET file: %s\n", PET_file_path);
	    }
	}

    }
  fprintf (stdout, "Done.\n\n");


  for (i=0; i < input_list->nitems; i++)
    {
      // copy the config file name and the grism file
      // name from the struct to loacal variables
      strcpy(conf_file, input_list->axe_items[i].config_file);
      strcpy(grism_file, input_list->axe_items[i].grism_file);

      build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);
      conf = get_aperture_descriptor (conf_file_path);

      build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);
      get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);

      if (bckmode){
	/* Build aperture file name */
	replace_file_extension (grism_file, aper_file, ".fits",
				".OAF", conf->science_numext);
	build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
      }
      else{
	if (backpet)
	  replace_file_extension (grism_file, aper_file, ".fits",
				  ".OAF", conf->science_numext);
	else
	  replace_file_extension (grism_file, aper_file, ".fits",
				  ".OAF", conf->science_numext);
	build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
      }

      fprintf(stdout,
	      "aXe_DRZPREP: Input image file list:    %s\n", list_file_path);
      fprintf(stdout,
	      "aXe_DRZPREP: Input grism image:        %s\n", grism_file_path);
      fprintf(stdout,
	      "aXe_DRZPREP: Input aperture file:      %s\n", aper_file_path);
      fprintf(stdout,
	      "aXe_DRZPREP: Input configuration file: %s\n", conf_file_path);

      obs = load_dummy_observation ();
      /* Loading the object list */
      fprintf (stdout, "aXe_DRZPREP: Loading object aperture list ... ");
      oblist = file_to_object_list_seq (aper_file_path, obs);

      if (oblist != NULL) {
	pixmax = get_npixels (grism_file_path, conf->science_numext);
	fprintf (stdout,"%d objects loaded.\n\n",object_list_size(oblist));
	nobjects2 = add_observation(grism_file_path, conf_file_path,
				    allobjects, nobjects2, oblist,
				    object_list_size(oblist), pixmax,
				    conf->science_numext);
	free_oblist(oblist);
      }
    }

  fprintf (stdout, "aXe_DRZPREP: %i objects in all observations.\n", nobjects2);
  fprintf (stdout, "aXe_DRZPREP: Starting to compute the mean values...");
  f_status =  make_refpoints(conf_file_path, grism_file_path,
			     pixmax, allobjects, nobjects2);
  fprintf (stdout, "Done.\n\n");

  for (in_index=0; in_index < input_list->nitems; in_index++)
    {
      // copy the config file name and the grism file
      // name from the struct to loacal variables
      strcpy(conf_file, input_list->axe_items[in_index].config_file);
      strcpy(grism_file, input_list->axe_items[in_index].grism_file);

      build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);
      conf = get_aperture_descriptor (conf_file_path);


      build_path (AXE_IMAGE_PATH, grism_file, grism_file_path);

      //
      // Start of working on one image
      //
      //
      //

      /* Determine where the various extensions are in the FITS file */
      get_extension_numbers(grism_file_path, conf,conf->optkey1,conf->optval1);
      sky_cps = (double)get_float_from_keyword(grism_file_path, conf->science_numext, "SKY_CPS");
      if (isnan(sky_cps))
	sky_cps = 0.0;
      //
      // try to get the descriptor 'sky_cps' from the 'sci'-extension
      //
      exptime = (double)get_float_from_keyword(grism_file_path, conf->science_numext, conf->exptimekey);
      if (isnan(exptime))
	exptime = (double)get_float_from_keyword(grism_file_path, 1, conf->exptimekey);
      if (isnan(exptime))
	exptime = 1.0;

      if (bckmode)
	{
	  /* Build aperture file name */
	  //	  replace_file_extension (grism_file, aper_file, ".fits",
	  //				  ".BAF", conf->science_numext);
	  replace_file_extension (grism_file, aper_file, ".fits",
				  ".OAF", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);

	  /* Build object PET file name */
	  replace_file_extension (grism_file, PET_file, ".fits",
				  ".BCK.PET.fits", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);

	  if (opt_extr)
	    {
	      /* Build second PET file name */
	      replace_file_extension (grism_file, SEC_file, ".fits",
				      ".PET.fits", conf->science_numext);
	      build_path (AXE_OUTPUT_PATH, SEC_file, SEC_file_path);
	    }
	}
      else
	{
	  /* Build aperture file name */
	  //	  replace_file_extension (grism_file, aper_file, ".fits",
	  //				  ".BAF", conf->science_numext);
	  replace_file_extension (grism_file, aper_file, ".fits",
				  ".OAF", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);

	  /* Build object PET file name */
	  replace_file_extension (grism_file, PET_file, ".fits",
				  ".PET.fits", conf->science_numext);
	  build_path (AXE_OUTPUT_PATH, PET_file, PET_file_path);
	}

      if ( (opt = get_online_option ("outputroot", argc, argv)) ){
	strcpy (outputroot, opt);
      }
      else{
	replace_file_extension (PET_file, outputroot, ".PET.fits", "",-1);
      }

      build_path (AXE_OUTPUT_PATH, outputroot, outputroot_path);
      sprintf(output_path,"%s.DPP.fits",outputroot_path);

      fprintf (stdout, "aXe_DRZPREP: Input PET file name:       %s\n",
	       PET_file_path);
      fprintf (stdout, "aXe_DRZPREP: Input Aperture file name:  %s\n",
	       aper_file_path);
      fprintf (stdout, "aXe_DRZPREP: Name of DPP file :         %s\n",
	       output_path);

      obs = load_dummy_observation ();
      /* Loading the object list */
      fprintf (stdout, "aXe_DRZPREP: Loading object aperture list...");
      oblist = file_to_object_list_seq (aper_file_path, obs);
      fprintf (stdout,"%d objects loaded.\n",object_list_size(oblist));

      //  Open the OPET file for reading
      fits_open_file (&PET_ptr, PET_file_path, READONLY, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_DRZPREP: Could not open file: %s\n",
		       PET_file_path);
	}

      if (!bckmode)
	// check whether there is quantitative contamination
	quant_cont = check_quantitative_contamination(PET_ptr);
      else
	quant_cont=1;


      if (bckmode && opt_extr)
	{
	  //  Open the OPET file for reading
	  fits_open_file (&SEC_ptr, SEC_file_path, READONLY, &f_status);
	  if (f_status)
	    {
	      ffrprt (stderr, f_status);
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "aXe_DRZPREP: Could not open file: %s\n",
			   SEC_file_path);
	    }
	}

      DPP_ptr = create_FITSimage_opened (output_path, 1);
      /* Copy the header info from the grism image */
      {
	FITScards *cards;
	if (bckmode && opt_extr)
	  cards = get_FITS_cards_opened (SEC_ptr);
	else
	  cards = get_FITS_cards_opened (PET_ptr);
	put_FITS_cards_opened(DPP_ptr,cards);
	free_FITScards(cards);
      }

      i = 0;
      if (oblist!=NULL)
	{
	  while (1)
	    {

	      /* Get the PET for this object */
	      pri_PET = get_ALL_from_next_in_PET(PET_ptr, &aperID, &beamID);
	      if ((aperID==-1) && (beamID==-1)) break;

	      if (bckmode && opt_extr)
		sec_PET = get_ALL_from_next_in_PET(SEC_ptr, &aperID, &beamID);
	      else
		sec_PET = NULL;

	      if (beamID >= 0 && beamID <= DRZMAX){
		//	      if (beamID == 0){
		fprintf (stdout, "aXe_DRZPREP: BEAM_%d%c.", aperID, BEAM(beamID));
		objindex =  find_object_in_object_list(oblist,aperID);

		sprintf (label, "%s.%d%c.ps/CPS", outputroot_path,
			 oblist[objindex]->ID, BEAM (oblist[objindex]->beams[beamID].ID));


		pixel.x = oblist[objindex]->beams[beamID].refpoint.x - conf->refx;
		pixel.y = oblist[objindex]->beams[beamID].refpoint.y - conf->refy;
		disp = get_dispstruct_at_pos(conf_file_path, 1,
					     oblist[objindex]->beams[beamID].ID,pixel);
		trace = oblist[objindex]->beams[beamID].spec_trace;

		outref = get_mean_refpoint(allobjects, nobjects2, oblist[objindex]->ID,
					   &boxwidth, &boxheight, &relx, &rely, &objwidth,
					   &orient, &cdref, &cdscale, &cdmeanscale,
					   &sprefreso, &spreso, &spmeanreso, &trlength);

		minxy_PET = get_minxy_from_PET(pri_PET);
		relx = oblist[objindex]->beams[beamID].refpoint.x - minxy_PET.x;
		rely = oblist[objindex]->beams[beamID].refpoint.y - minxy_PET.y;

		//*************************************************
		// patch to correct the reference point in case
		// that the the trace descritpion
		// does have a non negligeable first order term!
		gaga = oblist[objindex]->beams[beamID].spec_trace->data;
		rely = rely + gaga[1];
		//**************************************************

		outdisp = get_dispstruct_at_pos(conf_file_path, 1,
						  oblist[objindex]->beams[beamID].ID,outref);

		cdcorr = cdscale / cdref;
		spcorr = spreso  / sprefreso;

		// determine the trace length
		// NOTE: putting in this line is wrong, but opens to make
		//       drizzle results as in aXe-1.7
		//trlength = get_beam_trace_length(oblist[objindex]->beams[beamID]);

		drizzcoeffs = get_drizzle_coeffs( disp, trace, boxwidth, boxheight,
						  trlength, relx, rely, conf, orient,
						  outdisp, cdcorr, sprefreso, spmeanreso);

		drizzle_width = cdcorr * get_drizzle_width(oblist[objindex],beamID,trace);
		objwidth      = cdmeanscale / cdref * objwidth;

		refwave_pos = get_refwave_position( disp,  trace, pixel, conf);
		refwave_pos.x = refwave_pos.x - minxy_PET.x;
		refwave_pos.y = refwave_pos.y - minxy_PET.y;

		trlength = gsl_matrix_get(drizzcoeffs, 0,10);
		dimension =  get_drzprep_dim(pri_PET, oblist[objindex]->beams[beamID].width,
					     boxwidth, boxheight);
		{
		  px_point tmp_in;
		  px_point tmp_out;

		  d_point new_pos;

		  tmp_in.x = dimension.xsize;
		  tmp_in.y = dimension.ysize;

		  tmp_out.x = (int)trlength;
		  tmp_out.y = 2*(int)ceil(objwidth) + 10;
		  //		    gsl_matrix_fprintf (stdout, drizzcoeffs, "%f");

		  gsl_matrix_set(drizzcoeffs, 0,10,0.0);
		  new_pos = get_drz_position_free(refwave_pos, drizzcoeffs, tmp_in, tmp_out);
		}

		drzprep_stamps = stamp_img_drzprep(opt_extr, pri_PET, sec_PET,
						   oblist[objindex]->beams[beamID].width,
						   -1000000.0, quant_cont, dimension,
						   drizzcoeffs, exptime, sky_cps,
						   (double)conf->rdnoise, bckmode);

		// store the count-extension
		sprintf (label, "BEAM_%d%c",oblist[objindex]->ID,
			 BEAM (oblist[objindex]->beams[beamID].ID));
		gsl_to_FITSimage_opened (drzprep_stamps->counts, DPP_ptr ,0,label);
		cards = drzinfo_to_FITScards(oblist[objindex],beamID,
					     outref, conf, drizzcoeffs,
					     trlength, relx, rely,objwidth,
					     refwave_pos, sky_cps,
					     drizzle_width, cdref, spcorr);
		put_FITS_cards_opened(DPP_ptr,cards);

		// store the error-extension
		sprintf (label, "ERR_%d%c",oblist[objindex]->ID,
			 BEAM (oblist[objindex]->beams[beamID].ID));
		gsl_to_FITSimage_opened (drzprep_stamps->error, DPP_ptr ,0,label);
		put_FITS_cards_opened(DPP_ptr,cards);

		// store the contamination-extension
		sprintf (label, "CONT_%d%c",oblist[objindex]->ID,
			 BEAM (oblist[objindex]->beams[beamID].ID));
		gsl_to_FITSimage_opened (drzprep_stamps->cont, DPP_ptr ,0,label);
		put_FITS_cards_opened(DPP_ptr,cards);

		//		  if (opt_extr)
		if (drzprep_stamps->model)
		  {
		    // store the model-extension
		    sprintf (label, "MOD_%d%c",oblist[objindex]->ID,
			     BEAM (oblist[objindex]->beams[beamID].ID));
		    gsl_to_FITSimage_opened (drzprep_stamps->model, DPP_ptr ,0,label);
		    put_FITS_cards_opened(DPP_ptr,cards);
		  }

		if (drzprep_stamps->vari)
		  {
		    // store the variance extension
		    sprintf (label, "VAR_%d%c",oblist[objindex]->ID,
			     BEAM (oblist[objindex]->beams[beamID].ID));
		    gsl_to_FITSimage_opened (drzprep_stamps->vari, DPP_ptr ,0,label);
		    put_FITS_cards_opened(DPP_ptr,cards);
		  }

		// release memory
		free_drzprep(drzprep_stamps);
		free_FITScards(cards);
		free_dispstruct(disp);
		free_dispstruct(outdisp);
		gsl_matrix_free(drizzcoeffs);

		fprintf (stdout, " Done.\n");
		i++;
	      }

	      if (pri_PET!=NULL) free(pri_PET);
	      if (sec_PET!=NULL) free(sec_PET);
	    }
	}
      if (oblist!=NULL) free_oblist (oblist);

      fits_close_file (DPP_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_DRZPREP: Could not" " close DPP file: %s\n", output_path);
	}

      fits_close_file (PET_ptr, &f_status);
      if (f_status)
	{
	  ffrprt (stderr, f_status);
	  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		       "aXe_DRZPREP: Could not" " close PET file: %s\n", PET_file_path);
	}

      if (bckmode && opt_extr)
	{
	  fits_close_file (SEC_ptr, &f_status);
	  if (f_status)
	    {
	      ffrprt (stderr, f_status);
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "aXe_DRZPREP: Could not" " close PET file: %s\n", SEC_file_path);
	    }
	}

      fprintf (stdout, "aXe_DRZPREP: %s Done...\n\n", output_path);
      //
      // End of working on one image
      //
      //
      //

    }

  free_axe_inputs(input_list);
  free_objectobs(allobjects);
  exit (0);
}
