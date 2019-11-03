/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "inout_aper.h"
#include "aper_conf.h"
#include <gsl/gsl_multifit_nlin.h>
#include "trfit_utils.h"


#define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"


int
main(int argc, char *argv[])
{

  char           *opt;

  char            grism_image[MAXCHAR];
  char            grism_image_path[MAXCHAR];

  char            conf_file[MAXCHAR];
  char            conf_file_path[MAXCHAR];

  char            aper_file[MAXCHAR];
  char            aper_file_path[MAXCHAR];


  int             i, j;
  int             num;

  object        **oblist;
  observation    *obs;

  aperture_conf  *conf;

  double exptime;

  gsl_vector *fit_result;

  if ((argc < 1) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
	    "aXe_TFIT Version %s:\n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_TFIT: Starting...\n");

  // Get the data file name
  strcpy(grism_image, argv[1]);

  // Get the configuration file name
  strcpy(conf_file, argv[2]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);

  // Read the configuration file
  conf = get_aperture_descriptor(conf_file_path);

  // Determine where the various extensions are in the FITS file
  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);
  get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);

  // Get or set up the name of the output Aperture File
  if ((opt = get_online_option ("in_AF", argc, argv)))
    {
      /* get it */
      strcpy (aper_file, opt);
      strcpy (aper_file_path, opt);
    }
  else {
    // Build aperture file name
    replace_file_extension (grism_image, aper_file, ".fits",
			    ".OAF", conf->science_numext);
    build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
  }

  fprintf(stdout,
	  "aXe_TFIT: Input configuration file name:   %s\n",
	  conf_file_path);
  fprintf(stdout,
	  "aXe_TFIT: Input data file name:            %s\n",
	  grism_image_path);
  fprintf(stdout,
	  "aXe_TFIT: SCI extension number:            %d\n",
	  conf->science_numext);
  fprintf(stdout,
	  "aXe_TFIT: ERR extension number:            %d\n",
	  conf->errors_numext);
  fprintf(stdout,
	  "aXe_TFIT: DQ extension number:             %d\n",
	  conf->dq_numext);
  fprintf(stdout,
	  "aXe_TFIT: DQ mask:                         %d\n",
	  conf->dqmask);
  fprintf(stdout,
	  "aXe_TFIT: Input aperture file name:        %s\n",
	  aper_file_path);
  fprintf(stdout, "\n\n");

  /* Loading the observation data */
  fprintf(stdout, "aXe_TFIT: ");

  //
  // try to get the descriptor 'exptime' from the 'sci'-extension
  //
  exptime = (double)get_float_from_keyword(grism_image_path,
					   conf->science_numext,
					   conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_image_path, 1,
					     conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;

  // Load an image from AXE_DATA_PATH
  obs = load_image_t(grism_image_path, conf->science_numext,
		     conf->errors_numext, conf->dq_numext,
		     conf->dqmask, exptime, conf->rdnoise);

  // Loading the object list
  fprintf(stdout, "aXe_TFIT: Loading object list...");
  oblist = file_to_object_list_seq(aper_file_path, obs);
  fprintf(stdout, "%d objects loaded.\n", object_list_size(oblist));

  // initialize the counter
  i = 0;

  // check whether there are apertures
  if (oblist != NULL) {

    // go over all apertures
    while (oblist[i] != NULL) {

      // print the current aperture ID
      fprintf(stdout, "aXe_TFIT: Fitting object ID:%d", oblist[i]->ID);

      // go over all beams
      for (j = 0; j < oblist[i]->nbeams; j++) {

	// skip beam if ignore flag for thisbeam is set
	if (oblist[i]->beams[j].ignore !=0 || oblist[i]->beams[j].ID != 0)
	  {
	    /*fprintf(stdout,", %c ( Ignored )",
		    BEAM(oblist[i]->beams[j].ID));*/
	    continue;
	  }

	// print the current beam ID
	/*fprintf(stdout, ", %c (Fitted)", BEAM(oblist[i]->beams[j].ID));*/


	// fit the trace of the beam, modify the
	// beam definition with the new values
	// oblist[i]->beams[j].refpoint.y = fit_beamtrace(conf, obs, j, oblist[i]->beams[j]);
	fit_result = fit_beamtrace(conf, obs, j, oblist[i]->beams[j]);
	oblist[i]->beams[j].refpoint.x = gsl_vector_get(fit_result,0);
	oblist[i]->beams[j].refpoint.y = gsl_vector_get(fit_result,1);
	oblist[i]->beams[j].width = gsl_vector_get(fit_result,4);

	gsl_vector_free(fit_result);
      }

      // give feedback on the progress,
      // enhance the counter
      fprintf(stdout, " Done.\n");
      i++;
    }
  }

  // write the modified beam definitions
  // to an OAF file
  fprintf (stdout, "aXe_TFIT: Writing new aperture file...");fflush(stdout);
  num = object_list_to_file (oblist, aper_file_path, 1);
  fprintf (stdout, "%d beams written.\n", num);

  // release the memory
  free_observation(obs);
  if (oblist != NULL){
	free_oblist(oblist);}

  // print a last status message and exit
  fprintf(stdout, "aXe_TFIT: Done...\n");
  exit(0);
}
