/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "aXe_errors.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include "fringe_conf.h"
#include "fringe_model.h"


#define AXE_IMAGE_PATH  "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

int
main(int argc, char *argv[])
{
  char        *opt;

  char        fconf_file[MAXCHAR];
  char        fconf_file_path[MAXCHAR];

  char        fimage_file[MAXCHAR];
  char        fimage_file_path[MAXCHAR];

  fringe_conf  *fconf;
  //interpolator *filter_through;
  gsl_matrix   *fringe_image;

  if ((argc < 3) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
	    "aXe_FRIGEN Version %s:\n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_FRIGEN: Starting...\n");

  /* Get the name of the fringe configuration file*/
  strcpy(fconf_file, argv[1]);
  build_path(AXE_CONFIG_PATH, fconf_file, fconf_file_path);

  /* Get the configuration file name */
  strcpy(fimage_file, argv[2]);
  build_path(AXE_OUTPUT_PATH, fimage_file, fimage_file_path);

  /* report on the input and output that will be used */
  fprintf(stdout,
	  "aXe_FRIGEN: Input fringe configuration file:   %s\n",
	  fconf_file_path);
  fprintf(stdout,
	  "aXe_FRIGEN: Output fringe image name:          %s\n",
	  fimage_file_path);

  // load the fringe configuration file
  fconf = load_fringe_conf(fconf_file_path);

  // check whether all necessary information
  // is in place, provide defaults
  check_fringe_conf(fconf);

  // check whether a filter throughput file is set
  if (!fconf->filter_through)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_FRIGEN: No filter throughput table\n"
		 "is set in the fringe configuration file!\n");


  // comute the fringe amplitude
  fringe_image = compute_fringe_amplitude(fconf);

  // save the matrix to a fits-image
  gsl_to_FITSimage (fringe_image, fimage_file_path, 1, NULL);

  // release the space in the fringe image
  gsl_matrix_free(fringe_image);

  // release the allocated memory
  free_fringe_conf(fconf);

  // exit the program
  fprintf(stdout, "aXe_FRIGEN: Done...\n");
  exit(0);
}
