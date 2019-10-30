/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_sex.h"
#include "spc_utils.h"
#include "spc_CD.h"


#ifndef AXE_IMAGE_PATH
  #define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#endif
#ifndef AXE_OUTPUT_PATH
  #define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#endif
#ifndef AXE_CONFIG_PATH
  #define AXE_CONFIG_PATH "AXE_CONFIG_PATH"
#endif

int
main (int argc, char *argv[])
{
  char *opt;
  char direct_image[MAXCHAR];
  char direct_image_path[MAXCHAR];
  int direct_hdunum;

  char grism_image[MAXCHAR];
  char grism_image_path[MAXCHAR];
  int grism_hdunum;

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  char sex_catalog[MAXCHAR];
  char out_sex_catalog[MAXCHAR];
  char out_sex_catalog_path[MAXCHAR];

  aperture_conf *conf;

  struct WorldCoor *d_wcs;
  struct WorldCoor *g_wcs;

  int nodirim=0;
  int distortion=0;


  if ((argc<3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
	       "ST-ECF European Coordinating Facility\n"
	       "aXe_SEX2GOL Version %s: \n"
	       "             aXe task to create a grism/prism object list based either on\n"
	       "\n"
	       "             1. a direct image object list, one direct and one grism/prism\n"
	       "                image\n"
	       "             or \n"
	       "             2. an object list with world coordinates and one grism/prism image\n"
	       "\n"
	       "             All images must have valid WCS keywords to determine the \n"
	       "             position of the objects in image coordinates.\n"
	       "             The format of the object is that of one produced using Sextractor.\n"
	       "             Lines starting with a \";\" are ignored.\n"
	       "\n"
	       "             Input FITS mages are looked for   in $AXE_IMAGE_PATH.\n"
	       "             The aXe config file is looked for in $AXE_CONFIG_PATH.\n"
	       "             All output is written to             $AXE_OUTPUT_PATH.\n"
	       "\n"
	       "             For 1. the following columns must be present in the input catalog:\n"
	       "                NUMBER      Running object number\n"
	       "                X_IMAGE     Object position along x                    [pixel]\n"
	       "                Y_IMAGE     Object position along y                    [pixel]\n"
	       "                X_WORLD     Barycenter position along world x axis     [deg]\n"
	       "                Y_WORLD     Barycenter position along world y axis     [deg]\n"
	       "                A_IMAGE     Profile RMS along major axis               [pixel]\n"
	       "                B_IMAGE     Profile RMS along minor axis               [pixel]\n"
	       "                THETA_IMAGE Position angle (CCW/x)                     [deg]\n"
	       "                A_WORLD     Profile RMS along major axis (world units) [deg]\n"
	       "                B_WORLD     Profile RMS along minor axis (world units) [deg]\n"
	       "                THETA_WORLD Position angle (CCW/world-x)               [deg]\n"
	       "                MAG_AUTO    Kron-like elliptical aperture magnitude    [mag]\n"
	       "\n"
	       "             The IMAGE coordinates are recomputed using the world coordinates\n"
	       "             and the WCS descriptors\n"
	       "\n"
	       "             For 2. the following columns must be present in the input catalog:\n"
	       "                NUMBER      Running object number\n"
	       "                X_WORLD     Barycenter position along world x axis     [deg]\n"
	       "                Y_WORLD     Barycenter position along world y axis     [deg]\n"
	       "                A_WORLD     Profile RMS along major axis (world units) [deg]\n"
	       "                B_WORLD     Profile RMS along minor axis (world units) [deg]\n"
	       "                THETA_WORLD Position angle (CCW/world-x)               [deg]\n"
	       "             or THETA_SKY   Position angle (CCW, east of north)        [deg]\n"
	       "                MAG_AUTO    Kron-like elliptical aperture magnitude    [mag]\n"
	       "\n"
	       "             The IMAGE coordinates are only recomputed if the following columns\n"
	       "             are NOT present:\n"
	       "                A_IMAGE     Profile RMS along major axis               [pixel]\n"
	       "                B_IMAGE     Profile RMS along minor axis               [pixel]\n"
	       "                THETA_IMAGE Position angle (CCW/x)                     [deg]\n"
	       "                X_IMAGE     Object position along x                    [pixel]\n"
	       "                Y_IMAGE     Object position along y                    [pixel]\n"
	       "\n"
	       "             The flag '-no_direct_image' switches to option 2.\n"
	       "\n"
	       "Usage:\n"
	       "      aXe_SEX2GOL [direct image filename] [g/prism image filename] \n"
	       "                  [aXe config filename] [options]\n"
	       "\n"
	       "Options:\n"
	       "             -no_direct_image    - no direct image is given (as desribed in 2.)\n"
	       "             -dir_hdu=[integer]  - overwrites the default direct image\n"
	       "                                   extension to look for the CD matrix.\n"
	       "             -spec_hdu=[integer] - overwrites the default grism/prism image\n"
	       "                                     extension to look for the CD matrix.\n"
	       "             -in_SEX=[string]    - overwrites the default input object\n"
	       "                                   catalog name\n"
	       "             -out_SEX=[string]   - overwrites the default output object \n"
	       "                                   catalog name\n "
	       "\n",RELEASE);

      exit (1);
    }

  fprintf (stdout, "aXe_SEX2GOL: Starting...\n");
  if ( (opt = get_online_option ("no_direct_image", argc, argv)) )
    nodirim = 1;
  if ( (opt = get_online_option ("distortion", argc, argv)) )
    distortion = 1;

  if (nodirim){
    strcpy (grism_image, argv[1]);
    strcpy (conf_file, argv[2]);
  } else{
    strcpy (direct_image, argv[1]);
    strcpy (grism_image, argv[2]);
    strcpy (conf_file, argv[3]);
  }
  if (!nodirim)
    build_path (AXE_IMAGE_PATH, direct_image, direct_image_path);
  build_path (AXE_IMAGE_PATH, grism_image, grism_image_path);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

     /* Read the configuration file */
  conf = get_aperture_descriptor (conf_file_path);

  if (!nodirim){
    if ((opt = get_online_option ("dir_hdu", argc, argv)))
      {
	char str[MAXCHAR];
	strcpy(str,opt);
	sscanf(str,"%d",&direct_hdunum);
      } else {
	get_extension_numbers(direct_image_path, conf,conf->optkey1,conf->optval1);
	direct_hdunum = conf->science_numext; /* Default */
      }
  }

  if ((opt = get_online_option ("spec_hdu", argc, argv)))
    {
      char str[MAXCHAR];
      strcpy(str,opt);
      sscanf(str,"%d",&grism_hdunum);
    } else {
      get_extension_numbers(grism_image_path, conf,conf->optkey1,conf->optval1);
      grism_hdunum = conf->science_numext; /* Default */
    }


  if ((opt = get_online_option ("in_SEX", argc, argv)))
    {
      strcpy (sex_catalog, opt);
    }
  else
    {
      /* Generate a default catalog name */
      if (nodirim)
	replace_file_extension (grism_image_path, sex_catalog, ".fits",
				"_in.cat", grism_hdunum);
      else
	replace_file_extension (direct_image_path, sex_catalog, ".fits",
				".cat", direct_hdunum);
    }

  if ((opt = get_online_option ("out_SEX", argc, argv)))
    {
      strcpy (out_sex_catalog, opt);
      strcpy (out_sex_catalog_path, out_sex_catalog);
    }
  else
    {
      /* Generate a default output catalog name */
      replace_file_extension (grism_image, out_sex_catalog, ".fits",
			      ".cat", grism_hdunum);
      build_path (AXE_OUTPUT_PATH, out_sex_catalog,
		  out_sex_catalog_path);
    }


  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf (stdout,
	   "aXe_SEX2GOL: Main configuration file name:   %s\n",
	   conf_file_path);
  if (!nodirim)
    {
      fprintf (stdout,
	       "aXe_SEX2GOL: Input direct image file name:   %s\n",
	       direct_image_path);
      fprintf (stdout,
	       "aXe_SEX2GOL: Direct image extension number:  %d\n",
	       direct_hdunum);
      if (distortion)
	fprintf (stdout,"aXe_SEX2GOL: Image distortions are taken into account.\n");
      else
	fprintf (stdout,"aXe_SEX2GOL: Image distortions are NOT taken into account.\n");

    }
  fprintf (stdout,
	   "aXe_SEX2GOL: Input g/prism image file name:  %s\n",
	   grism_image_path);
  fprintf (stdout,
	   "aXe_SEX2GOL: g/prism image extension number: %d\n",
	   grism_hdunum);
  fprintf (stdout,
	   "aXe_SEX2GOL: Input catalog name:             %s\n",
	   sex_catalog);
  fprintf (stdout,
	   "aXe_SEX2GOL: Output catalog name:            %s\n",
	   out_sex_catalog_path);


  // Load the CD matrices
  if (!nodirim)
    d_wcs = get_wcs (direct_image_path, direct_hdunum);
  g_wcs = get_wcs (grism_image_path, grism_hdunum);
  if ( (!nodirim && d_wcs==NULL) || (g_wcs==NULL) )
    {
      fprintf(stdout,
	      "aXe_SEX2GOL: No valid WCS information from FITS headers.\n");
    }

  /* We allow for the input WCS coordinates to be recomputed using */
  /* the d_wcs, then all missing img coordinates are also          */
  /* recomputed to keep things consistant. The new WCS coordinates */
  /* in the g_wcs are then computed and the information is updated */
  if (!nodirim)
     catalog_to_wcs (grism_image_path, grism_hdunum, sex_catalog,
		    out_sex_catalog_path, d_wcs, g_wcs,distortion,1,1);
  else
    catalog_to_wcs_nodim (sex_catalog, out_sex_catalog_path, g_wcs,1,1);

  fprintf (stdout, "aXe_SEX2GOL: Done.\n");
  exit (0);
}
