/*
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_sex.h"
#include "spc_utils.h"
#include "inout_aper.h"

#ifndef AXE_IMAGE_PATH
  #define AXE_IMAGE_PATH "AXE_IMAGE_PATH"
#endif
#ifndef AXE_OUTPUT_PATH
  #define AXE_OUTPUT_PATH "AXE_OUTPUT_PATH"
#endif
#ifndef AXE_CONFIG_PATH
  #define AXE_CONFIG_PATH "AXE_CONFIG_PATH"
#endif


observation *
load_empty_image (const char *const fname, int hdunum)
{
  observation *obs = malloc (sizeof (observation));

  obs->grism = FITSnaxes_to_gsl (fname, hdunum);
  return obs;
}

int
main (int argc, char *argv[])
{
  char *opt;
  char grism_image[MAXCHAR];
  char grism_image_path[MAXCHAR];

  char sex_catalog[MAXCHAR];
  char sex_catalog_path[MAXCHAR];

  char aper_file[MAXCHAR];
  char aper_file_path[MAXCHAR];

  char conf_file[MAXCHAR];
  char conf_file_path[MAXCHAR];

  float mfwhm, dmag;

  int leaveout_ignored = 1;
  int auto_reorient = 1;
  //int no_orient = 0;
  //int slitless_geom=0;
  int num;
  char ext[MAXCHAR];
  double lambda_mark=0.0;
  int bck_mode=0;

  aperture_conf *conf;
  SexObject **os;
  object **oblist;
  observation *obs;

  if ((argc < 3) || (opt = get_online_option ("help", argc, argv)))
    {
      fprintf (stdout,
               "aXe_GOL2AF:\n"
               "            aXe task to create an Object Aperture File (OAF)\n"
               "            or a Background Aperture File (BAF) (-bck option)  using an\n"
               "            input Sextractor Grism Object List (GOL).\n"
               "            The AF files contain information about each object (aperture)\n"
               "            and information about the various beams (orders) in each of \n"
               "            these apertures. These are simple text files which can be\n"
               "            edited by hand. The IGNORE keywords can be set to 1 if a\n"
               "            particular beam (order) is to be ignored during the extraction\n"
               "            process. The MMAG_EXTRACT and MMAG_MARK keyword of each beam\n"
               "            listed in the aXe configuration file determine if a particluar\n "
               "            beam is flagged to be extracted or completely ignored.\n"
               "            The -dmag switch can be used to adjust these. See the aXe manual\n"
               "            for more details. An extraction width multiplicative factor can\n"
               "            be set with the -mfwhm option (e.g. to generate BAF containing\n"
               "            broader apertures).\n"
               "            By default, when the extraction angle is too close to the\n"
               "            dispersion direction, 90 degrees are added to the extraction\n"
               "            angle. This can be overwritten by using the -no_auto_orient\n"
               "            online parameter.\n"
               "\n"
               "            Input FITS mages are looked for in $AXE_IMAGE_PATH\n"
               "            aXe config file is looked for in $AXE_CONFIG_PATH\n"
               "            All outputs are writen to $AXE_OUTPUT_PATH\n"
               "\n"
               "Usage:\n"
               "        aXe_GOL2AF [g/prism image filename] [aXe filename] [options]\n"
               "\n"
               "Options:\n"
               "             -bck                - to generate a BAF instead of an OAF file\n"
               "             -SCI_hdu=[integer]  - overwrite the default from the \n"
               "                                   configuration file \n"
               "             -mfwhm=[float]      - an extraction width multiplicative factor.\n"
               "             -exclude_faint      - do not even list faint object in the result.\n"
               "             -out_AF=[string]    - overwrites the default output aper filename.\n"
               "             -in_GOL=[string]    - overwrites the default input catalog name.\n"
               "             -auto_orient=[0/1]  - set to 0 to disable the automatic extraction\n"
               "                                   orientation.\n"
               "             -no_orient          - disable tilted extraction (vertical\n"
               "                                   extraction only).\n"
               "             -dmag=[float]       - A number to add to the MMAG_EXTRACT and\n"
               "                                   MMAG_MARK values.\n"
               "\n"
               "Example: ./aXe_GOL2AF slim_grism.fits -bck -dmag=2 -mfwhm=2\n"
               "\n");
      exit (1);
    }

  fprintf (stdout, "aXe_GOL2AF: Starting...\n");

  /* Get the Grism/Prism image name */
  strcpy (grism_image, argv[1]);
  build_path (AXE_IMAGE_PATH, grism_image, grism_image_path);

  /* Get the name of the configuration file */
  strcpy (conf_file, argv[2]);
  build_path (AXE_CONFIG_PATH, conf_file, conf_file_path);

  /* Read the configuration file */
  conf = get_aperture_descriptor (conf_file_path);

  /* Determine where the various extensions are in the FITS file */
  get_extension_numbers(grism_image_path, conf,conf->optkey1,conf->optval1);

  if ((opt = get_online_option ("SCI_hdu", argc, argv)))
    {
      char str[MAXCHAR];
      strcpy(str,opt);
      conf->science_numext = atoi(str);
    }


  /* Get or set up the name of the output Aperture File */
  if ((opt = get_online_option ("out_AF", argc, argv)))
    {
      /* get it */
      strcpy (aper_file, opt);
      strcpy (aper_file_path, opt);
    }
  else
    {
      /* set the name automatically */
      /* Get the name of the aperture file type, OAF or BAF */
      if ((opt=get_online_option ("bck",argc,argv)))
        {
          strcpy (ext,".BAF");
          bck_mode=1;
        }
      else
        {
          strcpy (ext,".OAF");
        }
      strcpy (aper_file, argv[1]);
      replace_file_extension (grism_image, aper_file, ".fits", ext,
                              conf->science_numext);
      build_path (AXE_OUTPUT_PATH, aper_file, aper_file_path);
    }

     /* Set the option extraction width multiplicative factor */
  if ((opt = get_online_option ("mfwhm", argc, argv)))
    {
      {
        /*char str[MAXCHAR];
          strcpy (str, opt);
          sscanf (str, "%f", &mfwhm);*/
        mfwhm = atof(opt);
      }
    }
  else
    {
      mfwhm = 1.0;
    }

  if ((opt = get_online_option ("dmag", argc, argv)))
    {
      {
        /*char str[MAXCHAR];
          strcpy (str, opt);
          sscanf (str, "%f", &mfwhm);*/
        dmag = atof(opt);
      }
    }
  else
    {
      dmag = 0.0;
    }

  /* Set the option to include object that are too faint into the output
     aperture file.
  */
  if ((opt = get_online_option ("exclude_faint", argc, argv)))
    {
      leaveout_ignored = 1;
    }
  else
    {
      leaveout_ignored = 0;
    }

  /* Set the auto orientation mode for extraction on or off, default is on (=1)

  if ((opt = get_online_option ("auto_orient",argc,argv)))
    {
      if (!atof(opt))
        auto_reorient = 0;
    }
  else
    {
      auto_reorient = 1;
    }
  */

  // set the flag for using slitless geometry,
  // to say extraction parameters optimized for
  // slitless geometry
  if ((opt = get_online_option ("slitless_geom",argc,argv))){
      if (!atof(opt)) {
        //slitless_geom = 0;
        auto_reorient = 0;
      }
  } else {
    //slitless_geom = 1;
    auto_reorient = 1;
  }

  /* Set the option to disbale tilted extraction */
  if ((opt = get_online_option ("orient", argc, argv)))
    {
      if (!atof(opt))
        auto_reorient = 2;
    }

  /* Get or generate the name of the sextractor catalog to read */
  if ((opt = get_online_option ("in_GOL", argc, argv)))
    {
      strcpy (sex_catalog, opt);
      strcpy (sex_catalog_path, opt);
    }
  else
    {
      strcpy (sex_catalog, argv[1]);
      replace_file_extension (grism_image, sex_catalog, ".fits", ".cat",
                              conf->science_numext);
      build_path (AXE_OUTPUT_PATH, sex_catalog, sex_catalog_path);
    }

   if ((opt = get_online_option ("lambda_mark", argc, argv)))
    {
      lambda_mark = atof(opt);
    }
   else{
      lambda_mark = 800.0;
   }

   // check the configuration file for the smoothin keywords
   if (!check_conf_for_slitlessgeom(conf, auto_reorient))
   //if (!check_conf_for_slitlessgeom(conf, slitless_geom))
     aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
         "aXe_GOL2AF: Either the configuration file %s does not contain\n"
         "the necessary keywords for the slitless geometry: POBJSIZE,\n"
         "or this keyword has an unreasonable value < 0.0!\n",
         conf_file_path);

  fprintf (stdout,
           "aXe_GOL2AF: Main configuration file name:    %s\n",
           conf_file_path);
  fprintf (stdout,
           "aXe_GOL2AF: Input data file name:            %s\n",
           grism_image_path);
  fprintf (stdout,
           "aXe_GOL2AF: SCI extension number:            %d\n",
           conf->science_numext);
  fprintf (stdout,
           "aXe_GOL2AF: Input grism object list (GOL):   %s\n",
           sex_catalog_path);
  fprintf (stdout,
           "aXe_GOL2AF: Output aperture file name:       %s\n",
           aper_file_path);
  fprintf (stdout,
           "aXe_GOL2AF: FWHM multiplicative factor:      %f\n",
           mfwhm);
  fprintf (stdout,
           "aXe_GOL2AF: Dmag adjustment:      %f\n",
           dmag);
  if (leaveout_ignored)
    fprintf (stdout,
             "aXe_GOL2AF: Ignored object will not "
             "be included in the output aperture file\n");
  if (!leaveout_ignored)
    fprintf (stdout,
             "aXe_GOL2AF: Ignored object will be "
             "included in the output aperture file\n");
  if (auto_reorient != 2)
    fprintf (stdout,
             "aXe_GOL2AF: Tilted extraction will be performed.\n");
  if (auto_reorient == 0)
    fprintf (stdout,
             "aXe_GOL2AF: The extraction geometry follows the object\n");
  if (auto_reorient == 1)
    //if (slitless_geom)
    fprintf (stdout,
        "aXe_GOL2AF: Extraction geometry is optimized for slitless spectroscopy.\n");
    //else
    //  fprintf (stdout,
    //        "aXe_GOL2AF: The extraction geometry follows the object, but may switch.\n");
  if (auto_reorient==2)
    fprintf (stdout,
             "aXe_GOL2AF: NO tilted extraction will be performed (e.g. vertical only).\n");
  fprintf (stdout,
           "aXe_GOL2AF: Wavelength for object selection [nm]: %f\n",
           lambda_mark);
  fprintf(stdout,"\n\n");

  // extend the beams when creating
  // BAF's
  if (bck_mode)
    extend_config_beams(conf);

  fprintf (stdout, "aXe_GOL2AF: Loading object list...");
  os = get_SexObject_from_catalog (sex_catalog_path, lambda_mark);
  fprintf (stdout, "Done.\n");fflush(stdout);

  fprintf (stdout, "aXe_GOL2AF: Reading image...");
  obs = load_empty_image (grism_image_path, conf->science_numext);
  fprintf (stdout, "Done.\n");fflush(stdout);

  fprintf (stdout, "aXe_GOL2AF: Generating aperture list...");

  oblist = SexObjects_to_oblistII(os, obs, conf, conf_file_path, mfwhm, dmag, auto_reorient, bck_mode);

  /*fprintf (stdout, "Done.\n");fflush(stdout);*/

  fprintf (stdout, "aXe_GOL2AF: Writing aperture file...");fflush(stdout);
  num = object_list_to_file (oblist, aper_file_path, leaveout_ignored);
  fprintf (stdout, "%d beams written.\n", num);

  free_SexObjects (os);
  free_oblist (oblist);

  fprintf (stdout, "aXe_GOL2AF: Done...\n");

  exit (0);
}
