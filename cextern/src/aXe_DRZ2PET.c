/**
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_PET.h"
#include "inout_aper.h"
#include "trace_conf.h"
#include "aper_conf.h"
#include "spc_sex.h"
#include "disp_conf.h"
#include "spc_wl_calib.h"
#include "drz2pet_utils.h"

#include "fitsio.h"

#define AXE_IMAGE_PATH   "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH  "AXE_OUTPUT_PATH"
#define AXE_DRIZZLE_PATH "AXE_DRIZZLE_PATH"
#define AXE_CONFIG_PATH  "AXE_CONFIG_PATH"

int
main(int argc, char *argv[])
{

  //char           *opt;
  char            option_value[MAXCHAR];
  char            grism_image[MAXCHAR];
  char            grism_image_path[MAXCHAR];

  char            conf_file[MAXCHAR];
  char            conf_file_path[MAXCHAR];

  char            aper_file[MAXCHAR];
  char            aper_file_path[MAXCHAR];

  char            PET_file[MAXCHAR];
  char            PET_file_path[MAXCHAR];
  //char            label[MAXCHAR];
  char            hdu_name[MAXCHAR];
  char            keyword[FLEN_COMMENT];

  int             i;
  //int j, flags;

  object        **oblist;
  observation    *obs=NULL,*wobs=NULL;

  ap_pixel       *result = NULL;
  //d_point         pixel;
  aperture_conf  *conf;

  //dispstruct     *disp;
  //calib_function *wl_calibration;

  //FITScards      *cards;

  int             f_status = 0;

  int             bckmode = 0;
  int             opt_extr=0;
  fitsfile       *PET_fitsptr;
  //fitsfile     ME_fitsptr;

  int             index=0;
  int             extver=0;
  int             exp_ext=0;
  int             con_ext=0;
  int             wht_ext=0;
  int             mod_ext=0;
  int             var_ext=0;

  char list_file[MAXCHAR];
  char list_file_path[MAXCHAR];
  FILE *Filelist;
  static char Buffer[LINE_LEN_MAX];
  char *WorkPtr;
  char *WorkPtr2;
  //double refpntx, refpnty, xoffs;
  double lambda0, dlambda;
  char extname[FLEN_COMMENT];



  if ((argc < 3) || (get_online_option2("help", option_value, argc, argv))) {
    fprintf(stdout,
	    "aXe_DRZ2PET Version %s:\n"
	    "           aXe task that produces one Object or Background  Pixel Extraction\n"
	    "           Table (-bck option) from a set of images created with aXe_DRIZZLE.\n"
	    "           The image list as well as input aperture file (AF) and the aXe\n"
	    "           config are all automatically created by the aXe_DRIZZLE task.\n"
	    "           The format and the usage of the PET's created with aXe_DRZ2PET\n"
	    "           is identical to the PET's created by aXe_AF2PET.\n"
	    "\n"
	    "           The image list is looked for in the current directory.\n"
	    "           The images listed in the image list, the aXe config file and the\n"
	    "           aperture file are looked for in$AXE_DRIZZLE_PATH.\n"
	    "           All outputs are writen to $AXE_DRIZZLE_PATH\n"
	    "\n"
	    "Example:\n"
	    "        aXe_DRZ2PET [image list filename] [aXe config filename] [options]\n"
	    "\n"
	    "Options:\n"
	    "           -bck          - generating a PET for the drizzled background images\n"
	    "                           created by aXE_DRIZZLE using a BAF file instead\n"
	    "                           of a OAF file.\n"
	    "       -in_AF=[string]   - overwrite the automatically generated name\n"
	    "                           of the input aperture file\n"
	    "       -out_PET=[string] - overwrite the automatically generated name\n"
	    "                           of the output PET\n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_DRZ2PET: Starting...\n");

  index = 0;

  strcpy (list_file_path, argv[++index]);
  strcpy (list_file, list_file_path);

  /* Get the configuration file name */
  strcpy(conf_file, argv[++index]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);
  conf = get_aperture_descriptor(conf_file_path);

  // check the background flagg
  if (get_online_option2("bck", option_value, argc, argv))
    bckmode = 1;

  // check the background flagg
  if (get_online_option2("opt_extr", option_value, argc, argv))
    opt_extr = 1;

  /* Fix the PET file name */
  if ((get_online_option2("out_PET", option_value, argc, argv)))
    {
      strcpy (PET_file, option_value);
      strcpy (PET_file_path, option_value);
    }
  else
    {
      strcpy(PET_file, list_file);
      if (bckmode) {
	replace_file_extension(list_file, PET_file, ".fits",
			       ".BCK.PET.fits", -1);
      } else {
	replace_file_extension(list_file, PET_file, ".fits",
			       ".PET.fits", -1);
      }
      build_path(AXE_DRIZZLE_PATH, PET_file, PET_file_path);
    }

  /* Fix the aperture file name */
  if ((get_online_option2("in_AF", option_value, argc, argv))){
    strcpy (aper_file, option_value);
    strcpy (aper_file_path, option_value);
  }
  else {
    if (bckmode){
      replace_file_extension(list_file, aper_file, ".list",
			     ".BAF", -1);}
    else{
      replace_file_extension(list_file, aper_file, ".list",
			     ".OAF", -1);}
    build_path(AXE_DRIZZLE_PATH, aper_file, aper_file_path);
  }

  // give feedback onto the screen:
  // report on input and output
  // and also on specific parameter
  fprintf(stdout,
	  "aXe_DRZ2PET: Input Image List:          %s\n",
	  list_file_path);
  fprintf(stdout,
	  "aXe_DRZ2PET: Input configuration file name: %s\n",
	  conf_file_path);
  fprintf(stdout,
	  "aXe_DRZ2PET: Input aperture file name:      %s\n",
	  aper_file_path);
  fprintf(stdout,
	  "aXe_DRZ2PET: Output Pixel Extraction Table \n"
	  "            (PET) file name:                %s\n",
	  PET_file_path);
  if (opt_extr)
    fprintf(stdout,"aXe_DRZ2PET: Computing optimal weights.\n\n");


  /* Loading the object list */
  fprintf(stdout, "aXe_DRZ2PET: Loading object list...");
  oblist = file_to_object_list_seq(aper_file_path, obs);
  fprintf(stdout, "%d objects loaded.\n\n", object_list_size(oblist));

  /* Create a new empty PET file */
  create_PET(PET_file_path, 1);
  PET_fitsptr = create_PET_opened(PET_file_path, 1);



  fprintf(stdout, "aXe_DRZ2PET: Checking files ... ");
  // check whether all necessary extensions do exist
  i = 0;
  Filelist = fopen (list_file_path, "r");
  if (NULL == Filelist){
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s\n",
		 list_file_path);
  }
  while (NULL != fgets (Buffer, LINE_LEN_MAX, Filelist)){

    // cut leading white spaces:
    WorkPtr = Buffer;
    while (isspace ((int) *WorkPtr))
      WorkPtr++;
    //    if (0 == strlen (WorkPtr))
    if (0 == strlen (WorkPtr) || Buffer[0] == '#')
      continue;

    // cut after first item:
    WorkPtr2 = WorkPtr;
    while (!isspace ((int) *WorkPtr2))
      WorkPtr2++;
    *WorkPtr2 = '\0';
    WorkPtr2 = NULL;

    // detemrine the full name of the drizzled grism image
    strcpy (grism_image,WorkPtr);
    build_path (AXE_DRIZZLE_PATH, grism_image, grism_image_path);

    // find the extension numbers for SCI, ERR, CON and load them
    // into the structure 'obs'
    get_extension_numbers(grism_image_path, conf,conf->optkey1,conf->optval1);
    sprintf (hdu_name, "CON");
    con_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
				      conf->optkey1,conf->optval1,extver);
    sprintf (hdu_name, "EXPT");
    exp_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
    				     conf->optkey1,conf->optval1,extver);

  if (conf->science_numext < 0 || conf->errors_numext < 0
      ||  con_ext < 0 || exp_ext < 0)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_DRZ2PET: Grism file: %s does not have all necessary extensions!\n",
		 grism_image_path);

  if (opt_extr)
    {
      sprintf (hdu_name, "MOD");
      mod_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
					conf->optkey1,conf->optval1,extver);
      //      sprintf (hdu_name, "VAR");
      //     var_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
      //conf->optkey1,conf->optval1,extver);
      var_ext=-1;

      //      if (mod_ext < 0 || var_ext < 0)
      //	aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
      //		     "aXe_DRZ2PET: Grism file: %s does not have necessary extensions for optimale extraction!\n",
      //		     grism_image_path);

    }
  }

  fclose(Filelist);
  fprintf(stdout, "Done\n\n");




  i = 0;
  Filelist = fopen (list_file_path, "r");
  if (NULL == Filelist){
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s\n",
		 list_file_path);
  }
  while (NULL != fgets (Buffer, LINE_LEN_MAX, Filelist)){

    // cut leading white spaces:
    WorkPtr = Buffer;
    while (isspace ((int) *WorkPtr))
      WorkPtr++;
    //    if (0 == strlen (WorkPtr))
    if (0 == strlen (WorkPtr) || Buffer[0] == '#')
      continue;

    // cut after first item:
    WorkPtr2 = WorkPtr;
    while (!isspace ((int) *WorkPtr2))
      WorkPtr2++;
    *WorkPtr2 = '\0';
    WorkPtr2 = NULL;

    // detemrine the full name of the drizzled grism image
    strcpy (grism_image,WorkPtr);
    build_path (AXE_DRIZZLE_PATH, grism_image, grism_image_path);
    fprintf(stdout,
	    "aXe_DRZ2PET: Input data file name:      %s\n", grism_image_path);

    // find the extension numbers for SCI, ERR, CON and load them
    // into the structure 'obs'
    get_extension_numbers(grism_image_path, conf,conf->optkey1,conf->optval1);
    sprintf (hdu_name, "CON");
    con_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
    conf->optkey1,conf->optval1,extver);
    obs = load_image_t(grism_image_path, conf->science_numext,
		     conf->errors_numext, con_ext, conf->dqmask, 1.0, 0.0);

    // find the extension numbers for EXPT and load it
    // into the structure 'wobs'
    sprintf (hdu_name, "EXPT");
    exp_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
    				     conf->optkey1,conf->optval1,extver);
    sprintf (hdu_name, "MOD");
    mod_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
    				     conf->optkey1,conf->optval1,extver);
    sprintf (hdu_name, "VAR");
    var_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
				      conf->optkey1,conf->optval1,extver);

    // preliminar for inverse variance weighting!!
    //    var_ext=conf->errors_numext;

    fprintf(stdout, "aXe_DRZ2PET: Using exts %i,%i,%i as SCI,ERR,DQ\n", exp_ext,mod_ext,var_ext);
    wobs = load_image_t(grism_image_path,exp_ext,mod_ext,
			var_ext , -1, 1.0, 0.0);

    // normalize the exposure times to
    // generate  the weights
    normalize_weight(wobs, oblist[i], opt_extr);

    // get initial wavelength and dispersion
    // either from the configuration file an
    // according descriptor
    if (conf->drz_resol > 0.0){
      dlambda = conf->drz_resol;
    }
    else{
      sprintf (keyword, "DLAMBDA");
      dlambda = (double)get_float_from_keyword(grism_image_path, 2, keyword);
    }
    if (conf->drz_lamb0 > 0.0){
      lambda0 = conf->drz_lamb0;
    }
    else{
      sprintf (keyword, "LAMBDA0");
      lambda0 = (double)get_float_from_keyword(grism_image_path, 2, keyword);
    }

    // built the new vector with PET pixels
    result = (ap_pixel*) make_spc_drztable(obs, wobs, oblist[i], lambda0, dlambda);

    // compose the name of the new extension and generate this extension
    get_ID_num(grism_image,extname);
    strcat(extname,"A");
    add_ALL_to_PET(result, extname, PET_fitsptr, 0);

    // copy the keywords OBJECTID and BEAMID to the PET extension
    drzprep_fitstrans(grism_image_path, conf->science_numext, PET_fitsptr);
    //    cards = get_FITS_cards(grism_image_path, conf->science_numext);
    //    put_FITS_cards_opened(PET_fitsptr, cards);
    //    free_FITScards(cards);

    // determine he extension number of the weight
    // image and store the new weights there
    sprintf (hdu_name, "WHT");
    wht_ext = get_hdunum_from_hduname(grism_image_path, hdu_name,
				     conf->optkey1,conf->optval1,extver);
    if (wht_ext < 0){
      wht_ext = gsl_to_FITSimage (wobs->grism, grism_image_path, 0, hdu_name);
      fprintf(stdout, "aXe_DRZ2PET: Weights written to %s %i\n", grism_image_path, wht_ext);
    }
    else{
      wht_ext = gsl_to_FITSimageHDU (wobs->grism, grism_image_path, 0, hdu_name, wht_ext);
      fprintf(stdout, "aXe_DRZ2PET: Weights replaced in %s %i\n", grism_image_path, wht_ext);
    }

    // free the memory
    if (result!=NULL)
      {
	free(result);
	result = NULL;
      }
    free_observation(obs);
    free_observation(wobs);
    i++;
  }

  if (oblist != NULL)
    free_oblist(oblist);

  fclose(Filelist);
  fits_close_file(PET_fitsptr, &f_status);
  if (f_status) {
    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
		"aXe_DRZ2PET: " "Error closing PET: %s ",
		PET_file_path);
  } else {
    transport_cont_modname(grism_image_path, PET_file_path);
    fprintf(stdout, "aXe_DRZ2PET: PET file %s closed.\n", PET_file_path);
  }

  fprintf(stdout, "aXe_DRZ2PET: Done...\n");
  exit(0);
}
