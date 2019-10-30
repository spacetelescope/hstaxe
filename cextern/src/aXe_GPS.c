/*
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

  char            search_beam[MAXCHAR];
  //char            find_beam[MAXCHAR];

  int             i, j, flags;
  int             for_grism;
  object        **oblist;
  observation    *obs;

  ap_pixel       *result = NULL;
  d_point         pixel;
  //tracestruct    *trace;
  aperture_conf  *conf;

  dispstruct     *disp;
  calib_function *wl_calibration;

  //FITScards      *cards;

  //int             f_status = 0;

  //int             bckmode = 0;

  int             found = 0;

  int             xval = 0;
  int             yval = 0;
  beam           *curbeam;

  //fitsfile       *PET_fitsptr;

  double exptime;

  if ((argc < 6) || (opt = get_online_option("help", argc, argv))) {
    fprintf(stdout,
	    "aXe_GPS Version %s:\n"
	    "         aXe task to determine the beam properties of a single pixel on a\n"
	    "         grism image. The tasks lists:\n"
	    "\n"
	    "         1. the wavelength at pixel center\n"
	    "         2. the dispersion at pixel center\n"
	    "         3. the trace distance of the section point\n"
	    "         4. the distance of the pixel center to the section point\n"
	    "         5. the data value of the pixel.\n"
	    "\n"
	    "         aXe_GPS works on the .OAF file. The corresponding OAF file\n"
	    "         must therefore exist (run aXe_SEX2GOL and aXe_GOL2AF) before\n"
	    "         aXe_GPS gives a result.\n"
	    "         For numerical reasons a solution can only be guaranteed within\n"
	    "         the bounding box of the specified beam. The extraction width as\n"
	    "         specified with the parameter '-mfwhm=' in aXe_GOL2AF has an\n"
	    "         influence on the bounding box. In case that you do not get the\n"
	    "         desired information for the pixel of your interest you may repeat\n"
	    "         aXe_GOL2AF with a larger value of '-mfwhm=' to make the bounding\n"
	    "         box sufficiently large.\n"
	    "         The corner points which define the bounding box of the beam are\n"
	    "         listed in the output such that the user can understand why the\n"
	    "         pixel information could not be computed.\n"
	    "\n"
	    "Usage:\n"
	    "      aXe_GPS [g/prism image filename] [aXe config filename] beam_reference xval yval\n"
	    "\n"
	    "\n"
	    "Example:\n"
	    "        aXe_GPS j8qq11k0q_flt.fits ACS.WFC.CHIP2.Cycle11.4.conf 1083A 2400 1080\n"
	    "\n", RELEASE);
    exit(1);
  }
  fprintf(stdout, "aXe_GPS: Starting...\n");

  /* Get the data file name */
  strcpy(grism_image, argv[1]);

  /* Get the configuration file name */
  strcpy(conf_file, argv[2]);
  build_path(AXE_CONFIG_PATH, conf_file, conf_file_path);

  strcpy(search_beam, argv[3]);
  /* Read the configuration file */
  conf = get_aperture_descriptor(conf_file_path);

  xval = atoi(argv[4]);
  yval = atoi(argv[5]);

  if ((opt = get_online_option("xval", argc, argv)))
    xval = atoi(opt);
  if ((opt = get_online_option("yval", argc, argv)))
    yval = atoi(opt);


  /* Determine where the various extensions are in the FITS file */
  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);
  get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);


  build_path(AXE_IMAGE_PATH, grism_image, grism_image_path);
  get_extension_numbers(grism_image_path, conf, conf->optkey1, conf->optval1);  /* Fix the aperture file name */
  if ((opt = get_online_option("in_AF", argc, argv))) {
    strcpy(aper_file, opt);
    strcpy(aper_file_path, opt);
  } else {
  replace_file_extension(grism_image, aper_file, ".fits",
			 ".OAF", conf->science_numext);
  build_path(AXE_OUTPUT_PATH, aper_file, aper_file_path);

  }


  fprintf(stdout,
	  "aXe_GPS: Input configuration file name:   %s\n",
	  conf_file_path);
  fprintf(stdout,
	  "aXe_GPS: Input data file name:            %s\n",
	  grism_image_path);
  fprintf(stdout,
	  "aXe_GPS: SCI extension number:            %d\n",
	  conf->science_numext);
  fprintf(stdout,
	  "aXe_GPS: ERR extension number:            %d\n",
	  conf->errors_numext);
  fprintf(stdout,
	  "aXe_GPS: DQ extension number:             %d\n",
	  conf->dq_numext);
  fprintf(stdout,
	  "aXe_GPS: DQ mask:                         %d\n",
	  conf->dqmask);
  fprintf(stdout,
	  "aXe_GPS: Input aperture file name:        %s\n",
	  aper_file_path);
  fprintf(stdout,
	  "aXe_GPS: Reference beam:                  %s\n",
	  search_beam);
  fprintf(stdout,
	  "aXe_GPS: Run GPS on pixel:          (%i,%i)\n",
	  xval, yval);
  fprintf(stdout, "\n\n");

  //
  // try to get the descriptor 'exptime' from the 'sci'-extension
  //
  exptime = (double)get_float_from_keyword(grism_image_path, conf->science_numext, conf->exptimekey);
  if (isnan(exptime))
    exptime = (double)get_float_from_keyword(grism_image_path, 1, conf->exptimekey);
  if (isnan(exptime))
    exptime = 1.0;
  obs = load_image_t(grism_image_path, conf->science_numext,
		     conf->errors_numext, conf->dq_numext,
		     conf->dqmask, exptime, conf->rdnoise);



  /* Loading the object list */
  fprintf(stdout, "aXe_GPS: Loading object list...");
  oblist = file_to_object_list_seq(aper_file_path, obs);
  fprintf(stdout, "%d objects loaded.\n", object_list_size(oblist));

  i = 0;
  if (oblist != NULL){
    while (oblist[i] != NULL && found ==0) {
      for (j = 0; j < oblist[i]->nbeams; j++) {
	char            ID[60];
	if (found==1)
	  continue;
	sprintf(ID, "%d%c", oblist[i]->ID, BEAM(oblist[i]->beams[j].ID));
	fprintf(stdout, "aXe_GPS: comparing ID: %s\n",ID);
	/*
	 * skip beam if ignore flag for thisbeam is
	 * set
	 */
	if (!strcmp(ID,search_beam)){
	  fprintf(stdout, "--> Found BEAM_%s\n\n",search_beam);
	  found=1;

	  quad_to_bbox(oblist[i]->beams[j].corners,
		       oblist[i]->beams[j].bbox,
		       oblist[i]->beams[j].bbox + 1);
	  curbeam = oblist[i]->beams + j;
	  result = make_gps_table(oblist[i], j, &flags, xval-1, yval-1);

	  /************************/
	  /* Wavelength calibrate */
	  /************************/

	  /*
	   * check whether it is grism (for_grism=1)
	   * or prism (for_grism=0) data
	   */
	  for_grism = check_for_grism (conf_file_path,oblist[i]->beams[j].ID);

	  /*
	   * get the wavelength dispersion relation at
	   * position "refpoint". conf->refx and conf->refy
	   * are used at this point to allow for a non (0,0) centered
	   * 2D field dependence.
	   */
	  pixel.x = oblist[i]->beams[j].refpoint.x - conf->refx;
	  pixel.y = oblist[i]->beams[j].refpoint.y - conf->refy;
	  disp = get_dispstruct_at_pos(conf_file_path, for_grism,
				       oblist[i]->beams[j].ID,pixel);

	  wl_calibration = create_calib_from_gsl_vector(for_grism, disp->pol);
	  /*
	   * Apply the wavelength calibration to the
	   * PET
	   */
	  wl_calib(result, wl_calibration);
	  free_calib(wl_calibration);

	  if (result[0].p_x == -1 && result[0].p_y == -1){
	    fprintf(stdout, "aXe_GPS: Too far away from the reference point (%7.2f,%7.2f)\n",
		    (oblist[i]->beams[j].refpoint.x+1), (oblist[i]->beams[j].refpoint.y+1));
	    fprintf(stdout,"aXe_GPS: Corners of bounding box for BEAM_%s: (%i, %i), (%i, %i)\n",
		    ID,curbeam->bbox[0].x+1, curbeam->bbox[0].y+1,
		    curbeam->bbox[1].x+1, curbeam->bbox[1].y+1);
	  }
	  else{
	    fprintf(stdout, "aXe_GPS: Grism image: %s  BEAM_%s\n",grism_image_path, ID);
	    fprintf(stdout, "aXe_GPS: SCI extension number:            %d\n",conf->science_numext);
	    fprintf(stdout, "aXe_GPS: Beam reference point: (%7.2f,%7.2f)\n",
		    (oblist[i]->beams[j].refpoint.x+1), (oblist[i]->beams[j].refpoint.y+1));
	    fprintf(stdout,"aXe_GPS: Corners of beam bounding box: (%i, %i), (%i, %i)\n\n",
		    curbeam->bbox[0].x+1, curbeam->bbox[0].y+1,
		    curbeam->bbox[1].x+1, curbeam->bbox[1].y+1);
	    fprintf(stdout, "aXe_GPS: Information for pixel (%i,%i):\n",xval,yval);
	    fprintf(stdout, "-------------------------------------------\n");
	    fprintf(stdout, "aXe_GPS:                    lambda: %8.2f [AA],\n",result[0].lambda);
	    fprintf(stdout, "aXe_GPS:                dispersion: %8.2f [AA/px],\n",result[0].dlambda);
	    fprintf(stdout, "aXe_GPS: trace length of sect. pt.: %8.2f [px],\n",result[0].xi);
	    fprintf(stdout, "aXe_GPS: distance to section point: %8.2f [px],\n",result[0].dist);
	    fprintf(stdout, "aXe_GPS: data value               : %10.3e [cps],\n\n",result[0].count);
	  }

	  free_dispstruct(disp);
	  if (result!=NULL)
	    {
	      free(result);
	      result = NULL;
	    }
	}
      }
      i++;
    }
    if (found == 0)
      fprintf(stdout, "\naXe_GPS: BEAM_%s not found on image %s!\n\n", search_beam, grism_image_path);
  }
  free_observation(obs);
  if (oblist != NULL)
    free_oblist(oblist);
  fprintf(stdout, "aXe_GPS: Done...\n");

  exit(0);
}
