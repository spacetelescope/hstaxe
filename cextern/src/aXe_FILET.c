#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "aXe_grism.h"
#include "aXe_errors.h"
#include "aXe_utils.h"

#define AXE_IMAGE_PATH   "AXE_IMAGE_PATH"
#define AXE_OUTPUT_PATH  "AXE_OUTPUT_PATH"
#define AXE_CONFIG_PATH  "AXE_CONFIG_PATH"
#define AXE_DRIZZLE_PATH "AXE_DRIZZLE_PATH"

int
get_file_name (char *filename, int opt_extr, char *fits_val, int backmode, char *out_image);

char *
get_file_name2 (char *filename, char *fits_val, int backmode, char *out_image);

char *
get_ID_number (char *fits_val, char *ID_num);


int
main(int argc, char *argv[])
  {
    char           *opt;
    char            fits_image[MAXCHAR];
    char            fits_image_path[MAXCHAR];

    char            out_image[MAXCHAR];
    char            out_image_path[MAXCHAR];

    char            out_image_dir[MAXCHAR];

    char            fits_comm[FLEN_COMMENT];
    char            fits_key[FLEN_COMMENT];
    char            fits_val[FLEN_COMMENT];

    fitsfile       *input;
    fitsfile       *ouput;
    FILE           *in_file;


    //int             j, flags;
    int             i;
    int             f_status=0;
    int             hdunum, hdutype;
    int             backmode=0;
    int             opt_extr=0;
    int             copy_ext=0;

    float           exptime = 0.0;

    //char           *outim;

    if ((argc < 2) || (opt = get_online_option("help", argc, argv))) {
      fprintf(stdout,
          "aXe_FILET Version %s:\n"
          "\n", RELEASE);
      exit(1);
    }
    fprintf(stdout, "aXe_FILET: Starting...\n");

    // get the data file name
    strcpy(fits_image, argv[1]);

    // the full filename to the FITS file
    // print the full pathname to the screen
    build_path(AXE_OUTPUT_PATH, fits_image, fits_image_path);
    fprintf(stdout, "aXe_FILET: Input data file name:  %s\n", fits_image_path);

    if ( (opt = get_online_option("opt_extr", argc, argv)) )
      opt_extr=1;
    else
      opt_extr=0;

    // get the name of the output directory
    if ( (opt = get_online_option ("drztmp", argc, argv)) )
      strcpy (out_image_dir, opt);
    else
      strcpy (out_image_dir, getenv (AXE_DRIZZLE_PATH));

    // open the FITS-file
    fits_open_file (&input, fits_image_path, READONLY, &f_status);
    if (f_status)
      {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "aXe_FILET: " " Could not open file: %s\n",
            fits_image_path);
      }

    // read the extension keyword in
    sprintf (fits_key, "EXPTIME");
    fits_read_keyword(input, fits_key, fits_val, fits_comm, &f_status);
    if (f_status)
      {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "Could NOT read exposure time of file: %s\n", fits_image_path);
      }
    exptime = atof(fits_val);
    fprintf(stdout, "aXe_FILET: Exposure time:            %.1f\n", exptime);


    // determine the number of HDU's, print this
    // number ot the screen
    fits_get_num_hdus (input, &hdunum, &f_status);
    if (f_status)
      {
        ffrprt (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "Could get number of extension in file: %s\n", fits_image_path);
      }
    fprintf(stdout, "aXe_FILET: Number of extensions:            %i\n", hdunum-1);

    // determnie the backmode
    if (strstr(fits_image,".BCK.\0") != NULL)
      backmode=1;

    // go over each extention
    i=0;
    while (++i < hdunum){
      // move forward
      fits_movrel_hdu (input, 1, &hdutype, &f_status);

      // read the extension keyword in
      sprintf (fits_key, "EXTNAME");
      fits_read_keyword(input, fits_key, fits_val, fits_comm, &f_status);

      // get the name for the new fits file
      copy_ext = get_file_name(fits_image, opt_extr, fits_val, backmode, out_image);

      if (!copy_ext)
        continue;

      // create the full pathname to the new file
      strcpy (out_image_path, out_image_dir);
      strcat (out_image_path, "/");
      strcat (out_image_path, out_image);

      // open the new fits file
      in_file = fopen (out_image_path, "r");

      // delete an already existing file with the same name
      if (in_file != NULL)
        {
          fclose (in_file);
          unlink (out_image_path);
        }

      // create the fits file new
      fits_create_file (&ouput, out_image_path, &f_status);
      fprintf(stdout, "aXe_FILET: %s[%s] ---> %s\n",fits_image_path, fits_val, out_image_path);
      if (f_status)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
              "aXe_FILET: Could not open" " file: %s\n", out_image_path);
        }

      // copy the current extention to the new fits file
      fits_copy_hdu(input, ouput, 0, &f_status);
      if (f_status)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
              "aXe_FILET: " "Could not copy file: %s\n",
              out_image_path);
        }

      fits_write_key(ouput, TFLOAT, "EXPTIME", &exptime, "exposure time", &f_status);
      if (f_status)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
              "aXe_FILET: " "Could not write exposure time %.1f into file: %s\n",
              exptime, out_image_path);
        }

      // close the new fits file
      fits_close_file (ouput, &f_status);
      if (f_status)
        {
          fits_report_error (stderr, f_status);
          aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
              "aXe_FILET: " "Could not close file: %s\n",
              out_image_path);
        }
    }

    // close the DPP file
    fits_close_file (input, &f_status);
    if (f_status)
      {
        fits_report_error (stderr, f_status);
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
            "aXe_FILET: " "Could not close file: %s\n",
            fits_image_path);
      }

    // give the last feedback
    fprintf(stdout,"aXe_FILET: Done...\n");
    exit(0);
  }

/*
 * Function: get_file_name
 * The function determines the name of the new fits
 * file from the name of the DPP file and the
 * name of the current extension.
 * Example:
 * j8m820leq_flt_2.DPP.fits['CONT_3A'] --> j8m820leq_flt_2_con_ID3.fits
 *
 * Parameters:
 * @param filename  - the name of the DPP file
 * @param opt_extr  - indicates optimal extraction
 * @param fits_val  - the name of the current extension
 * @param backmode  - flagg to indicate background DPP as input
 * @param out_image - the name of the new fits file
 *
 * @return ret - the name of the new fits file
 */
int
get_file_name (char *filename, int opt_extr, char *fits_val, int backmode, char *out_image)
{

  char *ID;
  char *found;
  char in_image[MAXCHAR];
  char tmp_image[MAXCHAR];
  char ID_num[FLEN_COMMENT];

  char beam[5] = "BEAM\0";
  char err[4]  = "ERR\0";
  char cont[5] = "CONT\0";
  char mod[4]  = "MOD\0";
  char var[4]  = "VAR\0";
  //char under[2] = "_\0";
  char dot[2]   = ".\0";

  int ret=0;

  // copy the DPP name into a local variable
  strcpy(in_image,filename);

  // determine the object number from the extension name
  ID = get_ID_number(fits_val, ID_num);

  // cut the name after the '.'
  found = strstr(in_image, dot);
  *found = '\0';

  // check whether the extension name contains 'BEAM'
  if (strstr(fits_val, beam) != NULL)
    {
      // start composing the name for beams
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_flt_ID");
      ret = 1;
    }
  // check whether the extension name contains 'ERR'
  else if (strstr(fits_val, err) != NULL)
    {
      // start composing the name for errors
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_err_ID");
      ret = 1;
    }
  // check whether the extension name contains 'CONT'
  else if (strstr(fits_val, cont) != NULL)
    {
      // start composing the name for contamination
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_con_ID");
      ret = 1;
    }
  // check whether the extension name contains 'MOD'
  else if (strstr(fits_val, mod) != NULL)
    {
      // start composing the name for model
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_mod_ID");
      if (opt_extr)
	ret=1;
      else
	ret=0;
    }
  // check whether the extension name contains 'VAR'
  else if (strstr(fits_val, var) != NULL)
    {
      // start composing the name for model
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_var_ID");
      if (opt_extr)
	ret=1;
      else
	ret=0;
    }
  else
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_FILET: No idea what to do with\nextension %s of image %s!\n",
		   fits_val,filename);
    }

  // append the '.fits'
  strcat(tmp_image,ID_num);
  if (backmode)
    strcat(tmp_image,".BCK.fits");
  else
    strcat(tmp_image,".fits");

  // copy the result to both output channels
  strcpy(out_image, tmp_image);
  //  ret = tmp_image;

  return ret;
}

/*
 * Function: get_file_name2
 * The odl version of the function 'get_file_name'
 *
 */
char *
get_file_name2 (char *filename, char *fits_val, int backmode, char *out_image)
{

  char *ret;
  char *ID;
  char *found;
  char in_image[MAXCHAR];
  char tmp_image[MAXCHAR];
  char ID_num[FLEN_COMMENT];

  char beam[5] = "BEAM\0";
  char err[4]  = "ERR\0";
  char cont[5] = "CONT\0";
  //char mod[4]  = "MOD\0";
  //char under[2] = "_\0";
  char dot[2]   = ".\0";

  strcpy(in_image,filename);


  ID = get_ID_number(fits_val, ID_num);

  found = strstr(in_image, dot);
  *found = '\0';

  found = strstr(fits_val, beam);
  if (found != NULL)
    {
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_flt_ID");
      strcat(tmp_image,ID_num);
      if (backmode)
	strcat(tmp_image,".BCK.fits");
      else
	strcat(tmp_image,".fits");
    }

  found = strstr(fits_val, err);
  if (found != NULL)
    {
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_err_ID");
      strcat(tmp_image,ID_num);
      if (backmode)
	strcat(tmp_image,".BCK.fits");
      else
	strcat(tmp_image,".fits");
    }

  found = strstr(fits_val, cont);
  if (found != NULL)
    {
      strcpy(tmp_image,in_image);
      strcat(tmp_image,"_con_ID");
      strcat(tmp_image,ID_num);
      if (backmode)
	strcat(tmp_image,".BCK.fits");
      else
	strcat(tmp_image,".fits");
    }

  strcpy(out_image, tmp_image);
  ret = tmp_image;
  return ret;
}

/*
 * Function: get_ID_number
 * The function extracts the object number from
 * the name of an extension in a DPP file.
 * Example:
 * extension name 'BEAM_129A' --> object number '129'
 *
 * Parameters:
 * @param fits_val - the extension name
 * @param ID_num   - the number extracted from the extension name
 *
 * Returns:
 * @return WorkPtr2 - the number extracted from the extension name
 */
char *
get_ID_number (char *fits_val, char *ID_num){

  char *WorkPtr;
  char *WorkPtr2;
  char **t_err=NULL;
  char tmp[FLEN_COMMENT];

  int len;
  int i;

  // copy the extension name into the temp-array
  strcpy(tmp, fits_val);

  // allocate an error array
  t_err = (char **) malloc(sizeof(char *)*1);

  // cut whitespaces plus an additional
  // digit from behind
  WorkPtr = tmp+ strlen (tmp);
  *WorkPtr-- = '\0';
  *WorkPtr-- = '\0';
  while (isspace ((int) *WorkPtr))
    *WorkPtr-- = '\0';
  *WorkPtr-- = '\0';
  WorkPtr = NULL;

  // cut one digit from front unti
  // you find an expression which
  // can be converted to a number
  WorkPtr2 = tmp;
  len = strlen(WorkPtr2);
  i = 0;
  while (++i < len){
    strtol(WorkPtr2,t_err,10);
    if (strcmp(WorkPtr2, t_err[0]) != 0)
      break;
    WorkPtr2++;
  }

  // copy the result to
  strcpy(ID_num, WorkPtr2);

  // free allocate dspace
  free(t_err);

  // return the result
  return WorkPtr2;
}
