/* 
 * Interface for aXe_utils.c
 *
 */

#ifndef aXe_UTILS_H
#define aXe_UTILS_H

#include <unistd.h>
#include "aXe_grism.h"
#include "spc_FITScards.h"

#define MAXNCONF 10

struct Col_Descr
{
     char *ttype;
     char *tform;
     char *tunit;
};


extern gsl_matrix *
simulate_errors(gsl_matrix * img, float exptime, float gain);

extern gsl_matrix *
simulate_errors_t(gsl_matrix * img, const double exptime,
             const double rdnoise);

extern void
apply_DQ(observation * obs, int dqmask);

extern observation *
load_image(const char *const fname, int hdunum_data,
           int hdunum_err, int hdunum_dq, int dqmask,
           float exptime, float gain);

extern observation *
load_image_t(const char *const fname, int hdunum_data,
             int hdunum_err, int hdunum_dq, int dqmask,
             const double exptime, const double rdnoise);

extern observation *
load_sci_image(const char *const fname, int hdunum_data);

extern int
FITSextnum(const char *const fname);

extern gsl_matrix *
FITSnaxes_to_gsl (const char *const fname, int hdunum);

extern gsl_matrix *
FITSimage_to_gsl (const char *const fname, int hdunum, int fatal);

extern px_point
get_npixels (const char *const fname, int hdunum);

extern int
gsl_to_FITSimage (gsl_matrix * data, char filename[],
                  int overwrite, char ID[]);

extern int
gsl_to_FITSimageHDU (gsl_matrix * data, char filename[], int overwrite,
                     char ID[], int hdu_num);

extern void
quad_to_bbox (const px_point * const corners, px_point * const ll,
              px_point * const ur);

extern int
gsl_matrix_to_pgm (const char *const fname, const gsl_matrix * const m);

extern void
build_path (char envpathname[], const char file[], char *filename);

extern char **
alloc_char_arr(int nrows, int ncols);

extern void
free_char_arr(char **char_arr, int nrows);

extern int
build_config_files(char list_file[], char **list_files, int maxnrows);

extern void
replace_file_extension (char infile[], char outfile[], char from_ext[],
                        char to_ext[], int hdu);

extern void
free_oblist (object ** oblist);

extern void
fprintf_object (FILE * output, object * o);

extern char *
get_online_option (char option_name[], int argc, char *argv[]);

extern int
get_online_option2 (char option_name[], char option_value[],
                    int argc, char *argv[]);

extern void
free_observation (observation * obs);

extern int 
get_hdunum_from_hduname(char filename[], char extname[],
                        char keyword[], char keyval[], int extver);

extern FITScards *
get_FITS_cards (char filename[], int hdu);

extern FITScards *
get_FITS_cards_opened (fitsfile *input);

extern void
put_FITS_cards (char filename[], int hdu, FITScards *cards);

extern void
put_FITS_cards_opened (fitsfile *output, FITScards *cards);

extern void
gsl_to_FITSimage_opened (gsl_matrix * data, fitsfile *ouput,
                         int overwrite, char ID[]);

extern fitsfile *
create_FITSimage_opened (char filename[], int overwrite);

extern float
get_float_from_keyword(char *filename, int hdunum, char *keyword);

extern void
create_FITSimage(char *filename,int overwrite);

extern int
drzprep_fitstrans(char filename[], int hdunum, fitsfile *PET_ptr);

extern px_point
get_npixel(const observation *obs);
#endif
