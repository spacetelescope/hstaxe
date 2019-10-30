#ifndef SPC_UTILS_H
#define SPC_UTILS_H

void init_search (const char *string);  /* Pbmsrch.C      */
char *strsearch (const char *string);   /* Pbmsrch.C      */


#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include <gsl/gsl_version.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_vector.h>

#include "aXe_errors.h"
#include "aXe_grism.h"

#define LAMBDACENTRAL 8000.0
#define COLNAMELENGTH 12
#define MAXCOLS       200

/**
 * Structure: sexcol
 * Structure for a column header
 */
typedef struct sexcol
{
  char name[COLNAMELENGTH];
  int  number;
}
sexcol;

/**
 * Structure: colinfo
 * Structure for the table header
 */
typedef struct colinfo
{
  int numcols;
  sexcol columns[MAXCOLS];
}
colinfo;

extern int
is_valid_entry(double mag);

extern int
get_valid_entries(const gsl_vector *magnitudes);

extern int
check_worldcoo_input(const colinfo * actcatinfo, const int thsky);

extern int
check_imagecoo_input(const colinfo * actcatinfo);

extern void
make_GOL_header(FILE *fout, const colinfo * actcatinfo,
                const gsl_vector * waves, const gsl_vector * cnums,
                const px_point backwin_cols, const px_point modinfo_cols);

extern colinfo *
get_sex_col_descr (char *filename);

extern double
get_col_value (const colinfo * actcatinfo, const char key[],
               gsl_vector * v, int fatal);

extern double
get_col_value2 (const colinfo * actcatinfo, const char key[],
                gsl_vector * v, int fatal);

extern int
line_is_valid (const colinfo * actcatinfo, char line[]);

extern int
get_magauto_col(const gsl_vector *wavelength, const gsl_vector *colnums,
                const double lambda_mark);

extern int
has_magnitudes(const colinfo * actcatinfo);

extern px_point
has_backwindow(const colinfo * actcatinfo);

extern px_point
has_modelinfo(const colinfo * actcatinfo);

extern int 
get_magcols(const colinfo * actcatinfo, gsl_vector *wavelength,
            gsl_vector *colnums);

extern int
resolve_colname(const char colname[]);

extern void
get_columname(const colinfo * actcatinfo, const int colnum, char colname[]);

extern int
get_columnumber(const char cname[], const colinfo * actcatinfo);


extern char *
rmlead (char *str);

extern char *
stptok (const char *s, char *tok, size_t toklen, char *brk);

#define      strMove(d,s) memmove(d,s,strlen(s)+1)

extern void
lv1ws (char *str);

extern int
isnum2(char *string);

extern gsl_vector *
string_to_gsl_array (char *str);

extern void
check_libraries (void);

extern observation *
load_dummy_observation (void);
#endif
