/**
 * Subroutines used mostly for working on catalogs
 * in the tasks SEX2GOL and GOL2AF.  
 *
 */
#include "aper_conf.h"
#include "spc_sex.h"
#include "spc_utils.h"

/**
 *  Function: is_valid_entry
 *  The routine checks whether a magnitude is valid
 *  or not. Returns true/false in case the value
 *  is valid/not valid.
 *
 * Parameters:
 *  @param  mag   - the magnitude to be checked
 *
 * Returns:
 *  @return valid - 1/0 for valid/non valid magnitude, respectively
 */
int
is_valid_entry(double mag){

  int valid = 0;

  // Magnitudes must be in the range 0.0 < mag < 90.0.
  // Sextractor usually gives the value 99.0 if the flux
  // is below zero.
  if (mag < 0.0 || mag > 90.0)
    valid=0;
  else
    valid=1;

  return valid;
}

/**
 * Function: get_valid_entries
 *  The routine determines the number of valid
 *  magnitudes in a magnitude array.
 *
 * Parameters:
 *  @param  magnitudes - the vector with magnitude to be checked
 *
 * Returns:
 *  @return nvalid     - the number of valid magnitudes in the vector
 */
int 
get_valid_entries(const gsl_vector *magnitudes){

  int i=0;
  int nvalid = 0;
  
  // go over each entry in the vector
  for (i=0; i < (int)magnitudes->size; i++)
    {
      // check whether an entry is valid,
      // count up 'nvalid' if the entry is valid
      if (is_valid_entry(gsl_vector_get(magnitudes, i)))
        nvalid++;
    }

  // return the number of valid entries
  return nvalid;
}

/**
 * Function: check_worldcoo_input
 * Checks whether the columns for the world coordinates
 * are present in a catalog header. Returns the number
 * of missing columns
 *
 * Parameters:
 * @param  actcatinfo - the catalog header ot be examined
 *
 * Returns:
 * @return checksum   - the number of missing columns
 */
int
check_worldcoo_input(const colinfo * actcatinfo, const int thsky)
{
  int checksum=0;

  // check for each column
  if (!get_columnumber("NUMBER", actcatinfo))
    checksum++;
  if (!get_columnumber("X_WORLD", actcatinfo))
    checksum++;
  if (!get_columnumber("Y_WORLD", actcatinfo))
    checksum++;
  if (!get_columnumber("A_WORLD", actcatinfo))
    checksum++;
  if (!get_columnumber("B_WORLD", actcatinfo))
    checksum++;

  //  if (thsky)
  //    {
  //      if (!get_columnumber("THETA_SKY", actcatinfo))
  //    checksum++;
  //    }
  //  else
  //    {
  if (!get_columnumber("THETA_WORLD", actcatinfo))
    checksum++;
        //    }

  // return the number of missing columns
  return checksum;
}

/**
 * Function: check_imagecoo_input
 * Checks whether the columns for the image coordinates
 * are present in a catalog header. Returns the number
 * of missing columns
 *
 * Parameters:
 * @param  actcatinfo - the catalog header ot be examined
 *
 * Returns:
 * @return checksum   - the number of missing columns
 */
int
check_imagecoo_input(const colinfo * actcatinfo){

  int checksum=0;
  
  // check for each column 
  if (!get_columnumber("X_IMAGE", actcatinfo))
    checksum++;
  if (!get_columnumber("Y_IMAGE", actcatinfo))
    checksum++;
  if (!get_columnumber("A_IMAGE", actcatinfo))
    checksum++;
  if (!get_columnumber("B_IMAGE", actcatinfo))
    checksum++;
  if (!get_columnumber("THETA_IMAGE", actcatinfo))
    checksum++;

  // return the number of missing columns
  return checksum;
}

/**
 * Function: make_GOL_header
 * The function prints the header of a GOL file to
 * an output stream. Besides the standard columns 
 * it composes the names for the magnitude columns
 * and also makes the header entries for those columns
 *
 * Parameters:
 * @param fout       - the output stream
 * @param actcatinfo - the structure with the header information
 *                     of the SExtractor catalog
 * @param waves      - vector with the wavelengths of magnitude columns
 * @param cnums      - vector with the column numbers of magnitude columns
 * @param backwin_cols - vector with the column numbers of background window columns
 * @param modinfo_cols - vector with the column numbers of modelling index numbers
 *
 */
void 
make_GOL_header(FILE *fout, const colinfo * actcatinfo,
                const gsl_vector * waves, const gsl_vector * cnums,
                const px_point backwin_cols, const px_point modinfo_cols)
{

  char str[CATBUFFERSIZE];
  char cname[COLNAMELENGTH];
  int cbase=1;
  int i=0;

  sprintf (str, "#%3i NUMBER       Running object number\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i X_WORLD      Barycenter position along world x axis     [deg]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i Y_WORLD      Barycenter position along world Y axis     [deg]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i A_WORLD      Profile RMS along major axis (world units) [deg]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i B_WORLD      Profile RMS along minor axis (world units) [deg]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i THETA_WORLD  Position angle (CCW/world-x)               [deg]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i X_IMAGE      Object position along x                    [pixel]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i Y_IMAGE      Object position along x                    [pixel]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i A_IMAGE      Profile RMS along major axis               [pixel]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i B_IMAGE      Profile RMS along minor axis               [pixel]\n", cbase++);
  fputs (str, fout);
  sprintf (str, "#%3i THETA_IMAGE  Position angle (CCW/x)                     [deg]\n", cbase++);
  fputs (str, fout);

  if (backwin_cols.x != -1)
    {
      sprintf (str, "#%3i BACKWINDOW  Upper window for background                [pixel]\n", cbase++);
      fputs (str, fout);
      // old code for FORS2-MXU
      // sprintf (str, "#%3i BACKWIN_LOW  Lower window for background                [pixel]\n", cbase++);
      // fputs (str, fout);
    }

  if (modinfo_cols.x != -1)
    {
      sprintf (str, "#%3i MODSPEC      Index for model spectrum                  [pixel]\n", cbase++);
      fputs (str, fout);
    }
  if (modinfo_cols.y != -1)
    {
      sprintf (str, "#%3i MODIMAGE     Index for the object shape                [pixel]\n", cbase++);
      fputs (str, fout);
    }

  if ((int)waves->size == 1 && (int)gsl_vector_get (waves, 0) == 0)
    {
      sprintf (str, "#%3i MAG_AUTO     Kron-like elliptical aperture magnitude    [mag]\n", cbase++);
      fputs (str, fout);
    }
  else
    {
      for (i=0; i < (int)waves->size; i++){
        get_columname(actcatinfo, (int)gsl_vector_get (cnums, i), cname);
        sprintf (str, "#%3i %-13sAB-magnitude at wavelength %9i nm    [mag]\n", cbase+i, cname, (int)gsl_vector_get (waves, i));
        fputs (str, fout);
        
      }
    }
}

/**
 * Function: get_sex_col_descr
 *  The function parses through the catalog and extracts
 *  the header information , which mainly are the columnnames
 *  plus the associated column numbers. The information is stored
 *  and returned in a structure of type 'colinfo'.
 *
 * Parameters:
 *  @param  filename - the catalog name
 *
 * Returns:
 *  @return actinfo  - the header structure with the catalog information
 */
colinfo *
get_sex_col_descr (char *filename)
{
  FILE    *input;
  char     Buffer[MAXCHAR];
  int      num=0;
  int      count=0;
  char    *key;
  colinfo *actinfo;
     

  // allocate memory for the header structure
  actinfo = (colinfo *) malloc (sizeof (colinfo));
  if (actinfo == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "get_sex_col_descr: Could not allocate memory.");

  //  allocate memory for the dummy 'key
  key = (char *) malloc (MAXCHAR * sizeof (char));
  if (key == NULL)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "get_sex_col_descr: Could not allocate memory.");

  // open the catalog file
  input = fopen (filename, "r");
  if (!input)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_sex_col_descr: Could not open %s", filename);
    }
  
  //  go through the catalog file
  //  This is expected one column name in each line
  //  The way source extractor writes headers
  while (NULL != fgets (Buffer, MAXCHAR, input))
    {
      // see whether the current line is part of the header
      if (Buffer[0] == '#')
        {
          // read the column number and column name from the header line          
          sscanf (Buffer, "# %d %21s", &num, key);
          // store column number and column name in the header structure
          strcpy(actinfo->columns[count].name, key);
          actinfo->columns[count++].number = num;
        }
    }
  // fill in the number of columns
  actinfo->numcols = count;

  // free memory
  free (key);
  key = NULL;

  // return the header structure
  return actinfo;
}


/**
 * Function: get_col_value
 * The function identifies the column number for a given column name,
 * and extracts the data value for this column from the vector
 * representation of a catalog line.
 *
 * Parameters:
 *  @param  actcatinfo - the header of the catalog
 *  @param  key        - the column name
 *  @param  v          - a gsl_vector containing the data of a catalog line
 *  @param  fatal      - controls the severity of the error if
 *                       the column does not exist
 *
 * Returns:
 *  @return res        - the data value in the column
 */
double
get_col_value (const colinfo * actcatinfo, const char key[], gsl_vector * v, int fatal)
{
  double res;
  int    col;
  
  // determine the column number for the column name 
  col = get_columnumber(key, actcatinfo);
  // depending on 'fatal' exit or give back NaN
  // if the column does not exist
  if ( (fatal!=0) && (col == 0) )
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "%s not found in object catalog", key);
    }

  if ( (fatal==0) && (col == 0) )
    return GSL_NAN;

  // get the data value from the vector
  res = gsl_vector_get (v, col - 1);

  //  return the data value
  return res;
}

/*
 * Function: get_col_value2
 */
double
get_col_value2 (const colinfo * actcatinfo, const char key[],
                gsl_vector * v, int fatal)
{
  double res;
  int    col;
  
  // determine the column number for the column name 
  col = get_columnumber(key, actcatinfo);
  // depending on 'fatal' exit or give back NaN
  // if the column does not exist
  if ( (fatal!=0) && (col == 0) )
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "%s not found in object catalog", key);
    }

  if ( (fatal==0) && (col == 0) )
    return -1000000;

  // get the data value from the vector
  res = gsl_vector_get (v, col - 1);

  //  return the data value
  return res;
}

/**
 * Function: line_is_valid
 *  The function checks whether a catalog line is
 *  valid. Valid lines have a data entry for each
 *  catalog column
 *
 * Parameters:
 *  @param  actcatinfo - the header of the catalog
 *  @param  line       - line to be examined
 *
 *  
 * Returns:
 *  @return 0/1        - 0 for not valid lines, 1 for a valid line
 */
int
line_is_valid (const colinfo * actcatinfo, char line[])
{

  gsl_vector *v;

  // transform the line into a vector with the data values
  lv1ws (line);
  v = string_to_gsl_array (line);

  // an empty vector indicates a non valid line
  if (v==NULL) return 0;

  // compare the vector length to the number of columns
  // and return '0' for a non valid line if the numbers differ
  if ((int)v->size != actcatinfo->numcols)
    return 0;

  // return '1' for a valid line
  return 1;
}

/**
 * Function: resolve_colname
 *  The function extracts the wavelength from column names.
 *  The column names must have the format "MAG_CNNNNNN...C..."
 *  with C as character and N as a number. 
 *
 * Parameters:
 *  @param  colname - the column name
 *
 * Returns:
 *  @return iwave   - the wavelength extracted from the column name
 */
int
resolve_colname(const char colname[]){

  char tmp[COLNAMELENGTH];
  char **t_err=NULL;
  char *WorkPtr;
  char *WorkPtr2;

  int i, len;
  //int nmagcols = 0;
  int iwave=0;

  // copy everything after 'MAG_C' to tmp
  strcpy(tmp, &colname[5]);

  // place pointer at both ends of tmp
  WorkPtr = tmp;
  WorkPtr2 = tmp+ strlen (tmp);
  *WorkPtr2-- = '\0';
  len = strlen(tmp);

  // check whether the string is a number,
  // and cut away a digit at the end if not
  i = 0;
  while (i++ < len){
    iwave = (int)strtol(WorkPtr,t_err,10);
    if (iwave != 0)
      break;
    *WorkPtr2-- = '\0';
  }

  // return the wavelength
  return iwave;
}

/**
 * Function: has_magnitudes
 *  The functions searches for columns with magnitudes in the
 *  set of column names within the catalog header structure
 *
 * Parameters:
 *  @param  actcatinfo      - the structure with the column names 
 *
 * Returns:
 *  @return hasmags+magauto - the number of columns with magnitudes
 */
int
has_magnitudes(const colinfo * actcatinfo){

  int hasmags = 0, magauto=0;
  int i= 0;
  //int iwlength=0;

  // go over all columns
  for (i = 0; i < actcatinfo->numcols; i++)
    {
      // look whether the column name is suspicious
      if (strstr(actcatinfo->columns[i].name,"MAG_") != NULL)
        {
          // look whether the column name is 'MAG_AUTO'
          if (!strcmp(actcatinfo->columns[i].name,"MAG_AUTO"))
            {
              // increment the counter for mag_auto
              magauto++;
            }
          else
            {
              // look whether it is possible to extract a wavelength
              // from the column name
              if (resolve_colname(actcatinfo->columns[i].name))
                hasmags++;
              // increment the counter for columns with wavelengths coded
            }
        }
    }
  // exit if there is mag_auto plus other magnitude columns
  if (magauto && hasmags)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                 "There exist MAG_AUTO plus additional columns with magnitudes!\n");

  // return the number of magnitude columns
  return (hasmags+magauto);
}

/**
 * Function: has_backwindow
 * The functions searches for columns which describe the apertures
 * around the object to compute the background. This information
 * is stored in the columns :BACKWIN_UPP and :BACKWIN_LOW
 *
 * Parameters:
 *  @param  actcatinfo      - the structure with the column names 
 *
 * Returns:
 *  @return hasmags+magauto - the number of columns with magnitudes
 */
px_point
has_backwindow(const colinfo * actcatinfo)
{

  px_point ret;
  
  ret.x = -1;
  ret.y = -1;

  // check the backwindow column
  ret.x = get_columnumber("BACKWINDOW", actcatinfo);
  ret.y = -1;

  if (!ret.x)
    {
      ret.x=-1;
      ret.y=-1;
    }

  // -----------------------------
  //
  // code for FORS2-MXU
  // ret.x = get_columnumber("BACKWIN_UPP", actcatinfo);
  // ret.y = get_columnumber("BACKWIN_LOW", actcatinfo);
  // if (!ret.x || !ret.y)
  //  {
  //    ret.x=-1;
  //    ret.y=-1;
  //  }
   
  return ret;
}

/**
 * Function: has_modelinfo
 * The function searches for the columns 'MODSPEC' and 'MODIMAGE'.
 * The column numbers are returned. if not found, the column
 * index -1 is returned.
 *
 * Parameters:
 *  @param  actcatinfo      - the structure with the column names 
 *
 * Returns:
 *  @return ret - column numbers for object shape and spectral model
 */
px_point
has_modelinfo(const colinfo * actcatinfo)
{

  px_point ret;
  
  ret.x = -1;
  ret.y = -1;

  // check for each column
  ret.x = get_columnumber("MODSPEC", actcatinfo);
  ret.y = get_columnumber("MODIMAGE", actcatinfo);

  if (!ret.x)
    ret.x=-1;
  if (!ret.y)
    ret.y=-1;

  // return the column numbers
  return ret;
}

/**
 * Function: get_magcols
 *  The function searches for magnitude coulumns in the header structure,
 *  then it extracts the wavelength from column name and stores
 *  column number and wavelength in gsl_vectors.
 *
 * Parameters:
 *  @param  actcatinfo  - the header structure with the column names 
 *  @param  wavelength  - vector with the wavelengths
 *  @param  wavelength  - vector with the column numbers
 *
 * Returns:
 *  @return mags        - the number of columns with magnitudes
 */
int
get_magcols(const colinfo * actcatinfo,
            gsl_vector *wavelength, gsl_vector *colnums){

  int mags = 0;
  int i= 0, j=0;
  int iwlength=0;
  int entry = 0;
  int idest;


  for (i = 0; i < actcatinfo->numcols; i++)
    {
      // look whether the column name is suspicious
      if (strstr(actcatinfo->columns[i].name,"MAG_") != NULL)
        {
          // look whether the column name is "MAG_AUTO"
          if (!strcmp(actcatinfo->columns[i].name,"MAG_AUTO"))
            {
              // store the column number, store wavelength 0.0 to indicate
              // that it is column MAG_AUTO
              gsl_vector_set (colnums, entry, (double)actcatinfo->columns[i].number);
              printf("column number for wave is %i\n", actcatinfo->columns[i].number);
              gsl_vector_set (wavelength, entry++, 0.0);
              // increment the number of magnitude columns
              mags++;
            }
          else
            {
              // look whether it is possible to extract a wavelength
              // from the column name
              iwlength = resolve_colname(actcatinfo->columns[i].name);
              if (iwlength)
                {
                  // search for the correct position of the actual column
                  // in the set of magnitude columns. the magnitude
                  // columns are stored in ascending order 
                  idest=0;
                  while (iwlength > gsl_vector_get(wavelength, idest) && idest < entry)
                    idest++;
                  
                  // check whether the correct position is at the end of the
                  // magnitude column list
                  if (idest == entry)
                    {

                      // append the column to the list
                      gsl_vector_set (colnums, entry, (double)actcatinfo->columns[i].number);
                      gsl_vector_set (wavelength, entry++, (double)iwlength);
                    }

                  // in case that the correct position is somewhere in the middle
                  else
                    {


                      // move the other columns, starting from the last in the
                      // list until the column which is at the correct
                      // position of the new column
                      for (j=entry; j > idest; j--)
                        {
                          gsl_vector_set (colnums, j, gsl_vector_get(colnums, j-1));
                          gsl_vector_set (wavelength, j, gsl_vector_get(wavelength, j-1));
                        } 

                      // put the new column at the correct position
                      gsl_vector_set (colnums, idest, (double)actcatinfo->columns[i].number);
                      gsl_vector_set (wavelength, idest, (double)iwlength);                                   
                      entry++;
                    }

                  // increment the number of magnitude columns
                  mags++;
                }
            }
        }
    }
  // return the number of magnitude columns
  return mags;
}

/**
 * Function: get_columname
 *  The function extracts a specific column name
 *  from the header structure.
 *
 * Parameters:
 *  @param  actcatinfo - the header structure with the column names 
 *  @param  colnum     - the column number
 *  @param  cname      - string to store the column name
 *
 */
void
get_columname(const colinfo * actcatinfo, const int colnum, char cname[])
{
  // copy the column name into the buffer string
  strcpy(cname, actcatinfo->columns[colnum-1].name);
}

/**
 * Function: get_columnumber
 *  The function searches in the catalog header structure 
 *  for a column name and returns its  column number.
 *
 * Parameters: 
 *  @param  cname      - the column name 
 *  @param  actcatinfo - the header structure with the column names 
 *
 * Returns:
 *  @return cnum       - the number of columns with magnitudes
 */
int
get_columnumber(const char cname[], const colinfo * actcatinfo)
{
  int cnum = 0;
  int i;

  // go through all columns
  for (i = 0; i < actcatinfo->numcols; i++)
    {
      // compare each column name against the name 
      if (!strcmp(actcatinfo->columns[i].name,cname))
        // store the column number
        cnum = actcatinfo->columns[i].number;
    }
  // return the column number
  return cnum;
}

/**
 * Function: get_magauto_col
 *  The functions selects the catalog column with
 *  the closest wavelength to LAMBDACENTRAL.
 *  The column number is returned.
 *
 * Parameters:
 *  @param  wavelength - vector with the wavelengths of all magnitude columns
 *  @param  colnums    - vector with the column numbers of of all magnitude columns
 *
 * Return:
 *  @return select     - column number of the magnitude closest to lambda_mark
 */
int
get_magauto_col(const gsl_vector *wavelength, const gsl_vector *colnums,
                const double lambda_mark)
{
  int i =0, select=0;
  double minimum = 10E+06;

  // if there is only one magnitude column
  // gsl bugfix: updated to catch 0 size
  if (wavelength->size == 1 ){
    // return the column number of the magnitude column 
    return ((int)gsl_vector_get(colnums, 0));
  }
  else{
    // in case there are more magnitude columns, loop over all
    for (i=0; i < (int)wavelength->size; i++){
      // look whether the difference to LAMBDACENTRAL is 
      // smaller then the actual minimum
      if (fabs(lambda_mark - gsl_vector_get(wavelength,i)) < minimum)
        {
          // store the new minimum and its column number
          minimum = fabs(lambda_mark - gsl_vector_get(wavelength,i));
          select = (int)gsl_vector_get(colnums,i);
        }
    }
  }
  // return the column number
  // printf("column number for lambda mark is %i\n", lambda_mark);
  return select;
}


/**
 * Function: check_libraries
 * A helper function which checks the versin numbers of a few
 * runtime libraries to ensure that the proper versions of libraries
 * are being found.
 */
void
check_libraries (void)
{
  int n;
  int major = 0, minor = 0, micro = 0;
  int gsl_major_version = 0, gsl_minor_version = 0, gsl_micro_version = 0;
  
#define NEEDED_GSL_VERSION "0.9.1"
  
  n = sscanf (GSL_VERSION, "%d.%d.%d", &major, &minor, &micro);
  if (n != 2 && n != 3)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "check_libraries: %s, bad GSL version string\n",
                   GSL_VERSION);
    }
  n = sscanf (NEEDED_GSL_VERSION, "%d.%d.%d", &gsl_major_version,
              &gsl_minor_version, &gsl_micro_version);
  if (n != 2 && n != 3)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "check_libraries: %s, bad NEEEDED GSL version string\n",
                   GSL_VERSION);
    }
  if (((gsl_major_version > major)
       || ((gsl_major_version == major) && (gsl_minor_version > minor))
       || ((gsl_major_version == major) && (gsl_minor_version == minor)
           && (gsl_micro_version >= micro))))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "check_libraries: Need GSL version %s, loaded version is %s.",
                   NEEDED_GSL_VERSION, GSL_VERSION);
    }
}

/*
**        A Pratt-Boyer-Moore string search, written by Jerry Coffin
**  sometime or other in 1991.  Removed from original program, and
**  (incorrectly) rewritten for separate, generic use in early 1992.
**  Corrected with help from Thad Smith, late March and early
**  April 1992...hopefully it's correct this time. Revised by Bob Stout.
**
**  This is hereby placed in the Public Domain by its author.
**
**  10/21/93 rdg  Fixed bug found by Jeff Dunlop
*/
static size_t table[UCHAR_MAX + 1];
static size_t len;
static char *findme;

/*
 * Function: init_search
 * Call this with the string to locate to initialize the table
 *
 * Parameters:
 * @param string - no idea
 */
void
init_search (const char *string)
{
  size_t i;
  
  len = strlen (string);
  for (i = 0; i <= UCHAR_MAX; i++)      /* rdg 10/93 */
    table[i] = len;
  for (i = 0; i < len; i++)
    table[(unsigned char) string[i]] = len - i - 1;
  findme = (char *) string;
}

/*
 * Function: strsearch
 * Call this with a buffer to search
 *
 * Parameters:
 * @param string - no idea
 * 
 * Returns:
 * @return here - no idea
 */
char *
strsearch (const char *string)
{
  register size_t shift = 0;
  register size_t pos = len - 1;
  char *here;
  size_t limit = strlen (string);

  while (pos < limit)
    {
      while (pos < limit
             && (shift = table[(unsigned char) string[pos]]) > 0)
        {
          pos += shift;
        }
      if (0 == shift)
        {
          if (0 ==
              strncmp (findme, here =
                       (char *) &string[pos - len + 1], len))
            {
              return (here);
            }
          else
            pos++;
        }
    }
  return NULL;
}


/** 
 * Function: isnum2
 * determine if a string is a number 
 *
 * Parameters:
 * @param str a pointer to a string to test
 *
 * Retruns:
 * @return 1 if string contains a numeric entry, 0 otherwise
 */
int
isnum2 (char *string)
{
  int lstr, i, nd;
  char cstr, cstr1;
  int fpcode;
  
  /* Return 0 if string is NULL */
  if (string == NULL)
    return (0);
  
  lstr = strlen (string);
  nd = 0;
  fpcode = 1;

  /* Return 0 if string starts with a D or E */
  cstr = string[0];
  if (cstr == 'D' || cstr == 'd' || cstr == 'E' || cstr == 'e')
    {
      return (0);
    }

  /* Numeric strings contain 0123456789-+ and d or e for exponents */
  for (i = 0; i < lstr; i++)
    {
      cstr = string[i];
      if (cstr == '\n')
        break;
      if (cstr == ' ' && nd == 0)
        continue;
      if ((cstr < 48 || cstr > 57) && cstr != '+' && cstr != '-'
          && cstr != 'D' && cstr != 'd' && cstr != 'E' && cstr != 'e'
          && cstr != '.')
        return (0);
      else if (cstr == '+' || cstr == '-')
        {
          if (string[i + 1] == '-' || string[i + 1] == '+')
            return (0);
          else if (i > 0)
            {
              cstr1 = string[i - 1];
              if (cstr1 != 'D' && cstr1 != 'd' && cstr1 != 'E'
                  && cstr1 != 'e' && cstr1 != ' ')
                return (0);
            }
        }
      else if (cstr >= 47 && cstr <= 57)
        nd++;
      if (cstr == '.' || cstr == 'd' || cstr == 'e' || cstr == 'd'
          || cstr == 'e')
        fpcode = 2;
    }
  if (nd > 0)
    return (fpcode);
  else
    return (0);
}




/**
 * Function: string_to_gsl_array
 * Parse the content of a  string into a float GSL array.
 * Checks each entry to ensure that it is numeric. skip over
 * the non numerical ones. Entries containing ###,NaN, -NaN, +NaN are
 * set to GSL_NAN.
 *
 * Parameters:
 * @param str the pointer to a string containing floats to parse.
 *
 * Returns:
 * @return v - the gslvector with the table row
 */
gsl_vector *
string_to_gsl_array (char *str)
{
  char *ptr, buf[256];
  gsl_vector *v;
  int i, n = 0;
  

  if (str == NULL)
    return NULL;
  
  ptr = rmlead (str);
  lv1ws (ptr);
  
  /* Find the number of elements */
  do
    {
      ptr = stptok (ptr, buf, sizeof (buf), " ");
      if (isnum2 (buf))
        n++;
      if (!strcmp (buf, "###"))
        n++;            /* Allow entries to be ### */
      if (!strcmp (buf, "NaN"))
        n++;            /* Allow entries to be ### */
      if (!strcmp (buf, "nan"))
        n++;            /* Allow entries to be ### */
      if (!strcmp (buf, "-NaN"))
        n++;            /* Allow entries to be ### */
      if (!strcmp (buf, "+NaN"))
        n++;            /* Allow entries to be ### */
      
    }
  while (ptr && *ptr);
  
  if (n<1) return NULL;
  
  /* Allocate a vector */
  v = gsl_vector_alloc (n);
  //     gsl_vector_set_zero (v);
  /* Fill the vector */
  ptr = str;
  lv1ws (ptr);
  
  i = 0;
  do
    {
      ptr = stptok (ptr, buf, sizeof (buf), " ");
      if (isnum2 (buf))
        {
          gsl_vector_set (v, i, atof (buf));
          i++;
        }
      if ((!strcmp (buf, "###")) || (!strcmp (buf, "NaN"))
          || (!strcmp (buf, "-NaN")) || (!strcmp (buf, "+NaN"))
          || (!strcmp (buf, "nan")))
        {
          gsl_vector_set (v, i, GSL_NAN);
          i++;
        }
    }
  while (ptr && *ptr);

  return v;
}

/**
 * Function: stptok
 * A function to recursively tokenize an input string.
 *
 * Parameters:
 * @param s a pointer to the input string.
 * @param tok a pointer to the current token in s.
 * @param toklen the length of the tok.
 * @param brk the separator string to use.
 *
 * Returns:
 * @return a new position in the input string.
 */
char *
stptok (const char *s, char *tok, size_t toklen, char *brk)
{
  char *lim, *b;

  if (!*s)
    return NULL;

  lim = tok + toklen - 1;
  while (*s && tok < lim)
    {
      for (b = brk; *b; b++)
        {
          if (*s == *b)
            {
              *tok = 0;
              for (++s, b = brk; *s && *b; ++b)
                {
                  if (*s == *b)
                    {
                      ++s;
                      b = brk;
                    }
                }
              return (char *) s;
            }
        }
      *tok++ = *s++;
    }
  *tok = 0;
  return (char *) s;
}

/**
 * Function: lv1ws
 * Convert in-place  all blanks in a string into single space.
 *
 * Parameters:
 * @param str - a pointer to a string.
 */
void
lv1ws (char *str)
{
  char *ibuf, *obuf;
  int i, cnt;
  
  if (str)
    {
      ibuf = obuf = str;
      i = cnt = 0;
      
      while (*ibuf)
        {
          if (isspace ((int) *ibuf) && cnt)
            ibuf++;
          else
            {
              if (!isspace ((int) *ibuf))
                cnt = 0;
              else
                {
                  *ibuf = ' ';
                  cnt = 1;
                }
              obuf[i++] = *ibuf++;
            }
        }
      obuf[i] = '\0';
    }
}

/*
 * Function: rmlead
 *  Originally published as part of the MicroFirm Function Library
 *
 *  Copyright 1986, S.E. Margison
 *  Copyright 1989, Robert B.Stout
 *
 *  The user is granted a free limited license to use this source file
 *  to create royalty-free programs, subject to the terms of the
 *  license restrictions specified in the LICENSE.MFL file.
 *
 *  remove leading whitespace from a string
 *
 * Parameters:
 * @param str - pointer to the input string
 *
 * Returns:
 * @param str - pointer to the stripped string
 */
char *
rmlead (char *str)
{
  char *obuf;
  
  if (str)
    {
      for (obuf = str; *obuf && isspace ((int) *obuf); ++obuf);
      if (str != obuf)
        strMove (str, obuf);
    }
  return str;
}


/**
 * Function: load_dummy_observation
 * A function that returns a dummy allocated abservation
 * pointer.
 *
 Returns:
 * @return obs - pointer to the dummy observation
 */
observation *
load_dummy_observation ()
{
     observation *obs;

     obs = malloc (sizeof (observation));
     obs->grism = gsl_matrix_alloc (2, 2);
     obs->pixerrs = gsl_matrix_alloc (2, 2);

     return obs;
}
