/**
 */
#include "disp_conf.h"
#include <string.h>
#include <math.h>
#include "aXe_errors.h"
#include "spc_cfg.h"

/**
 * Routine to read in the MMAG_EXTRACT keyword for a specific beam
 * from the configuration file
 *
 * @param  filename     filename of the configuration file
 * @param  beamID       beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return mmag_extract mag_extract-value for that beam
 */
float
get_beam_mmag_extract (char *filename, int beamID)
{
  char beam[MAXCHAR];
  float mmag_extract = -1;
  int found = 0;
  int i=0;

  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;

  sprintf (beam, "MMAG_EXTRACT_%c", BEAM (beamID));
  CfgRead (filename, DispConfig);
  if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
    {
      mmag_extract = atof (DispConfig[0].data);
      found = 1;
    }
  if (!found)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_beam_mmag_extract: %s was not found in %s.\n",
                   beam, filename);
    }

  // release memory
  i=0;
  while(DispConfig[i].name!=NULL)
    free(DispConfig[i++].data);

  // return the magnitude
  return mmag_extract;
}

/**
 * Routine to read in the MMAG_MARK keyword for a specific beam
 * from the configuration file
 *
 * @param  filename     filename of the configuration file
 * @param  beamID       beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return mmag_mark mag_mark-value for that beam
 */
float
get_beam_mmag_mark (char *filename, int beamID)
{
  char beam[MAXCHAR];
  float mmag_mark = -1;
  int found = 0;
  int i=0;

  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;

  sprintf (beam, "MMAG_MARK_%c", BEAM (beamID));
  CfgRead (filename, DispConfig);
  if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
    {
      mmag_mark = atof (DispConfig[0].data);
      found = 1;
    }
  if (!found)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_beam_mmag_mark: %s was not found in %s.\n",
                   beam, filename);
    }

  // release memory
  i=0;
  while(DispConfig[i].name!=NULL)
    free(DispConfig[i++].data);

  // return the value
  return mmag_mark;
}

/**
 * Parses a configuration file and searches for DISP_ORDER_%c
 * keywords. Returns the order of the field dependent
 * dispersion relation for that order.
 *
 * @param  filename  filename of the configuration file
 * @param  beamID    beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return disporder the order of the dispersion relation for the beam
 *                   of interest
 */
int
get_beam_disp_norder (char *filename, int beamID)
{
  char beam[MAXCHAR];
  int disporder = -1, found = 0, i;

  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;

  sprintf (beam, "DISP_ORDER_%c", BEAM (beamID));
  CfgRead (filename, DispConfig);
  if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
    {
      disporder = atoi (DispConfig[0].data);
      found = 1;
    }
  if (!found)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_beam_disp_norder: %s was not found in %s.\n",
                   beam, filename);
    }
  i = 0;
  while(DispConfig[i].name!=NULL){
    free(DispConfig[i++].data);
  }

  //free(DispConfig[0].data);
  return disporder;

}

/**
 * Parses a configuration file and returns a gsl vector
 * containing the 2D polynomial coefficient allowing one to compute
 * a given dispersion relation coefficient anywhere on the detector.
 *
 * @param  filename filename of the configuration file
 * @param  beamID   beamID of the beam of interest (A=0, B=1, etc...)
 * @param  order    order of the coefficient for which the 2D field dependent
 *                  coefficients are queried.
 *
 * @return v        a gsl_vector() containing all of the 2D field
 *                  dependent polynomial coeffiecient
 */
gsl_vector *
get_beam_disp_order (char *filename, const int for_grism, int beamID,
                     int order)
{
  char beam[MAXCHAR];

  int norder, found = 0;
  int i;
  gsl_vector *v = NULL;
  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;
  norder = get_beam_disp_norder (filename, beamID);
  if (for_grism == 1)
    {
      sprintf (beam, "DLDP_%c_%d", BEAM (beamID), order);
    }
  else
    {
      sprintf (beam, "DLD1P_%c_%d", BEAM (beamID), order);
    }
  CfgRead (filename, DispConfig);
  if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
    {
      v = string_to_gsl_array (DispConfig[0].data);
      found = 1;
    }

  if (!found)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_beam_disp_order: %s was not found in %s.\n",
                   beam, filename);
    }
  i = 0;
  while(DispConfig[i].name!=NULL){
    free(DispConfig[i++].data);
  }
  return v;
}

/**
 * Parses a configuration file, computes and return the value
 * of the order^th dispersion polynomial coefficient at a given
 * detector pixel.
 *
 * @param filename filename of the configuration file
 * @param beamID   beamID of the beam of interest (A=0, B=1, etc...)
 * @param order    order of the coefficient for which the 2D field dependent
 * @param p        a d_point object containing the position to compute the
 *                 dispersion coefficient at.
 *
 * @return res     the value of the order6th coefficient of beamID computed at
 *                 detector pixel p.
 */
float
get_disp_coeff_at_pos (char *filename, const int for_grism, int beamID,
                       int order, d_point p)
{
  gsl_vector *v;
  float n, c, res;
  int i, j, k;

  v = get_beam_disp_order (filename, for_grism, beamID, order);
  n = 0.5 * (-1.0 + sqrt (1 + 8 * v->size));
  if ((floor (n) - n) != 0)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_disp_coeff_at_pos: "
                   "Order %d of Beam %d in %s does not contain a correct number of entries (i.e. "
                   "1,3,6,10,15...,n^2/2+n/2", order, beamID, filename);
    }

  i = 0;
  res = 0;
  for (j = 0; j < n; j++)
    {
      for (k = 0; k < (j + 1); k++)
        {
          c = gsl_vector_get (v, i);
          res = res + c * pow (p.x, (j - k)) * pow (p.y, k);
          i++;
        }
    }
  gsl_vector_free (v);
  return res;
}



/**
 * Parse a configuration file and return a gsl array containing all the
 * dispersion relation polynomial coefficient computed at a pixel p.
 *
 * @param filename filename of the configuration file
 * @param beamID   beamID of the beam of interest (A=0, B=1, etc...)
 * @param p        a d_point object containing the position to compute the
 *                 dispersion coefficients at.
 *
 * @return v       a gsl_vector containing all the polynomial dispersion
 *                 relation coefficients computed at detector pixel p.
 */
gsl_vector *
get_disp_coeffs_at_pos (char *filename, const int for_grism, int beamID,
                        d_point p)
{
  int norder, i;
  gsl_vector *v;
  float c;

  norder = get_beam_disp_norder (filename, beamID);

  v = gsl_vector_alloc (norder + 1);
  for (i = 0; i < norder + 1; i++)
    {
      c = get_disp_coeff_at_pos (filename, for_grism, beamID, i, p);
      gsl_vector_set (v, i, c);
    }

  return v;
}


/**
 * Generate and returns a dispersion structure that is
 * computed at the given pixel p detector location.
 *
 * @param filename filename of the configuration file
 * @param beamID   beamID of the beam of interest (A=0, B=1, etc...)
 * @param p        a d_point object containing the position to compute
 *                 the dispersion coefficients at.
 *
 * @return res     a dispersion structure computed for beam beamID at
 *                 detector pixel p.
 */
dispstruct *
get_dispstruct_at_pos (char *filename, const int for_grism, int beamID,
                       d_point p)
{
  dispstruct *res;
  gsl_vector *pol2;

  res = malloc (sizeof (dispstruct));


  res->pol = get_disp_coeffs_at_pos (filename, for_grism, beamID, p);
  pol2 = res->pol;

  res->ID = beamID;
  res->cpoint.x = p.x;
  res->cpoint.y = p.y;
  sprintf (res->file, "%s", filename);
  res->for_grism = for_grism;
  return res;
}

/**
 * Free the memory allocated to a dispstruct structure
 *
 * @param p  pointer to a dispstruct structure
 */
void free_dispstruct(dispstruct *p)
{
  if (p!=NULL)  gsl_vector_free(p->pol);
  free(p);
  p=NULL;
}

/**
 * The program checks if the dispersion relation of a certain beam
 * is given in "grism" form (polynomial) or in "prism" form
 * (inverse polynomial)
 *
 * @param  filename filename of the configuration file
 * @param  beamID   beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return c        integer to declare grism (=1) or prism (=0)
 */
int
check_for_grism (char *filename, int beamID)
{
  int norder, i;
  int c=0;

  // get the order of the dispersion solution
  norder = get_beam_disp_norder (filename, beamID);

  // go over all orders and check the assumption
  // that it is a grism
  for (i = 0; i < norder + 1; i++)
    c += check_disp_coeff (filename, 1, beamID, i);

  // if the assumption is right,
  // set the return value for grism
  if (c == 0)
    {
      c = 1;
    }
  // if the assumption is wrong
  else
    {
      // reset the variable
      c = 0;

      // check the assumption that it is a grism
      for (i = 0; i < norder + 1; i++)
        c += check_disp_coeff (filename, 0, beamID, i);

      // if also this assumption is wrong,
      // raise an error that the configuration file
      // is inconsistent
      if (c != 0)
        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                     "check_for_grism: "
                     "The dispersion coefficients for beam %c in file",
                     "%s are inconsistent.", beamID, filename);
    }

  // return the result
  return c;
}

/**
 * The program checks an order of a dispersion solution for the
 * "suspicion" to be either grismatic or prismatic.
 *
 * @param  filename  filename of the configuration file
 * @param  for_grism the suspicion (=1 for grism, =0 for prism)
 * @param  beamID    beamID of the beam of interest (A=0, B=1, etc...)
 * @param  order     the order of the dispersion solution to be invetigated
 *
 * @return found     integer to confirm (=0) or falsily (=1) the suspicion
 */
int
check_disp_coeff (char *filename,const int for_grism,int beamID,int order)
{
  char beam[MAXCHAR];

  //int norder;
  int found = 1;
  int i;

  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;
  if (for_grism == 1)
    {
      sprintf (beam, "DLDP_%c_%d", BEAM (beamID), order);
    }
  else
    {
      sprintf (beam, "DLD1P_%c_%d", BEAM (beamID), order);
    }

  CfgRead (filename, DispConfig);
  if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
    {
      found = 0;
    }
  i = 0;
  while(DispConfig[i].name!=NULL){
    free(DispConfig[i++].data);
  }
  return found;
}

/**
 * Routine to read in the DLD1P_?_PRANGE keyword ? for a specific beam
 * from the configuration file. This is only for a prismatic dispersion
 * solution
 *
 * @param  filename filename of the configuration file
 * @param  beamID   beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return v        gsl-vector [pmin, pmax]
 */
gsl_vector *
get_prange (char *filename, int beamID)
{
  char beam[MAXCHAR];

  gsl_vector *v = NULL;

  struct CfgStrings DispConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };
  DispConfig[0].name = beam;

  sprintf (beam, "DLD1P_%c_PRANGE", BEAM (beamID));

  CfgRead (filename, DispConfig);

  if (DispConfig[0].data != NULL){
    v = string_to_gsl_array (DispConfig[0].data);
    if (v->size != 2)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_prange: two items, PMIN and PMAX must be given, not %i\n", v->size);
    if (gsl_vector_get(v,0) > gsl_vector_get(v,1))
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_prange: PMIN must be smaller than PMAX, however PMIN: %f, PMAX %f\n",
                   gsl_vector_get(v,0), gsl_vector_get(v,1));

  }
  return v;
}

global_disp *
get_global_disp(char *filename, const int for_grism, int beamID)
{
  global_disp *gdisp;

  gsl_vector *tmp;

  int order;


  gdisp =  (global_disp *)malloc (sizeof(global_disp));

  gdisp->n_order = get_beam_disp_norder(filename, beamID);

  //  fprintf (stdout, "Order of Polynomial: %i\n", gdisp->n_order);

  gdisp->all_coeffs = (gsl_vector **)malloc((gdisp->n_order+1)*sizeof(gsl_vector *));

  for (order=0; order <  gdisp->n_order+1; order++)
    {

      tmp = get_beam_disp_order (filename, for_grism, beamID, order);

      if (check_disp_order(tmp, beamID, order))
        gdisp->all_coeffs[order] = tmp;

    }


  gdisp->ID = beamID;
  sprintf (gdisp->file, "%s", filename);
  gdisp->for_grism = for_grism;

  return gdisp;
}

gsl_vector *
get_calvector_from_gdisp(const global_disp *gdisp, const d_point p)
{

  gsl_vector *calvector;

  double value;

  int order;


  calvector = gsl_vector_alloc(gdisp->n_order+1);

  for (order=0; order <  gdisp->n_order+1; order++)
    {
      value = get_2Ddep_value(gdisp->all_coeffs[order], p);
      gsl_vector_set(calvector, order, value);
      //      fprintf(stdout, "coefficient: %f\n", value);
    }

  return calvector;
}


double
get_2Ddep_value(const gsl_vector *coeffs, const d_point p)
{
  double value=0.0;
  double c=0.0;

  int i, j, k, n;

  i=0;
  n = 0.5 * (-1.0 + sqrt (1 + 8 * coeffs->size));
  for (j = 0; j < n; j++)
    {
      for (k = 0; k < (j + 1); k++)
        {
          c = gsl_vector_get (coeffs, i);
          value = value + c * pow (p.x, (j - k)) * pow (p.y, k);
          i++;
        }
    }
  return value;
}

void
print_global_disp(const global_disp *gdisp)
{
  int order;
  int i;

  fprintf(stdout, "Filename: %s\n", gdisp->file);
  fprintf(stdout, "BEAM: %c\n", BEAM(gdisp->ID));
  fprintf(stdout, "Grism: %i\n", gdisp->for_grism);
  fprintf(stdout, "Number of orders: %i\n", gdisp->n_order);

  for (order=0; order <  gdisp->n_order+1; order++)
    {
      fprintf(stdout, "Dimension in order %i: %zi", order, gdisp->all_coeffs[order]->size);
      for (i=0; i < (int)gdisp->all_coeffs[order]->size; i++)
        fprintf(stdout, " %e ", gsl_vector_get(gdisp->all_coeffs[order], i));
      fprintf(stdout, "\n");
    }
}

int
check_disp_order(const gsl_vector *tmp, const int beamID, const int order)
{

  int n;

  n = 0.5 * (-1.0 + sqrt (1 + 8 * tmp->size));
  if ((floor (n) - n) != 0)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "get_disp_coeff_at_pos: "
                   "Order %d of Beam %d does not contain a correct number of entries (i.e. "
                   "1,3,6,10,15...,n^2/2+n/2", order, beamID);
    }

  return 1;
}

void
free_global_disp(global_disp *gdisp)
{
  int order;

  for (order=0; order <  gdisp->n_order+1; order++)
    gsl_vector_free(gdisp->all_coeffs[order]);

  free(gdisp->all_coeffs);
  free(gdisp);
  gdisp = NULL;
}
