/**
 * File: trace_conf.c
 */

#include "trace_conf.h"


/**
 * Parses a configuration file and searches for DYDX_ORDER_%c
 * keywords. Returns the order of the field dependent
 * trace relation for that order.
 * @param filename filename of the configuration file
 * @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return disporder the order of the dispersion relation for the beam of interest
 */
int get_beam_trace_norder(char *filename, int beamID)
  {
    char beam[MAXCHAR];
    int disporder = -1;
    int found = 0;
    int i=0;

    struct CfgStrings DispConfig[] =
      {
        { NULL, NULL },
        { NULL, NULL } };
    DispConfig[0].name = beam;
    sprintf(beam, "DYDX_ORDER_%c", BEAM (beamID));
    CfgRead(filename, DispConfig);
    if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
      {
        disporder = atoi(DispConfig[0].data);
        found = 1;
      }
    if (!found)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_beam_trace_norder: %s was not found in %s.\n",
        beam, filename);
      }
    if (disporder<0)
      aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
      "Trace polynomial must at least have one coefficients (i.e. zeroth order).\n");

    // release memory
    i=0;
    while (DispConfig[i].name!=NULL)
      free(DispConfig[i++].data);

    return disporder;
  }

/**
 *  Parses a configuration file and returns a gsl vector
 *  containing the 2D polynomial coefficient allowing one to compute
 *  a given trace coefficient anywhere on the detector.
 *  @param filename filename of the configuration file
 *  @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *  @param order order of the coefficient for which the 2D field dependent
 *  coefficients are queried.
 *
 * @return v A gsl_vector() containing all of the 2D field dependent polynomial coeffiecient
 */
gsl_vector * get_beam_trace_order(char *filename, int beamID, int order)
  {
    char beam[MAXCHAR];

    int norder;
    int found = 0;
    int i=0;

    gsl_vector *v= NULL;
    struct CfgStrings DispConfig[] =
      {
        { NULL, NULL },
        { NULL, NULL } };
    DispConfig[0].name = beam;
    norder = get_beam_disp_norder(filename, beamID);
    sprintf(beam, "DYDX_%c_%d", BEAM (beamID), order);
    CfgRead(filename, DispConfig);
    if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
      {
        v = string_to_gsl_array(DispConfig[0].data);
        found = 1;
      }

    if (!found)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_beam_trace_order: %s was not found in %s.\n",
        beam, filename);
      }

    // release memory
    i=0;
    while (DispConfig[i].name!=NULL)
      free(DispConfig[i++].data);

    return v;
  }

/**
 *   Parses a configuration file, computes and return the value
 *   of the order^th trace polynomial coefficient at a given
 *   detector pixel.
 *
 *   @param filename filename of the configuration file
 *   @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *   @param order order of the coefficient for which the 2D field dependent
 *   @param p A d_point object containing the position to compute the trace coefficient at.
 *
 *   @return res The value of the order6th coefficient of beamID computed at detector pixel p.
 */
float get_trace_coeff_at_pos(char *filename, int beamID, int order, d_point p)
  {
    gsl_vector *v;
    float n, c, res;
    int i, j, k;

    v = get_beam_trace_order(filename, beamID, order);
    n = 0.5 * (-1.0 + sqrt(1 + 8 * v->size));
    if ((floor(n) - n) != 0)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_trace_coeff_at_pos: "
        "Order %d of Beam %d in %s does not contain a correct number of entries (i.e. "
        "1,3,6,10,15...,n^2/2+n/2", order, beamID, filename);
      }

    i = 0;
    res = 0;
    for (j = 0; j < n; j++)
      {
        for (k = 0; k < (j + 1); k++)
          {
            c = gsl_vector_get(v, i);
            res = res + c * pow(p.x, (j - k)) * pow(p.y, k);
            i++;
          }
      }

    return res;
  }

/**
 * Parses a configuration file and returns a gsl vector
 * containing the 2D polynomial coefficient allowing one to compute
 * direct/grism images x-offset trace corection (caused by optics) anywhere on the detector.
 * @param filename filename of the configuration file
 * @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *
 *  @return v A gsl_vector() containing all of the 2D field dependent polynomial coeffiecient
 */
gsl_vector * get_beam_trace_xoff(char *filename, int beamID)
  {
    char beam[MAXCHAR];

    int norder;
    int found = 0;
    int i=0;

    gsl_vector *v= NULL;
    struct CfgStrings DispConfig[] =
      {
        { NULL, NULL },
        { NULL, NULL } };
    DispConfig[0].name = beam;
    norder = get_beam_disp_norder(filename, beamID);
    sprintf(beam, "XOFF_%c", BEAM (beamID));
    CfgRead(filename, DispConfig);
    if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
      {
        v = string_to_gsl_array(DispConfig[0].data);
        found = 1;
      }

    if (!found)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_beam_trace_xoff: %s was not found in %s.\n",
        beam, filename);
      }

    // release memory
    i=0;
    while (DispConfig[i].name!=NULL)
      free(DispConfig[i++].data);

    return v;
  }

/**
 *  Parses a configuration file and returns a gsl vector
 *  containing the 2D polynomial coefficient allowing one to compute
 *  direct/grism images y-offset corection (caused by optics) anywhere on the detector.
 *  @param filename filename of the configuration file
 *  @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *
 *  @return v A gsl_vector() containing all of the 2D field dependent polynomial coeffiecient
 */
gsl_vector * get_beam_trace_yoff(char *filename, int beamID)
  {
    char beam[MAXCHAR];

    int norder;
    int found = 0;
    int i=0;
    gsl_vector *v= NULL;
    struct CfgStrings DispConfig[] =
      {
        { NULL, NULL },
        { NULL, NULL } };
    DispConfig[0].name = beam;
    norder = get_beam_trace_norder(filename, beamID);
    sprintf(beam, "YOFF_%c", BEAM (beamID));
    CfgRead(filename, DispConfig);
    if ((DispConfig[0].name == beam) && (DispConfig[0].data != NULL))
      {
        v = string_to_gsl_array(DispConfig[0].data);
        found = 1;
      }

    if (!found)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_beam_trace_yoff: %s was not found in %s.\n",
        beam, filename);
      }

    // release memory
    i=0;
    while (DispConfig[i].name!=NULL)
      free(DispConfig[i++].data);

    return v;
  }

/**
 *  The function evaluates the a polynomial given in the input
 *  at a specific point, which also is given in the input.
 *  The polynomial follows the description of 2D vatriable
 *  quantities as given in the aXe manual.
 *
 *  @param coeffs - the coefficients of the polynomial
 *  @param p      - the point to evaluate the polynomial
 *  @param beamID - beamID of the beam of interest (A=0, B=1, etc...)
 *
 *  @return res   - of the polynomial at point p
 *
 */
float eval_trace_off_at_pos(gsl_vector *coeffs, d_point p, int beamID)
  {
    float n, c, res;
    int i, j, k;

    n = 0.5 * (-1.0 + sqrt(1 + 8 * coeffs->size));
    if ((floor(n) - n) != 0)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_trace_xoff_at_pos: "
        "Beam %d in configuration file does not contain a correct number of entries (i.e. "
        "1,3,6,10,15...,n^2/2+n/2", beamID);
      }

    i = 0;
    res = 0;
    for (j = 0; j < n; j++)
      {
        for (k = 0; k < (j + 1); k++)
          {
            c = gsl_vector_get(coeffs, i);
            res = res + c * pow(p.x, (j - k)) * pow(p.y, k);
            i++;
          }
      }

    return res;
  }

/**
 * @param filename filename of the configuration file
 * @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *
 * @return res The value of the x-offset of beamID computed at detector pixel p.
 */
float get_trace_xoff_at_pos(char *filename, int beamID, d_point p)
  {
    gsl_vector *v;
    float n, c, res;
    int i, j, k;

    v = get_beam_trace_xoff(filename, beamID);
    n = 0.5 * (-1.0 + sqrt(1 + 8 * v->size));
    if ((floor(n) - n) != 0)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_trace_xoff_at_pos: "
        "Beam %d in %s does not contain a correct number of entries (i.e. "
        "1,3,6,10,15...,n^2/2+n/2", beamID, filename);
      }

    i = 0;
    res = 0;
    for (j = 0; j < n; j++)
      {
        for (k = 0; k < (j + 1); k++)
          {
            c = gsl_vector_get(v, i);
            res = res + c * pow(p.x, (j - k)) * pow(p.y, k);
            i++;
          }
      }

    return res;
  }



/**
 * @param filename filename of the configuration file
 * @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 *
 *   @return res The value of the y-offset of beamID computed at detector pixel p.
 */
float get_trace_yoff_at_pos(char *filename, int beamID, d_point p)
  {
    gsl_vector *v;
    float n, c, res;
    int i, j, k;

    v = get_beam_trace_yoff(filename, beamID);
    n = 0.5 * (-1.0 + sqrt(1 + 8 * v->size));
    if ((floor(n) - n) != 0)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_trace_yoff_at_pos: "
        "Beam %d in %s does not contain a correct number of entries (i.e. "
        "1,3,6,10,15...,n^2/2+n/2", beamID, filename);
      }

    i = 0;
    res = 0;
    for (j = 0; j < n; j++)
      {
        for (k = 0; k < (j + 1); k++)
          {
            c = gsl_vector_get(v, i);
            res = res + c * pow(p.x, (j - k)) * pow(p.y, k);
            i++;
          }
      }

    return res;
  }



/**
 * Parse a configuration file and return a gsl array containing all the
 * trace polynomial coefficient computed at a pixel p.
 *
 * @param filename filename of the configuration file
 * @param beamID beamID of the beam of interest (A=0, B=1, etc...)
 * @param p A d_point object containing the position to compute the dispersion coefficients at.
 *
 *  @return v A gsl_vector containing all the polynomial trace coefficients
 * computed at detector pixel p.
 */
gsl_vector * get_trace_coeffs_at_pos(char *filename, int beamID, d_point p)
  {
    int norder, i;
    gsl_vector *v;
    float c;

    norder = get_beam_trace_norder(filename, beamID);

    v = gsl_vector_alloc(norder + 1);
    for (i = 0; i < norder + 1; i++)
      {
        c = get_trace_coeff_at_pos(filename, beamID, i, p);
        gsl_vector_set(v, i, c);
      }

    return v;
  }


/**
   Generate and returns a trace structure that is
   computed at the given pixel p detector location.
   @param filename filename of the configuration file
   @param beamID beamID of the beam of interest (A=0, B=1, etc...)
   @param p A d_point object containing the position to compute the dispersion coefficients at.

   @return res A dispersion structure computed for beam beamID at detector pixel p.
*/
tracestruct * get_tracestruct_at_pos(char *filename, int beamID, d_point p)
  {
    tracestruct *res;

    res = malloc(sizeof(tracestruct));
    if (res == NULL)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "get_tracestruct_at_pos:"
        "Could not allocate tracestruct");
      }
    res->pol = get_trace_coeffs_at_pos(filename, beamID, p);
    res->offset.x = get_trace_xoff_at_pos(filename, beamID, p);
    res->offset.y = get_trace_yoff_at_pos(filename, beamID, p);
    res->ID = beamID;
    res->cpoint.x = p.x;
    res->cpoint.y = p.y;
    sprintf(res->file, "%s", filename);
    return res;
  }

void free_tracestruct(tracestruct * trace)
  {
    gsl_vector_free(trace->pol);
    if (trace->pol == NULL)
      {
        fprintf(stderr, "properly set to NULL!");
        exit(0);
      }
    trace = NULL;
  }

/**
 *Displays the content of a tracestruct structure to *FILE
 */
void tracestruct_fprintf(FILE * file, tracestruct * trace)
  {
    fprintf(file, "Trace ID: %d\n", trace->ID);
    fprintf(file, "Trace detector location: %f %f\n", trace->cpoint.x,
        trace->cpoint.y);
    gsl_vector_fprintf(file, trace->pol, "Trace pol: %f");
    fprintf(file, "Trace xoff: %f\n", trace->offset.x);
    fprintf(file, "Trace yoff: %f\n", trace->offset.y);
    fprintf(file, "Trace file: %s\n", trace->file);
  }
