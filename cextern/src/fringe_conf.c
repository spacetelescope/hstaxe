/**
 */
#include <math.h>
#include "fitsio.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include "spc_cfg.h"
#include "aXe_utils.h"
#include "aXe_errors.h"
#include "fringe_conf.h"

#define AXE_CONFIG_PATH "AXE_CONFIG_PATH"

/**
 * Function: load_CCD_layer
 * The function creates and returns a structure for an optical
 * CCD layer. The input parameters are interpreted, and
 * the data which defines the CCD layer is loaded.
 *
 * Parameters:
 * @param refr_table - the name of the refractive index table
 * @param thickness  - the value for the thicknes keyword
                       (number or image name)
 *
 * Returns:
 * @return opt_layer  - the optical layer created
 */
ccd_layer *
load_CCD_layer(const char refr_table[], const char thickness[])
{
  ccd_layer *opt_layer;

  float tvalue;

  char **t_err=NULL;

  char refr_table_path[MAXCHAR];
  char thickness_path[MAXCHAR];

  // build the full pathname to the refraction table
  build_path (AXE_CONFIG_PATH, refr_table, refr_table_path);


  // allocate an error array
  t_err = (char **) malloc(sizeof(char *)*1);

  // allocate space for the return structure;
  // complain if this fails
  opt_layer = (ccd_layer *)malloc(sizeof(ccd_layer));

  // initialize both, the single value
  // as well as the matrix to NULL;
  opt_layer->thickness   = 0.0;
  opt_layer->thickness2D = NULL;

  // convert the keyvalue to a float
  // tvalue = atof(thickness);
  tvalue = strtod(thickness,t_err);

  // in case the float value is NULL = 0.0
  // the conversion failed, and it must be a string
  //  if (tvalue)
  if (strcmp(thickness, t_err[0]))
    {
#ifdef DEBUGFCONF
      fprintf(stderr, "Setting the layer thickness to value: %f mum\n",tvalue);
#endif
      // set the fixed thickness value
      opt_layer->thickness = tvalue;
    }
  else
    {
#ifdef DEBUGFCONF
      fprintf(stderr, "Loading layer thickness image: %s\n", thickness);
#endif

      // build the full pathname to the thickness image
      build_path (AXE_CONFIG_PATH, thickness, thickness_path);

      // load the 2D image for the thickness
      opt_layer->thickness2D = FITSimage_to_gsl(thickness_path, 1, 1);
    }

#ifdef DEBUGFCONF
      fprintf(stderr, "Loading refractive index table: %s\n", refr_table);
#endif
  // load the real part of the refraction index
  opt_layer->re_refraction =
    create_interp_ftable(refr_table_path, 2, "WAVELENGTH", "N",
			 REFRAC_INTERP_TYPE);

  // load the real part of the refraction index
  opt_layer->im_refraction =
    create_interp_ftable(refr_table_path, 2, "WAVELENGTH", "K",
			 REFRAC_INTERP_TYPE);

  // return the optical layer
  return opt_layer;
}


/**
 * Function: free_CCD_layer
 * The function releases all the memory allocated
 * in a structure for a CCD layer.
 *
 * Parameters:
 * @param opt_layer - structure for CCD layer
 *
 * Returns:
 * @return -
 */
void
free_CCD_layer(ccd_layer *opt_layer)
{
  // free the thickness matrix
  if (opt_layer->thickness2D)
    gsl_matrix_free(opt_layer->thickness2D);

  // free both interpolators
  free_interp(opt_layer->re_refraction);
  free_interp(opt_layer->im_refraction);

  // free the rest
  free(opt_layer);

  // set the structure to NULL
  opt_layer = NULL;
}


/**
 * Function: load_CCD_layers
 * The function extracts from a fringe configuration file
 * all information on individual CCD layers. It creates
 * the appropriate structure for every CCD layer,
 * and creates and returns a structure which completely
 * describes all layers in a CCD.
 *
 * Parameters:
 * @param fring_conf_path - the name of a fringe configuration file
 *
 * Returns:
 * @return opt_layers - the structure for the optical layers created
 */
ccd_layers *
load_CCD_layers(char fring_conf_path[])
{
  ccd_layers *opt_layers;

  char layer[MAXCHAR];
  char refr_table[MAXCHAR];

  //int nlayers;
  int i=0;

  struct CfgStrings LayerConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };

  LayerConfig[0].name = layer;

  // allocate space for the return structure;
  // complain if this fails
  opt_layers = malloc (sizeof (ccd_layers));

  // initialize the substrate
  opt_layers->substrate = NULL;

  // detemrine the number of layers
  // which are described in the
  //configuration file
  opt_layers->num_layers = 0;
  for (i = 0; i < MAX_BEAMS; i++)
    {
      sprintf (layer, "REFR_INDEX_%c", BEAM (i));
      CfgRead (fring_conf_path, LayerConfig);
      if (LayerConfig[0].data != NULL)
	{
	  opt_layers->num_layers += 1;
	  LayerConfig[0].data = NULL;
	}
    }

  //  if (opt_layers->num_layers < NLAYERS_MIN)
  //    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
  //		 "The number of layers in the configuration\n file is %i. "
  //		 "There must be at least: %i!\n",
  //		 opt_layers->num_layers, NLAYERS_MIN);


  // allocate space for the array to
  // the layers described in the
  // configuration file
  opt_layers->opt_layer =
    (ccd_layer **)malloc (opt_layers->num_layers*sizeof(ccd_layers));

  // successively find the
  // description and finally
  // load all layers
  for (i = 0; i < MAX_BEAMS; i++)
    {
      // chekc for the first relevant keyword
      sprintf (layer, "REFR_INDEX_%c", BEAM (i));
      CfgRead (fring_conf_path, LayerConfig);
      if (LayerConfig[0].data != NULL)
	{
	  sprintf (refr_table, "%s", LayerConfig[0].data);
	  LayerConfig[0].data = NULL;

	  // check for the second relevant keyword
	  sprintf (layer, "THICKNESS_%c", BEAM (i));
	  CfgRead (fring_conf_path, LayerConfig);
	  if (LayerConfig[0].data != NULL)
	    {
	      // load the layer
	      opt_layers->opt_layer[i] =
		load_CCD_layer(refr_table, LayerConfig[0].data);

	      LayerConfig[0].data = NULL;
	    }
	  else
	    {
	      // complain that a keyword is missing
	      // for a layer
	      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			   "Could load index table %s,\n"
			   "but not thickness information "
			   "for beam %c\n", refr_table, BEAM (i));
	    }
	}
    }

  // return the whole structure
  return opt_layers;
}


/**
 * Function: free_CCD_layers
 * The function releases all memory allocated in a
 * CCD layers structure.
 *
 * Parameters:
 * @param opt_layers - the CCD layers structure
 *
 * Returns:
 * @return -
 */
void
free_CCD_layers(ccd_layers *opt_layers)
{
  int index=0;

  // free each layer individually
  for (index=0; index < opt_layers->num_layers; index++)
      free_CCD_layer(opt_layers->opt_layer[index]);

  // free the interpolator for the substrate
  //  free_linint(opt_layers->substrate);
  free_interp(opt_layers->substrate);

  // free everything
  free(opt_layers);

  // set the structure to NULL
  opt_layers = NULL;
}


/**
 * Function: load_fringe_conf
 * The function loads the fringe configuration file given as parameter.
 * All relevant data marked with keywords is extracted. A fringe
 * configuration structure is built up and returned from the
 * data given in the keyvalues.
 *
 * Parameters:
 * @param fring_conf_path - the complete path-name to a fringe configuration file
 *
 * Returns:
 * @return fconf - the fringe configuration structure created
 */
fringe_conf *
load_fringe_conf(char fring_conf_path[])
{
  char beam[MAXCHAR];

  char ffile_name[MAXCHAR];
  char ffile_name_path[MAXCHAR];

  int index=0;

  fringe_conf *fconf=NULL;

  gsl_vector *v=NULL;

  struct CfgStrings FringeConfig[] =
    {
      {"FRINGE_AMPLITUDE", NULL},
      {"FRINGE_PHASE", NULL},
      {"FRINGE_RANGE", NULL},
      {"MAX_DISPERSION", NULL},
      {"FILTER_NAME", NULL},
      {"FRINGE_STEP",NULL},
      {"NUM_STEPS",NULL},
      {"SUBSTR_TRANS",NULL},

      {NULL, NULL},
      {NULL, NULL}		/* array terminator. REQUIRED !!! */
    };


  FringeConfig[8].name = beam;

  // allocate space for the return structure;
  // complain if this fails
  fconf = malloc (sizeof (fringe_conf));
  if (fconf == NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Could not allocate memory for fringe configuration");
    }

  // load the optical layers directly from the file
  fconf->opt_layers = load_CCD_layers(fring_conf_path);

  // read in the file
  CfgRead (fring_conf_path, FringeConfig);

  // make some intitializations
  fconf->fringe_amp     = 0.0;
  fconf->fringe_phase   = 1.0e+32;
  fconf->fringe_step    = 0.0;
  fconf->num_steps      = 0;
  fconf->fringe_range   = NULL;
  fconf->max_dispersion = 0.0;
  fconf->filter_through = NULL;


  for (index = 0; index < 9; index++)
    {

      // read in the fringe amplitude
      if (!strcmp (FringeConfig[index].name, "FRINGE_AMPLITUDE"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      fconf->fringe_amp = atof(FringeConfig[index].data);
	    }
	}

      // read in the fringe phase
      if (!strcmp (FringeConfig[index].name, "FRINGE_PHASE"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      fconf->fringe_phase = atof(FringeConfig[index].data);
	    }
	}


      // read in the fringe step
      if (!strcmp (FringeConfig[index].name, "FRINGE_STEP"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      fconf->fringe_step = atof(FringeConfig[index].data);
	    }
	}

      // read in the number of steps
      if (!strcmp (FringeConfig[index].name, "NUM_STEPS"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      fconf->num_steps = atoi(FringeConfig[index].data);
	    }
	}

      // read in the fringe range
      if (!strcmp (FringeConfig[index].name, "FRINGE_RANGE"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      v = string_to_gsl_array (FringeConfig[index].data);
	      fconf->fringe_range = v;
	    }
	}

      // read in the fringe amplitude
      if (!strcmp (FringeConfig[index].name, "MAX_DISPERSION"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      fconf->max_dispersion = atof(FringeConfig[index].data);
	    }
	}

      // read in the filter name
      if (!strcmp (FringeConfig[index].name, "FILTER_NAME"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      sprintf (ffile_name, "%s", FringeConfig[index].data);

	      // build the full pathname to the refraction table
	      build_path (AXE_CONFIG_PATH, ffile_name, ffile_name_path);

	      fconf->filter_through =
		create_interp_ftable(ffile_name_path, 2, "WAVELENGTH",
				     "THROUGHPUT", FILTER_INTERP_TYPE);
	    }
	}

      // read in the substrate transmission
      if (!strcmp (FringeConfig[index].name, "SUBSTR_TRANS"))
	{
	  if (FringeConfig[index].data != NULL)
	    {
	      sprintf (ffile_name, "%s", FringeConfig[index].data);

	      // build the full pathname to the refraction table
	      build_path (AXE_CONFIG_PATH, ffile_name, ffile_name_path);

	      fconf->opt_layers->substrate =
		create_interp_ftable(ffile_name_path, 2, "WAVELENGTH",
				     "TRANSMISSION", REFRAC_INTERP_TYPE);
	    }
	}
    }

  return fconf;
}

/**
 * Function: check_fringe_conf
 *
 * Parameters:
 * @param fconf  - the fringe configuration structure
 *
 * Returns:
 * @return
 */
void
check_fringe_conf(fringe_conf *fconf)
{
  // make  sure the fringe amplitude is set
  if (!fconf->fringe_amp)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "fringe_conf: The fringing amplitude must\n"
		 "be set to be able determining pixel fringing!\n");

  // make  sure the fringe phase is set
  if (fconf->fringe_phase > 1.0e+31)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "fringe_conf: The fringing phase must\n"
		 "be set to be able determining pixel fringing!\n");

  // make  sure the transmission of the substrate
  // layer is given
  if (!fconf->opt_layers->substrate)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "fringe_conf: The substrate transmission\n"
		 "must be given to be able determining pixel fringing!\n");

  // make sure that there is a minimum
  // number of layers.
  if (fconf->opt_layers->num_layers < NLAYERS_MIN)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "The number of layers in the configuration\n file is %i. "
		 "There must be at least: %i!\n",
		 fconf->opt_layers->num_layers, NLAYERS_MIN);

  // if not defined, give a meaningfull
  // default for the fringing range
   if (!fconf->fringe_range)
    {
      fconf->fringe_range = gsl_vector_alloc(2);
      gsl_vector_set(fconf->fringe_range, 0, DEFAULT_MIN_LAMBDA);
      gsl_vector_set(fconf->fringe_range, 1, DEFAULT_MAX_LAMBDA);
      fprintf(stdout, "Setting the fringing range to defaults: [ %f, %f].\n",
	      DEFAULT_MIN_LAMBDA, DEFAULT_MAX_LAMBDA);
    }

   // give reasonable default for
   // the fringe step
   if (!fconf->fringe_step)
     {
       fconf->fringe_step = DEFAULT_FRINGE_STEP;
       fprintf(stdout, "Setting the fringing step to default: %f.\n",
	       DEFAULT_FRINGE_STEP);
     }

   // give a reasonable default for the
   // number of steps
   if (!fconf->num_steps)
     {
       fconf->num_steps = DEFAULT_NUM_STEPS;
       fprintf(stdout, "Setting the number of steps to default: %i.\n",
	       DEFAULT_NUM_STEPS);
     }

   // give a reasonable default for the
   // maximum dispersion to correct for
   if (!fconf->max_dispersion)
     {
       fconf->max_dispersion = DEFAULT_MAX_DISPERSION;
       fprintf(stdout, "Setting the maximal dispersion to default: %f.\n",
	       DEFAULT_MAX_DISPERSION);
     }
}

/**
 * Function: free_fringe_conf
 * The function deallocates all memory in a
 * fringe configuration structure.
 *
 * Parameters:
 * @param fconf  - the fringe configuration structure
 *
 * Returns:
 * @return -
 */
void
free_fringe_conf(fringe_conf *fconf)
{


  // free the filter throughput
  if (fconf->filter_through)
    free_interp(fconf->filter_through);

  // free the optical layers
  free_CCD_layers(fconf->opt_layers);

  // free the fringe range
  gsl_vector_free(fconf->fringe_range);
  fconf->fringe_range = NULL;

  // free everything
  free(fconf);
}


/**
 * Function: create_interp_ftable
 * This function cretaes an interpolator from data values stored
 * in a fits table. An interpolator is a structure to hold all
 * data to compute the interpolated data values for different
 * interpolation methods. Core of this structure is the interpolator
 * types offered in the gsl-library.
 * This function is a tool to extract the relevant data points from
 * a fits table and to set up and return the interpolator.
 *
 *
 * Parameters:
 * @param table_name  - the name of the fits table
 * @param hdunum      - the extension number with the data
 * @param xcol        - the column name with the indepenent values
 * @param ycol        - the column name with the dependent values
 * @param interp_type - the interpolation type
 *
 * Returns:
 * @return interp - the interpolation structure created
 */
interpolator *
create_interp_ftable(const char table_name[], const int hdunum,
		     char xcol[], char ycol[],
		     const gsl_interp_type *interp_type)
{

  interpolator *interp;

  double *x;
  double *y;
  int    f_status = 0;
  int hdutype, anynul;
  long nrows=0;
  int colnum=0;

  fitsfile *input;

  // allocate space for the return structure;
  // complain if this fails
  //    interp = (interpolator *)malloc (sizeof (interpolator));
  //   if (interp == NULL)
  //     {
  //        aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
  //  		   "Could not allocate memory for interpolator");
  //      }

#ifdef DEBUGFCONF
  fprintf(stderr, "Loading columns %s and %s of fitstable: %s\n", xcol, ycol, table_name);
#endif

  //  Open the file for reading
  fits_open_file (&input, table_name, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_response_function_fromFITS: "
		   "Could not open" " file: %s",
		   table_name);
    }

  /* Move to the correct hdu */
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable:"
		   "Could not read extention %d from file: %s",
		   hdunum, table_name);
    }

  /* Get number of rows */
  fits_get_num_rows (input, &nrows, &f_status);
  if (f_status) {
    ffrprt (stderr, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "create_interp_ftable: "
		 "Could not determine the number of rows in"
		 " table %s",table_name);
  }

  /* Allocate temporary memory space */
  x = (double *) malloc(nrows*sizeof(double));
  if (!x) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }
  y = (double *) malloc(nrows*sizeof(double));
  if (!y) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  /**************************/
  /* Read the X-column      */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, xcol, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not determine column %s in "
		   " table %s", xcol, table_name);
    }

  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, nrows, NULL, x,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not read content of WAVELENGTH column "
		   " from BINARY table %s",table_name);
    }

  /**************************/
  /* Read the y-column      */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, ycol, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not determine column %s in "
		   " table %s", ycol, table_name);
    }

  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, nrows, NULL, y,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not read column %s"
		   " from BINARY table %s", ycol, table_name);
    }

  fits_close_file(input,&f_status);
  if (f_status) {
      aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		   "Could not close %s", table_name);
  }

  // create the interpolator
  interp = create_interp(nrows, interp_type, x, y);

  // return the interpolator
  return interp;
}


/**
 * Function: create_interp
 * The function creates, and intializes and returns an interpolator
 * from the basic input, which is the data, the number of data
 * items and the interpolator type requested.
 *
 * Parameters:
 * @param nvals       - the number of data values in the arrays
 * @param interp_type - the interpolation type
 * @param xvals       - the array with the independent data values
 * @param yvals       - the array with the dependent data values
 *
 * Returns:
 * @return interp - the interpolation structure created
 */
interpolator *
create_interp(const int nvals, const gsl_interp_type *interp_type,
	      double *xvals, double *yvals)
{
  interpolator *interp;

  // allocate space for the return structure;
  // complain if this fails
  interp = (interpolator *)malloc (sizeof (interpolator));
  if (interp == NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Could not allocate memory for interpolator");
    }

  // store the min and max values
  interp->xmin    = xvals[0];
  interp->xmax    = xvals[nvals-1];

  // check for ascending order in the
  // independent values
  if (interp->xmax < interp->xmin)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "The independent data values to be stored\n"
		 " in an interplator must be in INCREASING order!\n");

  // store the number of data values
  interp->nvals   = nvals;

  // store the data arrays
  interp->xvals   = xvals;
  interp->yvals   = yvals;

  // create and intitialize the gsl inteprolator
  interp->acc     = gsl_interp_accel_alloc();
  interp->interp  = gsl_interp_alloc(interp_type, interp->nvals);
  gsl_interp_init(interp->interp, interp->xvals, interp->yvals, interp->nvals);

  // return the new structure
  return interp;
}

/**
 * Function: print_interp
 *
 * Parameters:
 * @param interp - the interpolator structure
 *
 * Returns:
 * @return -
 */
void
print_interp(interpolator *interp)
{
  int index;
  double x_typical;

  x_typical = interp->xmin+(interp->xmax-interp->xmin)/2.0;

  fprintf(stdout, "xmin: %e, xmax: %e\n",interp->xmin, interp->xmax);
  fprintf(stdout, "number of data values: %i\n",interp->nvals);
  fprintf(stdout, "Interpolation type: %s\n", gsl_interp_name(interp->interp));
  fprintf(stdout, "Characteristic value pair: (x,y) = (%e, %e)\n",
	  x_typical, gsl_interp_eval(interp->interp, interp->xvals,
				     interp->yvals, x_typical, interp->acc));

  fprintf(stdout, "Alternative value pair: (x,y) = (%e, %e)\n",
	  x_typical, eval_interp(interp, x_typical));
  for (index=0; index < interp->nvals; index++)
    fprintf(stdout, "xvalue: %e, yvalue: %e\n",interp->xvals[index], interp->yvals[index]);
  fprintf(stdout, "\n");
}


/**
 * Function: eval_interp
 * The function computes and returns the interpolated value
 * at a given position for an interpolater.
 *
 * Parameters:
 * @param interp - the interpolator
 * @param xval   - the position to evaluate the interpolator
 *
 * Returns:
 * @return (value) - the interpolated data value
 */
double
eval_interp(interpolator *interp, const double xval)
{
  // check whether the x-value is within
  // the range spanned by the data;
  // complain if the x-value is outside
  if (xval < interp->xmin || xval > interp->xmax)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "independent interpolation value %f "
		   "is outside interval (%f, %f)\n", xval,
		   interp->xmin, interp->xmax);

  // evaluate and return the interpolated value
  // in on the spot
  return  gsl_interp_eval(interp->interp, interp->xvals, interp->yvals,
			  xval, interp->acc);
}


/**
 * Function: free_interp
 * the function frees all memory allocated in
 * an interpolator structure.
 *
 * Parameters:
 * @param interp - the interpolator structure
 *
 * Returns:
 * @return -
 */
void
free_interp(interpolator *interp)
{
  // free the data vectors
  free(interp->xvals);
  free(interp->yvals);

  // free the two gsl structures
  gsl_interp_accel_free(interp->acc);
  gsl_interp_free (interp->interp);

  // free the rest
  free(interp);

  // set it to NULL
  interp = NULL;
}


/**
 * Function: create_linint_ftable
 * This function creates a linear interpolator from data values stored
 * in a fits table. An linear interpolator is a structure to hold all
 * data to compute the linear interpolated data values at any
 * point bracketed by the data.
 * This function is a tool to extract the relevant data points from
 * a fits table and to set up and return the linear interpolator.
 *
 * Parameters:
 * @param table_name  - the name of the fits table
 * @param hdunum      - the extension number with the data
 * @param xcol        - the column name with the indepenent values
 * @param ycol        - the column name with the dependent values
 * @param interp_type - the interpolation type
 *
 * Returns:
 * @return lin_int - the linear interplator structure created
 */
linear_interp *
create_linint_ftable(const char table_name[], const int hdunum,
		     char xcol[], char ycol[])
{

  linear_interp *lin_int;

  double *x;
  double *y;
  int    f_status = 0;
  int hdutype, anynul;
  long nrows=0;
  int colnum;
  int index=0;

  gsl_vector *xvals;
  gsl_vector *yvals;

  fitsfile *input;

  // allocate space for the return structure;
  // complain if this fails
  lin_int = (linear_interp *)malloc (sizeof (linear_interp));
  if (lin_int == NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Could not allocate memory for linear interpolator");
    }

#ifdef DEBUGFCONF
  fprintf(stderr, "Loading columns %s and %s of fitstable: %s\n", xcol, ycol, table_name);
#endif

  //  Open the file for reading
  fits_open_file (&input, table_name, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "get_response_function_fromFITS: "
		   "Could not open" " file: %s",
		   table_name);
    }

  /* Move to the correct hdu */
  fits_movabs_hdu (input, hdunum, &hdutype, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable:"
		   "Could not read extention %d from file: %s",
		   hdunum, table_name);
    }

  /* Get number of rows */
  fits_get_num_rows (input, &nrows, &f_status);
  if (f_status) {
    ffrprt (stderr, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "create_interp_ftable: "
		 "Could not determine the number of rows in"
		 " table %s",table_name);
  }

  /* Allocate temporary memory space */
  x = (double *) malloc(nrows*sizeof(double));
  if (!x) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }
  y = (double *) malloc(nrows*sizeof(double));
  if (!y) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  /**************************/
  /* Read the X-column      */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, xcol, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not determine column %s in "
		   " table %s", xcol, table_name);
    }

  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, nrows, NULL, x,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not read content of WAVELENGTH column "
		   " from BINARY table %s",table_name);
    }

  /**************************/
  /* Read the y-column      */
  /**************************/
  /* Get column number */
  fits_get_colnum (input, CASEINSEN, ycol, &colnum, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not determine column %s in "
		   " table %s", ycol, table_name);
    }

  /* Read the data */
  fits_read_col (input, TDOUBLE, colnum, 1, 1, nrows, NULL, y,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not read column %s"
		   " from BINARY table %s", ycol, table_name);
    }

  fits_close_file(input,&f_status);
  if (f_status) {
      aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		   "Could not close %s", table_name);
  }

#ifdef DEBUGFCONF
  fprintf(stderr, "Number of rows: %i\n", nrows);
#endif

  // allocate space for the two vectors
  xvals = gsl_vector_alloc(nrows);
  yvals = gsl_vector_alloc(nrows);

  // transport the values from the
  // arrays to the gsl-arrays
  for (index=0; index < nrows; index++)
    {
      gsl_vector_set(xvals, index, x[index]);
      gsl_vector_set(yvals, index, y[index]);
    }

  // reverse both vectors if necessary
  if (gsl_vector_get(xvals, 0) > gsl_vector_get(xvals, nrows-1))
    {
      gsl_vector_reverse(xvals);
      gsl_vector_reverse(yvals);
    }

  // set the two helper elements
  lin_int->act_index = 0;
  lin_int->num_elem  = nrows;

  // set the minimum and maximum
  lin_int->xmin = gsl_vector_get(xvals, 0);
  lin_int->xmax = gsl_vector_get(xvals, nrows-1);

  // set the two vectors
  lin_int->xvals = xvals;
  lin_int->yvals = yvals;

  // release memory in the temp variables
  free(x);
  free(y);

  // return the new structure
  return lin_int;
}


/**
 * Function: free_linint
 * The function frees all memory allocated
 * in a linear interpolator structure.
 *
 * Parameters:
 * @param lin_int - the linear inteprolator
 *
 * Returns:
 * @return -
 */
void
free_linint(linear_interp *lin_int)
{
  // free the two gsl-vectors
  gsl_vector_free(lin_int->xvals);
  gsl_vector_free(lin_int->yvals);

  // free the whole structure
  free(lin_int);

  // set the structure to NULL
  lin_int = NULL;
}
