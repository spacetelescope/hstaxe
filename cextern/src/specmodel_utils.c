/**
 * File: specmodel_utils.c
 *
 */
#include "specmodel_utils.h"
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spc_wl_calib.h"
#include "aXe_errors.h"
#include "model_utils.h"


/**
 * Function: load_object_models
 * The function loads all direct emission models stored in a mulit-extension
 * fits image into a corresponding object model structure.
 *
 * Parameters:
 * @param object_models_file - the file with the object model data
 *
 * Returns:
 * @return objmodels - the object model structure
 */
object_models *
load_object_models(const char object_models_file[])
 {
  object_models *objmodels;

  int n_extent=0;

  // give feedback to the user
  fprintf (stdout, "aXe_PETCONT: Loading object models: %s", object_models_file);

  // allocate space for the direct emission structure
  objmodels = (object_models *)malloc(sizeof(object_models));

  // get the number of models
  n_extent = get_num_extensions(object_models_file);

  // load the list with the individual direct emission
  objmodels->obj_list = load_diremission_list(object_models_file, n_extent);

  // store the number of direct emission images
  objmodels->n_models =  n_extent;

  // report the number of models loaded
  fprintf (stdout,"%d models loaded.\n", objmodels->n_models);

  // return the list structure
  return objmodels;
}


/**
 * Function: load_diremission_list
 * The function loads all the direct emission models given in a
 * multi-extension fits file and creates a list of direct emission
 * models. This list is returned.
 *
 * Parameters:
 * @param object_models_file - the file with the object model data
 * @param n_models           - the number of object models in the file
 *
 * Returns:
 * @return obj_list - the list of direct emission models
 */
dirim_emission **
load_diremission_list(const char object_models_file[], const int n_models)
{
  int i=0;

  dirim_emission **obj_list;

  // allocate space for the SED list
  obj_list = (dirim_emission **)malloc(n_models * sizeof(dirim_emission *));

  // go over the nuymber of models
  for (i=0; i < n_models; i++)
    {
      // load one model
      // the 'i + 2' is necessary due to the
      // different counting in "fitsio.h"
      obj_list[i] =  load_diremission(object_models_file, i+2);
    }

  // return the loaded structure
  return obj_list;
}


/**
 * Function: load_diremission
 * The function loads one extension of a multi-extension fits file
 * with emission model and creates a single direct emission object.
 * This object is eventually returned.
 *
 * Parameters:
 * @param object_models_file - the file with the object model data
 * @param ext_num            - the extension number to load
 *
 * Returns:
 * @return diremission - a direct emission model
 */
dirim_emission *
load_diremission(const char object_models_file[], const int ext_num)
{
  dirim_emission *diremission;

  // allocate memory for the structure
  diremission = (dirim_emission *) malloc(sizeof(dirim_emission));

  // load the image in the gsl
  diremission->modimage = FITSimage_to_gsl(object_models_file, ext_num, 1);

  // transfer the image dimension
  diremission->dim_x = (int)diremission->modimage->size1;
  diremission->dim_y = (int)diremission->modimage->size2;

  // compute and store the mean image coordinates
  diremission->xmean = (float)(diremission->dim_x-1) / 2.0;
  diremission->ymean = (float)(diremission->dim_y-1) / 2.0;

  // return the structure
  return diremission;
}

/**
 * Function: free_object_models
 * The function frees all memory allocated in a object model structure.
 *
 * Parameters:
 * @param objmodels - the object model structure to free
 *
 * Returns:
 * @return -
 */
void
free_object_models(object_models *objmodels)
{
  int iii;

  // go over each model
  for (iii=0; iii < objmodels->n_models; iii++)
    {
      // free every model separately
      free_dirim_emission(objmodels->obj_list[iii]);
    }

  // free the rest
  free(objmodels->obj_list);
  free(objmodels);

  // set the structure to null
  objmodels = NULL;
}

/**
 * Function: free_dirim_emission
 * The function frees all memory allocated in a direct
 * emission object.
 *
 * Parameters:
 * @param diremission - the direct emission model to free
 *
 * Returns:
 * @return -
 */
void
free_dirim_emission(dirim_emission *diremission)
{
  // free the matrix with the direct emission
  gsl_matrix_free(diremission->modimage);

  // free the whole structure
  free(diremission);

  // set the struct to NULL
  diremission = NULL;
}

/**
 * Function: print_object_models
 * The function prints a direct object model structure.
 *
 * Parameters:
 * @param objmodels - the object model structure to 'print'
 *
 * Returns:
 * @return -
 */
void
print_object_models(const object_models *objmodels)
{
  int i;

  // go over each model
  for (i=0; i < objmodels->n_models; i++)
    {
      // print each model
      fprintf(stdout, "Model No.: %i\n", i);
      print_dirim_emission(objmodels->obj_list[i]);
    }
}

/**
 * Function: print_dirim_emission
 * The function prints a direct emission object.
 *
 * Parameters:
 * @param diremission - the direct emission model to print
 *
 * Returns:
 * @return -
 */
void
print_dirim_emission(const dirim_emission *diremission)
{
  fprintf(stdout, "dim_x: %i, dim_y: %i\n", diremission->dim_x, diremission->dim_y);
  fprintf(stdout, "xmean: %f, ymean: %f\n", diremission->xmean, diremission->ymean);
}

/**
 * Function: get_dirim_emission
 * The function returns from a direct object model the direct emission
 * object addressed by the index. The function translates from the
 * object index, which is in [1, n_model] to the structure index
 * in [0, n_model-1]
 *
 * Parameters:
 * @param objmodels - the object model structure
 * @param objspec   - the 'index' of the desired direct emission model
 *
 * Returns:
 * @return dirim - the desired direct emission model
 */
dirim_emission *
get_dirim_emission(const object_models *objmodels, const int objspec)
{
  dirim_emission *dirim;

  if (objspec >  objmodels->n_models || objspec < 1)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_PETCONT: " "Direct emission model number %i does NOT exist!",
		 objspec);

  // extract the desired SED
  dirim = objmodels->obj_list[objspec-1];

  // return the SED
  return dirim;
}

/**
 * Function: has_aperture_dirim
 * The function checks whether an object has an associated direct emission
 * structure or not. The corresponding boolean integer value is returned.
 *
 * Parameters:
 * @param objmodels - the object models structure
 * @param actobject - the object structure
 *
 * Returns:
 * @return has_dirim - boolean whether direct emission exists (=1) or not (=0)
 */
int
has_aperture_dirim(const object_models *objmodels, const object *actobject)
{
  int has_dirim=0;

  // check whether the direct emission model
  // structure is defined
  if (objmodels)
    // check whether the object has beams
    if (actobject->nbeams > 0)
      // check whether the direct emission number makes sense
      if (actobject->beams[0].modimage > 0) //&& actobject->beams[0].modimage <= objmodels->n_models)
	// report an existing direct emission
	has_dirim = 1;

  // return the result
  return has_dirim;
}

/**
 * Function: get_diremission_value
 * The function computes the direct emission value for a given direct
 * emission object and coordinates. The coordinates are given in x- and
 * y-position with respect to the mean position of the direct emission
 * object.
 *
 * Parameters:
 * @param diremission - the direct emission model
 * @param xpos        - the relative x-position
 * @param ypos        - the relative y-position
 *
 * Returns:
 * @return value - the emission value
 */
double
get_diremission_value(const dirim_emission *diremission,
		      const double xpos, const double ypos)
{
  double x_abs;
  double y_abs;

  double value=0.0;

  // compute the absolute positions
  // in the coordinate system
  // of the matrix
  x_abs = xpos + diremission->xmean;
  y_abs = ypos + diremission->ymean;

  //fprintf(stderr, " from x,y= % e,% e ",  x_abs, y_abs);

  // check whether the position
  // is outside the matrix
  if (x_abs < 0.0 || y_abs < 0.0
      || x_abs >= (float)(diremission->dim_x-1) || y_abs >= (float)(diremission->dim_y-1))
    // assign 0.0
    {
      //fprintf(stderr, "FIX %i,%i  % e,% e ", (int)x_abs, (int)y_abs, x_abs, y_abs);
      value = 0.0;
      //fprintf(stdout, "NI \n",  x_abs, y_abs);
    }
  else
    {  //  fprintf(stdout, "from x,y= %e,%e",  x_abs, y_abs);
    // assign the bilinear interpolated value
      value = bilin_interp_matrix(diremission->modimage, x_abs, y_abs);
      //fprintf(stdout, "IN \n",  x_abs, y_abs);
    }
  // return the value
  return value;
}

/**
 * Function: bilin_interp_matrix
 * The function computes and returns the bilinear interpolated value
 * of a given matrix at a given coordinate position. The algorithm
 * was taken from the numerical recipes.
 *
 * Parameters:
 * @param modimage - the direct emission model
 * @param x        - the x-position
 * @param y        - the y-position
 *
 * Returns:
 * @return value - the interpolated value
 */
double
bilin_interp_matrix(const gsl_matrix *modimage, const double x, const double y)
{
  int x_lower;
  int y_lower;

  double v1, v2, v3, v4;
  double t, u;

  double value;

  // get the indices of the 'lower'
  // pixel corner
  x_lower = (int) x;
  y_lower = (int) y;

  // check whether a pixel center is directly hit
  if (fabs((float)x_lower-x) < 1.0e-06 && fabs((float)y_lower-y) < 1.0e-06)
    {
      // give the full pixel value
      value =  gsl_matrix_get(modimage, x_lower  , y_lower  );
    }
  else
    {
      // extract the matrix values at the corners surrounding
      // the requested interpolation position
      v1 = gsl_matrix_get(modimage, x_lower  , y_lower  );
      v2 = gsl_matrix_get(modimage, x_lower+1, y_lower  );
      v3 = gsl_matrix_get(modimage, x_lower+1, y_lower+1);
      v4 = gsl_matrix_get(modimage, x_lower  , y_lower+1);

      // compute the pixel offset
      // from the lower corner
      t = x - (float)x_lower;
      u = y - (float)y_lower;

      // compute the bi-linear interpolated value
      // (algorithm taken from 'Numerical Recipes'
      value = (1-t)*(1-u)*v1 + t*(1-u)*v2 + t*u*v3 + (1-t)*u*v4;
    }

  // return the interpolated value
  return value;
}


/**
 * Function: load_spectral_models
 * The function file a spectral models structure with models extracted
 * from a fits file.
 *
 * Parameters:
 * @param  spectral_models_file - pathname to the spectral models file
 *
 * Returns:
 * @return smodels - the tracedata structure
 */
spectral_models *
load_spectral_models(const char spectral_models_file[])
{
  int n_models=0;
  spectral_models *smodels;

  // give feedback to the user
  fprintf (stdout, "aXe_PETCONT: Loading spectral models: %s", spectral_models_file);

  // allocate space for the spectral models structure
  smodels = (spectral_models *)malloc(sizeof(spectral_models));

  // get the number of models
  n_models = get_num_extensions(spectral_models_file);

  // allocate space for the models
  //smodels->SEDlist = (energy_distrib **)malloc(n_models * sizeof(energy_distrib *));
  smodels->SEDlist = load_SED_list(spectral_models_file, n_models);
  smodels->n_models =  n_models;

  // report the number of models loaded
  fprintf (stdout,"%d models loaded.\n", smodels->n_models);

  // return the structure
  return smodels;
}

/**
 * Function: print_spectral_models
 * The function prints the content of a spectral models
 * structure onto the screen
 *
 * Parameters:
 * @param  smodels - the spectral models structure
 *
 * Returns:
 * @return -
 */
void
print_spectral_models(const spectral_models *smodels)
{
  int i=0;
  int n=0;

  for (i=0; i < smodels->n_models; i++)
    {
      // print the model number
      fprintf(stdout,"Model number: %i\n", i);

      // print the model data
      for (n=0; n < smodels->SEDlist[i]->npoints; n++)
	fprintf(stdout,"%i Wavelength: %e, Flux: %e\n", n, smodels->SEDlist[i]->wavelength[n], smodels->SEDlist[i]->flux[n]);
    }
}

/**
 * Function: load_SED_list
 * The function creates a list of energy distributions from
 * data stored in several extensions of a fits file.
 *
 * Parameters:
 * @param  spectral_models_file - pathname to the spectral models file
 * @param  n_models             - the number of spectral models in the file
 *
 * Returns:
 * @return SEDlist              - the list of energy distributions
 */
energy_distrib **
load_SED_list(const char spectral_models_file[], const int n_models)
{
  int i;
  int f_status=0;
  int hdutype;
  energy_distrib **SEDlist;
  fitsfile *s_models;

  // allocate space for the SED list
  SEDlist = (energy_distrib **)malloc(n_models * sizeof(energy_distrib *));

  // open the fits file
  // report any error
  fits_open_file (&s_models, spectral_models_file, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETCONT: " "Could not open file: %s",
		   spectral_models_file);
    }

  // go over all extensions
  for (i=0; i < n_models; i++)
    {
      // move to the right extension
      fits_movabs_hdu (s_models, i+2, &hdutype, &f_status);

      // load the energy distribution
      SEDlist[i] = load_SED_from_fitsext(spectral_models_file, s_models);
    }

  // close the fits file
  fits_close_file (s_models, &f_status);

  // return the SED list
  return SEDlist;
}


/**
 * Function: load_SED_from_fitsext
 * The function creates a energy distribution from the data stored
 * in a fits file extension. The data must be stored in the columns
 * "wavelength" and "flux".
 *
 * Parameters:
 * @param  spectral_models_file - pathname to the spectral models file
 * @param  s_models             - pointer to the fits file extension
 *
 * Returns:
 * @return sed - the energy distribution created
 */
energy_distrib *
load_SED_from_fitsext(const char spectral_models_file[], fitsfile *s_models)
{
  int f_status=0;
  int anynul;
  long nrows=0;
  int colnum1;
  int colnum2;

  energy_distrib *sed;
  double *sed_wavs;
  double *sed_flux;

  // allocate memory for the energy distribution
  sed = (energy_distrib *) malloc(sizeof(energy_distrib));

  // get number of rows
  fits_get_num_rows (s_models, &nrows, &f_status);
  if (f_status) {
    ffrprt (stderr, f_status);
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "load_SED_from_fitsext: "
		 "Could not determine the number of rows in"
		 " table %s",spectral_models_file);
  }

  // allocate memory for the data
  sed_wavs = (double *) malloc(nrows*sizeof(double));
  if (!sed_wavs) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }
  sed_flux = (double *) malloc(nrows*sizeof(double));
  if (!sed_flux) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // get the column number for the wavelength
  fits_get_colnum (s_models, CASEINSEN, "WAV_NM", &colnum1, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_interp_ftable: "
		   "Could not determine column %s in "
		   " table %s", "WAV_NM", spectral_models_file);
    }
  // read the wavelength
  fits_read_col (s_models, TDOUBLE, colnum1, 1, 1, nrows, NULL, sed_wavs,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "load_SED_from_fitsext: "
		   "Could not read content of WAVELENGTH column "
		   " from BINARY table %s", spectral_models_file);
    }

  // get the column number for the flux
  fits_get_colnum (s_models, CASEINSEN, "FLUX", &colnum2, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
  		   "create_interp_ftable: "
  		   "Could not determine column %s in "
  		   " table %s", "FLUX", spectral_models_file);
    }
  // read the flux column
  fits_read_col (s_models, TDOUBLE, colnum2, 1, 1, nrows, NULL, sed_flux,
                    &anynul, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "load_SED_from_fitsext: "
		   "Could not read content of FLUX column "
		   " from BINARY table %s", spectral_models_file);
    }

  // transfer the vector and length information
  // to the SED object
  sed->npoints    = nrows;
  sed->wavelength = sed_wavs ;
  sed->flux       = sed_flux;

  // allocate some structures for the interpolator
  sed->interp     = gsl_interp_alloc (SMODEL_INTERP_TYPE, (size_t)sed->npoints );
  sed->accel      = gsl_interp_accel_alloc ();

  // initialize the iterpolator
  gsl_interp_init (sed->interp, sed->wavelength, sed->flux, (size_t)sed->npoints);

  // return the energy distribution
  return sed;
}


/**
 * Function: get_num_extensions
 * The function determines the number of extensions in
 * the fits file. Here the primary extension does not count.
 *
 * Parameters:
 * @param  spectral_models_file - pathname to the spectral models file
 *
 * Returns:
 * @return n_models              - the number of spectral models
 */
int
get_num_extensions(const char spectral_models_file[])
{
  int n_ext=0;
  int n_models=0;
  //int i;
  int f_status=0;


  fitsfile *s_models;

  // open the fits file
  // report any error
  fits_open_file (&s_models, spectral_models_file, READONLY, &f_status);
  if (f_status)
    {
      ffrprt (stderr, f_status);
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "aXe_PETCONT: " "Could not open file: %s",
		   spectral_models_file);
    }

  // do until a fits error occurs
  while (!f_status)
    {
      // move one extension forward
      fits_movrel_hdu (s_models, 1, NULL, &f_status);

      // count up
      n_ext++;
    }

  // close the fits file
  fits_close_file (s_models, &f_status);

  // the zeroth extension
  // does nnot count!!
  n_models = n_ext - 1;

  // return the
  return n_models;
}


/**
 * Function: free_spectral_models
 * Releases the memory allocated in a spectral
 * models structure
 *
 * Parameters:
 * @param  smodels - the spectral models structure
 *
 * Returns:
 * @return -
 */
void
free_spectral_models(spectral_models *smodels)
{
  int i=0;

  // go over all energy distributions
  for (i=0; i < smodels->n_models; i++)
    {
      // free the energy distribution
      free_enerdist(smodels->SEDlist[i]);
      smodels->SEDlist[i] = NULL;
    }

  // free the SED list
  free(smodels->SEDlist);
  smodels->SEDlist = NULL;

  // free the structure
  free(smodels);
  smodels = NULL;
}


/**
 * Function: get_model_sed
 * Extracts and returns the desired model spectrum from th model structure.
 * Translates from the object index [1, n_models] to the spectral model
 * structure index [0, n_model-1].
 *
 * Parameters:
 * @param  spectral_models - the spectral models structure
 * @param  modspec         - the index of the desired SED
 *
 * Returns:
 * @return -
 */
energy_distrib *
get_model_sed(const spectral_models *spec_mod, const int modspec)
{
  energy_distrib *sed;

  if (modspec >  spec_mod->n_models)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "aXe_PETCONT: " "Model number %i does NOT exist!",
		 modspec);

  // extract the desired SED
  sed = spec_mod->SEDlist[modspec-1];

  // return the SED
  return sed;
}
