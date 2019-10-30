/**
 */
#include "aXe_grism.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "fringe_conf.h"
#include "fringe_model.h"
#include "fringe_utils.h"

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


/**
 * Function: fringe_correct_PET
 * Performs fringe corrections to PET pixel lists
 *
 * The main function to compute and perform the fringe
 * correction to a list of PET pixels. for every pixel in
 * the PET list the pixel throughput function is 
 * computed, and the relevant information from
 * the fringing model is extracted. Then the fringe
 * amplitude is computed and the count value of the pixel
 * is corrected. Also the PET pixel list is simultaneously
 * parsed and corrected for the fringing.
 *
 * Parameters:
 * @param fconf    - the fringe configuration file
 * @param act_beam - the beam aperture
 * @param obj_pet  - the PET pixel
 * @param bck_pet  - the PET pixel
 */
void
fringe_correct_PET(const fringe_conf *fconf, const beam act_beam,
		   ap_pixel *obj_pet, ap_pixel *bck_pet)
{
  ap_pixel *obj_pix=NULL;
  ap_pixel *bck_pix=NULL;

  //d_point angles;
  pixel_tput *p_tput ;
  gsl_vector **tput_vectors=NULL;

  //int ii=0;
  //int jj=0;
  int index;
  int pix_index;

  double lambda_mean;
  double pixel_ampl;
  double fringe_factor;
  double fringe_tot;

  optical_property *optprops;

  // allocate the pixel throughput structure
  p_tput = alloc_pixel_tput();

  // allocate memory for the optical property structure  
  optprops = alloc_optprops_list(fconf);

  // compute the mean wavelength
  lambda_mean = (gsl_vector_get(fconf->fringe_range, 0) +
    (gsl_vector_get(fconf->fringe_range, 1)
     - gsl_vector_get(fconf->fringe_range, 0))/2.0)/1.0e+04;

  // initialize some values in the optical property list
  init_optprops_list(fconf, lambda_mean, optprops);

  // initialize the PET pixel storages
  obj_pix = obj_pet;
  bck_pix = bck_pet;

  fringe_tot = 0.0;
  pix_index = 0;
  // go over the whole PET list
  while (obj_pix->p_x != -1)
    {
      // check whether the pixel has to be fringe corrected
      // at all, based on boundaries given in the configuration
      if ( obj_pix->lambda < gsl_vector_get(fconf->fringe_range, 0)
      	   || obj_pix->lambda > gsl_vector_get(fconf->fringe_range, 1)
      	   || obj_pix->dlambda > fconf->max_dispersion)
      	{

      	  // count up the object PET 
      	  obj_pix++;

      	  // count up the background PET
      	  if (bck_pet)
      	    bck_pix++;
      
      	  // jump to the next PET pixel
      	  continue;
      	}

      // determine the filter throughput values
      //      tput_vectors = evaluate_pixel_throughput(fconf, act_beam, obj_pix);
      // for the ISR with the WFC value
      tput_vectors = get_gauss_throughput(fconf, obj_pix, 80.0);

      // fill the optical thickness of the layers
      // into the structure
      fill_optprops_thickness(fconf->opt_layers, obj_pix->p_x,
			      obj_pix->p_y, optprops);

      pixel_ampl = 0.0;
      for (index=0; index < (int)tput_vectors[0]->size; index++)
	{
	  // fill all information in the optical
	  // property list
	  fill_optprops_all(fconf->opt_layers,
			    gsl_vector_get(tput_vectors[0],index),
			    optprops);
	  
	  // compute and add the contribution at a wavelength
	  pixel_ampl += gsl_vector_get(tput_vectors[1],index)*
	    fringe_contrib_single(optprops, fconf);
	}

      //put together the fringe correction factor
      fringe_factor = fconf->fringe_amp * pixel_ampl + 1.0;
#ifdef DEBUGFCONF
      fprintf(stderr, "Wavelength: %f, dispersion: %f, fringe factor: %f, (x,y): (%i, %i)\n", 
	      obj_pix->lambda, obj_pix->dlambda, fringe_factor, obj_pix->p_x, obj_pix->p_y);
#endif

      // correct the object PET
      obj_pix->count /= fringe_factor;
      //obj_pix->count = fringe_factor;
      //fprintf(stdout, "%e ",obj_pix->count); 

      fringe_tot += fringe_factor;
      pix_index++;


      // treat the background PET pixel
      if (bck_pet)
	{
	  if (bck_pix->p_x != obj_pix->p_x || bck_pix->p_y != obj_pix->p_y)
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "aXe_FRINGECORR: The sequence of pixels\n"
			 "in the object PET and background PET differs.\n"
			 "Something must be wrong!\n");

	  // correct also the background PET
	  bck_pix->count /= fringe_factor;

	  // count up the background PET
	  bck_pix++;
	}

      // count up the object PET
      obj_pix++;
      //      obj_pix->p_x = -1;
    }

  //  if (pix_index)
  //    fprintf(stdout,"Mean fringe factor: %f\n", fringe_tot/(double)pix_index);

  if (tput_vectors)
    {
      // release the memory in the vectors
      gsl_vector_free(tput_vectors[0]);
      gsl_vector_free(tput_vectors[1]);
      free(tput_vectors);
    }

  // free the optical property structure
  free_optprops_list(optprops);

  // free the memory
  free_pixel_tput(p_tput);
}


/**
 * Function: evaluate_pixel_throughput
 * Computes the wavelength and filter throughputs for a PET pixel.
 *
 * The function computes and returns two vector, one with wavelength
 * steps and the second with the pixel throughputs at those wavelengths.
 * This data is computed to a single PET pixel according to the settings
 * in the fringe configuration file.
 * 
 * Parameters:
 * @param fconf    - the fringe configuration file
 * @param act_beam - the beam aperture
 * @param obj_pix  - the PET pixel
 *
 * Returns:
 * @return double_vector - the wavelength and filter throughput steps
 */
gsl_vector **
evaluate_pixel_throughput(const fringe_conf *fconf,
			  const beam act_beam, const ap_pixel *obj_pix)
{
  pixel_tput *p_tput ;

  gsl_vector **double_vector;
  gsl_vector *lambda_values;
  gsl_vector *through_values;

  interpolator *interp;

  double stepsize;
  double lambda_act;
  double through_act;
  double through_tot;

  double abs_min;
  double abs_max;

  int index;

  // allocate the pixel throughput structure
  p_tput = alloc_pixel_tput();

  // allocate the space for the vectors
  lambda_values  = gsl_vector_alloc(fconf->num_steps + 1);
  through_values = gsl_vector_alloc(fconf->num_steps + 1);

  // allocate space for the return vector
  double_vector  = (gsl_vector **)malloc(2*sizeof (gsl_vector *));

  compute_pixel_tput(act_beam, obj_pix, p_tput);
#ifdef DEBUGFCONF
  print_pixel_tput(p_tput);
#endif
  // convert the pixel throughput information
  // into an interplator 
  interp = create_interp_tput(p_tput);

  abs_min = MAX(interp->xmin, gsl_vector_get(fconf->fringe_range, 0));
  abs_max = MIN(interp->xmax, gsl_vector_get(fconf->fringe_range, 1));
  
  // compute the steps in wavelength
  stepsize = (abs_max - abs_min)/(double)fconf->num_steps;

  // set the lower wavelength end
  lambda_act = abs_min;

  // loop over the wavelength
  through_tot = 0.0;
  for (index=0; index < fconf->num_steps; index++)
    {
      // set the wavelength value
      gsl_vector_set(lambda_values, index, lambda_act);

      // compute the throughput
      through_act = eval_interp(interp, lambda_act);

      // compute the total throughput
      through_tot += through_act;

      // set the thoughput value
      gsl_vector_set(through_values, index, through_act);

      lambda_act +=  stepsize;
    }

  // treat the last step separately, avoiding
  // to leave the allowed interpolation range
  // due to rounding errors 
  gsl_vector_set(lambda_values, fconf->num_steps, abs_max);
  through_act = eval_interp(interp, abs_max);
  through_tot += through_act;
  gsl_vector_set(through_values, fconf->num_steps, through_act);

  for (index=0; index < (int)through_values->size; index++)
    {
      // normalize the filter throughput values
      gsl_vector_set(through_values, index,
		     gsl_vector_get(through_values, index)/through_tot);

      // convert the wavelength from AA to micron
      gsl_vector_set(lambda_values, index,
		     gsl_vector_get(lambda_values, index)*1.0e-04);
    }
#ifdef DEBUGFCONF
  through_tot = 0.0;
  for (index=0; index < through_values->size; index++)
    {
      fprintf(stderr, "Wavelength: %f, Throughput: %f\n",
	      gsl_vector_get(lambda_values, index),
	      gsl_vector_get(through_values, index));
      through_tot += gsl_vector_get(through_values, index);
    }
  fprintf(stderr, "Total throughput: %f\n", through_tot);
#endif
  
  // build up the output array
  double_vector[0] = lambda_values;
  double_vector[1] = through_values;

  
  // free the memory of the interpolator
  free_interp(interp);

  // free the memory of the pixel throughput
  free_pixel_tput(p_tput);

  // return the two gsl vectors
  return double_vector;
}


/**
 * Function: get_gauss_throughput
 *
 * Parameters:
 * @param fconf    - the fringe configuration file
 * @param act_beam - the beam aperture
 * @param obj_pix  - the PET pixel
 *
 * Returns:
 * @return double_vector - the wavelength and filter throughput steps
 */
gsl_vector **
get_gauss_throughput(const fringe_conf *fconf, const ap_pixel *obj_pix,
		     const double fwhm)
{
  pixel_tput *p_tput ;

  gsl_vector **double_vector;
  gsl_vector *lambda_values;
  gsl_vector *through_values;

  //interpolator *interp;

  double stepsize;
  double lambda_act;
  double through_act;
  double through_tot;

  double abs_min;
  double abs_max;

  double sigma;
  double arg;
  
  int index;

  // convert the FWHM into sigma
  sigma = fwhm / 2.35482;

  // allocate the pixel throughput structure
  p_tput = alloc_pixel_tput();

  // allocate the space for the vectors
  lambda_values  = gsl_vector_alloc(fconf->num_steps + 1);
  through_values = gsl_vector_alloc(fconf->num_steps + 1);

  // allocate space for the return vector
  double_vector  = (gsl_vector **)malloc(2*sizeof (gsl_vector *));

  // compute the minimum and maximum wavelength
  abs_min = MAX(obj_pix->lambda - 2.5 * fwhm,
		gsl_vector_get(fconf->fringe_range, 0));
  abs_max = MIN(obj_pix->lambda + 2.5 * fwhm,
		gsl_vector_get(fconf->fringe_range, 1));

  // compute the steps in wavelength
  stepsize = (abs_max - abs_min)/(double)fconf->num_steps;

  // set the lower wavelength end
  lambda_act = abs_min;

  // loop over the wavelength
  through_tot = 0.0;
  for (index=0; index < fconf->num_steps; index++)
    {
      // set the wavelength value
      gsl_vector_set(lambda_values, index, lambda_act);

      // z=(coord-g(2))/g(3)
      // sumgauss=sumgauss+DBLE(g(1)*EXP(-0.5*z**2))
      arg = (lambda_act - obj_pix->lambda) / sigma;

      // compute the throughput
      through_act = exp(-0.5*arg*arg);

      // compute the total throughput
      through_tot += through_act;

      // set the thoughput value
      gsl_vector_set(through_values, index, through_act);

      // increment the wavelength
      lambda_act += stepsize;
    }

  // treat the last step separately, avoiding
  // to leave the allowed interpolation range
  // due to rounding errors 
  gsl_vector_set(lambda_values, fconf->num_steps, abs_max);
  arg = (abs_max - obj_pix->lambda) / sigma;
  through_act = exp(-0.5*arg*arg);
  through_tot += through_act;
  gsl_vector_set(through_values, fconf->num_steps, through_act);


  for (index=0; index < (int)through_values->size; index++)
    {
      // normalize the filter throughput values
      gsl_vector_set(through_values, index,
		     gsl_vector_get(through_values, index)/through_tot);

      // convert the wavelength from AA to micron
      gsl_vector_set(lambda_values, index,
		     gsl_vector_get(lambda_values, index)*1.0e-04);
    }

#ifdef DEBUGFCONF
  through_tot = 0.0;
  for (index=0; index < through_values->size; index++)
    {
      fprintf(stderr, "Wavelength: %f, Throughput: %f\n",
	      gsl_vector_get(lambda_values, index),
	      gsl_vector_get(through_values, index));
      through_tot += gsl_vector_get(through_values, index);
    }
  fprintf(stderr, "Total throughput: %f\n", through_tot);
#endif

  // build up the output array
  double_vector[0] = lambda_values;
  double_vector[1] = through_values;

  return double_vector;
}


/**
 * Function: create_interp_tput
 * The function creates and returns an interpolator
 * for a pixel throughput function. The interpolator
 * structure is allocated, data is filled in and
 * finally initialized.
 *
 * Parameters:
 * @param p_tput   - the pixel throughput structure
 * 
 * Returns:
 * @return interp - the new interpolator
 */
interpolator *
create_interp_tput(const pixel_tput *p_tput)
{
  interpolator *interp;

  double *xvals;
  double *yvals;

  // allocate memory for the independent values 
  xvals = (double *) malloc(4 * sizeof(double));
  if (!xvals) { 
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // allocate memory for the dependent values 
  yvals = (double *) malloc(4 * sizeof(double));
  if (!yvals) {
    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
		 "Memory allocation failed");
  }

  // set the independent vector values 
  xvals[0] = gsl_vector_get(p_tput->lambda_values, 0);
  xvals[1] = gsl_vector_get(p_tput->lambda_values, 1);
  xvals[2] = gsl_vector_get(p_tput->lambda_values, 2);
  xvals[3] = gsl_vector_get(p_tput->lambda_values, 3);

  // set the dependent vector values
  yvals[0] = 0.0;
  yvals[1] = p_tput->tp_max;
  yvals[2] = p_tput->tp_max;
  yvals[3] = 0.0;


  // create the interpolator
  interp = create_interp(4, FILTER_INTERP_TYPE, xvals, yvals);

  // return the interpolator
  return interp;
}


/**
 * Function: compute_pixel_tput
 * The function computes all relevant data to characterize
 * the pixel throughput function for an individual PET pixel.
 * The information is written into a pixel throughput structure.
 *
 * Parameters:
 * @param act_beam - the beam aperture
 * @param obj_pix  - the PET pixel
 * @param p_tput   - the pixel throughput structure
 */
void
compute_pixel_tput(const beam act_beam, const ap_pixel *obj_pix,
		   pixel_tput *p_tput)
{
  d_point angles;

  double theta_2a;
  double theta_3a;

  double theta_1b;
  double theta_2b;

  double sin_theta_2a;
  double sin_theta_3a;

  double sin_theta_1b;
  double cos_theta_1b;
  double sin_theta_2b;

  double dist1;
  double dist2;

  // compute the anngles alpha and beta
  angles = compute_tput_angles(act_beam, obj_pix);

  // the computation of p_min-p3 in the
  // design document
  theta_1b = M_PI_2 + angles.x - angles.y;
  theta_2b = angles.y;

  // compute the sin's and cos's neede
  sin_theta_1b = sin(theta_1b);
  cos_theta_1b = cos(theta_1b);
  sin_theta_2b = sin(theta_2b);

  // compute one delta
  dist2 = obj_pix->dlambda * sin_theta_1b / sin_theta_2b;

  // the computation of p_h-p0 in the
  // design document
  theta_2a = M_PI - angles.y;
  theta_3a = angles.y - angles.x - M_PI_4;

  // next line due to sin(angle) = sin(180-angle)
  sin_theta_2a = sin_theta_2b;
  sin_theta_3a = sin(theta_3a);

  // compute the second delta 
  dist1 = obj_pix->dlambda * M_SQRT1_2 * sin_theta_3a / sin_theta_2a;

  // determine the maximum throughput
  p_tput->tp_max = 1.0/MAX(fabs(sin_theta_1b), fabs(cos_theta_1b));

  // put together the distances, using
  // the deltas and the dispersion
  gsl_vector_set(p_tput->lambda_values, 0, obj_pix->lambda - dist1 - dist2);
  gsl_vector_set(p_tput->lambda_values, 1, obj_pix->lambda - dist1);
  gsl_vector_set(p_tput->lambda_values, 2, obj_pix->lambda + dist1);
  gsl_vector_set(p_tput->lambda_values, 3, obj_pix->lambda + dist1 + dist2);

  // sort the distances
  gsl_sort_vector(p_tput->lambda_values);
}


/**
 * Function: compute_tput_angles
 * The function computes the angles which characterize
 * a pixel in the lambda-cross dispersion plane.
 * In detail this is the angle between the lambda-axis
 * and the x-axis, and the angle between the lambda abd
 * the cross dispersion axis.
 * The two angles are computed from the basic pixel
 * and beam data.
 *
 *
 * Parameters:
 * @param act_beam - the beam aperture
 * @param act_pet  - the PET pixel
 *
 * Returns:
 * @return angles - the relevant pixel angles
 */
d_point
compute_tput_angles(const beam act_beam, const ap_pixel *act_pet)
{
  d_point angles;

  // compute the angle 'alpha' between lambda-axis and
  // x-axis
  angles.x = -atan2(act_beam.spec_trace->deriv(act_pet->xi,
					       act_beam.spec_trace->data),1.0);
  // compute the angle beta between lambda and cross dispersion axis
  angles.y = act_beam.orient + angles.x;

  // return the two angles
  return angles;
}


/**
 * Function: alloc_pixel_tput
 * The function allocates and returns memory
 * for a pixel throughput function
 *
 * Returns:
 * @return p_tput -the new allocated pixel throughput structure
 */
pixel_tput *
alloc_pixel_tput()
{
  pixel_tput *p_tput;

  // allcoate the pointer
  p_tput = (pixel_tput *) malloc(sizeof(pixel_tput));

  // allocate the vector in the structure
  p_tput->lambda_values = gsl_vector_alloc(4);

  // return the structure
  return p_tput;
}

/**
 * Function: print_pixel_tput
 * The function prints the content of a pixel throughput
 * structure onto the screen. 
 *
 * Parameters:
 * @param p_tput - the pixel throughput structure
 */
void
print_pixel_tput(pixel_tput *p_tput)
{
  fprintf(stdout, "l_min: %f, l_1: %f, l_2: %f, l_max: %f; tp_max: %f\n\n",
	  gsl_vector_get(p_tput->lambda_values,0),
	  gsl_vector_get(p_tput->lambda_values,1),
	  gsl_vector_get(p_tput->lambda_values,2),
	  gsl_vector_get(p_tput->lambda_values,3),
	  p_tput->tp_max);
}


/**
 * Function: free_pixel_tput
 * The function releases the memory allocated
 * in a pixel throughput function.
 *
 * Parameters:
 * @param p_tput - the pixel throughput structure
 */
void
free_pixel_tput(pixel_tput *p_tput)
{
  // free the gsl-vector
  gsl_vector_free(p_tput->lambda_values);

  // free everything
  free(p_tput);

  // set it to NULL
  p_tput = NULL;
}

