/**
 */
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include "aXe_utils.h"
#include "aXe_errors.h" 
#include "fringe_conf.h"
#include "fringe_model.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

/**
 * Function: compute_fringe_amplitude
 * The function computes the fringe image for a CCD setup stored
 * in a fringe configuration structure. This fringe configuration
 * structure completely describes the problem, and this function
 * executes the loops for every pixel over the wavelength
 * range spanned by the filter.
 *
 * Parameters
 * @param fconf - the fringe configuration structure
 *
 * Returns:
 * @return fringe_image - the image with the computed fringe amplitudes
 */
gsl_matrix  *
compute_fringe_amplitude(fringe_conf *fconf)
{
  gsl_matrix *fringe_image;
  
  gsl_vector **filter_vectors;

  int index=0;
  int ii=0;
  int jj=0;

  double lambda_mean;
  double pixel_ampl;
  //double phase_number;

  optical_property *optprops;

  // allocate the fringe image
  fringe_image = alloc_fringe_image(fconf->opt_layers);

  // allocate memory for the optical property structure  
  optprops = alloc_optprops_list(fconf);

  // find the exact wavelength range,
  // and define the wavelength values and
  // normalized fiter throughput there
  filter_vectors = evaluate_wavelength_steps(fconf);
  //filter_vectors = get_PET_calibration_data();

  // compute the mean wavelength
  lambda_mean = gsl_vector_get(filter_vectors[0],filter_vectors[0]->size-1)/2.0
		 + gsl_vector_get(filter_vectors[0],0)/2.0;

  // initialize some values in the optical property list
  init_optprops_list(fconf, lambda_mean, optprops);

    for (ii=0; ii < (int)fringe_image->size1; ii++)
      //    for (ii=0; ii < 2; ii++)
    {
      fprintf(stderr, "Computing row No.: %i\n", ii);
            for (jj=0; jj < (int)fringe_image->size2; jj++)
	      //for (jj=0; jj < 2; jj++)
	{
	  // fill the optical thickness of the layers
	  // into the structure
	  fill_optprops_thickness(fconf->opt_layers, ii, jj, optprops);

	  pixel_ampl = 0.0;
	  for (index=0; index < (int)filter_vectors[0]->size; index++)
	    {
	      // fill all information in the optical
	      // property list
	      fill_optprops_all(fconf->opt_layers,
				gsl_vector_get(filter_vectors[0],index),
				optprops);

	      // compute and add the contribution at a wavelength
	      pixel_ampl += gsl_vector_get(filter_vectors[1],index)*
		fringe_contrib_single(optprops, fconf);
	    }

	  // finally set the pixel value
	  // in the output image
	  gsl_matrix_set(fringe_image, ii, jj,
			 fconf->fringe_amp * pixel_ampl + 1.0);
	}
    }

  // release the memory in the vectors
  gsl_vector_free(filter_vectors[0]);
  gsl_vector_free(filter_vectors[1]);
  free(filter_vectors);

  // free the optical property structure
  free_optprops_list(optprops);

  // return the fringe image
  return fringe_image;
}


/**
 * Function: alloc_fringe_image
 * The function browses through the the structure for the CCD layers
 * and extracts all information on image sizes stored there in
 * one or several thickness images. This information is checked
 * for consistency.
 * Finally, a matrix for a fringe image is allocated, initialized
 * and returned.
 *
 * Parameters:
 * @param opt_layers - the optical layers in the CCD
 *
 * Returns:
 * @return fringe_image - the allocated matrix for the fringe image
 */
gsl_matrix  *
alloc_fringe_image(const ccd_layers *opt_layers)
{
  int n1         = 0;
  int n2         = 0;
  int n1_new     = 0;
  int n2_new     = 0;
  int index      = 0;
  int is_defined = 0;
  
  gsl_matrix *fringe_image=NULL;
  // go through all CCD layers
  for (index=0; index < opt_layers->num_layers; index++)
    {
      // check whether the thickness information
      // is represented by an image
      if (opt_layers->opt_layer[index]->thickness2D != NULL)
	{

	  // mark that an image was found
	  is_defined=1;

	  // get the dimension of the image
	  n1_new = opt_layers->opt_layer[index]->thickness2D->size1;
	  n2_new = opt_layers->opt_layer[index]->thickness2D->size2;
	
	  // check the frst axix value against
	  // previous values, if possible
	  if (n1 && n1_new != n1)
	    // give an error if the new value is different
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "Thickness image of layer %i "
			 "has %i pixels in first axis."
			 "This differs from the previous value %i\n",
			 index, n1_new, n1);
	  else
	    // store the new value
	    n1 = n1_new;
	  
	  // check the frst axix value against
	  // previous values, if possible
	  if (n2 && n2_new != n2)
	    // give an error if the new value is different
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "Thickness image of layer %i "
			 "has %i pixels in second axis."
			 "This differs from the previous value %i\n",
			 index, n2_new, n2);
	  else
	    // store the new value
	    n2 = n2_new;
	}
    }

#ifdef DEBUGFCONF
  fprintf(stdout, "Allocating an image with size: (%i, %i)\n", n1, n2);
#endif

  // check whether at least one image was
  // found as thickness information
  if (is_defined)
    {
      // allocate the image
      fringe_image =  gsl_matrix_alloc(n1, n2);
      
      // set all vallues to zero
      gsl_matrix_set_all(fringe_image, 0.0);
    }
  else
    {
      // no image found; report the error
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "The size of the fringe image is unknown\n"
		   "since none of the layers has an image "
		   "to specify its thickness!\n");
    }

  // return the image
  return fringe_image;
}


/**
 * Function: get_calibration_data()
 * This function is a data storage for filter throughput values
 * as computed by the orioginal program 'final.f' of J.Walsh.
 * The cubic spline routine used there turned out to give
 * (slightly) different results than the gsl-routine used
 * here. To be able making a in-depth comparison the throughuts
 * were computed for the filter 'ccd_ref.fits' and stored here.
 *
 * Parameters:
 *
 * Returns:
 * @return double_vector - interpolated standard filter throughputs
 */
gsl_vector **
get_calibration_data()
{
  gsl_vector *lambda_values;
  gsl_vector *through_values;

  gsl_vector **double_vector;


  // allocate space for the return vector
  double_vector  = (gsl_vector **)malloc(2*sizeof (gsl_vector *));

  // allocate the space for the vectors
  lambda_values  = gsl_vector_alloc(43);
  through_values = gsl_vector_alloc(43);


  // what follows is the wavelength steps as well as the normalized
  // filter throughputs which were computed in the original
  // program 'final.f' of J.Walsh.
  // The gsl-cubic spline routine yields slightly different values.
  // The values fixed here are good for a detailed cross-check
  // of the two programs.
  gsl_vector_set(through_values, 0, 0.0000000000000000);gsl_vector_set(lambda_values, 0,   0.9179);
  gsl_vector_set(through_values, 1, 0.0095066184313808);gsl_vector_set(lambda_values, 1,   0.9180);
  gsl_vector_set(through_values, 2, 0.0176006475945001);gsl_vector_set(lambda_values, 2,   0.9181);
  gsl_vector_set(through_values, 3, 0.0240259215044616);gsl_vector_set(lambda_values, 3,   0.9182);
  gsl_vector_set(through_values, 4, 0.0289185490230508);gsl_vector_set(lambda_values, 4,   0.9183);
  gsl_vector_set(through_values, 5, 0.0324146730103703);gsl_vector_set(lambda_values, 5,   0.9184);
  gsl_vector_set(through_values, 6, 0.0346504181940869);gsl_vector_set(lambda_values, 6,   0.9185);
  gsl_vector_set(through_values, 7, 0.0357619115684218);gsl_vector_set(lambda_values, 7,   0.9186);
  gsl_vector_set(through_values, 8, 0.0358852710613785);gsl_vector_set(lambda_values, 8,   0.9187);
  gsl_vector_set(through_values, 9, 0.0351566282002872);gsl_vector_set(lambda_values, 9,   0.9188);
  gsl_vector_set(through_values,10, 0.0337121099793692);gsl_vector_set(lambda_values,10,   0.9189);
  gsl_vector_set(through_values,11, 0.0316878433928458);gsl_vector_set(lambda_values,11,   0.9190);
  gsl_vector_set(through_values,12, 0.0292199509018293);gsl_vector_set(lambda_values,12,   0.9191);
  gsl_vector_set(through_values,13, 0.0264445572339866);gsl_vector_set(lambda_values,13,   0.9192);
  gsl_vector_set(through_values,14, 0.0234977916500934);gsl_vector_set(lambda_values,14,   0.9193);
  gsl_vector_set(through_values,15, 0.0205157811443709);gsl_vector_set(lambda_values,15,   0.9194);
  gsl_vector_set(through_values,16, 0.0176346504444861);gsl_vector_set(lambda_values,16,   0.9195);
  gsl_vector_set(through_values,17, 0.0149905288112146);gsl_vector_set(lambda_values,17,   0.9196);
  gsl_vector_set(through_values,18, 0.0127195375723916);gsl_vector_set(lambda_values,18,   0.9197);
  gsl_vector_set(through_values,19, 0.0109578003224067);gsl_vector_set(lambda_values,19,   0.9198);
  gsl_vector_set(through_values,20, 0.0098414542549762);gsl_vector_set(lambda_values,20,   0.9199);
  gsl_vector_set(through_values,21, 0.0095066184313808);gsl_vector_set(lambda_values,21,   0.9200);
  gsl_vector_set(through_values,22, 0.0100435357173845);gsl_vector_set(lambda_values,22,   0.9201);
  gsl_vector_set(through_values,23, 0.0113589226644173);gsl_vector_set(lambda_values,23,   0.9202);
  gsl_vector_set(through_values,24, 0.0133136264280558);gsl_vector_set(lambda_values,24,   0.9203);
  gsl_vector_set(through_values,25, 0.0157684646986687);gsl_vector_set(lambda_values,25,   0.9204);
  gsl_vector_set(through_values,26, 0.0185842800987238);gsl_vector_set(lambda_values,26,   0.9205);
  gsl_vector_set(through_values,27, 0.0216219095843027);gsl_vector_set(lambda_values,27,   0.9206);
  gsl_vector_set(through_values,28, 0.0247421844451007);gsl_vector_set(lambda_values,28,   0.9207);
  gsl_vector_set(through_values,29, 0.0278059450370313);gsl_vector_set(lambda_values,29,   0.9208);
  gsl_vector_set(through_values,30, 0.0306740203832353);gsl_vector_set(lambda_values,30,   0.9209);
  gsl_vector_set(through_values,31, 0.0332072463065172);gsl_vector_set(lambda_values,31,   0.9210);
  gsl_vector_set(through_values,32, 0.0352664608962358);gsl_vector_set(lambda_values,32,   0.9211);
  gsl_vector_set(through_values,33, 0.0367124999751954);gsl_vector_set(lambda_values,33,   0.9212);
  gsl_vector_set(through_values,34, 0.0374061970996460);gsl_vector_set(lambda_values,34,   0.9213);
  gsl_vector_set(through_values,35, 0.0372083835592830);gsl_vector_set(lambda_values,35,   0.9214);
  gsl_vector_set(through_values,36, 0.0359798997100196);gsl_vector_set(lambda_values,36,   0.9215);
  gsl_vector_set(through_values,37, 0.0335815745749970);gsl_vector_set(lambda_values,37,   0.9216);
  gsl_vector_set(through_values,38, 0.0298742553097917);gsl_vector_set(lambda_values,38,   0.9217);
  gsl_vector_set(through_values,39, 0.0247187709375447);gsl_vector_set(lambda_values,39,   0.9218);
  gsl_vector_set(through_values,40, 0.0179759414151792);gsl_vector_set(lambda_values,40,   0.9219);
  gsl_vector_set(through_values,41, 0.0095066184313808);gsl_vector_set(lambda_values,41,   0.9220);
  gsl_vector_set(through_values,42, 0.0000000000000000);gsl_vector_set(lambda_values,42,   0.9221);
  
  // build up the output array
  double_vector[0] = lambda_values;
  double_vector[1] = through_values;

  // return the two vectors
  return double_vector;
}


/**
 * Function: get_PET_calibration_data()
 */
gsl_vector **
get_PET_calibration_data()
{
  gsl_vector *lambda_values;
  gsl_vector *through_values;

  gsl_vector **double_vector;


  // allocate space for the return vector
  double_vector  = (gsl_vector **)malloc(2*sizeof (gsl_vector *));

  // allocate the space for the vectors
  lambda_values  = gsl_vector_alloc(21);
  through_values = gsl_vector_alloc(21);

  gsl_vector_set(lambda_values, 0,1.069460);  gsl_vector_set(through_values, 0, 0.000000);
  gsl_vector_set(lambda_values, 1,1.069653);  gsl_vector_set(through_values, 1, 0.012630);
  gsl_vector_set(lambda_values, 2,1.069846);  gsl_vector_set(through_values, 2, 0.025261);
  gsl_vector_set(lambda_values, 3,1.070039);  gsl_vector_set(through_values, 3, 0.037891);
  gsl_vector_set(lambda_values, 4,1.070232);  gsl_vector_set(through_values, 4,   0.050522);
  gsl_vector_set(lambda_values, 5,1.070424);  gsl_vector_set(through_values, 5,   0.063152);
  gsl_vector_set(lambda_values, 6,1.070617);  gsl_vector_set(through_values, 6,   0.069010);
  gsl_vector_set(lambda_values, 7,1.070810);  gsl_vector_set(through_values, 7,   0.069010);
  gsl_vector_set(lambda_values, 8,1.071003);  gsl_vector_set(through_values, 8,   0.069010);
  gsl_vector_set(lambda_values, 9,1.071196);  gsl_vector_set(through_values, 9,   0.069010);
  gsl_vector_set(lambda_values, 10,1.071389); gsl_vector_set(through_values, 10,   0.069010);
  gsl_vector_set(lambda_values, 11,1.071582); gsl_vector_set(through_values, 11,   0.069010);
  gsl_vector_set(lambda_values, 12,1.071775); gsl_vector_set(through_values, 12,   0.069010);
  gsl_vector_set(lambda_values, 13,1.071967); gsl_vector_set(through_values, 13,   0.069010);
  gsl_vector_set(lambda_values, 14,1.072160); gsl_vector_set(through_values, 14,   0.069010);
  gsl_vector_set(lambda_values, 15,1.072353); gsl_vector_set(through_values, 15,   0.063152);
  gsl_vector_set(lambda_values, 16,1.072546); gsl_vector_set(through_values, 16,   0.050522);
  gsl_vector_set(lambda_values, 17,1.072739); gsl_vector_set(through_values, 17,   0.037891);
  gsl_vector_set(lambda_values, 18,1.072932); gsl_vector_set(through_values, 18,   0.025261);
  gsl_vector_set(lambda_values, 19,1.073125); gsl_vector_set(through_values, 19,  0.012630);
  gsl_vector_set(lambda_values, 20,1.073318); gsl_vector_set(through_values, 20,   0.000000);

  //Wavelength: 10713.888672, fringe factor: 1.097945, (x,y): (18, 745)

  
  // build up the output array
  double_vector[0] = lambda_values;
  double_vector[1] = through_values;

  // return the two vectors
  return double_vector;
}

/**
 * Function: evaluate_wavelength_steps
 * The function defines and computes the wavelength steps and the
 * filter throughputs which are used to compute the fringe amplitude
 * for all pixels. The basis for the wavelength data are the
 * filter data and the general are where fringing is significant
 * and computed, which is defined in the fringe
 * configuration structure.
 * First all throughput values below a certain value are discarded
 * at the short and long wavelength edge. Then the final wavelength
 * range is determined using the fringing range. Then all wavelength
 * values and the filter throughputs are determined and returned
 * as two gsl-vectors.
 *
 * Parameters:
 * @param fconf - the fringe configuration structure
 *
 * Returns:
 * @return double_vector - a set of two gsl-vectors with wavlength
                           and throughput
 */
gsl_vector **
evaluate_wavelength_steps(fringe_conf *fconf)
{
  int index=0;
  int lower=0;
  int upper=0;

  int nsteps;

  double lambda_min;
  double lambda_max;
  double lambda_act;
  double through_act=0.0;
  double through_tot=0.0;

  gsl_vector *lambda_values;
  gsl_vector *through_values;

  gsl_vector **double_vector;


  // allocate space for the return vector
  double_vector  = (gsl_vector **)malloc(2*sizeof (gsl_vector *));

  // search through the filter data
  // to narrow the usable range from below
  lower=0;
  while (lower < fconf->filter_through->nvals
	 && fconf->filter_through->yvals[lower] < FILTER_THRESHOLD)
    lower++;
  lower = MAX(0,lower-1);

  // search through the filter data
  // to narrow the usable range from above
  upper=fconf->filter_through->nvals - 1;
  while (upper > -1
	 && fconf->filter_through->yvals[upper] < FILTER_THRESHOLD)
    upper--;
  upper = MIN(fconf->filter_through->nvals - 1, upper+1);

#ifdef DEBUGFCONF
  fprintf(stderr, "New lower range index: %i, new upper: %i, old number: %i\n",
	  lower, upper, fconf->filter_through->nvals);
#endif

  // replace the old filter interpolator
  // with a new one, if necessary
  if (lower != 0 || upper != (fconf->filter_through->nvals - 1))
    fconf->filter_through = redefine_filter_throughput(lower, upper, fconf);

  // compute the lower wavelength range by comparing
  // the filter throughput with the range important
  // for fringing; the latter is given in AA!!
  lambda_min = MAX(gsl_vector_get(fconf->fringe_range, 0),
		   fconf->filter_through->xmin);

  // compute the upper wavelength range by comparing
  // the filter throughput with the range important
  // for fringing; the latter is given in AA!!
  lambda_max = MIN(gsl_vector_get(fconf->fringe_range, 1),
		   fconf->filter_through->xmax);

  // find the number of steps to cover the wavelength range
  nsteps = (int)ceil((lambda_max-lambda_min)/(fconf->fringe_step)) + 1;

  // allocate the space for the vectors
  lambda_values  = gsl_vector_alloc(nsteps);
  through_values = gsl_vector_alloc(nsteps);

  // fill the vectors with the wavelength steps and the 
  // filter throughput at those wavelengths.
  lambda_act = lambda_min;
  for (index=0; index < nsteps-1; index+=1)
    {
      gsl_vector_set(lambda_values, index, lambda_act);
      through_act = eval_interp(fconf->filter_through, lambda_act);
      through_tot += through_act;
      gsl_vector_set(through_values, index, through_act);
      lambda_act += fconf->fringe_step;
    }

  // add also the values at lambda_amax
  gsl_vector_set(lambda_values, nsteps-1, lambda_max);
  through_act = eval_interp(fconf->filter_through, lambda_max);
  through_tot += through_act;
  gsl_vector_set(through_values, nsteps-1, through_act);

 
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
 * Function: redefine_filter_throughput
 * The function re-defines an existing interpolator using only the
 * existing data values in an index range given as parameters.
 *
 * Parameters:
 * @param lower - index with the lowest data point
 * @param upper - index with the highest data point
 *
 * Returns:
 * @return new_interp - the new interpolator
 */
interpolator *
redefine_filter_throughput(const int lower, const int upper,
			   fringe_conf *fconf)
{
  int new_nvals;
  int index;

  double *new_x;
  double *new_y;

  interpolator *new_interp;
 
   // allocate space for the return structure;
  // complain if this fails
  new_interp = (interpolator *)malloc (sizeof (interpolator));
  if (new_interp == NULL)
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Could not allocate memory for interpolator");

  // compute the new length of the data arrays
  new_nvals = upper-lower+1;

  // allocate space for the new data arrays
  new_x = (double *)malloc(new_nvals*sizeof(double));
  new_y = (double *)malloc(new_nvals*sizeof(double));

  // transfer the data
  for (index=lower; index < upper+1; index++)
    {
      new_x[index-lower] = fconf->filter_through->xvals[index];
      new_y[index-lower] = fconf->filter_through->yvals[index];
    }

  // free the old structure
  free_interp(fconf->filter_through);

  // build up the interpolator structure
  new_interp->xmin    = new_x[0];
  new_interp->xmax    = new_x[new_nvals-1];
  new_interp->nvals   = new_nvals;
  new_interp->xvals   = new_x;
  new_interp->yvals   = new_y;
  new_interp->acc     = gsl_interp_accel_alloc();
  new_interp->interp  = gsl_interp_alloc(FILTER_INTERP_TYPE,
					 new_interp->nvals);

  // initialize the cubic spline
  gsl_interp_init(new_interp->interp, new_interp->xvals, new_interp->yvals, 
		  new_interp->nvals);
  
  // return the new structure
  return new_interp;
}


/**
 * Function: eval_linear_interp
 * The function computes and returns the interpolated value
 * at a given position for a linear interpolator.
 *
 * Parameters:
 * @param lin_int - the linear interpolator
 * @param xval    - the position to evaluate the linear interpolator
 *
 * Returns:
 * @return (value) - the interpolated data value
 */
double
eval_linear_interp(linear_interp *lin_int, const double xval)
{

  double factor=0.0;

  // check whether the x-value is within 
  // the range spanned by the data;
  // complain if the x-value is outside
  if (xval < lin_int->xmin || xval > lin_int->xmax)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		 "independent linear interpolation value %f "
		 "is outside interval (%f, %f)\n", xval,
		 lin_int->xmin, lin_int->xmax);

  // check whether you have to search upwards or downwards
  if (xval >= gsl_vector_get(lin_int->xvals, lin_int->act_index))
    {

      // in case that you search upwards, go up
      // the independent values until you find the right interval
      lin_int->act_index += 1;
      while(xval > gsl_vector_get(lin_int->xvals, lin_int->act_index))
	lin_int->act_index++;
    }
  else
    {

      // in case that you search downwards, go down
      // the independent values  until you find the right interval
      //      while(wavelength < resp->spec[nact-1].lambda_mean)
      while(xval < gsl_vector_get(lin_int->xvals, lin_int->act_index-1))
	lin_int->act_index--;
    }

  // interpolate within the interval to calculate the 
  // sensitivity
  factor = (xval - gsl_vector_get(lin_int->xvals, lin_int->act_index-1))/
    (gsl_vector_get(lin_int->xvals, lin_int->act_index) - gsl_vector_get(lin_int->xvals, lin_int->act_index-1));

  // compute and terutn the interpolated value
    return (gsl_vector_get(lin_int->yvals, lin_int->act_index-1)
    + factor * (gsl_vector_get(lin_int->yvals, lin_int->act_index)-gsl_vector_get(lin_int->yvals, lin_int->act_index-1)));
}



/**
 * Function: get_layer_thickness
 * The function determines the thickness of an individual layer
 * at a given pixel position.
 *
 * Parameters:
 * @param opt_layer - the optical layer
 * @param ii        - the first pixel index
 * @param jj        - the second pixel index
 *
 * Returns:
 * @return (value) - the thickness of the layer at the requested position
 */
double
get_layer_thickness(const ccd_layer *opt_layer, const int ii, const int jj)
{
  // check whether the thickness is 2D
  if (opt_layer->thickness2D)
    // return the 2D value
    return gsl_matrix_get(opt_layer->thickness2D, ii, jj);
  else
    // return the constant value
    return opt_layer->thickness;
}


/**
 * Function: get_complex_refindex
 * The function computes and returns the complex refraction index of
 * a given CCD layer at a given wavelength.
 *
 * Parameters:
 * @param opt_layer - the optical layer
 * @param lambda    - the wavelength
 *
 * Returns:
 * @return (value) - the complex refraction index
 */
gsl_complex
get_complex_refindex(const ccd_layer *opt_layer, const double lambda)
{
  double re;
  double im;

  // compute the real part
  re = eval_interp(opt_layer->re_refraction, lambda);

  // comute the imaginary part
  im = eval_interp(opt_layer->im_refraction, lambda);

  // return the complex number
  return gsl_complex_rect(re, im);
}


/**
 * Function: compute_reflection
 * The function computes and returns the reflection index (in %)
 * for two complex refraction indices  at a given wavelength.
 *
 * Parameters:
 * @param refract_l1 - refraction index of one layer 
 * @param refract_l2 - refraction index of the second layer 
 *
 * Returns:
 * @return (value) - the reflection index
 */
double
compute_reflection(const gsl_complex refract_l1,
		   const gsl_complex refract_l2)
{
  gsl_complex conj_l1;
  gsl_complex conj_l2;

  gsl_complex refl_compl;

  // get the comlex conjugated values of the input
  conj_l1 = gsl_complex_conjugate(refract_l1);
  conj_l2 = gsl_complex_conjugate(refract_l2);

  // compute the complex reflection
  refl_compl = gsl_complex_div(gsl_complex_sub(conj_l1, conj_l2),
			       gsl_complex_add(conj_l1, conj_l2));

  // return the magnitude of the complex reflection
  return gsl_complex_abs2(refl_compl);
}


/**
 * Function: compute_transmission
 * The function computes and returns the transmission index (in %)
 * for two complex refraction indices at a given wavelength.
 *
 * Parameters:
 * @param refract_l1 - refraction index of one layer 
 * @param refract_l2 - refraction index of the second layer 
 *
 * Returns:
 * @return (value) - the trnasmission index
 */
double
compute_transmission(const gsl_complex refract_l1,
		     const gsl_complex refract_l2)
{
  gsl_complex conj_l1;
  gsl_complex conj_l2;

  gsl_complex trans_compl;

  // get the comlex conjugated values of the input
  conj_l1 = gsl_complex_conjugate(refract_l1);
  conj_l2 = gsl_complex_conjugate(refract_l2);

  // compute the complex transmission
  trans_compl =
    gsl_complex_div(gsl_complex_sqrt(gsl_complex_mul(conj_l1, conj_l2)),
		    gsl_complex_add(conj_l1, conj_l2));
  
  // finish the computation of the complex transmission
  trans_compl = gsl_complex_mul_real (trans_compl, 2.0);

  // return the magnitude of the complex transmission
  return gsl_complex_abs2(trans_compl);
}


/**
 * Function: compute_attenuation
 * The function computes the attenuation/damping of a plane wave
 * with a given inverse wavelelength (phase number) traversing
 * TWO TIMES a layer of a given thickness and complex refraction
 * index.
 *
 * Parameters:
 * @param refract      - the complex refraction index
 * @param thickness    - the thickness of the layer
 * @param phase_number - the inverse wavelength
 *
 * Returns:
 * @return (value) - the attenuation
 */
double
compute_attenuation(const gsl_complex refract, const double thickness,
		    const double phase_number)
{
  // check if there is something to compute
  if (GSL_IMAG(refract)) 
    // compute and return the value
    return exp(-2.0*GSL_IMAG(refract)*thickness*phase_number);
  else
    // return the value for exp(0.0)
    return 1.0;
}


/**
 * Function: compute_pshift
 * The function computes the phase shift of a plane wave
 * with a given inverse wavelelength (phase number) traversing
 * TWO TIMES a layer of a given thickness and complex refraction
 * index.
 *
 * Parameters:
 * @param refract      - the complex refraction index
 * @param thickness    - the thickness of the layer
 * @param phase_number - the inverse wavelength
 *
 * Returns:
 * @return (value) - the complex phase shift
 */
gsl_complex
compute_pshift(const gsl_complex refract, const double thickness,
	       const double phase_number)
{
  // all in one line:
  return gsl_complex_exp(gsl_complex_rect(0.0,2.0*GSL_REAL(refract)*thickness*phase_number));
}


/**
 * Function: alloc_optprops_list
 * The function allocates space for a list of optical property
 * elements. the size of the list is derived from a fringe
 * configuration structure given in the input.
 *
 * Parameters:
 * @param fconf - the fringe configuration file
 *
 * Returns:
 * @return optprops - the list of optical properties
 */
optical_property *
alloc_optprops_list(const fringe_conf *fconf)
{
  optical_property *optprops;

  // allocate large enough memory 
  optprops =
    (optical_property *)malloc(fconf->opt_layers->num_layers*sizeof(ap_pixel));

  // return the pointer
  return optprops;
}


/**
 * Function: print_optprops_list
 * The function prints the content of an optical property
 * list onto the screen.
 *
 * Parameters:
 * @param optprops - the optical property list
 *
 * Returns:
 * @return -
 */
void
print_optprops_list(const optical_property *optprops, const int num_entries)
{
  int index;
  
  for (index = 0; index < num_entries; index++)
    {
      fprintf(stdout, "#%i: %f %f %f %f %f %f \n",index,
	      optprops[index].trans_upper, optprops[index].reflect_upper,
	      GSL_REAL(optprops[index].double_phshift),
	      GSL_IMAG(optprops[index].double_phshift),
	      optprops[index].trans_lower, optprops[index].reflect_lower);
    }
}

/**
 * Function: free_optprops_list
 * The function releases the space allocated in an
 * optical property list.
 *
 * Parameters:
 * @param optprops - the optical property list
 *
 * Returns:
 * @return -
 */
void
free_optprops_list(optical_property *optprops)
{
  if (optprops != NULL)
    {
      // free the structure
      free(optprops);

      // set it to NULL
      optprops = NULL;
    }
}


/*
 * Function: init_optprops_list
 * The function initializes an optical property list. This means
 * it computes and stores values for the list elements. Here
 * only elements are fixed, which will not change during the
 * whole fringe computation and remain constant for all pixels
 * at all wavelengths.
 *
 * Parameters:
 * @params fconf       - the fringe configuration structure
 * @params lambda_mean - the mean wavelength of analyzed waveband
 * @params optprops    - the list of optical property elements
 *
 * Returns:
 * @return -
 */
void
init_optprops_list(const fringe_conf *fconf, const double lambda_mean,
		   optical_property *optprops)
{
  int num_layers;
  int index;

  double n1;
  double n2;

  // determine the number of layers
  num_layers = fconf->opt_layers->num_layers;

  if (num_layers > 0)
    {
      // set the transmission and reflection
      // vacuum -- first CCD layer
      // that's a bit heuristic, since
      // the true value using the true
      // refration index for vacuum would
      // give a different result
      optprops[0].trans_upper   = 1.0;
      optprops[0].reflect_upper = 0.0;      

      // also that one is a bit heuristic
      // but assuming n(vacuum)=1.0 its also true 
      optprops[0].sign_upper = 1.0;

      // set the sign of the lowest layer -- substrate
      // this is purely heuristic and works 
      // only for the HRC/WFC fringing models!! 
      optprops[num_layers-1].sign_lower = -1.0;
    }

  for (index=0; index < num_layers-1; index++)
    {
      // get the real part of the refraction index
      // of the upper layer
      n1 = GSL_REAL(get_complex_refindex(fconf->opt_layers->opt_layer[index],
					 lambda_mean));

      // get the real part of the refraction index
      // of the lower layer
      n2 = GSL_REAL(get_complex_refindex(fconf->opt_layers->opt_layer[index+1],
					 lambda_mean));

      // distribute the signs
      // for the upper and lower layer
      if (n1 < n2)
	{
	  optprops[index].sign_lower   =  1.0;
	  optprops[index+1].sign_upper = -1.0;
	}
      else
	{
	  optprops[index].sign_lower   = -1.0;
	  optprops[index+1].sign_upper =  1.0;
	}
    }
}


/*
 * Function: fill_optprops_thickness
 * The function fills the thickness of the optival layers
 * at a given pixel into a list of optical properties.
 *
 * Parameters:
 * @params opt_layers - the structure for the optical layers
 * @params ii         - the mean wavelength of analyzed waveband
 * @params jj         - the list of optical property elements
 * @params optprops   - the list of optical property elements
 *
 * Returns:
 * @return -
 */
void
fill_optprops_thickness(const ccd_layers *opt_layers, const int ii,
			const int jj, optical_property *optprops)
{
  int index;

  // go over all layers
  for (index=0; index < opt_layers->num_layers; index++)
    // get and store the thickness
    optprops[index].thickness =
      get_layer_thickness(opt_layers->opt_layer[index], ii, jj);
}


/*
 * Function: fill_optprops_all
 * This function computes and fills in all pixel AND wavelength
 * dependent values into a list of optical properties. The completely
 * filled list can then be used to determine the fringe contribution
 * at that wavelength according to a certain model which defines
 * the various light-paths.
 *
 * Parameters:
 * @params opt_layers - the structure for the optical layers
 * @params lambda     - the actual wavelength
 * @params optprops   - the list of optical property elements
 *
 * Returns:
 * @return -
 */
void
fill_optprops_all(const ccd_layers *opt_layers, const double lambda,
		  optical_property *optprops)
{
  int index;

  double phase_number;

  gsl_complex compl_upp;
  gsl_complex compl_low;

  int num_layers;

  // store the number of layers
  num_layers = opt_layers->num_layers;

  // compute the phase number
  phase_number = 2.0*M_PI / lambda;

  // get the complex refractive index for the first layer
  compl_upp =
    get_complex_refindex(opt_layers->opt_layer[0],lambda);

  for (index=0; index < num_layers-1; index++)
    {

      // get the refractive index for the next layer
      compl_low =
      	get_complex_refindex(opt_layers->opt_layer[index+1],lambda);

      // compute the reflection and transmission
      // to the lower layers
      optprops[index].reflect_lower =
	compute_reflection(compl_upp, compl_low);
      optprops[index].trans_lower   = 1.0 - optprops[index].reflect_lower;

      // compute the reflection and transmission
      // to the upper layers
      optprops[index+1].reflect_upper = optprops[index].reflect_lower;
      optprops[index+1].trans_upper   = 1.0 - optprops[index+1].reflect_upper;

      // compute the attenuation
      optprops[index].double_attenuation =
      	compute_attenuation(compl_upp,optprops[index].thickness,phase_number);

      // compute the phase shift
      optprops[index].double_phshift =
	compute_pshift(compl_upp, optprops[index].thickness, phase_number);


      // prepare the next iteration:
      // the now upper layer will then be the lower
      compl_upp = compl_low;
    }
  // compute the attenuation for the last layer
  optprops[num_layers-1].double_attenuation =
    compute_attenuation(compl_upp,optprops[num_layers-1].thickness,phase_number);
 
  // compute the phase shift for the last layer
  optprops[num_layers-1].double_phshift =
    compute_pshift(compl_upp, optprops[num_layers-1].thickness, phase_number);

  // compute the transmission of the last layer
  // towards the substrate
  optprops[num_layers-1].trans_lower = 
    eval_interp(opt_layers->substrate, lambda);

  // compute the reflection of the last layer
  // towards the substrate
  optprops[num_layers-1].reflect_lower =
    1.0 - optprops[num_layers-1].trans_lower;

}


/*
 * Function: fill_contrib_single
 * The function computes the contribution to the fringe amplitude
 * on the basis of an appropriately filled list of optical properties.
 * In the properties the ingredients such as phase shifts and
 * attenuation and transmissions and reflections are defined.
 * This function just combines the ingredients using a model
 * which includes all single reflected beam combinations
 * as contributors to the fringing.
 *
 * Parameters:
 * @params optprops   - the list of optical property elements
 * @params fconf      - the fringe configuration structure
 *
 * Returns:
 * @return (value)    - the fringe contribution 
 */
double
fringe_contrib_single(const optical_property *optprops,
		      const fringe_conf *fconf)
{
  int index=2;
 
  gsl_complex amp_tmp;
  gsl_complex amp_tot;


  // compute the initial phase at the top of the
  // first layer
  amp_tot =
    gsl_complex_mul_real(gsl_complex_exp(gsl_complex_rect(0.0,fconf->fringe_phase*M_PI/180.0)), 1.0);
  

  // transfer the intitial value 
  // to the running variable
  amp_tmp = amp_tot; 
   
  // pile up the amplitude,
  // also putting in the transmission
  amp_tmp =
    gsl_complex_mul(amp_tmp,
		    gsl_complex_mul_real(optprops[0].double_phshift,
					 1.0*1.0*optprops[0].double_attenuation));

  // compute the intermediate result 
  // for the first layer
  amp_tot = 
    gsl_complex_add(amp_tot, 
		    gsl_complex_mul_real(amp_tmp,
					 optprops[0].sign_lower*optprops[0].reflect_lower));

  // go over layer 2 to the end
  for (index=1; index < fconf->opt_layers->num_layers; index++)
    {
      // pile up the amplitude 
      // to that layer
      amp_tmp =
	gsl_complex_mul(amp_tmp,
			gsl_complex_mul_real(optprops[index].double_phshift,
					   optprops[index-1].trans_lower*optprops[index].trans_upper*optprops[index].double_attenuation));

      // compute the intermediate result 
      // for the that layer
      amp_tot = 
	gsl_complex_add(amp_tot, 
			gsl_complex_mul_real(amp_tmp,
					     optprops[index].sign_lower*optprops[index].reflect_lower));
      //      if (index==3)
      //	fprintf(stdout,"ccc:%f, %f\n", GSL_REAL(amp_tmp),GSL_IMAG(amp_tmp) );

    }

  // return the final value
  return (GSL_REAL(amp_tot) - 1.0);
}
