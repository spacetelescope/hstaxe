/**
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "aXe_grism.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "drz2pet_utils.h"


/**
 * Function: make_spc_drztable
 * The function generates and fills an ap_pixel vector for one
 * drizzled 2D grism image.
 *
 * Parameters:
 * @param  obs     - pointer to observation structure with image values
 * @param  obs     - pointer to observation structure with image values
 * @param  ob      - pointer to the description of the current beam
 * @param  lambda0 - the starting wavelength in this drizzled image
 * @param  dlambda - the dispersion in this drizzled image
 *
 * Returns:
 * @return table   - pointer to an ap_pixel vector
 */
ap_pixel *
make_spc_drztable(observation * const obs, observation * const wobs,
                  const object *ob, const double lambda0,
                  const double dlambda){

  ap_pixel *table, *cur_ap;

  int npix1, npix2;
  int nsub1=1, nsub2=1;
  int i, j;
  int beamInt = 0;

  npix1 = obs->grism->size1;
  npix2 = obs->grism->size2;


  if (!(table = malloc ((npix1 * npix2 + 1) * sizeof (ap_pixel) * nsub1 * nsub2)))
    return NULL;

  cur_ap = table;

  for (i = 0; i < npix1; i++){
    for (j = 0; j < npix2; j++){

      //* check whether the pixel is within the extraction box
      if (fabs((double)j - ob->beams[beamInt].refpoint.y) > ob->beams[beamInt].width + .5)
        continue;

      //* fill the ap_pixel
      cur_ap->p_x  = i;
      cur_ap->p_y  = j;

      cur_ap->x       = (double)i - ob->beams[beamInt].refpoint.x;
      cur_ap->y       = (double)j - ob->beams[beamInt].refpoint.y;
      cur_ap->dist    = (double)j - ob->beams[beamInt].refpoint.y;
      cur_ap->xs      = (double)i - ob->beams[beamInt].refpoint.x;
      cur_ap->ys      = 0.0;
      cur_ap->dxs     = 0.0;  // be reminded: dxs is the LOCAL TRACEANGLE, not the pixels size..
      cur_ap->xi      = (double)i - ob->beams[beamInt].refpoint.x;
      cur_ap->lambda  = lambda0 + dlambda * ((double)i - ob->beams[beamInt].refpoint.x);
      cur_ap->dlambda = dlambda;

      // transfer the cps values and the weight values
      cur_ap->count  = gsl_matrix_get(obs->grism, i, j);
      cur_ap->weight = gsl_matrix_get(wobs->grism, i, j);

      cur_ap->error  = gsl_matrix_get(obs->pixerrs, i, j);
      cur_ap->contam = gsl_matrix_get(obs->dq, i, j);

      // dq is always set to zero, bad (=unexposed) pixels
      // have exptime=weight=0.0
      cur_ap->dq     = 0;

      // go to the next
      // PET pixel
      cur_ap++;
    }
  }

  cur_ap->p_x = -1;
  cur_ap->p_y = -1;
  cur_ap->count = -1;

  return table;
}

/**
 * Function: get_ID_num
 * The function generates the extension name used in the PET
 * for a 2D drizzled grism image
 *
 * Parameters:
 * @param  grism_image - the name of the 2D drizzled grism image
 * @param  extname     - The ID number of the image
 *
 * Returns:
 * @return WorkPtr     - Pointer to the ID-string
 */
char
*get_ID_num (char *grism_image, char *extname){

  char *WorkPtr;
  char *WorkPtr2;
  char **t_err=NULL;
  char tmp[MAXCHAR];
  char ID[4];

  int len;
  int i=0;

  strcpy(ID,"_ID");
  strcpy(tmp, grism_image);

  t_err = (char **) malloc(sizeof(char *)*1);

  WorkPtr = strstr(tmp, ID)+3;
  WorkPtr2 = WorkPtr + strlen (WorkPtr);
  while (isspace ((int) *WorkPtr2))
    *WorkPtr2-- = '\0';
  WorkPtr2 = WorkPtr2 - 1;
  len = strlen(WorkPtr);
  while  (i++ < len){
    strtol(WorkPtr,t_err,10);
    if (strlen(t_err[0])==0)
      break;
    *WorkPtr2-- = '\0';
  }
  strcpy(extname,WorkPtr);
  free(t_err);
  return WorkPtr;
}

/**
 * Function: normalize_weight
 * The function computes the weights used in the extraction
 * from the exposure times given for each pixel in a column.
 *
 * Parameters:
 * @param  wobs - pointer to the structure which carries the pixel data
 * @param  ob   - pointer to the beam description
 */
void
normalize_weight(observation *wobs, object *ob,
                 const int opt_extr)
{
  gsl_matrix *weight;
  //int i, j;
 // // check for optimal extraction
  if (opt_extr)
    {
      // let the optimal weights compute
      fprintf(stdout, "aXe_DRZ2PET: Use otpimal weighting.\n");
      weight = comp_opt_weight(wobs->pixerrs, wobs->dq, ob);
    }
  else
    {
      // let the exposure time weights compute
      fprintf(stdout, "aXe_DRZ2PET: Use normal, flat  weighting.\n");
      weight = comp_equ_weight(wobs->grism, ob);

      // code lines for exposure time weighting
      //fprintf(stdout, "use exposure time  weighting\n");
      //weight = comp_exp_weight(wobs->grism, ob);
    }

  // replace the weight extension with the new weights
  gsl_matrix_free(wobs->grism);
  wobs->grism = weight;
}

/*
 * Function: comp_equ_weight
 * The function creates and fills a weight map with flat entries.
 * Pixels with non-zero exposure time which which contribute to the
 * extracted spectrum all receive the weight 1.0, all other
 * get the weight 0.0.
 *
 * Parameters:
 * @param exp_map - the matrix with the exposure times
 * @param ob      - the object structure
 *
 * Returns:
 * @return weight - the matrix with the exposure time weights
 */
gsl_matrix *
comp_equ_weight(gsl_matrix *exp_map, const object *ob)
{
  gsl_matrix *weight;
  int i, j;
  int beamInt = 0;

  // allocate the matrix and set it to the default
  weight = gsl_matrix_alloc(exp_map->size1, exp_map->size2);
  gsl_matrix_set_all(weight, 0.0);

   // go over all columns
  for (i=0; i < (int)exp_map->size1; i++)
    // go over all rows
    for (j=0; j < (int)exp_map->size2; j++)
      // check whether the pixel is inside the extraction
      // area and has a positive exposure time
      if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5 &&
          gsl_matrix_get(exp_map, i, j) > 0.0 )
        // set the weight to one
        gsl_matrix_set(weight, i,j,1.0);

  // return the resulting weight map
  return weight;
}


/*
 * Function: comp_exp_weight
 * The function uses the exposure time map and computes
 * the weights which are proportional to the relative
 * exposure time of the pixel. Pixels outside of the extraction
 * area are set to 1000.0 Inside the default is 1.0 in unexposed
 * columns and 0.0 within exposed columns.
 *
 * Parameters:
 * @param exp_map - the matrix with the exposure times
 * @param ob      - the object structure
 *
 * Returns:
 * @return weight - the matrix with the exposure time weights
 */
gsl_matrix *
comp_exp_weight(gsl_matrix *exp_map, const object *ob)
{
  gsl_matrix *weight;
  int i, j;
  double sum, contr, norm, allweight;
  int beamInt = 0;

  // allocate the matrix and set it to the default
  weight = gsl_matrix_alloc(exp_map->size1, exp_map->size2);
  gsl_matrix_set_all(weight, 1000.0);

    //* go over all columns
  for (i=0; i < (int)exp_map->size1; i++)
    {
      sum = 0.0;
      contr = 0.0;
      allweight = 0.0;

      // determine for each column the total exposure time
      // and the number of pixels with non-zero exposure time
      for (j=0; j < (int)exp_map->size2; j++)
        {
          // check whether the pixel is inside the extraction
          // area and has a positive weight
          if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5 &&
              gsl_matrix_get(exp_map, i, j) > 0.0 )
            {
              // sum up the weigth and the number of pixels
              sum = sum + gsl_matrix_get(exp_map, i, j);
              contr = contr + 1.0;
            }
        }

      // check whether the row has exposure time at all
      if (sum > 0.0)
        {
          // if the row is exposed:
          // determine the mean exposure time
          norm = sum / contr;

          // go over each pixel
          for (j=0; j < (int)exp_map->size2; j++)
            {
              // check whether the pixel is inside the extraction are
              if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5)
                {
                  // compute and store the relative exposure time (==weight)
                  // for positive exposure time
                  if (gsl_matrix_get(exp_map, i, j) > 0.0)
                    {
                      gsl_matrix_set(weight, i,j,gsl_matrix_get(exp_map, i, j)/norm);
                      allweight = allweight + gsl_matrix_get(weight, i, j);
                    }
                  // make a default to avoid values <0.0
                  else
                    {
                      gsl_matrix_set(weight, i,j,0.0);
                    }
                }
            }
        }
      // treat unexposed columns differently:
      else
        {

          // go over each row
          for (j=0; j < (int)exp_map->size2; j++)
            {
              // check whether the pixel is inside the extraction are
              if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5)
                {
                  // set inside pixels to the defaulkt 1.0
                  gsl_matrix_set(weight, i,j,1.0);
                  allweight = allweight + gsl_matrix_get(weight, i, j);
                }
            }
        }
    }

  // return the resulting weight map
  return weight;
}

/*
 * Function: comp_opt_weight
 * The function uses the model map and and the variance
 * map and computes optimale weigths according to the
 * Hoorne method.
 * Pixels outside of the extraction
 * area are set to 1000.0
 *
 * Parameters:
 * @param mod_map - the matrix with the model values
 * @param var_map - the matrix with the variance values
 * @param ob      - the object structure
 *
 * Returns:
 * @return weight - the matrix with the exposure time weights
 */
gsl_matrix *
comp_opt_weight(gsl_matrix *mod_map,
                gsl_matrix *var_map, const object *ob)
{
  gsl_matrix *weight;
  int i, j;
  double mod_sum, weight_sum, contr, norm, allweight;
  int beamInt = 0;
  double mod_val;
  //double var_val;

  // allocate the weight matrix and set the default
  weight = gsl_matrix_alloc(mod_map->size1, mod_map->size2);
  gsl_matrix_set_all(weight, 0.0);

  //* go over all columns
  for (i=0; i < (int)mod_map->size1; i++)
    {
      mod_sum = 0.0;
      contr = 0.0;
      allweight = 0.0;
      weight_sum=0.0;

      // determine for each column the total model counts
      // and the number of pixels with non-zero model_counts
      for (j=0; j < (int)mod_map->size2; j++)
        {
          if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5 &&
              gsl_matrix_get(mod_map, i, j) > 0.0 )
            {
              mod_sum = mod_sum + gsl_matrix_get(mod_map, i, j);
              contr = contr + 1.0;
            }
        }

      // check whether the column has model values.
      // normalize the model values and compute
      // optimal weights if yes.
      if (mod_sum > 0.0)
        {
          //* determine the mean model value
          //      norm = mod_sum / contr;
          norm = mod_sum;

          // go over each row
          for (j=0; j < (int)mod_map->size2; j++)
            {
              // check whether the pixel is inside the extraction area
              if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5 &&
                  gsl_matrix_get(mod_map, i, j) > 0.0 )
                {
                  // normalize the model counts
                  mod_val = gsl_matrix_get(mod_map, i, j)/norm;

                  // store the normalized model counts
                  gsl_matrix_set(mod_map, i, j, mod_val);

                  // add up the normalization value for the optimal weights
                  weight_sum = weight_sum + mod_val*mod_val/gsl_matrix_get(var_map, i, j);
                }
            }

          // finally compute and write the weights:
          // go over each pixel
          for (j=0; j < (int)mod_map->size2; j++)
            {
              // check whether the pixel is inside the extraction area
              if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5)
                {
                  if (gsl_matrix_get(mod_map, i, j) > 0.0 )
                    {
                      // compute and set the individual pixel weight
                      gsl_matrix_set(weight, i, j, gsl_matrix_get(mod_map, i, j)/gsl_matrix_get(var_map, i, j)/weight_sum);
                      allweight = allweight + gsl_matrix_get(weight, i, j);
                    }
                  else
                    {
                      // set the default value if no model value is zero
                      gsl_matrix_set(weight, i, j,0.0);
                    }
                }
            }
        }
      // if the column does not have model values at all:
      else
        {
          // go over each row
          for (j=0; j < (int)mod_map->size2; j++)
            {
              // check whether the pixel is within the extraction area
              if (fabs(ob->beams[beamInt].refpoint.y-(double)j) <= ob->beams[beamInt].width+0.5)
                {
                  // set the inside default value
                  gsl_matrix_set(weight, i,j,1.0);
                  allweight = allweight + gsl_matrix_get(mod_map, i, j);
                }
            }
        }
    }

  // return the weight matrix
  return weight;
}
