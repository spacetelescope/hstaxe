/**
 * APER_CHECK - create a FITS image showing the locations of
 *               apertures.
 */
#include "aper_check.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "spce_output.h"


void mark_trace_in_aXe_mask(ap_pixel * ap_p, aXe_mask *mask)
{
  gsl_vector_int *trace_inds;
  int i, ind;

  trace_inds = (gsl_vector_int *) get_trace_inds (ap_p);


  for (i = 0; i < (int)trace_inds->size; i++)
    {
      ind = gsl_vector_int_get (trace_inds, i);
      fprintf(stderr,"%d %d\n",ap_p[ind].p_x,ap_p[ind].p_y);
    }

  gsl_vector_int_free(trace_inds);
}

/**
   Function to create and return a gsl array of the same size as the grism
   part of n observation structure.

   @param ob a pointer to an observation structure
   @return a pointer to a newly allocated aXe_mask structure

*/
aXe_mask *
aXe_mask_init (observation * ob)
{
  aXe_mask *mask;
  long i,j;
  double v;

     /* Allocated space for the mask array */
  mask = (aXe_mask *) malloc (sizeof (aXe_mask));
  if (!mask)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory.");
    }

  /* Create the GSL image array */
  mask->img = gsl_matrix_alloc (ob->grism->size1, ob->grism->size2);
  //  fprintf(stderr,"%ld %ld\n",ob->grism->size1,ob->grism->size2);

  /* Initialise the mask array */
  for (i=0;i<(int)ob->grism->size1;i++) {
    for (j=0;j<(int)ob->grism->size2;j++) {
      v = gsl_matrix_get(ob->grism,i,j);
      gsl_matrix_set(mask->img, i,j, v);

    }
  }

  //gsl_matrix_set_all (mask->img, 1.0);

  return mask;
}

/**
 * This function sets all the pixels listed in an ap_pixel list
 * in an aXe_mask to 0.0

    @param a pointer to a ap_pixel array
    @param mask an existing aXe_mask structure

*/
void
add_ap_p_to_aXe_mask (ap_pixel * ap_p, aXe_mask * mask)
{
  ap_pixel *cur_p;

  if (ap_p==NULL) return;

  for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
    {

      /* Loop over pixels which are valid */
      if (!isnan (cur_p->count))
        {
          //   fprintf(stderr,"%d %d\n",cur_p->p_x, cur_p->p_y);

          gsl_matrix_set (mask->img, cur_p->p_x, cur_p->p_y, 0.0);
        }
    }
}
