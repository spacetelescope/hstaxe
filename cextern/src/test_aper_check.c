//Verify that aper_check is functioning correctly with basic data
#include "aper_check.h"
#include <stdio.h>


int
main (void)
{
  int i, j;
  observation *obs;
  aXe_mask *mask;
  ap_pixel *pixel;

  obs = malloc (sizeof(observation));

  //Make 10x10 unitary matrices of fake grism data
  obs->grism = gsl_matrix_alloc (10, 10);
  obs->pixerrs = gsl_matrix_alloc (10, 10);
  obs->dq = gsl_matrix_alloc (10, 10);

  //populate the images
  for (i = 0; i < 10; i++)
    for (j = 0; j < 10; j++)
      gsl_matrix_set (obs->grism, i, j, 1.);
      gsl_matrix_set (obs->pixerrs, i, j, 0.1*i*j);
      gsl_matrix_set (obs->dq, i, j, 0.);

  //set a couple bad pixels in the mask
  gsl_matrix_set (obs->dq, 1, 1, 1);
  gsl_matrix_set (obs->dq, 3, 4, 16);

  //make an aXe_mask for the observation
  mask= malloc (sizeof(aXe_mask));
  mask = aXe_mask_init (obs);

  //setup a pixel object
  pixel = malloc (sizeof(ap_pixel));
  pixel->p_x = 5;
  pixel->p_y = 5;
  pixel->x = 6;
  pixel->y = 6;
  pixel->dist = 0.5;
  pixel->xs = 0.5;
  pixel->ys = 0.5;
  pixel->dxs = 0.1;
  pixel->xi = 0.25;
  pixel->lambda = 0.25;
  pixel->dlambda = 0.25;
  pixel->count = 100;
  pixel->weight = 1.0;
  pixel->error=1.0;
  pixel->contam = 1.0;
  pixel->model = 99.;
  pixel->dq=0;

  add_ap_p_to_aXe_mask (pixel, mask);
  printf("add_ap_p_to_aXe_mask successfull\n");

  mark_trace_in_aXe_mask(pixel, mask);
  printf("mark_trace_in_aXe_mask successfull\n");

  gsl_matrix_free (obs->grism);
  gsl_matrix_free (obs->pixerrs);
  gsl_matrix_free (obs->dq);
  free(obs);
  gsl_matrix_free (mask->img);
  free(mask);
  free(pixel);

  return 0;
}
