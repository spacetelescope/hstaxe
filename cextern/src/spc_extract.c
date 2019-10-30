/**
 * The driver function and some auxillaries
 * to compute the spectrum from an image and
 * information on the location and distortion of the spectra.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_sect.h"
#include "spce_is_in.h"
#include "spce_pathlength.h"

#define DEBUG_ME 0x10


/**
  fill in the path length field (xi in ap_pixel) from the abscissas (xs in
  ap_pixel) (helper function for spc_extract).

  @param func the spectrum trace to use for the transformation from xs to xi
  @param table a pointer to the table of aperture pixels, terminated with
  x=-1
*/
static int
transform_to_pathlen (const trace_func * const func, ap_pixel * const table)
{
  int i, len;
  ap_pixel *cur_p = table;
  gsl_vector *section_points;

  while (cur_p->p_x != -1)
    cur_p++;
  len = cur_p - table;
  if (len==0) return -1; // exit now if the table is empty

  section_points = gsl_vector_alloc (len);
  for (i = 0; i < len; i++)
    {
      gsl_vector_set (section_points, i, table[i].xs);
    }

  if (abscissa_to_pathlength (func, section_points))
    {
      gsl_vector_free (section_points);
      return -1;
    }

  for (i = 0; i < len; i++)
    {
      table[i].xi = gsl_vector_get (section_points, i);
    }

  gsl_vector_free (section_points);

  return 0;
}



/**
   creates an ap_pixel.

   @param x x coordinate of pixel relative to beam's reference point
   @param y y coordinate of pixel relative to beam's reference point
   @param px absolute x coordinate of the pixel
   @param py absolute y coordinate of the pixel
   @param ob the observation to get the pixels from
   @param sf pointer to sectionfun structure for this beam
   @param cur_ap a pointer to the next free aperture pixel
   @return a pointer to the next free aperture pixel after the new pixels
    have been added.

*/
static ap_pixel *
handle_one_pixel (const double x, const double y, const int px, const int py,
		  const observation * const obs, sectionfun * const sf,
		  const trace_func * const tracefun, ap_pixel * cur_ap)
{
  double res;
  double sect_y;
  double tmp;
  double phi_trace;

  if (find_section_point (sf, x, y, &res))
    {
      return cur_ap;	/* FIXME: Issue warning here */
    }
  /* set the extraction weight to 1. */
  cur_ap->weight = 1.;
  cur_ap->xs = res;
  phi_trace = atan (tracefun->deriv (cur_ap->xs, tracefun->data));
  cur_ap->dxs = phi_trace;
  cur_ap->ys = tracefun->func (cur_ap->xs, tracefun->data);
  sect_y = tracefun->func (res, tracefun->data);
  tmp = tracefun->func (x, tracefun->data);

  cur_ap->contam = -1.0;
  cur_ap->model = 0.0;
  cur_ap->p_x = px;
  cur_ap->p_y = py;
  cur_ap->x = x;
  cur_ap->y = y;
  cur_ap->dist =
    sqrt ((sect_y - y) * (sect_y - y) +
	  (cur_ap->xs - x) * (cur_ap->xs - x));

  if ( y < tmp ) {
    cur_ap->dist = cur_ap->dist * (-1.); /* If pixel is bwlow the trace, dist is neg. */
  }
  cur_ap->count = gsl_matrix_get (obs->grism, px, py);
  cur_ap->error = gsl_matrix_get (obs->pixerrs, px, py);
  if (obs->dq != NULL) cur_ap->dq = (long) gsl_matrix_get (obs->dq, px, py);
  else cur_ap->dq = 0;

  cur_ap++;
  return cur_ap;
}


/**
   Does some sanity checks on make_spc_table's input.

   @param ob the object to check
   @return 0 if everything is ok, -1 otherwise
*/
static int
sanitycheck (object * const ob)
{
  int i;

  for (i = 0; i < ob->nbeams; i++)
    {
      while (ob->beams[i].orient < 0)
	{
	  ob->beams[i].orient += M_PI;
	}
      while (ob->beams[i].orient > M_PI)
	{
	  ob->beams[i].orient -= M_PI;
	}
    }
  return 0;
}


#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
/**
   computes a table of aperture pixels, i.e. of tuples containing the
   source coordinates, the distance to the spectrum trace, the path
   length along the trace, and the intensity.  The end of the table
   is marked with an ap_pixel->x==-1.  This routine can do subsampling,
   whereby each pixel is divided into n_sub*n_sub smaller pixels.
   If n_sub==1, a special (faster) handling is enabled.

   @param ob the object struct to examine
   @param beamorder the order of the spectrum to examine
   @param flags warning flags like in the warning field of the
     spectrum structure.  This has to be passed in initialized.
*/
ap_pixel *
make_spc_table (object * const ob, const int beamorder,
		int *const flags)
{
  int bb_x, bb_y, bb_w, bb_h;
  int x, y;
  beam *curbeam = ob->beams + beamorder;
  double dx, dy;
  ap_pixel *table, *cur_ap;
  is_in_descriptor iid;
  sectionfun sf;
  trace_func *tracefun = curbeam->spec_trace;

  quad_to_bbox (curbeam->corners, curbeam->bbox, curbeam->bbox + 1);
  bb_x = curbeam->bbox[0].x;
  bb_y = curbeam->bbox[0].y;
  bb_w = curbeam->bbox[1].x - curbeam->bbox[0].x + 1;
  bb_h = curbeam->bbox[1].y - curbeam->bbox[0].y + 1;

  if (sanitycheck (ob))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Input data failed sanity check");
      return NULL;
    }

  if (fill_is_in_descriptor (&iid, curbeam->corners))
    return NULL;

  if (fill_in_sectionfun (&sf, curbeam->orient, curbeam))
    return NULL;

  if (!(table =
	malloc ((bb_w * bb_h + 1) * sizeof (ap_pixel))))
    return NULL;


  /* We have a little coordinate confusion here.  There are three
     systems:
     (a) the absolute one, rooted in (0,0) of the image.
     (b) the relative one, rooted in the reference point and used
     when evaluating the spectrum trace.  This one is used
     for finding the section points, etc.
     (c) one relative to the bounding box; x and y below live in this
     system and have to be converted into whatever system is
     required.  dx, dy are used for the conversion to (b)
  */
  dx = bb_x - curbeam->refpoint.x;
  dy = bb_y - curbeam->refpoint.y;
  cur_ap = table;
  for (y = 0; y <= bb_h; y++)
    {
      for (x = 0; x < bb_w; x++)
	{
	  if ((bb_x + x < 0) || (bb_y + y < 0)
	      || (bb_x + x >= (int)ob->grism_obs->grism->size1)
	      || (bb_y + y >= (int)ob->grism_obs->grism->size2))
	    {
	      *flags |= SPC_W_BORDER;
	      continue;
	    }


	  // check whether the trace description has
	  // an order higher than linear
	  if (curbeam->spec_trace->type > 1)
	    {
	      // new criteria based on the true trace distance
	      // which means the true distance from the section point
	      if (!tracedist_criteria(x + dx, y + dy, &sf, tracefun, curbeam->width+2.0))
		{
		  continue;
		}
	    }
	  else
	    {
	      // old criteria based on the box model
	      if (!is_in (bb_x + x, bb_y + y, &iid))
		{
		  continue;
		}
	    }

	  if (isnan
	      (gsl_matrix_get
	       (ob->grism_obs->grism, bb_x + x, bb_y + y)))
	    {
	      continue;
	    }
	  //	  if (ob->ID == 11 && x + bb_x == 59)
	    //	    fprintf(stdout, "xx: %i, yy: %i: %i\n", x + bb_x, y + bb_y, is_in (x + bb_x, y + bb_y, &iid));
	  cur_ap =
	    handle_one_pixel (x + dx, y + dy, x + bb_x,
			      y + bb_y, ob->grism_obs, &sf,
			      tracefun, cur_ap);

	}
    }

  cur_ap->p_x = -1;
  cur_ap->p_y = -1;
  cur_ap->count = -1;

  free_sectionfun (&sf);


  if (transform_to_pathlen (curbeam->spec_trace, table))
    {
      //   free (table);
      //   table = NULL;
      //   return NULL;
    }
  return table;
}

/**
   Similar to make_spc_table, it just works only
   on one spot given in the parameters.
*/
ap_pixel *
make_gps_table (object * const ob, const int beamorder,
		int *const flags, int xval, int yval)
{
  int bb_x, bb_y, bb_w, bb_h;
  int x, y;
  beam *curbeam = ob->beams + beamorder;
  double dx, dy;
  ap_pixel *table, *cur_ap;
  is_in_descriptor iid;
  sectionfun sf;
  trace_func *tracefun = curbeam->spec_trace;

  quad_to_bbox (curbeam->corners, curbeam->bbox, curbeam->bbox + 1);
  bb_x = curbeam->bbox[0].x;
  bb_y = curbeam->bbox[0].y;
  bb_w = curbeam->bbox[1].x - curbeam->bbox[0].x + 1;
  bb_h = curbeam->bbox[1].y - curbeam->bbox[0].y + 1;

  if (sanitycheck (ob))
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "Input data failed sanity check");
      return NULL;
    }

  if (fill_is_in_descriptor (&iid, curbeam->corners))
    return NULL;

  if (fill_in_sectionfun (&sf, curbeam->orient, curbeam))
    return NULL;

  if (!
      (table =
       malloc (2 * sizeof (ap_pixel))))
    return NULL;

  /* We have a little coordinate confusion here.  There are three
     systems:
     (a) the absolute one, rooted in (0,0) of the image.
     (b) the relative one, rooted in the reference point and used
     when evaluating the spectrum trace.  This one is used
     for finding the section points, etc.
     (c) one relative to the bounding box; x and y below live in this
     system and have to be converted into whatever system is
     required.  dx, dy are used for the conversion to (b)
  */
  dx = bb_x - curbeam->refpoint.x;
  dy = bb_y - curbeam->refpoint.y;
  x  = xval - bb_x;
  y  = yval - bb_y;

  cur_ap = table;

  cur_ap =
    handle_one_pixel (x + dx, y + dy, x + bb_x,
		      y + bb_y, ob->grism_obs, &sf,
		      tracefun, cur_ap);

  cur_ap->p_x = -1;
  cur_ap->p_y = -1;
  cur_ap->count = -1;

  free_sectionfun (&sf);

  transform_to_pathlen (curbeam->spec_trace, table);

  return table;
}



void
print_ap_pixel_table (const ap_pixel * ap_p)
{
  printf ("# x y pathlen distance lambda count error\n");
  while (ap_p->p_x != -1)
    {
      printf ("%d %d %f %f %f %f %f %f %f %ld\n", ap_p->p_x, ap_p->p_y,
	      ap_p->x, ap_p->y, ap_p->xi, ap_p->dist, ap_p->lambda,
	      ap_p->count, ap_p->error, ap_p->dq);
      ap_p++;
    }
}

/*
int
tracedist_criteria(const double x, const double y, sectionfun * const sf,
		   const trace_func *tracefun, const double width)
{
  double x_sect, y_sect;
  double dist, max_dist;

  int ireturn=0;


  find_section_point (sf, x, y, &x_sect);
  y_sect = tracefun->func (x_sect, tracefun->data);

  dist =  sqrt((y_sect-y)*(y_sect-y)+(x_sect-x)*(x_sect-x));

  if (dist < width)
    ireturn=1;

  return ireturn;
}
*/
