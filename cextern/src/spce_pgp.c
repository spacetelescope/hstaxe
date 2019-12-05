/**
 * File: spce_pgp.c
 * Various PGplot routines amd helper functions.
 */

#include "spce_pgp.h"

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define NCOL 64

#ifdef HAVE_PGPLOT

#include <cpgplot.h>


/**
    A helper function to convert a gsl_matrix into a linear
    integer array suitble to plot with pgplot's routines

    @param ar 2D gal_matrix
    @param r a pointer to an already allocated amount of memeory
*/
void
gsl_to_pgarr (gsl_matrix * ar, int *r)
{
     int i, j, k;

     for (i = 0; i < ar->size1; i++)
       {
	    for (j = 0; j < ar->size2; j++)
	      {
		   k = j * ar->size1 + i;
		   r[k] = (int) gsl_matrix_float_get (ar, i, j);
	      }
       }
}


/**
    Outputs a stamp image using pgplot and using the information for an
    ap_pixel structure.

    @param ap_p an existing ap_pixel structure
    @param n_sub subsampling used (usualy 1.)
    @param ob a pointer to the object correspoding to ap_p
    @param beamnum the numeric ID of the beam corresponding to ap_p
    @param filename a pointer to a char array containing the name of the file to output the
    postscript stamp image to

*/
void
pgp_stamp_image (const ap_pixel * const ap_p, const int n_sub,
		 const object * ob, int beamnum, char filename[],
		 const int negative)
{
     int ncol = NCOL;
     double r, g, b;
     beam *curbeam = ob->beams + beamnum;
     gsl_matrix *fpix;
     const ap_pixel *cur_p;
     double pxmin = ap_p->count, pxmax = ap_p->count;
     int i;
     int n;
     int m;
     float *x, *y;
     int p_x, p_y;
     double minx = ap_p->x, miny = ap_p->y;
     //gsl_vector_int * trace_inds = get_trace_inds(ap_p);
     gsl_vector_int *trace_inds;
     int *k;
     int ind;
     int icilo, icihi;
     char label[255];


     quad_to_bbox (curbeam->corners, curbeam->bbox, curbeam->bbox + 1);
     n = (curbeam->bbox[1].x - curbeam->bbox[0].x + 1) * n_sub;
     m = (curbeam->bbox[1].y - curbeam->bbox[0].y + 1) * n_sub;

     fpix =
	  gsl_matrix_alloc ((curbeam->bbox[1].x - curbeam->bbox[0].x +
			     1) * n_sub,
			    (curbeam->bbox[1].y - curbeam->bbox[0].y +
			     1) * n_sub);

     trace_inds = (gsl_vector_int *) get_trace_inds (ap_p);

     x = malloc (trace_inds->size * sizeof (float));
     y = malloc (trace_inds->size * sizeof (float));
     if ((x == NULL) || (y == NULL))
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "pgp_stamp_image: Could not allocate temporary memory space");
	    return;
       }
     if (!fpix)
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "Could not allocate GSL array");
	    return;
       }
     sprintf (label, "Object %d Beam %c at (%4.2f,%4.2f)", ob->ID,
	      BEAM (beamnum), curbeam->refpoint.x, curbeam->refpoint.y);

     for (i = 0; i < trace_inds->size; i++)
       {
	    if (gsl_vector_int_get (trace_inds, i) != -1)
	      {
		   break;
	      }
       }
     if (i == trace_inds->size)
       {
	    aXe_message (aXe_M_ERROR, __FILE__, __LINE__,
			 "No valid trace found");
	    return;
       }

     for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
       {
	    if (!isnan (cur_p->count))
	      {
		   pxmin = MIN (pxmin, cur_p->count);
		   pxmax = MAX (pxmax, cur_p->count);
	      }
	    minx = MIN (minx, cur_p->x);
	    miny = MIN (miny, cur_p->y);
       }
     //cpgopen("plot.ps/CPS");
     cpgopen (filename);

     cpgpage ();
     cpgscr (0, 0, 0.3, 0.2);
     cpgsvp (0.05, 0.95, 0.05, 0.95);
     cpgwnad (0., n, 0., m);
     cpgqcir (&icilo, &icihi);

     for (i = 0; i < trace_inds->size; i++)
       {
	    ind = gsl_vector_int_get (trace_inds, i);
	    x[i] = floor (ap_p[ind].x - minx) + .5;
	    y[i] = floor (ap_p[ind].y - miny) + .5;
       }

     ncol = icihi - icilo;
     /* Create a vector fot pgplot to plot */
     k = malloc (n * m * sizeof (int));
     if (k == NULL)
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Out of memory.");
       }
     gsl_matrix_set_all (fpix, negative ? 0 : pxmax);

     for (cur_p = ap_p; cur_p->p_x != -1; cur_p++)
       {
	    p_x = (int) floor ((cur_p->x - minx) * n_sub);
	    p_y = (int) floor ((cur_p->y - miny) * n_sub);
	    gsl_matrix_set (fpix, p_x, p_y,
			    (cur_p->count - pxmin) * ncol / (pxmax - pxmin) +
			    icilo + 1);
       }

     gsl_to_pgarr (fpix, k);

     /* Set up a uniform BW color table */
     for (i = 0; i < ncol; i++)
       {
	    r = (i - 1.) / (ncol - 1.) * 0.8 + 0.2;
	    g = MAX (0.0, 2. * (i - 1. - ncol / 2.0) / (ncol - 1.));
	    b = 0.2 + 0.4 * (ncol - i) / ncol;
	    cpgscr (i + 15, r, g, b);
	    //cpgscr(i + 15, 1 - i * 1. / ncol, 1 - i * 1. / ncol, 1 - i * 1. / ncol);
	    //cpgscr(i + 15, i * 1. / ncol, i * 1. / ncol, i * 1. / ncol);
       }
     //cpgpixl(k, 10, 10, 1, 10, 1, 10, 0.0, 10, 0., 10);

     cpgpixl (k, n, m, 1, n, 1, m, 0.0, n, 0., m);
     cpgline (trace_inds->size, x, y);

     cpgsci (1);
     cpgmtxt ("t", 1., 0.0, 0., label);
     cpgbox ("bcnts", 0., 0, "bcnts", 0., 0);
     cpgend ();

     free (k);
     k = NULL;
     free (x);
     x = NULL;
     free (y);
     y = NULL;
}
#endif
#ifndef HAVE_PGPLOT
void
pgp_stamp_image (const ap_pixel * const ap_p, const int n_sub,
		 const object * ob, int beamnum, char filename[],
		 const int negative)
{
     aXe_message (aXe_M_WARN4, __FILE__, __LINE__,
		  "pgp_stamp_image: PGplot support was not compiled in at compile time. ");
     return;

}
#endif
