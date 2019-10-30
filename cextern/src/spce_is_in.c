/**
 * A function to decide whether a point is in a given quadrangle
 * (and its helpers)
 *
 * Usage:
 *   is_in_descriptor iid;
 *   fill_is_in_descriptor(&iid, corners);
 *   for (i=0; i<30; i++) if (is_in(i,i,&iid)) printf("Yes");
 */

#include <math.h>
#include "aXe_grism.h"
#include "aXe_utils.h"
#include "spce_sect.h"
#include "spce_is_in.h"


#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


/**
 * Fills a triangle structure and makes sure the angle at points[0] 
 * is <=90 deg.
 *
 * @param trig the struct to fill.
 * @param points the three point at the corners of the triangle.
 */
void
fill_triag_struct (triangle * const trig, px_point points[3])
{
  trig->root = points[0];
  trig->s.x = points[1].x - points[0].x;
  trig->s.y = points[1].y - points[0].y;
  trig->t.x = points[2].x - points[0].x;
  trig->t.y = points[2].y - points[0].y;
  trig->sabs = sqrt (SQR (trig->s.x) + SQR (trig->s.y));
  trig->tabs = sqrt (SQR (trig->t.x) + SQR (trig->t.y));
  trig->s_t = trig->s.x * trig->t.x + trig->s.y * trig->t.y;
  
  if (trig->s_t / trig->sabs / trig->tabs < 0)
    {
      px_point tmp;
      tmp = points[0];
      points[0] = points[1];
      points[1] = tmp;
      fill_triag_struct (trig, points);
    }
}


/**
 * pre-computes some constant items for deciding whether a point is in
 * a given convex quadrangle. 
 * 
 * @param iid the descriptor to fill
 * @param corners the corners of the quadrangle
 * @return 0 for success, -1 for error
 */
int
fill_is_in_descriptor (is_in_descriptor * const iid,
		       const px_point * const corners)
{
  px_point points[3];

  points[0] = corners[0];
  points[1] = corners[1];
  points[2] = corners[3];
  iid->mini =
    MIN (MIN (MIN (corners[0].x, corners[1].x), corners[2].x),
	 corners[3].x);
  iid->maxi =
    MAX (MAX (MAX (corners[0].x, corners[1].x), corners[2].x),
	 corners[3].x);
  iid->minj =
    MIN (MIN (MIN (corners[0].y, corners[1].y), corners[2].y),
	 corners[3].y);
  iid->maxj =
    MAX (MAX (MAX (corners[0].y, corners[1].y), corners[2].y),
	 corners[3].y);

  fill_triag_struct (&(iid->trig1), points);
  points[0] = corners[2];
  points[1] = corners[1];
  points[2] = corners[3];
  fill_triag_struct (&(iid->trig2), points);
  return 0;
}


/**
 * decides if a point is within a triangle.
 *
 * @param pix_x x-coordinate of the point to check.
 * @param pix_y y-coordinate of the point to check.
 * @param trig a pointer to the triangle structure.
 * @return 1 if point is in triangle, 0 otherwise.
 */
/* This works by computing normalized covariant coordinates */
static int
is_in_triag (const int pix_x, const int pix_y, const triangle * const trig)
{
  int x, y;
  double p_s, p_t, pabs, xsqrtarg, ysqrtarg, xp, yp;
  
  x = pix_x - trig->root.x;
  y = pix_y - trig->root.y;
  p_s = trig->s.x * x + trig->s.y * y;
  p_t = trig->t.x * x + trig->t.y * y;
  pabs = sqrt (SQR (x) + SQR (y));
  xsqrtarg =
    (1 - SQR (p_s) / SQR (pabs) / SQR (trig->sabs)) / (1 - SQR (trig->s_t) / SQR (trig->sabs) / SQR (trig->tabs));
  ysqrtarg =
    (1 - SQR (p_t) / SQR (pabs) / SQR (trig->tabs)) / (1 - SQR (trig->s_t) / SQR (trig->sabs) / SQR (trig->tabs));

  if (xsqrtarg < 0)
    {
      if (xsqrtarg > -1e-8)
	{
	  xsqrtarg = 0;
	}
      else
	{
	  printf ("is_in_triag: %d %d x%g\n", x, y, xsqrtarg);
	  xsqrtarg = 0;
	  
	  //abort ();
	}
    }
  if (ysqrtarg < 0)
    {
      if (ysqrtarg > -1e-8)
	{
	  ysqrtarg = 0;
	}
      else
	{
	  printf ("is_in_triag: y%g\n", ysqrtarg);
	  ysqrtarg = 0;
	  
	  //abort ();
	}
    }
  
  xp = 1. / SQR (trig->sabs) * (p_s -
				trig->s_t * pabs / trig->tabs *
				sqrt (xsqrtarg));
  yp = 1. / SQR (trig->tabs) * (p_t -
				trig->s_t * pabs / trig->sabs *
				sqrt (ysqrtarg));

  
  /* The 1.0xx1 should be harmless for 
     quadrangles that are not too extended since these points are at
     the point where the triangles are stitched together. */
  if ((-1e-10 <= xp) && (-1e-10 <= yp) && ((xp + yp) < (1.0 + 1e-10)))
    return 1;
  
  return 0;
}


/**
  Decides if a point is within a quadrangle.

  @param x The x-coordinate of the point.
  @param y The y-coordinate of the point.
  @param iid The descriptor of the quadrangle.
  @return 1 if point is in quadrangle, 0 otherwise.
  @see fill_is_in_descriptor
*/
int
is_in (const int x, const int y, const is_in_descriptor * const iid)
{
  if ((x < iid->mini) || (x > iid->maxi) || (y < iid->minj)
      || (y > iid->maxj))
    return 0;
  
  if (is_in_triag (x, y, &(iid->trig1)))
    return 1;
  if (is_in_triag (x, y, &(iid->trig2)))
    return 1;
  return 0;
}

int
tracedist_criteria(const double x, const double y, sectionfun *sf,
		   const trace_func *tracefun, const double width)
{
  double x_sect, y_sect;
  double dist;

  int ireturn=0;


  find_section_point (sf, x, y, &x_sect);
  y_sect = tracefun->func (x_sect, tracefun->data);
 
  dist =  sqrt((y_sect-y)*(y_sect-y)+(x_sect-x)*(x_sect-x));

  if (dist < width)
    ireturn=1;

  return ireturn;
}
