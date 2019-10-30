#ifndef _SPC_IS_IN_H
#define _SPC_IS_IN_H


/**
  A triangle together with some pre-computed properties, consisting of
  the coordinates  of the root vertex, the coordinates of the ends
  of the sides s, t relative to the root and some pre-computed
  properties.
*/
typedef struct
{
  px_point root;		/* pixel coordinate of the angle point */
  px_point s, t;		/* the  end points of the sides relative to root */
  double sabs, tabs;		/* length of s and t */
  double s_t;	                /* scalar product of s and t */
}
triangle;

/**
  An opaque structure describing a quadrangle through two triangles
  @see fill_is_in_descriptor
*/
typedef struct
{
  triangle trig1, trig2;	/* The quadrangle is split into two triangles */
  int mini, maxi;		/* The mininum and maximum row-span of the
				   aperture */
  int minj, maxj;		/* The mininum and maximum col-span of the
				   aperture */
}
is_in_descriptor;


/* public */

extern int
is_in (const int x, const int y, const is_in_descriptor * const iid);

extern int 
fill_is_in_descriptor (is_in_descriptor * const iid,
		       const px_point * const corners);

extern int
tracedist_criteria(const double x, const double y, sectionfun *sf,
		   const trace_func *tracefun, const double width);
#endif

