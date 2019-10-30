#include "spc_driz.h"

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


/**
Calculate the area common to input clockwise polygon x(n), y(n) with
square (is, js) to (is+1, js+1).
This version is for a quadrilateral.
W.B. Sparks STScI 2-June-1990.
Phil Hodge        20-Nov-1990  Change calling sequence; single precision.
Richard Hook ECF  24-Apr-1996  Change coordinate origin
                               so that centre of pixel has integer position
                  03-Jan-2001  Removed accuracy check
Nor Pirzkal       27-Jan-2003  C version
**/
double boxer(int is,int js,double *xx,double *yy)
{
  double px[4],py[4], sum;
  int i;
	
  /*
    Set up coords relative to unit square at origin
    Note that the +0.5s were added when this code was
    included in DRIZZLE
  */	
  for(i=0;i<4;i++) {
    px[i] = xx[i] - is + 0.5;
    py[i] = yy[i] - js + 0.5;
    // fprintf(stderr,"%d %d %f %f\n",is,js,px[i],py[i]);
  }
  
  /*
    For each line in the polygon (or at this stage, input quadrilateral)
    calculate the area common to the unit square (allow negative area for
    subsequent `vector' addition of subareas).
  */
  sum = 0.;
  for (i=0;i<3;i++) {
    sum += sgarea(px[i],py[i],px[i+1],py[i+1],is,js);
  }
  sum += sgarea(px[3],py[3],px[0],py[0],is,js);
  
  return sum;
  
}

double sgarea(double x1, double y1, double x2, double y2, int is, int js)
{
  double m,c,dx;
  double xlo,xhi,ylo,yhi,xtop;
  int negdx=0;
  double res;
  
  dx = x2-x1;
  
  /* 
     Trap vertical line
  */
  if (dx==0.0) {
    res = 0.0;
    return res;
  }
  
  /*
    Order the two input points in x
  */
  if (x1<x2) {
    xlo = x1;
    xhi = x2;
  } else {
    xlo = x2;
    xhi = x1;
  }

  /* 
     And determine the bounds ignoring y for now
  */
  if (xlo>=1.0) {
    res = 0.0;
    return res;
  }
  
  if (xhi<=0.0) {
    res = 0.0;
    return res;
  }
  
  xlo = MAX(xlo,0.0);
  xhi = MIN(xhi,1.0);
  
  /*
    Now look at y basic info about the line y = mx+c
  */
  negdx = (dx<0.0);
  m = (y2-y1)/dx;
  c = y1-m*x1;
  ylo = m*xlo+c;
  yhi = m*xhi+c;
  
  /*
    Trap segment entirely below axis
  */
  if ( (ylo<=0.0)&&(yhi<=0.0) ) {
    res = 0.0;
    return res;
  }
  
  /*
    Adjust bounds if segment crosses axis (to exclude anything below axis)
  */
  if (ylo<0.0) {
    ylo = 0.0;
    xlo = -c/m;
  }
  if (yhi<0.0) {
    yhi = 0.0;
    xhi = -c/m;
  }
  
  /* 
     There are four possibilities: 
     both y below 1
     both y above 1
     and one of each
  */
	

  if ( (ylo>=1.0) && (yhi>=1.0) ) {
    /* Line segment is entirrely above square */
    if (negdx) {
      res = xlo-xhi;
    } else {
      res = xhi-xlo;
    }
    return res;
  }
  
  if ( (ylo<=1.0) && (yhi<=1.0) ) {
    /* segment is entirely within square */
    if (negdx) {
      res = 0.5*(xlo-xhi)*(yhi+ylo);
    } else {
      res = 0.5*(xhi-xlo)*(yhi+ylo);
    }
    return res;
  }
  
  
  /* Otherwise it must cross the top square */
  xtop = (1.0 -c)/m;
  
  if (ylo<1.) {
    if (negdx) {
      res = -(0.5*(xtop-xlo)*(1.0+ylo)+xhi-xtop);
    } else {
      res = 0.5*(xtop-xlo)*(1.0+ylo)+xhi-xtop;
    }
    return res;
  }
  
  if (negdx) {
    res = -(0.5*(xhi-xtop)*(1.0+yhi)+xtop-xlo);
  } else {
    res = 0.5*(xhi-xtop)*(1.0+yhi)+xtop-xlo;
  }	
  
  return res;
}
