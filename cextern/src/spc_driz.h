#ifndef _SPC_DRIZ_H
#define _SPC_DRIZ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>

double sgarea(double x1, double y1, double x2, double y2, int is, int js);
double boxer(int is,int js,double *x,double *y);


#endif
