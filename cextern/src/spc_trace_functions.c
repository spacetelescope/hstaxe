/**
 * Some functions to parametrize spectrum traces in aXe grism exposures
 */

#include "spc_trace_functions.h"

/*
A trace function ideally consists of the function itself, the derivative,
and a path length.  The path length function may be NULL, and in the
future the derivative, too.

A new type of trace function needs to provide an alloc and a free function
each, after the model of create_poly2 and free_poly2.

*/

/**
  A polynom of degree two.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer four doubles, the length (3) and the coefficients of the polynom.
  @returns the value of the polynom at x.
*/
double
poly2 (const double x, const void *const pars)
{
     const double *const coeffs = pars;

     return coeffs[1] + (coeffs[2] + (coeffs[3] * x)) * x;
}

/**
  A polynom of degree N.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer N+1 doubles, the length (N) and the coefficients of the polynom.
  @returns the value of the polynom at x.
*/
double
polyN (const double x, const void *const pars)
{
     const double *const coeffs = pars;
	 int i,N;
	 double p=0.0;
	 N = (int) coeffs[0];
	//fprintf(stderr,"---->length:%d %d %d\n",sizeof(*coeffs),sizeof(double),N);

	 for (i=0;i<N;i++) {
	 	p += coeffs[i+1]*pow(x,i);
	 }
     return p;
}

/**
  The derivative of a polynom of degree two.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer four doubles, the length (3) and the coefficients of the polynom.
  @returns the derivative of the polynom at x.
*/
double
poly2_deriv (const double x, const void *const pars)
{
     const double *const coeffs = pars;

     return coeffs[2] + 2 * coeffs[3] * x;
}

/**
  The derivative of a polynom of degree N.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer N+1 doubles, the length (N) and the coefficients of the polynom.
  @returns the derivative of the polynom at x.
*/
double
polyN_deriv (const double x, const void *const pars)
{
     const double *const coeffs = pars;
	 int i, N;
	 double p=0.0;
	 
	 N = (int) coeffs[0];
	 for (i=1;i<N;i++) {
		//fprintf(stderr,"%d %f %d %f %f %f\n",i,coeffs[i+1],i-1,x,pow(x,i-1),p);
	 	p += i*coeffs[i+1]*pow(x,i-1);
		 }
		 
	 //fprintf(stderr,"poly2deriv:%f polyNderiv:%f\n",poly2_deriv(x,pars),p);
     return p;
}

/**
  sqrt(1+(dy/dx)^2) of an N^th order polynomial

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer N+1 doubles, the length (N) and the coefficients of the polynom.
  @returns the derivative of the polynom at x.
*/
double
polyN_ds (double x, void *pars)
{	
	 return sqrt(1.+pow(polyN_deriv(x,pars),2.));
}


/**
  The path length of a polynom of degree two, relative to the point poly4(0).
  For higher polynoms, probably no such closed expression exists.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer four doubles, the length (3) and the coefficients of the polynom.
  @returns the path lenght for the polynom at x.
*/
double
poly2_pathlen (const double x, const void *const pars)
{
     double a1 = ((double *) pars)[2], a2 = ((double *) pars)[3];

     if (fabs (a2) <= 1e-20)
	  return sqrt (1 + a1 * a1) * x;


     return 1 / 4. / a2 * (2.0 *
			   sqrt (1.0 + a1 * a1 + 4.0 * a1 * a2 * x +
				 4.0 * a2 * a2 * x * x) * a2 * x +
			   log (2.0 * a2 * x + a1 +
				sqrt (1.0 + a1 * a1 + 4.0 * a1 * a2 * x +
				      4.0 * a2 * a2 * x * x)) + sqrt (1.0 +
								      a1 *
								      a1 +
								      4.0 *
								      a1 *
								      a2 * x +
								      4.0 *
								      a2 *
								      a2 * x *
								      x) *
			   a1 - sqrt (1.0 + a1 * a1) * a1 - log (a1 +
								 sqrt (1.0 +
								       a1 *
								       a1)));
}

/**
  The path length of a polynom of degree two, relative to the point poly4(0).
  For higher polynoms, probably no such closed expression exists.

  @param x the abscissa, usually relative to the beam's reference point.
  @param pars a pointer N+1 doubles, the length (N) and the coefficients of the polynom.
  @returns the path lenght for the polynom at x.
*/
double 
polyN_pathlen (const double x, const void *const pars)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_function F;
	F.function = &polyN_ds;

	// new to avoid warning during compilation
	F.params = (void *)pars;
	// old, with warning during compilation
	//	F.params = pars;
	//fprintf(stderr,"HEREEEE 3\n");
	//gsl_integration_qags (&F, 0, x, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_qag (&F, 0, x, 0, 1e-7, 1000, 6, w, &result, &error);
	gsl_integration_workspace_free(w);
	//result=x;
	//fprintf(stderr,"===> %f %f\n",result,poly2_pathlen(x,pars));
	return result;

}
/**
  creates a 2nd order polynomial trace function
  
  @param a0 coefficient for x^0
  @param a1 coefficient for x^1
  @param a2 coefficient for x^2
  @returns allocated trace_func structure or NULL for failure
*/
trace_func *
create_poly2 (const double a0, const double a1, const double a2)
{
     trace_func *func = malloc (sizeof (trace_func));
     double *cf;

     if (!func)
       {
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "create_poly2: Out of memory");
       }
     if (!(cf = malloc (3 * sizeof (double))))
       {
	    free (func);
	    func = NULL;
	    aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			 "create_poly2: Out of memory");
       }
     func->func = poly2;
     func->deriv = poly2_deriv;
     func->path_len = poly2_pathlen;
     func->data = cf;
     func->type = 2;
     cf[0] = a0;
     cf[1] = a1;
     cf[2] = a2;

     return func;
}

/**
  creates a Nth order polynomial trace function

  @param a pointer to a gsl_vector containing the coeeficients
  @returns allocated trace_func structure or NULL for failure
*/
trace_func *
create_polyN (gsl_vector *v)
{
  trace_func *func = malloc (sizeof (trace_func));
  double *cf;
  int i;

  //fprintf(stderr,"CREATING %d\n",v->size);
  if (!func)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_polyN: Out of memory");
    }
  if (!(cf = malloc ((v->size +1) * sizeof (double))))
    {
      free (func);
      func = NULL;
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		   "create_polyN: Out of memory");
    }

  func->func = polyN;
  func->deriv = polyN_deriv;
  func->path_len = polyN_pathlen;
  func->data = cf;
  func->type = v->size - 1;
  cf[0] = v->size;
  for (i=0;i<(int)v->size;i++)
    {
      cf[i+1]=gsl_vector_get (v, i);
    }

  return func;
}


/**
  creates a 2nd order polynomial trace function
  using the coefficients contained in a gsl vector 
  of length 1,2, or 3.
  
  @param v, a gsl vector containing polynomial coefficients. 
  @return a pointer to new trace_func object which contains a second order polynomial
  trace description
*/
trace_func *
vector_to_trace_poly2 (gsl_vector *v)
{
	    	
	float a0, a1, a2;

	if (v->size<1)
		aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
		"Trace polynomial must at least have one coefficients (i.e. zeroth order).\n"
		"Only %d found.",v->size);
	a0 = gsl_vector_get (v, 0);
	if (v->size>1)
		a1 = gsl_vector_get (v, 1);
	else
		a1 = 0.0;
	if (v->size>2)
		a2 = gsl_vector_get (v, 2);
	else
		a2 = 0.0;
	if (v->size>3)
		aXe_message (aXe_M_WARN4, __FILE__, __LINE__,	
		"%d trace coefficients read but only 2nd order used!",v->size);
		
	return create_poly2 (a0, a1, a2);


}

/**
  creates a Nth order polynomial trace function
  using the coefficients contained in a gsl vector 
  of length N.
  
  @param v, a gsl vector containing polynomial coefficients. 
  @return a pointer to new trace_func object which contains an Nth order polynomial
  trace description
*/
trace_func *
vector_to_trace_polyN (gsl_vector *v)
{		
  return create_polyN (v);
}

/**
  Free the trace_func struct for a 2nd order polynomial.

  @param func the trace_func struct.
*/
void
free_poly2 (trace_func * func)
{
     free (func->data);
     func->data = NULL;
     free (func);
     func = NULL;
}

/**
  Free the trace_func struct for a Nth order polynomial.

  @param func the trace_func struct.
*/
void
free_polyN (trace_func * func)
{
  free (func->data);
  func->data = NULL;
  free (func);
  func = NULL;
}
