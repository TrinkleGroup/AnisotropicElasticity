#ifndef __INTEGRATE_H
#define __INTEGRATE_H

/*
  Program: integrate.H
  Author:  D. Trinkle
  Date:    August 15, 2003 (from Numerical Recipes)
  Purpose: Computation of abscissae and weights for Gauss-Legendre
           quadrature on the interval (-1,1).

	   Code from "Numerical Recipes in C," 2/e.

	   Converted to double and to use reasonable indices for arrays.

	   Also, constants for use in a bootstrapping integration method.
*/

#include <stdio.h>
#include <math.h>

// This scheme comes by using Simpson's extended rule on crack 
// (Num. Rec., p. 134)

// f_-2 = -2.*f_2 + 8.*f_1 - 5.*f_2
const double f_min_2[3] = {-2., 8., -5.};
// f_-1 =    -f_2 +    f_1 +    f_0
const double f_min_1[3] = { -1., 1., 1.};
// Now, as for our integration scheme, if we've evaluated the integral
// up to F_(N-1), then define f0 = f(n), f1 = f(n-1), f2 = f(n-2), f3 = f(n-3)
// Then, F_N - F_(N-1) = h*(9*f0 + 19*f1 - 5*f2 + f3)/24.:
const double int_weight[4] = {9./24., 19./24., -5./24., 1./24.};
// When "bootstrapped" with the previously defined values for f(-2) and
// f(-1), F(0)=0, then F(1) is the trapezoid rule, and F(2) is Simpson's rule.

const double EPS = 3.0e-13;

void gauleg(double x1, double x2, double x[], double w[], int n)
{
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);
    for (i=0; i<m; ++i) {
      z = cos(M_PI*(i+0.75)/(n+0.5));
      do {
	p1 = 1.;
	p2 = 0.;
	for (j=0; j<n; ++j) {
	  p3 = p2;
	  p2 = p1;
	  p1 = ((2.*j + 1.)*z*p2 - j*p3)/(j+1.);
	}
	pp = n*(z*p1-p2)/(z*z-1.);
	z1 = z;
	z = z1 - p1/pp;
      } while ( fabs(z-z1) > EPS);
      x[i] = xm - xl*z;
      x[n-1-i] = xm + xl*z;
      w[i] = 2.0*xl/((1.-z*z)*pp*pp);
      w[n-1-i] = w[i];
    }
}

#endif
