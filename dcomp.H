#ifndef __DCOMP_H
#define __DCOMP_H
#include <math.h>

/*
  Program: dcomp.H
  Author:  D. Trinkle
  Date:    2008 May 7
  Purpose: Define dcomp, insidecell
*/

// Some general functions, etc.

const double TOLER = 1e-7;

inline int dcomp (const double &a, const double &b) 
{
  if (a>b)
    return ( (a-b) < TOLER );
  else
    return ( (b-a) < TOLER );
}

inline double insidecell (const double &x) 
{
  double y;
  // Return a value less than 1, and greater than or equal to 0.
  for (y = x + 1 + 1.e-13 - (int)(x); y >= 1.0; y -= 1.0) ;
  y -= 1.e-13;
  if (y < 0.) return 0;
  return y;
}

#endif
