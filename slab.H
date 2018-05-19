#ifndef __SLAB_H
#define __SLAB_H

/*
  Program: slab.H
  Author:  D. Trinkle
  Date:    August 21, 2003
  Purpose: Construct a cylindrical slab.  Pretty icky stuff, but it's here.
*/

#include <stdio.h>
#include <math.h>
#include "dcomp.H"
#include "matrix.H"

// Makes a cylindrical slab, where:
// t is the vertical axis (will be the z axis)
// m is the cut axis (will be the x axis)
// n is the mutual perp (will be the y axis)
// Rcut is the radius in the xy plane
// c is the cart. coord. of the center of the dislocation
// |m|=|n| = 1, while |t| is the thickness of the slab
// [a] is the cartesian coordinates of the lattice
// u is the matrix of atoms, in unit cell coord.
// Natoms is the number of atoms in the unit cell
//   We output:
// Nslab: number of atoms in the slab
// xyz: xyz positions of atoms in the slab (xyz[n][1,2,3])


int construct_slab (double t[3], double m[3], double n[3], double c[3],
		    double Rcut,
		    double a[9], double** u, char** names, int Natoms, 
		    int& Nslab, double** &xyz, char** &types) 
{
  int ERROR = 0;
  // S = (m,n,t): transformation matrix
  // Sa = S^-1*[a]
  // aS = [a]^-1*S
  double S[9], Sa[9], aS[9];
  double ainv[9], Sinv[9];
  double deter;
  double Rcut2;
  double tmagn;

  Rcut2 = Rcut*Rcut;

  // Construct the rotation matrix
  S[0] = m[0];  S[1] = n[0];  S[2] = t[0];
  S[3] = m[1];  S[4] = n[1];  S[5] = t[1];
  S[6] = m[2];  S[7] = n[2];  S[8] = t[2];
  tmagn = det(S);

  deter = 1./inverse(S, Sinv);
  for (int d=0; d<9; ++d) Sinv[d] *= deter;
  mult(Sinv, a, Sa);

  deter = 1./inverse(a, ainv);
  for (int d=0; d<9; ++d) ainv[d] *= deter;
  mult(ainv, S, aS); // aS = [a]^-1 * S
  // We use this to make our maximums for i0, i1, i2:
  int imax[3];
  for (int d=0; d<3; ++d) {
    imax[d] = (int)(Rcut*(fabs(aS[3*d]) + fabs(aS[3*d+1])) + fabs(aS[3*d+2])
		    + 1.999);
  }
  
  int Napprox;

  // Guess how many atoms we'll end up with at the end of the day...
  Napprox = (int)(2*Natoms * M_PI * Rcut*Rcut * tmagn / det(a));
  if (Napprox < 128) Napprox = 128;
  double** s = new double*[3];
  for (int d=0; d<3; ++d) s[d] = new double[Napprox];
  int* atomtype = new int[Napprox];

  // Calculate the atom positions for the unit cell, with the shift
  // built in:
  double** s_atom;
  double cu[3]; // Center, in unit coord.

  mult_vect(ainv, c, cu);
  for (int d=0; d<3; ++d) cu[d] = insidecell(cu[d]);

  s_atom = new double*[Natoms];
  for (int j=0; j<Natoms; ++j) {
    s_atom[j] = new double[3];
    double ushift[3];
    // Start inside the unit cell:
    for (int d=0; d<3; ++d) ushift[d] = insidecell(u[j][d] - cu[d]);
    mult_vect(Sa, ushift, s_atom[j]);
    // Make sure we're still inside the cylinder in the z direction:
    s_atom[j][2] = insidecell(s_atom[j][2]);
  }

  // Now, try to put all the atoms in.
  int cell[3];
  double scell[3], stry[3];
  int N = 0;
  for (cell[0]=-imax[0]; cell[0]<=imax[0]; ++cell[0])
    for (cell[1]=-imax[1]; cell[1]<=imax[1]; ++cell[1])
      for (cell[2]=-imax[2]; cell[2]<=imax[2]; ++cell[2]) {
	mult_vect(Sa, cell, scell);
	for (int j=0; j<Natoms; ++j) {
	  for (int d=0; d<3; ++d) stry[d] = scell[d] + s_atom[j][d];
	  if ( (stry[0]*stry[0] + stry[1]*stry[1]) < Rcut2 ) {
	    // Attempt to check for "collisions"
	    int found = dcomp(stry[2], 0.) || dcomp(stry[2], 1.);
	    if (found) {
	      // Do the *real* collision detection
	      found = 0;
	      for (int k=0; (k<N) && (!found); ++k)
		found = dcomp(stry[0], s[0][k]) && dcomp(stry[1], s[1][k]);
	      stry[2] = insidecell(stry[2]);
	    }
	    else
	      // Else, only add if we're between 0 and 1 along z axis:
	      found = (stry[2] < 0.) || ( stry[2] > 1.);
	    // We only add if *not* found.
	    if (!found) {
	      if (N >= Napprox) {++N; break;}
	      for (int d=0; d<3; ++d) s[d][N] = stry[d];
	      atomtype[N] = j;	// which atom this was--used only for the names
	      ++N;
	    }
	  }
	}
      }
  if (N>Napprox) {
    fprintf(stderr, "Error constructing cylinder... something went horribly awry\n");
    fprintf(stderr, "Found %d atoms, when estimated %d atoms.\n", N, Napprox);
    fprintf(stderr, "Continuing, but atoms will be missing.  Change construct_slab() routine.\n");
    ERROR = -1;
    N = Napprox;
  }
  // Now, convert into xyz coordinates:
  Nslab = N;
  xyz = new double*[N];
  for (int n=0; n<N; ++n) {
    xyz[n] = new double[3];
    xyz[n][0] = s[0][n];
    xyz[n][1] = s[1][n];
    xyz[n][2] = s[2][n] * tmagn;
  }
  if (names != NULL) {
    types = new char*[N];
    for (int n=0; n<N; ++n) types[n] = names[atomtype[n]];
  }
  
  // Garbage collection
  for (int j=0; j<Natoms; ++j) delete[] s_atom[j];
  delete[] s_atom;
  
  for (int d=0; d<3; ++d) delete[] s[d];
  delete[] s;
  delete[] atomtype;

  return ERROR;
}

inline
int construct_slab (double t[3], double m[3], double n[3], double c[3],
		    double Rcut,
		    double a[9], double** u, int Natoms,
		    int& Nslab, double** &xyz) 
{ char** types=(char**)NULL;
  return construct_slab(t, m, n, c, 
			Rcut,
			a, u, (char**)NULL, Natoms, 
			Nslab, xyz, types); }



void free_slab(int Nslab, double** &xyz) 
{
  if (xyz == NULL) return;
  for (int n=0; n<Nslab; ++n) delete[] xyz[n];
  delete[] xyz;
  xyz = NULL;
}


#endif
