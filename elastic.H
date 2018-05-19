#ifndef __ELASTIC_H
#define __ELASTIC_H

/*
  Program: elastic.H
  Author:  D. Trinkle
  Date:    Oct. 2, 2002; modified August 14, 2003
  Purpose: Define the series of 11 crystal classes and their elastic
           constants for anisotropic elastic theory.

	   For any crystal, there are 21 elastic constants: C_ij,
	   i=1..6, j=1..6.  Where:

   	     1->11;    2->22;    3->33;
	     4->23,32; 5->13,31; 6->12,21

	   thanks to symmetry.

	   This code is basically a bookkeeping nightmare come
	   to life.  We have 11 different classes, and for each
	   class, we need to keep track of:

	   (a) how many unique C_ij there are
	   (b) how to map the C_ij into the unique list.

	   All notation, etc., comes from "Physical Properties
	   of Crystals" by J. F. Nye (1969); QD931 .N9 1969.  There
	   are no longer any 2d classes included.
*/

#include "matrix.H"

const int NCLASSES = 11; // 10 3D + isotropic + 5 2D classes

// 3D Crystal classes:
//  0: triclinic                         (a != b != c, alpha!=beta!=gamma)
//  1: monoclinic, diad || x_2           (a != b != c, alpha==gamma==90!=beta)
//  2: monoclinic, diad || x_3           (a != b != c, alpha==beta==90!=gamma)
//  3: orthorhombic                      (a != b != c, alpha==beta==gamma==90)
//  4: cubic          _                  (a == b == c, alpha==beta==gamma==90)
//  5: tetragonal (4  4  4|m)       _    (a == b != c, alpha==beta==gamma==90)
//  6: tetragonal (4mm _422  4|mmm  42m) (a == b != c, alpha==beta==gamma==90)
//  7: trigonal,  (3   3)  _        (a == b != c, alpha==beta==90, gamma==120)
//  8: trigonal,  (32  3m  3m)      (a == b != c, alpha==beta==90, gamma==120)
//  9: hexagonal                    (a == b != c, alpha==beta==90, gamma==120)
// 10: isotropic

const int CLASS_TRICLINIC =         0;
const int CLASS_MONOCLINIC_X2 =     1;
const int CLASS_MONOCLINIC_X3 =     2;
const int CLASS_ORTHORHOMBIC =      3;
const int CLASS_CUBIC =             4;
const int CLASS_TETRAGONAL_NO45 =   5;
const int CLASS_TETRAGONAL      =   6;
const int CLASS_TRIGONAL_NOMIRROR = 7;
const int CLASS_TRIGONAL =          8;
const int CLASS_HEXAGONAL  =        9;
const int CLASS_ISOTROPIC =        10;


// Point groups if we have inversion...
const char* CLASS_NAME_INVERSION[NCLASSES] = {
  "triclinic (a!=b!=c, alpha!=beta!=gamma): -1, order=2",
  "monoclinic (a!=b!=c, alpha==gamma==90): 2|m, order=4",
  "monoclinic (a!=b!=c, alpha==beta==90): 2|m, order=4",
  "orthorhombic (a!=b!=c, alpha=beta=gamma=90): mmm, order=8",
  "cubic (a=b=c, alpha=beta=gamma=90): m3m, order = 48  or  m3, order=24",
  "tetragonal--no diag mirror plane (a=b!=c, alpha=beta=gamma=90): 4|m, order=8",
  "tetragonal (a=b!=c, alpha=beta=gamma=90): 4|mmm, order=16",
  "trigonal--no diag mirror plane (a=b!=c, alpha=beta=90, gamma=120): -3, order=6",
  "trigonal (a=b!=c, alpha=beta=90, gamma=120): -3m, order=12",
  "hexagonal (a=b!=c, alpha=beta=90, gamma=120): 6|mmm, order=24  or  6|m, order=12",
  "isotropic"
};

// Point groups if we don't have any inversion... (2d is the same)
const char* CLASS_NAME_NOINVERSION[NCLASSES] = {
  "triclinic (a!=b!=c, alpha!=beta!=gamma): 1, order=1",
  "monoclinic (a!=b!=c, alpha==gamma==90): 2 m, order=2",
  "monoclinic (a!=b!=c, alpha==beta==90): 2 m, order=2",
  "orthorhombic (a!=b!=c, alpha=beta=gamma=90): 222 mm2, order=4",
  "cubic (a=b=c, alpha=beta=gamma=90): -43m 432, order = 24  or  23, order=12",
  "tetragonal--no diag mirror plane (a=b!=c, alpha=beta=gamma=90): 4 -4, order=4",
  "tetragonal (a=b!=c, alpha=beta=gamma=90): -42m 4mm 422, order=8",
  "trigonal--no diag mirror plane (a=b!=c, alpha=beta=90, gamma=120): 3, order=3",
  "trigonal (a=b!=c, alpha=beta=90, gamma=120): 3m 32, order=6",
  "hexagonal (a=b!=c, alpha=beta=90, gamma=120): -6m2 6mm 622, order=12  or  6 -6, order=6",
  "isotropic"
};

const char* CRYSTAL_CLASS = 
"  0: triclinic                         (a != b != c, alpha!=beta!=gamma)\n\
  1: monoclinic, diad || x_2           (a != b != c, alpha==gamma==90!=beta)\n\
  2: monoclinic, diad || x_3           (a != b != c, alpha==beta==90!=gamma)\n\
  3: orthorhombic                      (a != b != c, alpha==beta==gamma==90)\n\
  4: cubic          _                  (a == b == c, alpha==beta==gamma==90)\n\
  5: tetragonal (4  4  4|m)       _    (a == b != c, alpha==beta==gamma==90)\n\
  6: tetragonal (4mm _422  4|mmm  42m) (a == b != c, alpha==beta==gamma==90)\n\
  7: trigonal,  (3   3)  _        (a == b != c, alpha==beta==90, gamma==120)\n\
  8: trigonal,  (32  3m  3m)      (a == b != c, alpha==beta==90, gamma==120)\n\
  9: hexagonal                    (a == b != c, alpha==beta==90, gamma==120)\n\
 10: isotropic";


// This is basically just the elastic constant matrix from p.140-1
// from Nye's book.  My version of the notation:
//   * If an element is listed as "0", there is no elastic constant
//   * If an element has a positive value, it has a elastic constant
//   * Matching numbers have the same elastic constant
//   * If an element has a NEGATIVE value, it has the same elastic
//     constant as the positive valued C_ij, with opposite sign.
//   * Finally, element 99 == 0.5*(C11 - C12) (marked as "x" by Nye)
const int CIJ_matrix[NCLASSES][6][6] = {
  //  3D lattice classes:
  // 0: triclinic
  {{1, 2, 3, 4, 5, 6},
   {2, 7, 8, 9,10,11},
   {3, 8,12,13,14,15},
   {4, 9,13,16,17,18},
   {5,10,14,17,19,20},
   {6,11,15,18,20,21}},
  // 1: monoclinic, diad || x_2 (standard orientation)
  {{1, 2, 3, 0, 4, 0},
   {2, 5, 6, 0, 7, 0},
   {3, 6, 8, 0, 9, 0},
   {0, 0, 0,10, 0,11},
   {4, 7, 9, 0,12, 0},
   {0, 0, 0,11, 0,13}},
  // 2: monoclinic, diad || x_3
  {{1, 2, 3, 0, 0, 4},
   {2, 5, 6, 0, 0, 7},
   {3, 6, 8, 0, 0, 9},
   {0, 0, 0,10,11, 0},
   {0, 0, 0,11,12, 0},
   {4, 7, 9, 0, 0,13}},
  // 3: orthorhombic
  {{1, 2, 3, 0, 0, 0},
   {2, 4, 5, 0, 0, 0},
   {3, 5, 6, 0, 0, 0},
   {0, 0, 0, 7, 0, 0},
   {0, 0, 0, 0, 8, 0},
   {0, 0, 0, 0, 0, 9}},
  // 4: cubic
  {{1, 2, 2, 0, 0, 0},
   {2, 1, 2, 0, 0, 0},
   {2, 2, 1, 0, 0, 0},
   {0, 0, 0, 3, 0, 0},
   {0, 0, 0, 0, 3, 0},
   {0, 0, 0, 0, 0, 3}},
  // 5: tetragonal, classes 4  -4  4|m
  {{1, 2, 3, 0, 0, 4},
   {2, 1, 3, 0, 0,-4},
   {3, 3, 5, 0, 0, 0},
   {0, 0, 0, 6, 0, 0},
   {0, 0, 0, 0, 6, 0},
   {4,-4, 0, 0, 0, 7}},
  // 6: tetragonal, classes 4mm -422  4|mmm  -42m
  {{1, 2, 3, 0, 0, 0},
   {2, 1, 3, 0, 0, 0},
   {3, 3, 4, 0, 0, 0},
   {0, 0, 0, 5, 0, 0},
   {0, 0, 0, 0, 5, 0},
   {0, 0, 0, 0, 0, 6}},
  // 7: trigonal,   classes 3   3
  {{ 1, 2, 3, 4,-5, 0},
   { 2, 1, 3,-4, 5, 0},
   { 3, 3, 6, 0, 0, 0},
   { 4,-4, 0, 7, 0, 5},
   {-5, 5, 0, 0, 7,-4},
   { 0, 0, 0, 5,-4,99}},
  // 8: trigonal,   classes 32  3m  -3m
  {{ 1, 2, 3, 4, 0, 0},
   { 2, 1, 3,-4, 0, 0},
   { 3, 3, 5, 0, 0, 0},
   { 4,-4, 0, 6, 0, 0},
   { 0, 0, 0, 0, 6,-4},
   { 0, 0, 0, 0,-4,99}},
  // 9: hexagonal
  {{1, 2, 3, 0, 0, 0},
   {2, 1, 3, 0, 0, 0},
   {3, 3, 4, 0, 0, 0},
   {0, 0, 0, 5, 0, 0},
   {0, 0, 0, 0, 5, 0},
   {0, 0, 0, 0, 0,99}},
  // 10: isotropic
  {{1, 2, 2, 0, 0, 0},
   {2, 1, 2, 0, 0, 0},
   {2, 2, 1, 0, 0, 0},
   {0, 0, 0,99, 0, 0},
   {0, 0, 0, 0,99, 0},
   {0, 0, 0, 0, 0,99}}
};

// Class elastic constants (simply stored for human readability):
const int class_Cij[NCLASSES][21] = {
  // 0: triclinic
  {11,12,13,14,15,16,22,23,24,25,26,33,34,35,36,44,45,46,55,56,66},
  // 1: monoclinic, diad || x_2 (standard orientation)
  {11,12,13,15,22,23,25,33,35,44,46,55,66, 0, 0, 0, 0, 0, 0, 0, 0},
  // 2: monoclinic, diad || x_3
  {11,12,13,16,22,23,26,33,36,44,45,55,66, 0, 0, 0, 0, 0, 0, 0, 0},
  // 3: orthorhombic
  {11,12,13,22,23,33,44,55,66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 4: cubic
  {11,12,44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 5: tetragonal, classes 4  -4  4|m
  {11,12,13,16,33,44,66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 6: tetragonal, classes 4mm -422  4|mmm  -42m
  {11,12,13,33,44,66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 7: trigonal,   classes 3   3
  {11,12,13,14,25,33,44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 8: trigonal,   classes 32  3m  -3m
  {11,12,13,14,33,44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 9: hexagonal
  {11,12,13,33,44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  // 10: isotropic
  {11,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};

// Class length:
const int class_len[NCLASSES] = {21,13,13,9,3,7,6,7,6,5,2};


// Now, the routine to convert an input list of Cmn that is
// consistent with the above lists into a real, usable, Cijkl matrix.

// First, from ij to m:
const int ij2m[3][3] = 
  {{0, 5, 4}, 
   {5, 1, 3},
   {4, 3, 2}};

// This assumes that Cijkl has already been declared as
//   double Cijkl[3][3][3][3]
// and that Cmn_list is a single array of Cmn values, according
// to class_Cij above, so that if I want the correct value
// for a general Cmn for class c, I look at:
//   Cmn_list[ CIJ_matrix[c][m][n] - 1 ]
// being careful to (a) check for 0, (b) negative values, and (c) 99.
// The end result should be the correct Cijkl matrix... woohoo!

void make_Cijkl (int c, double* Cmn_list, double Cijkl[9][9]) 
{
  int i, j, k, l;
  double C99;
  int cij_in;

  // First, calculate "C99" = 0.5*(C11-C12)  (may not be used):
  C99 = 0.5 * ( Cmn_list[ CIJ_matrix[c][0][0] - 1 ] -
		Cmn_list[ CIJ_matrix[c][0][1] - 1 ] );

  for (i=0; i<3; ++i) 
    for (j=0; j<3; ++j) 
      for (k=0; k<3; ++k) 
	for (l=0; l<3; ++l) {
	  cij_in = CIJ_matrix[c][ij2m[i][j]][ij2m[k][l]];
	  if (cij_in == 0)
	    Cijkl[index(i,j)][index(k,l)] = 0.;
	  else
	    if (cij_in == 99)
	      Cijkl[index(i,j)][index(k,l)] = C99;
	    else
	      if (cij_in > 0)
		Cijkl[index(i,j)][index(k,l)] = Cmn_list[ cij_in - 1 ];
	      else
		Cijkl[index(i,j)][index(k,l)] = - Cmn_list[ (-cij_in) - 1 ];
	}
}

#endif
