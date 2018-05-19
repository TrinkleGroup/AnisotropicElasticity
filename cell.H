#ifndef __CELL_H
#define __CELL_H

/*
  Program: cell.H
  Author:  D. Trinkle
  Date:    Oct. 6, 2002, modified Aug. 21, 2003 for dislocations
  Purpose: This reads in the cell information.

	   We write our Cartesian vectors in column primary form:

	   | a1x a2x a3x |
	   | a1y a2y a3y |
	   | a1z a2z a3z |
	   
	   Our cell input file contains the following information:

	   ==== input file ====
	   a0                            # Scale factor for unit cell
	   a1.x a1.y a1.z                # Cartesian coord of unit cell
	   a2.x a2.y a2.z
	   a3.x a3.y a3.z
	   crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
	   Natoms                        # Number of atoms in first unit cell
           u1.1 u1.2 u1.3                # Atom locations, in direct coord.
           ...
           uN.1 uN.2 uN.3
	   ==== input file ====

	   Notes:
	   1. There can be comments on every line, though on the
	      crystal class C_ij line, the order of C_ij is important.
	   2. In theory, you can not read the atoms, by setting
	      Natoms = -1.

	   This header file contains the following subroutines:

	   read_cell (cell_file, cart, crystal_class, C_ij, u, Natoms)
	     --reads in the input file (formatted as above).  Assumes that
	       cell_file is an already open file (so that we can feed
	       stdin as an option).  Doesn't close the file either.
	       Allocates memory for each of the pointers, and fills it up.
	       Return value: ERROR code.

	   verbose_output_cell (cart, crystal_class, C_ij, u, Natoms) 
	     --outputs everything we read in from the input file.

           free_cell (C_ij, u, Natoms)
	     --frees up the associated memory (done for completeness only).
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "dcomp.H"
#include "io-short.H"
#include "matrix.H"
#include "elastic.H"

//****************************** ERROR FLAGS ***************************

// Note: we reuse error codes, since once we throw an error, everything
//       will stop anyway.

// errors for read_cell:
const int ERROR_NOFILE =      1; // General error if a file can't be opened.
const int ERROR_BADFILE =     2; // File has wrong format.
const int ERROR_ZEROVOL =     4; // Cartesian cell has zero volume
const int ERROR_LEFTHANDED =  8; // Cartesian cell has negative volume
const int ERROR_MEMORY =     16; // Not enough memory

// Flags for read_cell:
const int NO_ATOMS = -1;

//****************************** SUBROUTINES ***************************

int read_cell (FILE *cell_file, double cart[9], int &crystal_class, 
	       double* &Cmn_list, double** &u, int &Natoms);

void verbose_output_cell (double cart[9], int crystal_class, 
                          double* Cmn_list, double** u, int Natoms);

void free_cell (double* &Cmn_list, double** &u, int Natoms);


//****************************** read_cell *****************************
// We read in the input file (formatted as above).  Assumes that
// cell_file is an already open file (so that we can feed stdin as an
// option).  Doesn't close the file either.
// Allocates memory for each of the pointers, and fills it up.
// Note: if Nop = NO_GROUP then we don't read the group info.
//       if Natoms = NO_ATOMS then we don't read the atom info.
// Return value: ERROR code (see above)
// ==== input file ====
// a0                            # Scale factor for unit cell
// a1.x a1.y a1.z                # Cartesian coord of unit cell
// a2.x a2.y a2.z
// a3.x a3.y a3.z
// crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
// Natoms                        # Number of atoms in first unit cell
// u1.1 u1.2 u1.3                # Atom locations, in direct coord.
// ...
// uN.1 uN.2 uN.3
// ==== input file ====
int read_cell (FILE *cell_file, double cart[9], int &crystal_class, 
	       double* &Cmn_list, double** &u, int &Natoms) 
{
  double a0; // Scale factor.
  double det_cart;           // det(cart)
  char dump[512]; // For reading files (allows extra characters at the end
                  // of each line)
  int i, k;

  if (cell_file == NULL)
    return ERROR_NOFILE;

  // ========================== read input file ========================

  // a0                            # Scale factor for unit cell
  nextnoncomment(dump, sizeof(dump), cell_file);
  sscanf(dump, "%lf", &a0);

  // a1.x a1.y a1.z                # Cartesian coord of unit cell
  // a2.x a2.y a2.z
  // a3.x a3.y a3.z
  for (i=0; i<3; ++i) {
    nextnoncomment(dump, sizeof(dump), cell_file);
    sscanf(dump, "%lf %lf %lf", 
	   cart+index(0,i), cart+index(1,i), cart+index(2,i));
  }
  for (i=0; i<9; ++i)
    cart[i] *= a0;

  // TEST: Determine determinants:
  det_cart = det(cart);

  if (dcomp(det_cart, 0)) return ERROR_ZEROVOL;
  if (det_cart < 0)       return ERROR_LEFTHANDED;


  // crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
  nextnoncomment(dump, sizeof(dump), cell_file);
  sscanf(dump, "%d", &crystal_class);

  // TEST: Good crystal class?
  if ( (crystal_class < 0 ) || (crystal_class > NCLASSES) )
    return ERROR_BADFILE;

  // Now, let's read those elastic constants.
  // This is a little tricky; watch and learn:
  Cmn_list = new double [class_len[crystal_class]];
  char *startp, *endp;
  // We need to advance past the first number (the crystal class)
  startp = dump;
  strtod(startp, &endp);
  startp = endp;
  for (i=0; i<class_len[crystal_class]; ++i) {
    Cmn_list[i] = strtod(startp, &endp);
    // check for conversion:
    if (startp == endp) break;
    startp = endp;
  }
  // TEST: Did we read enough C_ij's?
  if (i != class_len[crystal_class]) {
    return ERROR_BADFILE;
  }


  if (Natoms == NO_ATOMS) {
    // Read just the number of atoms.

    // Natoms                        # Number of atoms in first unit cell
    nextnoncomment(dump, sizeof(dump), cell_file);
    sscanf(dump, "%d", &Natoms);
    u = NULL;
  }
  else {
    // Read in everything.

    // Natoms                        # Number of atoms in first unit cell
    nextnoncomment(dump, sizeof(dump), cell_file);
    sscanf(dump, "%d", &Natoms);
    // u1.1 u1.2 u1.3                # Atom locations, in direct coord.
    // ...
    // uN.1 uN.2 uN.3
    u = new double *[Natoms];
    for (i=0; i<Natoms; ++i) {
      u[i] = new double[3];
      nextnoncomment(dump, sizeof(dump), cell_file);
      sscanf(dump, "%lf %lf %lf", u[i], u[i]+1, u[i]+2);
      for (k=0; k<3; ++k)
	u[i][k] = insidecell(u[i][k]);
    }
  }

  // If we've gotten this far, it all went well.
  return 0;
}



// ************************* verbose_output_cell ***********************
// Little routine to output everything we read in from cell.

void verbose_output_cell (double cart[9], int crystal_class, 
                          double* Cmn_list, double** u, int Natoms)
{
  int i;

  printf("# Cartesian coordinates:\n");
  for (i=0; i<3; ++i)
    printf("# a%1d = %8.5lf %8.5lf %8.5lf\n", i+1, 
	   cart[i], cart[i+3], cart[i+6]);

  printf("# Atoms in cell (%d):\n", Natoms);
  if (u == NULL)
    printf("# atom positions not read.\n");
  else
    for (i=0; i<Natoms; ++i)
      printf("# u%1d = %8.5lf %8.5lf %8.5lf\n", i+1, 
	     u[i][0], u[i][1], u[i][2]);

  printf("# You've chosen the crystal class %d:\n# %s\n",
         crystal_class, CLASS_NAME_INVERSION[crystal_class]);
  printf("# which has the following %d elastic constants:\n#", 
	 class_len[crystal_class]);
  for (i=0; i<class_len[crystal_class]; ++i)
    printf(" C_%2d = %.3lf", class_Cij[crystal_class][i], Cmn_list[i]);
  printf("\n");
}



//****************************** free_cell *****************************
// Free up the group and atom info:

void free_cell (double* &Cmn_list, double** &u, int Natoms)
{
  int i;

  if (Cmn_list != NULL)
    delete[] Cmn_list;
  Cmn_list = NULL;

  if (u != NULL) {
    for (i=0; i<Natoms; ++i)
      delete[] u[i];
    delete[] u;
  }
  u = NULL;
}

#endif
