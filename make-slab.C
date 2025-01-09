/*
  Program: make-slab.C
  Author:  D. Trinkle
  Date:    2008 May 16
  Purpose: Construct the undislocated slab, with the appropriate center, as
           an XYZ file

  Param.:  cell infile Rcut
           cell:     cell file (see below for format)
           infile:   input file (see below for format)
	   Rcut:     cutoff for building slab

	   ==== cell ====
           a0                            # Scale factor for unit cell
           a1.x a1.y a1.z                # Cartesian coord of unit cell
           a2.x a2.y a2.z
           a3.x a3.y a3.z
           crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
           Natoms                        # Number of atoms in first unit cell
           u1.1 u1.2 u1.3 [name1]        # Atom locations, in direct coord.
           ...                           #   with optional name for each atom
           uN.1 uN.2 uN.3 [nameN]
	   ==== cell ====
	   
	   ==== infile ====
	   t1 t2 t3     # dislocation line direction (unit cell coord.)
	   b1 b2 b3 bd  # burgers vector (unit cell coord.)/bd
	   m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)
	   c1 c2 c3 cd  # center of dislocation in unit cell ([c1 c2 c3]/cd)
	   c1' c2' c3'  # center of dislocation (shifts are added)
	   ==== infile ====

  Flags:   VERBOSE: 
	   TESTING: 

  Algo.:   Call to construct_slab after some initial setup.

*/

// ************************** COMPILIATION OPTIONS ***********************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include "dcomp.H"
#include "io-short.H"
#include "matrix.H"
#include "elastic.H"
#include "cell.H"
#include "slab.H"

// This is the permutation matrix; eps[i][j][k] =
//  1: if ijk is an even permutation of (012)
// -1: if ijk is an odd permutation of (012)
//  0: otherwise
const int eps[3][3][3] = {
  {{0,0,0}, {0,0,1}, {0,-1,0}},
  {{0,0,-1}, {0,0,0}, {1,0,0}},
  {{0,1,0}, {-1,0,0}, {0,0,0}}
};

// ****************************** SUBROUTINES ****************************

inline double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "[-hvt] [-a atomname] [-e] cell infile Rcut";

const char* ARGEXPL = 
"  cell:     cell file (-h for format)\n\
  infile:   input file (-h for format)\n\
  Rcut:     cutoff for building slab\n\
\n\
  -a atomname  replace atomnames in cell file (needed if names missing)\n\
  -e           assume all atom positions in cell file equivalent (with -a)\n\
  -v           verbosity\n\
  -t           testing\n\
  -h           help";

const char* FILEEXPL =
"==== cell ====\n\
a0                            # Scale factor for unit cell\n\
a1.x a1.y a1.z                # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                        # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3 [name1]        # Atom locations, in direct coord.\n\
...                           #   with optional name for each atom\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== infile ====\n\
t1 t2 t3     # dislocation line direction (unit cell coord.)\n\
b1 b2 b3 bd  # burgers vector (unit cell coord.)/bd \n\
m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)\n\
c1 c2 c3 cd  # center of dislocation in unit cell ([c1 c2 c3]/cd)\n\
c1' c2' c3'  # center of dislocation (shifts are added)\n\
==== infile ====\n\
\n\
==== undisloc ====\n\
N               # standard xyz format\n\
comment         # this *should* be the threading length\n\
atomtype x y z\n\
...\n\
==== undisloc ====\n";

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)

  char ch;
  char* atomname=NULL;
  int EQUIV = 0;
  while ((ch = getopt(argc, argv, "vthea:")) != -1) {
    switch (ch) {
    case 'a':
      atomname = new char[strlen(optarg)+1];
      strncpy(atomname, optarg, sizeof(atomname));
      break;
    case 'e':
      EQUIV = 1;
      break;
    case 'v':
      VERBOSE = 1;
      break;
    case 't':
      TESTING = 1;
      VERBOSE = 1;
      break;
    case 'h':
    case '?':
    default:
      ERROR = 1;
    }
  }
     
  argc -= optind; if (argc<NUMARGS && !ERROR) ERROR = 2;
  argv += optind;

  // argument compatibility check
  if (EQUIV && (atomname==NULL)) ERROR = 1;

  // All hell broken loose yet?
  if (ERROR != 0) {
    fprintf(stderr, "%s %s\n%s\n", progname, ARGLIST, ARGEXPL);
    if (ERROR == 1) {
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "Crystal classes:\n%s\n", CRYSTAL_CLASS);
      fprintf(stderr, "\nElastic constants ordering:\n");
      for (int k=0; k<NCLASSES; ++k) {
	fprintf(stderr, "  Class %2d (%d):", k, class_len[k]);
	for (int i=0; i<class_len[k]; ++i)
	  fprintf(stderr, " C_%2d", class_Cij[k][i]);
	fprintf(stderr, "\n");
      }
    }
    exit(ERROR);
  }


  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char *cell_name = argv[0];
  char *infile_name = argv[1];
  double Rcut;
  sscanf(argv[2], "%lf", &Rcut);

  double a0, cart[9];

  // First, read in the cell.
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(1);
  }
  // hacked from read_cell...
  // a0                            # Scale factor for unit cell
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%lf", &a0);

  // a1.x a1.y a1.z                # Cartesian coord of unit cell
  // a2.x a2.y a2.z
  // a3.x a3.y a3.z
  for (int i=0; i<3; ++i) {
    nextnoncomment(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", 
           cart+index(0,i), cart+index(1,i), cart+index(2,i));
  }
  for (int d=0; d<9; ++d) cart[d] *= a0;
  
  {
    // TEST: Determine determinants:
    double det_cart = det(cart);
    if (dcomp(det_cart, 0)) ERROR = ERROR_ZEROVOL;
    if (det_cart < 0)       ERROR = ERROR_LEFTHANDED;
  }

  // crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
  nextnoncomment(dump, sizeof(dump), infile);

  // Natoms                        # Number of atoms in first unit cell
  int Natoms;
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Natoms);
  // u1.1 u1.2 u1.3 [name1]        # Atom locations, in direct coord.
  // ...
  // uN.1 uN.2 uN.3
  double** u = new double*[Natoms];
  char** name = new char*[Natoms];
  
  for (int n=0; n<Natoms; ++n) {
    u[n] = new double[3];
    name[n] = new char[512];
    nextnoncomment(dump, sizeof(dump), infile);
    if (atomname==NULL)
      sscanf(dump, "%lf %lf %lf %s", u[n], u[n]+1, u[n]+2, name[n]);
    else {
      sscanf(dump, "%lf %lf %lf", u[n], u[n]+1, u[n]+2);
      if (EQUIV) strncpy(name[n], atomname, sizeof(atomname)+1);
      else sprintf(name[n], "%s.%d", atomname, n);
    }
    for (int d=0; d<3; ++d) u[n][d] = insidecell(u[n][d]);
  }

  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }

  if (TESTING) {
    printf("# Cartesian coordinates:\n");
    for (int d=0; d<3; ++d)
      printf("# a%1d = %8.5lf %8.5lf %8.5lf\n", d+1, 
	     cart[d], cart[d+3], cart[d+6]);

    printf("# Atoms in cell (%d):\n", Natoms);
    for (int n=0; n<Natoms; ++n)
      printf("# %s u%1d = %8.5lf %8.5lf %8.5lf\n", name[n],
	     n+1, u[n][0], u[n][1], u[n][2]);
  }


  // disl. line, burgers vect, cut, center of dislocation (all in unit coord)
  int tu0[3], bu0[3], mu0[3], cu0[3];	// all in unit cell coord; must be int.
  int bu_denom, cu_denom;		// denominator for burgers vector (partials)
  double cint0[3], c0[3];		// c0 will be the *true* center; cint0
					// is for converting cu0
  double t0[3], b0[3], m0[3], n0[3];	// n0 = t0 x m0, in cart. coord. 

  // Now, read in the dislocation information
  infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", infile_name);
    exit(1);
  }

  // **** NOTE: all input in unit cell coord, so first three vect. are int.
  //  t1 t2 t3            # dislocation line
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d", &tu0[0], &tu0[1], &tu0[2]);

  //  b1 b2 b3            # burgers vector
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d", &bu0[0], &bu0[1], &bu0[2], &bu_denom);
  // For backwards compatibility...
  if (bu_denom == 0) bu_denom = 1;

  //  m1 m2 m3            # dislocation cut vector (perp. to t)
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d", &mu0[0], &mu0[1], &mu0[2]);

  //  c1 c2 c3 cd           # center of dislocation
  //  c1' c2' c3'
  char *p_fgets = fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d", &cu0[0], &cu0[1], &cu0[2], &cu_denom);
  if (cu_denom == 0) cu_denom = 1;
  p_fgets = fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf %lf %lf", &cint0[0], &cint0[1], &cint0[2]);
  for (int d=0; d<3; ++d) cint0[d] += ((double) cu0[d])/((double) cu_denom);

  myclose(infile);

  // Now, convert vectors from unit cell to cartesian coord.:
  mult_vect(cart, tu0, t0);
  mult_vect(cart, bu0, b0); for (int d=0; d<3; ++d) b0[d] *= 1./bu_denom;
  mult_vect(cart, mu0, m0);
  mult_vect(cart, cint0, c0);

  // Sanity check on vectors:
  if ( dot(t0, t0) < 1e-8 ) {
    fprintf(stderr, "Bad t vector.\n");
    ERROR = ERROR_BADFILE;
  }
  if ( dot(b0, b0) < 1e-8 ) {
    fprintf(stderr, "Bad b vector.\n");
    ERROR = ERROR_BADFILE;
  }
  // We also need to project out any t components of m, and place
  // it in the slip plane (provided t x b isn't 0):
  for (int d=0; d<3; ++d)
    m0[d] -= dot(m0, t0)/dot(t0,t0) * t0[d];
  // Now, calculate n0 (we'll recalc it later, correctly)
  for (int d=0; d<3; ++d) {
    n0[d] = 0.;
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
	n0[d] += eps[d][j][k]*t0[j]*b0[k];
  }
  if (! dcomp(dot(n0,n0), 0.) ) 
    // We have a non-screw dislocation...
    for (int d=0; d<3; ++d)
      m0[d] -= dot(m0, n0)/dot(n0,n0) * n0[d];

  if ( dcomp(dot(m0, m0), 0.) ) {
    fprintf(stderr, "Bad m0 vector (parallel to t or out of the t x b slip plane).\n");
    ERROR = ERROR_BADFILE;
  }

  // Now, normalize:
  double magn;
  magn = 1./sqrt(dot(m0,m0));
  for (int d=0; d<3; ++d) m0[d] *= magn;

  
  if (VERBOSE) {
    printf("# Run dislocation along (%.5lf %.5lf %.5lf)\n",t0[0],t0[1],t0[2]);
    printf("# Burgers vector        (%.5lf %.5lf %.5lf), magn = %.5lf\n", 
	   b0[0],b0[1],b0[2], sqrt(dot(b0,b0)));
    printf("# Cut direction         (%.5lf %.5lf %.5lf)\n",m0[0],m0[1],m0[2]);

  }

  // ***************************** ANALYSIS **************************

  if (VERBOSE) {
    double comp;
    comp = fabs(dot(b0,t0)/sqrt(dot(b0,b0)*dot(t0,t0)));
    printf("# Screw component: %5.2lf%%  Edge component: %5.2lf%%\n",
	   comp*100.0, (1.-comp)*100.0);
  }
  
  // Now, compute n0 = t0 x m0:
  for (int d=0; d<3; ++d) {
    n0[d] = 0;
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
	n0[d] += eps[d][j][k] * t0[j] * m0[k];
  }
  // Normalize:
  magn = 1./sqrt(dot(n0,n0));
  for (int d=0; d<3; ++d) n0[d] *= magn;
  
  if (TESTING) {
    printf("##\n## Normalized vectors:\n");
    printf("## Run dislocation along (%.5lf %.5lf %.5lf)\n", t0[0],t0[1],t0[2]);
    printf("## Cut direction         (%.5lf %.5lf %.5lf)\n", m0[0],m0[1],m0[2]);
    printf("## Perp direction        (%.5lf %.5lf %.5lf)\n", n0[0],n0[1],n0[2]);
    printf("## Dislocation center    (%.5lf %.5lf %.5lf)\n", c0[0],c0[1],c0[2]);
  }
  if (VERBOSE) {
     printf("# %17.12lf %17.12lf %17.12lf : normalized x axis\n", m0[0], m0[1], m0[2]);
     printf("# %17.12lf %17.12lf %17.12lf : normalized y axis\n", n0[0], n0[1], n0[2]);
     printf("# %17.12lf %17.12lf %17.12lf : normalized z axis\n",
            t0[0]/sqrt(dot(t0,t0)), t0[1]/sqrt(dot(t0,t0)), t0[2]/sqrt(dot(t0,t0))); 
  }

  int Nslab;
  double** xyz=NULL;
  char** types=NULL;
  ERROR = construct_slab(t0, m0, n0, c0, Rcut, cart, u, name, Natoms, Nslab, xyz, types);
  
  // ****************************** OUTPUT ***************************

  // Output XYZ file
  printf("%d\n", Nslab);
  printf("%.15lf = z: undislocated slab, t = [%d %d %d], b = [%d %d %d]",
	 sqrt(dot(t0,t0)),
	 tu0[0], tu0[1], tu0[2], 
	 bu0[0], bu0[1], bu0[2]);
  if (bu_denom != 1) fprintf(infile, "/%d", bu_denom);
  printf(" Rmax = %.3lf\n", Rcut);
  for (int n=0; n<Nslab; ++n)
    printf("%s %20.15lf %20.15lf %20.15lf\n", types[n],
	    xyz[n][0], xyz[n][1], xyz[n][2]);

  // ************************* GARBAGE COLLECTION ********************
  free_slab(Nslab, xyz);
  delete[] types;
  for (int n=0; n<Natoms; ++n) delete[] name[n];
  delete[] name;
  for (int n=0; n<Natoms; ++n) delete[] u[n];
  delete[] u;

  return 0;
}
